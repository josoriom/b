use std::{io::Cursor, str};

use miniz_oxide::inflate::decompress_to_vec_zlib;
use zstd::{bulk::decompress as zstd_decompress, stream::decode_all as zstd_decode_all};

use crate::utilities::{attr_meta::*, cv_table, mzml::*};

const HEADER_SIZE: usize = 192;
const INDEX_ENTRY_SIZE: usize = 32;
const BLOCK_DIR_ENTRY_SIZE: usize = 32;

const ACC_MZ_ARRAY: u32 = 1_000_514;
const ACC_INTENSITY_ARRAY: u32 = 1_000_515;
const ACC_TIME_ARRAY: u32 = 1_000_595;

const ACC_32BIT_FLOAT: u32 = 1_000_521;
const ACC_64BIT_FLOAT: u32 = 1_000_523;

const ACC_ZLIB_COMPRESSION: u32 = 1_000_574;
const ACC_NO_COMPRESSION: u32 = 1_000_576;

const HDR_CODEC_MASK: u8 = 0x0F;
const HDR_CODEC_ZLIB: u8 = 0;
const HDR_CODEC_ZSTD: u8 = 1;

const HDR_FLAG_SPEC_META_COMP: u8 = 1 << 4;
const HDR_FLAG_CHROM_META_COMP: u8 = 1 << 5;
const HDR_FLAG_GLOBAL_META_COMP: u8 = 1 << 6;

const HDR_ARRAY_FILTER_OFF: usize = 178;
const ARRAY_FILTER_NONE: u8 = 0;
const ARRAY_FILTER_BYTE_SHUFFLE: u8 = 1;

const ACC_ISO_TARGET_MZ: u32 = 1_000_827;
const ACC_ISO_LOWER_OFFSET: u32 = 1_000_828;
const ACC_ISO_UPPER_OFFSET: u32 = 1_000_829;

const ACC_MZ: u32 = 1_000_040;
const ACC_CHARGE_STATE: u32 = 1_000_041;

const ACC_COLLISION_INDUCED_DISSOCIATION: u32 = 1_000_133;
const ACC_IN_SOURCE_CID: u32 = 1_001_880;
const ACC_COLLISION_ENERGY: u32 = 1_000_045;

const ACC_CENTROID_SPECTRUM: u32 = 1_000_127;
const ACC_LOWEST_MZ: u32 = 1_000_528;
const ACC_HIGHEST_MZ: u32 = 1_000_527;
const ACC_BASE_PEAK_MZ: u32 = 1_000_504;
const ACC_BASE_PEAK_INTENSITY: u32 = 1_000_505;
const ACC_TOTAL_ION_CURRENT: u32 = 1_000_285;

const ACC_SCAN_TIME: u32 = 1_000_016;
const ACC_FILTER_STRING: u32 = 1_000_512;

const ACC_SCAN_MZ_LOWER_LIMIT: u32 = 1_000_501;
const ACC_SCAN_MZ_UPPER_LIMIT: u32 = 1_000_500;

const ACC_XCALIBUR_RAW_FILE: u32 = 1_000_563;
const ACC_SHA1: u32 = 1_000_569;

const ACC_CONTACT_NAME: u32 = 1_000_586;
const ACC_CONTACT_ADDRESS: u32 = 1_000_587;
const ACC_CONTACT_URL: u32 = 1_000_588;
const ACC_CONTACT_EMAIL: u32 = 1_000_589;

const ACC_POSITIVE_SCAN: u32 = 1_000_130;
const ACC_FULL_SCAN: u32 = 1_000_498;

const ACC_PRODUCT_ION_MZ: u32 = 1_001_225;
const ACC_DWELL_TIME: u32 = 1_000_502;
const ACC_COMPLETION_TIME: u32 = 1_000_747;

#[derive(Clone, Copy)]
struct BlockDirEntry {
    comp_off: u64,
    comp_size: u64,
    uncomp_bytes: u64,
}

struct Container<'a> {
    compressed_region: &'a [u8],
    dir: Vec<BlockDirEntry>,
    block_start_elems: Vec<u64>,
    cache: Vec<Option<Vec<u8>>>,
    codec: u8,
    compression_level: u8,
    elem_size: usize,
    array_filter: u8,
    scratch: Vec<u8>,
}

impl<'a> Container<'a> {
    fn empty() -> Self {
        Self {
            compressed_region: &[],
            dir: Vec::new(),
            block_start_elems: vec![0],
            cache: Vec::new(),
            codec: HDR_CODEC_ZLIB,
            compression_level: 0,
            elem_size: 1,
            array_filter: ARRAY_FILTER_NONE,
            scratch: Vec::new(),
        }
    }

    fn new(
        file: &'a [u8],
        off: usize,
        size: usize,
        block_count: u32,
        codec: u8,
        compression_level: u8,
        elem_size: usize,
        array_filter: u8,
    ) -> Result<Self, String> {
        if size == 0 || block_count == 0 {
            return Ok(Self::empty());
        }
        if elem_size == 0 {
            return Err("Invalid elem_size".to_string());
        }

        let container_bytes = read_slice(file, off, size)?;
        let block_count = block_count as usize;

        let dir_bytes = block_count
            .checked_mul(BLOCK_DIR_ENTRY_SIZE)
            .ok_or_else(|| "Block directory size overflow".to_string())?;
        if dir_bytes > container_bytes.len() {
            return Err("Container too small for block directory".to_string());
        }

        let mut dir = Vec::with_capacity(block_count);
        let mut base = 0usize;
        for _ in 0..block_count {
            dir.push(BlockDirEntry {
                comp_off: read_u64_at(container_bytes, base)?,
                comp_size: read_u64_at(container_bytes, base + 8)?,
                uncomp_bytes: read_u64_at(container_bytes, base + 16)?,
            });
            base += BLOCK_DIR_ENTRY_SIZE;
        }

        let compressed_region = &container_bytes[dir_bytes..];

        let mut block_start_elems = Vec::with_capacity(block_count + 1);
        block_start_elems.push(0);

        let elem_size_u64 = elem_size as u64;
        let mut acc = 0u64;
        for e in &dir {
            acc = acc.saturating_add(e.uncomp_bytes / elem_size_u64);
            block_start_elems.push(acc);
        }

        Ok(Self {
            compressed_region,
            dir,
            block_start_elems,
            cache: vec![None; block_count],
            codec,
            compression_level,
            elem_size,
            array_filter,
            scratch: Vec::new(),
        })
    }

    #[inline]
    fn block_count(&self) -> usize {
        self.dir.len()
    }

    fn block_bytes(&mut self, block_id: u32) -> Result<&[u8], String> {
        let id = block_id as usize;
        if id >= self.block_count() {
            return Err("Invalid block id".to_string());
        }

        let e = self.dir[id];
        let comp_off = e.comp_off as usize;
        let comp_size = e.comp_size as usize;
        let end = comp_off
            .checked_add(comp_size)
            .ok_or_else(|| "Block range overflow".to_string())?;

        let needs_owned = self.compression_level != 0
            || (self.array_filter == ARRAY_FILTER_BYTE_SHUFFLE && self.elem_size > 1);

        if !needs_owned {
            return self
                .compressed_region
                .get(comp_off..end)
                .ok_or_else(|| "EOF".to_string());
        }

        if self.cache[id].is_none() {
            let comp = self
                .compressed_region
                .get(comp_off..end)
                .ok_or_else(|| "EOF".to_string())?;

            let mut block = if self.compression_level == 0 {
                if e.uncomp_bytes != 0 && comp.len() != e.uncomp_bytes as usize {
                    return Err("Uncompressed block size mismatch".to_string());
                }
                comp.to_vec()
            } else {
                let inflated = match self.codec {
                    HDR_CODEC_ZLIB => decompress_to_vec_zlib(comp)
                        .map_err(|_| "Zlib decompression failed".to_string())?,
                    HDR_CODEC_ZSTD => zstd_decompress(comp, e.uncomp_bytes as usize)
                        .map_err(|_| "Zstd decompression failed".to_string())?,
                    _ => return Err("Unsupported container codec".to_string()),
                };

                if e.uncomp_bytes != 0 && inflated.len() != e.uncomp_bytes as usize {
                    return Err("Inflated block size mismatch".to_string());
                }

                inflated
            };

            if self.array_filter == ARRAY_FILTER_BYTE_SHUFFLE
                && self.elem_size > 1
                && !block.is_empty()
            {
                self.scratch.resize(block.len(), 0);
                unshuffle_into(&mut self.scratch, &block, self.elem_size)?;
                std::mem::swap(&mut block, &mut self.scratch);
            }

            self.cache[id] = Some(block);
        }

        Ok(self.cache[id].as_deref().unwrap_or(&[]))
    }

    fn slice_elems(
        &mut self,
        block_id: u32,
        global_elem_off: u64,
        elem_len: u32,
    ) -> Result<&[u8], String> {
        let id = block_id as usize;
        if id + 1 >= self.block_start_elems.len() {
            return Err("Invalid block id".to_string());
        }

        let block_start = self.block_start_elems[id];
        if global_elem_off < block_start {
            return Err("Element offset before block start".to_string());
        }

        let local_elems = (global_elem_off - block_start) as usize;

        let byte_off = local_elems
            .checked_mul(self.elem_size)
            .ok_or_else(|| "Byte offset overflow".to_string())?;
        let byte_len = (elem_len as usize)
            .checked_mul(self.elem_size)
            .ok_or_else(|| "Byte length overflow".to_string())?;
        let end = byte_off
            .checked_add(byte_len)
            .ok_or_else(|| "Slice range overflow".to_string())?;

        self.block_bytes(block_id)?
            .get(byte_off..end)
            .ok_or_else(|| "EOF".to_string())
    }
}

#[inline]
fn unshuffle_into(dst: &mut [u8], src: &[u8], elem_size: usize) -> Result<(), String> {
    if dst.len() != src.len() {
        return Err("unshuffle size mismatch".to_string());
    }
    if elem_size <= 1 {
        dst.copy_from_slice(src);
        return Ok(());
    }
    if src.len() % elem_size != 0 {
        return Err("unshuffle: invalid byte length".to_string());
    }

    let n = src.len() / elem_size;
    for b in 0..elem_size {
        let col = b
            .checked_mul(n)
            .ok_or_else(|| "unshuffle overflow".to_string())?;
        for i in 0..n {
            dst[i * elem_size + b] = src[col + i];
        }
    }
    Ok(())
}

enum BytesMaybeOwned<'a> {
    Borrowed(&'a [u8]),
    Owned(Vec<u8>),
}

impl<'a> BytesMaybeOwned<'a> {
    #[inline]
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Borrowed(b) => b,
            Self::Owned(v) => v.as_slice(),
        }
    }
}

#[inline]
fn decompress_zlib_allow_pad0(input: &[u8]) -> Result<Vec<u8>, String> {
    if let Ok(v) = decompress_to_vec_zlib(input) {
        return Ok(v);
    }
    let mut end = input.len();
    for _ in 0..7 {
        if end == 0 || input[end - 1] != 0 {
            break;
        }
        end -= 1;
        if let Ok(v) = decompress_to_vec_zlib(&input[..end]) {
            return Ok(v);
        }
    }
    Err("Zlib decompression failed".to_string())
}

#[inline]
fn decompress_zstd_allow_pad0(input: &[u8]) -> Result<Vec<u8>, String> {
    if let Ok(v) = zstd_decode_all(Cursor::new(input)) {
        return Ok(v);
    }
    let mut end = input.len();
    for _ in 0..7 {
        if end == 0 || input[end - 1] != 0 {
            break;
        }
        end -= 1;
        if let Ok(v) = zstd_decode_all(Cursor::new(&input[..end])) {
            return Ok(v);
        }
    }
    Err("Zstd decompression failed".to_string())
}

fn decompress_meta_if_needed<'a>(
    codec: u8,
    is_compressed: bool,
    bytes: &'a [u8],
) -> Result<BytesMaybeOwned<'a>, String> {
    if !is_compressed {
        return Ok(BytesMaybeOwned::Borrowed(bytes));
    }
    match codec {
        HDR_CODEC_ZLIB => Ok(BytesMaybeOwned::Owned(decompress_zlib_allow_pad0(bytes)?)),
        HDR_CODEC_ZSTD => Ok(BytesMaybeOwned::Owned(decompress_zstd_allow_pad0(bytes)?)),
        _ => Err("Unsupported meta codec".to_string()),
    }
}

#[inline]
fn is_isolation_window_tail(t: u32) -> bool {
    matches!(
        t,
        ACC_ISO_TARGET_MZ | ACC_ISO_LOWER_OFFSET | ACC_ISO_UPPER_OFFSET
    )
}

#[inline]
fn is_selected_ion_tail(t: u32) -> bool {
    matches!(t, ACC_MZ | ACC_CHARGE_STATE)
}

#[inline]
fn is_activation_tail(t: u32) -> bool {
    matches!(
        t,
        ACC_IN_SOURCE_CID | ACC_COLLISION_INDUCED_DISSOCIATION | ACC_COLLISION_ENERGY
    )
}

#[inline]
fn is_spectrum_description_tail(t: u32) -> bool {
    matches!(
        t,
        ACC_CENTROID_SPECTRUM
            | ACC_LOWEST_MZ
            | ACC_HIGHEST_MZ
            | ACC_BASE_PEAK_MZ
            | ACC_BASE_PEAK_INTENSITY
            | ACC_TOTAL_ION_CURRENT
    )
}

#[inline]
fn is_scan_tail(t: u32) -> bool {
    matches!(t, ACC_SCAN_TIME | ACC_FILTER_STRING)
}

#[inline]
fn is_scan_window_tail(t: u32) -> bool {
    matches!(t, ACC_SCAN_MZ_LOWER_LIMIT | ACC_SCAN_MZ_UPPER_LIMIT)
}

#[inline]
fn is_source_file_tail(t: u32) -> bool {
    matches!(t, ACC_XCALIBUR_RAW_FILE | ACC_SHA1)
}

#[inline]
fn is_contact_tail(t: u32) -> bool {
    matches!(
        t,
        ACC_CONTACT_NAME | ACC_CONTACT_ADDRESS | ACC_CONTACT_URL | ACC_CONTACT_EMAIL
    )
}

#[inline]
fn base64_encoded_len(byte_len: usize) -> usize {
    ((byte_len + 2) / 3) * 4
}

#[inline]
fn float_acc_from_fmt(fmt: u8) -> u32 {
    if fmt == 1 {
        ACC_32BIT_FLOAT
    } else {
        ACC_64BIT_FLOAT
    }
}

#[inline]
fn is_attr_param(p: &CvParam) -> bool {
    if p.cv_ref.as_deref() == Some(CV_REF_ATTR) {
        return true;
    }
    match p.accession.as_deref() {
        Some(a) => {
            let pref = CV_REF_ATTR.as_bytes();
            let b = a.as_bytes();
            b.starts_with(pref) && b.get(pref.len()) == Some(&b':')
        }
        None => false,
    }
}

#[inline]
fn attr_string_value(p: &CvParam) -> Option<String> {
    if let Some(v) = p.value.as_ref() {
        if !v.is_empty() {
            return Some(v.clone());
        }
    }
    (is_attr_param(p) && !p.name.is_empty()).then(|| p.name.clone())
}

/// <precursorList>
fn infer_precursor_list_from_spectrum_cv(params: &mut Vec<CvParam>) -> Option<PrecursorList> {
    let mut spectrum_ref: Option<String> = None;
    let mut iso = Vec::<CvParam>::new();
    let mut sel = Vec::<CvParam>::new();
    let mut act = Vec::<CvParam>::new();

    let mut rest = Vec::<CvParam>::with_capacity(params.len());
    for p in params.drain(..) {
        let tail = parse_acc_tail(p.accession.as_deref());

        if is_attr_param(&p) && tail == ACC_ATTR_SPECTRUM_REF {
            if spectrum_ref.is_none() {
                spectrum_ref = attr_string_value(&p);
            }
            continue;
        }

        if is_isolation_window_tail(tail) {
            iso.push(p);
        } else if is_selected_ion_tail(tail) {
            sel.push(p);
        } else if is_activation_tail(tail) {
            act.push(p);
        } else {
            rest.push(p);
        }
    }
    *params = rest;

    if iso.is_empty() && sel.is_empty() && act.is_empty() && spectrum_ref.is_none() {
        return None;
    }

    let isolation_window = (!iso.is_empty()).then(|| IsolationWindow {
        cv_params: iso,
        ..Default::default()
    });

    let selected_ion_list = (!sel.is_empty()).then(|| SelectedIonList {
        count: Some(1),
        selected_ions: vec![SelectedIon {
            cv_params: sel,
            ..Default::default()
        }],
    });

    let activation = (!act.is_empty()).then(|| Activation {
        cv_params: act,
        ..Default::default()
    });

    Some(PrecursorList {
        count: Some(1),
        precursors: vec![Precursor {
            spectrum_ref,
            isolation_window,
            selected_ion_list,
            activation,
            ..Default::default()
        }],
    })
}

/// <spectrumDescription>
fn infer_spectrum_description_from_spectrum_cv(
    params: &mut Vec<CvParam>,
) -> Option<SpectrumDescription> {
    let mut desc = Vec::<CvParam>::new();
    let mut scan = Vec::<CvParam>::new();
    let mut scan_window = Vec::<CvParam>::new();

    let mut rest = Vec::<CvParam>::with_capacity(params.len());
    for p in params.drain(..) {
        let tail = parse_acc_tail(p.accession.as_deref());
        if is_spectrum_description_tail(tail) {
            desc.push(p);
        } else if is_scan_tail(tail) {
            scan.push(p);
        } else if is_scan_window_tail(tail) {
            scan_window.push(p);
        } else {
            rest.push(p);
        }
    }
    *params = rest;

    if desc.is_empty() && scan.is_empty() && scan_window.is_empty() {
        return None;
    }

    let scan_window_list = (!scan_window.is_empty()).then(|| ScanWindowList {
        count: Some(1),
        scan_windows: vec![ScanWindow {
            cv_params: scan_window,
            ..Default::default()
        }],
    });

    let scan_obj = if scan.is_empty() && scan_window_list.is_none() {
        None
    } else {
        Some(Scan {
            cv_params: scan,
            scan_window_list,
            ..Default::default()
        })
    };

    Some(SpectrumDescription {
        cv_params: desc,
        scan_list: scan_obj.map(|s| ScanList {
            count: Some(1),
            scans: vec![s],
        }),
        ..Default::default()
    })
}

/// <fileDescription>
fn split_file_description_from_cv_params(mut params: Vec<CvParam>) -> FileDescription {
    let mut file_content_cv = Vec::<CvParam>::new();
    let mut contact_cv = Vec::<CvParam>::new();

    let mut source_files = Vec::<SourceFile>::new();
    let mut cur_id = String::new();
    let mut cur_name = String::new();
    let mut cur_location = String::new();
    let mut cur_cv = Vec::<CvParam>::new();

    #[inline]
    fn flush_source_file(
        out: &mut Vec<SourceFile>,
        id: &mut String,
        name: &mut String,
        location: &mut String,
        cv: &mut Vec<CvParam>,
    ) {
        if id.is_empty() && name.is_empty() && location.is_empty() && cv.is_empty() {
            return;
        }
        let mut sf = SourceFile::default();
        sf.id = std::mem::take(id);
        sf.name = std::mem::take(name);
        sf.location = std::mem::take(location);
        sf.cv_param = std::mem::take(cv);
        out.push(sf);
    }

    for p in params.drain(..) {
        let tail = parse_acc_tail(p.accession.as_deref());

        if is_attr_param(&p) {
            let v = attr_string_value(&p).unwrap_or_default();
            if v.is_empty() {
                continue;
            }

            if tail == ACC_ATTR_ID {
                if !cur_id.is_empty() && cur_id != v {
                    flush_source_file(
                        &mut source_files,
                        &mut cur_id,
                        &mut cur_name,
                        &mut cur_location,
                        &mut cur_cv,
                    );
                }
                cur_id = v;
            } else if tail == ACC_ATTR_NAME {
                cur_name = v;
            } else if tail == ACC_ATTR_LOCATION {
                cur_location = v;
            }
            continue;
        }

        if is_source_file_tail(tail) {
            cur_cv.push(p);
        } else if is_contact_tail(tail) {
            contact_cv.push(p);
        } else {
            file_content_cv.push(p);
        }
    }

    flush_source_file(
        &mut source_files,
        &mut cur_id,
        &mut cur_name,
        &mut cur_location,
        &mut cur_cv,
    );

    let mut fd = FileDescription::default();
    fd.file_content.cv_params = file_content_cv;

    if !source_files.is_empty() {
        fd.source_file_list.count = Some(source_files.len());
        fd.source_file_list.source_file = source_files;
    }

    if !contact_cv.is_empty() {
        let mut c = Contact::default();
        c.cv_params = contact_cv;
        fd.contacts.push(c);
    }

    fd
}

/// <spectrum>
fn split_spectrum_attrs(
    item_idx: usize,
    x_len: u32,
    params: Vec<CvParam>,
) -> (String, Option<u32>, Option<usize>, Vec<CvParam>) {
    split_item_attrs(item_idx, x_len, params, "spectrum")
}

/// <chromatogram>
fn split_chromatogram_attrs(
    item_idx: usize,
    x_len: u32,
    params: Vec<CvParam>,
) -> (String, Option<u32>, Option<usize>, Vec<CvParam>) {
    split_item_attrs(item_idx, x_len, params, "chromatogram")
}

#[inline]
fn split_item_attrs(
    item_idx: usize,
    x_len: u32,
    params: Vec<CvParam>,
    prefix: &str,
) -> (String, Option<u32>, Option<usize>, Vec<CvParam>) {
    let mut id: Option<String> = None;
    let mut index: Option<u32> = None;
    let mut default_array_length: Option<usize> = None;

    let mut out = Vec::with_capacity(params.len());
    for p in params {
        if is_attr_param(&p) {
            let tail = parse_acc_tail(p.accession.as_deref());
            if tail == ACC_ATTR_ID {
                id = attr_string_value(&p);
                continue;
            }
            if tail == ACC_ATTR_INDEX {
                index = Some(
                    attr_string_value(&p)
                        .and_then(|s| s.parse::<u32>().ok())
                        .unwrap_or(item_idx as u32),
                );
                continue;
            }
            if tail == ACC_ATTR_DEFAULT_ARRAY_LENGTH {
                default_array_length = Some(
                    attr_string_value(&p)
                        .and_then(|s| s.parse::<usize>().ok())
                        .unwrap_or(x_len as usize),
                );
                continue;
            }
        }
        out.push(p);
    }

    (
        id.unwrap_or_else(|| format!("{}_{}", prefix, item_idx)),
        Some(index.unwrap_or(item_idx as u32)),
        Some(default_array_length.unwrap_or(x_len as usize)),
        out,
    )
}

#[inline]
fn filter_spectrum_top_level_cv_params(params: &mut Vec<CvParam>) {
    params.retain(|p| {
        let tail = parse_acc_tail(p.accession.as_deref());
        if tail == ACC_POSITIVE_SCAN || tail == ACC_FULL_SCAN {
            return false;
        }
        !(p.cv_ref.as_deref() == Some("MS") && p.name.is_empty())
    });
}

#[inline]
fn make_binary_data_array(
    array_len: u32,
    elem_size: usize,
    fmt: u8,
    array_accession_tail: u32,
    decoded_f32: Vec<f32>,
    decoded_f64: Vec<f64>,
) -> BinaryDataArray {
    let byte_len = (array_len as usize) * elem_size;
    let enc_len = base64_encoded_len(byte_len);

    let mut ba = BinaryDataArray::default();
    ba.array_length = Some(array_len as usize);
    ba.encoded_length = Some(enc_len);
    ba.is_f32 = Some(fmt == 1);
    ba.is_f64 = Some(fmt == 2);
    ba.cv_params = vec![
        ms_cv_param(float_acc_from_fmt(fmt)),
        ms_cv_param(ACC_NO_COMPRESSION),
        ms_cv_param(array_accession_tail),
    ];
    ba.decoded_binary_f32 = decoded_f32;
    ba.decoded_binary_f64 = decoded_f64;
    ba
}

/// <mzML>
pub fn decode(bytes: &[u8]) -> Result<MzML, String> {
    if bytes.len() < HEADER_SIZE {
        return Err("Buffer too small for header".to_string());
    }
    let header = &bytes[..HEADER_SIZE];

    if &header[..4] != b"B000" {
        return Err("Invalid binary magic number".to_string());
    }
    if read_u8_at(header, 4)? != 0 {
        return Err("Unsupported endianness flag".to_string());
    }

    let off_spec_index = read_u64_at(header, 8)? as usize;
    let off_chrom_index = read_u64_at(header, 16)? as usize;
    let off_spec_meta = read_u64_at(header, 24)? as usize;
    let off_chrom_meta = read_u64_at(header, 32)? as usize;
    let off_global_meta = read_u64_at(header, 40)? as usize;

    let size_container_spect_x = read_u64_at(header, 48)? as usize;
    let off_container_spect_x = read_u64_at(header, 56)? as usize;
    let size_container_spect_y = read_u64_at(header, 64)? as usize;
    let off_container_spect_y = read_u64_at(header, 72)? as usize;
    let size_container_chrom_x = read_u64_at(header, 80)? as usize;
    let off_container_chrom_x = read_u64_at(header, 88)? as usize;
    let size_container_chrom_y = read_u64_at(header, 96)? as usize;
    let off_container_chrom_y = read_u64_at(header, 104)? as usize;

    let spectrum_count = read_u32_at(header, 112)?;
    let chrom_count = read_u32_at(header, 116)?;

    let spec_meta_count = read_u32_at(header, 120)?;
    let spec_num_count = read_u32_at(header, 124)?;
    let spec_str_count = read_u32_at(header, 128)?;

    let chrom_meta_count = read_u32_at(header, 132)?;
    let chrom_num_count = read_u32_at(header, 136)?;
    let chrom_str_count = read_u32_at(header, 140)?;

    let global_meta_count = read_u32_at(header, 144)?;
    let global_num_count = read_u32_at(header, 148)?;
    let global_str_count = read_u32_at(header, 152)?;

    let block_count_spect_x = read_u32_at(header, 156)?;
    let block_count_spect_y = read_u32_at(header, 160)?;
    let block_count_chrom_x = read_u32_at(header, 164)?;
    let block_count_chrom_y = read_u32_at(header, 168)?;

    let codec_flags = read_u8_at(header, 172)?;
    let codec = codec_flags & HDR_CODEC_MASK;

    let chrom_x_fmt = read_u8_at(header, 173)?;
    let chrom_y_fmt = read_u8_at(header, 174)?;
    let spect_x_fmt = read_u8_at(header, 175)?;
    let spect_y_fmt = read_u8_at(header, 176)?;
    let compression_level = read_u8_at(header, 177)?;
    let array_filter = read_u8_at(header, HDR_ARRAY_FILTER_OFF)?;

    let spect_x_elem_size = fmt_elem_size(spect_x_fmt)?;
    let spect_y_elem_size = fmt_elem_size(spect_y_fmt)?;
    let chrom_x_elem_size = fmt_elem_size(chrom_x_fmt)?;
    let chrom_y_elem_size = fmt_elem_size(chrom_y_fmt)?;

    let spectrum_index_len = (spectrum_count as usize)
        .checked_mul(INDEX_ENTRY_SIZE)
        .ok_or_else(|| "Index overflow".to_string())?;
    let chromatogram_index_len = (chrom_count as usize)
        .checked_mul(INDEX_ENTRY_SIZE)
        .ok_or_else(|| "Index overflow".to_string())?;

    let spectrum_index_bytes = read_slice(bytes, off_spec_index, spectrum_index_len)?;
    let chromatogram_index_bytes = read_slice(bytes, off_chrom_index, chromatogram_index_len)?;

    if off_chrom_meta < off_spec_meta || off_global_meta < off_chrom_meta {
        return Err("Invalid meta offsets".to_string());
    }

    let spec_meta_region = read_slice(bytes, off_spec_meta, off_chrom_meta - off_spec_meta)?;
    let chrom_meta_region = read_slice(bytes, off_chrom_meta, off_global_meta - off_chrom_meta)?;

    let first_container_off = min_nonzero_usize(&[
        off_container_spect_x,
        off_container_spect_y,
        off_container_chrom_x,
        off_container_chrom_y,
    ])
    .unwrap_or(bytes.len());

    if first_container_off < off_global_meta {
        return Err("Invalid global meta/container offsets".to_string());
    }

    let global_meta_region = read_slice(
        bytes,
        off_global_meta,
        first_container_off - off_global_meta,
    )?;

    let spec_meta_bytes = decompress_meta_if_needed(
        codec,
        (codec_flags & HDR_FLAG_SPEC_META_COMP) != 0,
        spec_meta_region,
    )?;
    let chrom_meta_bytes = decompress_meta_if_needed(
        codec,
        (codec_flags & HDR_FLAG_CHROM_META_COMP) != 0,
        chrom_meta_region,
    )?;
    let global_meta_bytes = decompress_meta_if_needed(
        codec,
        (codec_flags & HDR_FLAG_GLOBAL_META_COMP) != 0,
        global_meta_region,
    )?;

    let mut spect_x_container = Container::new(
        bytes,
        off_container_spect_x,
        size_container_spect_x,
        block_count_spect_x,
        codec,
        compression_level,
        spect_x_elem_size,
        array_filter,
    )?;
    let mut spect_y_container = Container::new(
        bytes,
        off_container_spect_y,
        size_container_spect_y,
        block_count_spect_y,
        codec,
        compression_level,
        spect_y_elem_size,
        array_filter,
    )?;
    let mut chrom_x_container = Container::new(
        bytes,
        off_container_chrom_x,
        size_container_chrom_x,
        block_count_chrom_x,
        codec,
        compression_level,
        chrom_x_elem_size,
        array_filter,
    )?;
    let mut chrom_y_container = Container::new(
        bytes,
        off_container_chrom_y,
        size_container_chrom_y,
        block_count_chrom_y,
        codec,
        compression_level,
        chrom_y_elem_size,
        array_filter,
    )?;

    let spec_meta_by_item = decode_meta_block(
        spec_meta_bytes.as_slice(),
        spectrum_count,
        spec_meta_count,
        spec_num_count,
        spec_str_count,
    )?;
    let chrom_meta_by_item = decode_meta_block(
        chrom_meta_bytes.as_slice(),
        chrom_count,
        chrom_meta_count,
        chrom_num_count,
        chrom_str_count,
    )?;

    let (
        cv_list,
        file_description,
        referenceable_param_group_list,
        sample_list,
        instrument_list,
        software_list,
        data_processing_list,
        scan_settings_list,
    ) = decode_global_meta_structs(
        global_meta_bytes.as_slice(),
        global_meta_count,
        global_num_count,
        global_str_count,
    )?;

    let mut spectra = Vec::with_capacity(spectrum_count as usize);
    for (i, item_params) in spec_meta_by_item.into_iter().enumerate() {
        let (x_off, y_off, x_len, y_len, x_block, y_block) =
            read_index_entry_with_blocks(spectrum_index_bytes, i)?;

        let mz_bytes = spect_x_container.slice_elems(x_block, x_off, x_len)?;
        let in_bytes = spect_y_container.slice_elems(y_block, y_off, y_len)?;

        let (mz_f32, mz_f64) = decode_array_by_fmt_from_bytes(mz_bytes, spect_x_fmt)?;
        let (in_f32, in_f64) = decode_array_by_fmt_from_bytes(in_bytes, spect_y_fmt)?;

        let mz_ba = make_binary_data_array(
            x_len,
            spect_x_elem_size,
            spect_x_fmt,
            ACC_MZ_ARRAY,
            mz_f32,
            mz_f64,
        );
        let inten_ba = make_binary_data_array(
            y_len,
            spect_y_elem_size,
            spect_y_fmt,
            ACC_INTENSITY_ARRAY,
            in_f32,
            in_f64,
        );

        let mut spectrum_params = item_params;
        strip_binary_array_cv_params(&mut spectrum_params);

        let (id, index, default_array_length, mut spectrum_params) =
            split_spectrum_attrs(i, x_len, spectrum_params);

        let precursor_list = infer_precursor_list_from_spectrum_cv(&mut spectrum_params);
        let spectrum_description =
            infer_spectrum_description_from_spectrum_cv(&mut spectrum_params);

        filter_spectrum_top_level_cv_params(&mut spectrum_params);

        spectra.push(Spectrum {
            id,
            index,
            default_array_length,
            cv_params: spectrum_params,
            spectrum_description,
            precursor_list,
            binary_data_array_list: Some(BinaryDataArrayList {
                count: Some(2),
                binary_data_arrays: vec![mz_ba, inten_ba],
            }),
            ..Default::default()
        });
    }

    let mut chromatograms = Vec::with_capacity(chrom_count as usize);
    for (j, item_params) in chrom_meta_by_item.into_iter().enumerate() {
        let (x_off, y_off, x_len, y_len, x_block, y_block) =
            read_index_entry_with_blocks(chromatogram_index_bytes, j)?;

        let t_bytes = chrom_x_container.slice_elems(x_block, x_off, x_len)?;
        let in_bytes = chrom_y_container.slice_elems(y_block, y_off, y_len)?;

        let (t_f32, t_f64) = decode_array_by_fmt_from_bytes(t_bytes, chrom_x_fmt)?;
        let (in_f32, in_f64) = decode_array_by_fmt_from_bytes(in_bytes, chrom_y_fmt)?;

        let time_ba = make_binary_data_array(
            x_len,
            chrom_x_elem_size,
            chrom_x_fmt,
            ACC_TIME_ARRAY,
            t_f32,
            t_f64,
        );
        let inten_ba = make_binary_data_array(
            y_len,
            chrom_y_elem_size,
            chrom_y_fmt,
            ACC_INTENSITY_ARRAY,
            in_f32,
            in_f64,
        );

        let mut chrom_params = item_params;
        strip_binary_array_cv_params(&mut chrom_params);

        let (id, index, default_array_length, chrom_params) =
            split_chromatogram_attrs(j, x_len, chrom_params);

        chromatograms.push(Chromatogram {
            id,
            index,
            default_array_length,
            cv_params: chrom_params,
            binary_data_array_list: Some(BinaryDataArrayList {
                count: Some(2),
                binary_data_arrays: vec![time_ba, inten_ba],
            }),
            ..Default::default()
        });
    }

    Ok(MzML {
        cv_list,
        file_description,
        referenceable_param_group_list,
        sample_list,
        instrument_list,
        software_list,
        data_processing_list,
        scan_settings_list,
        run: Run {
            id: "run".to_string(),
            spectrum_list: Some(SpectrumList {
                count: Some(spectrum_count as usize),
                spectra,
                ..Default::default()
            }),
            chromatogram_list: Some(ChromatogramList {
                count: Some(chrom_count as usize),
                chromatograms,
                ..Default::default()
            }),
            ..Default::default()
        },
    })
}

#[inline]
fn min_nonzero_usize(xs: &[usize]) -> Option<usize> {
    let mut m: Option<usize> = None;
    for &x in xs {
        if x == 0 {
            continue;
        }
        m = Some(m.map_or(x, |cur| cur.min(x)));
    }
    m
}

#[inline]
fn fmt_elem_size(fmt: u8) -> Result<usize, String> {
    match fmt {
        1 => Ok(4),
        2 => Ok(8),
        _ => Err("Invalid float format".to_string()),
    }
}

#[inline]
fn strip_binary_array_cv_params(params: &mut Vec<CvParam>) {
    params.retain(|p| {
        let tail = parse_acc_tail(p.accession.as_deref());
        !matches!(
            tail,
            ACC_MZ_ARRAY
                | ACC_INTENSITY_ARRAY
                | ACC_TIME_ARRAY
                | ACC_32BIT_FLOAT
                | ACC_64BIT_FLOAT
                | ACC_ZLIB_COMPRESSION
                | ACC_NO_COMPRESSION
        )
    });
}

fn read_index_entry_with_blocks(
    index_bytes: &[u8],
    item_idx: usize,
) -> Result<(u64, u64, u32, u32, u32, u32), String> {
    let base = item_idx
        .checked_mul(INDEX_ENTRY_SIZE)
        .ok_or_else(|| "Index overflow".to_string())?;
    let end = base
        .checked_add(INDEX_ENTRY_SIZE)
        .ok_or_else(|| "Index overflow".to_string())?;
    if end > index_bytes.len() {
        return Err("Index overflow".to_string());
    }

    let x_off = u64::from_le_bytes(index_bytes[base..base + 8].try_into().unwrap());
    let y_off = u64::from_le_bytes(index_bytes[base + 8..base + 16].try_into().unwrap());
    let x_len = u32::from_le_bytes(index_bytes[base + 16..base + 20].try_into().unwrap());
    let y_len = u32::from_le_bytes(index_bytes[base + 20..base + 24].try_into().unwrap());
    let x_block = u32::from_le_bytes(index_bytes[base + 24..base + 28].try_into().unwrap());
    let y_block = u32::from_le_bytes(index_bytes[base + 28..base + 32].try_into().unwrap());

    Ok((x_off, y_off, x_len, y_len, x_block, y_block))
}

#[inline]
fn decode_array_by_fmt_from_bytes(bytes: &[u8], fmt: u8) -> Result<(Vec<f32>, Vec<f64>), String> {
    match fmt {
        1 => Ok((bytes_to_f32_exact(bytes)?, Vec::new())),
        2 => Ok((Vec::new(), bytes_to_f64_exact(bytes)?)),
        _ => Err("Invalid float format".to_string()),
    }
}

fn bytes_to_f64_exact(bytes: &[u8]) -> Result<Vec<f64>, String> {
    if bytes.len() % 8 != 0 {
        return Err("Invalid f64 byte length".to_string());
    }
    let n = bytes.len() / 8;

    if cfg!(target_endian = "little") {
        let mut out: Vec<f64> = Vec::with_capacity(n);
        unsafe {
            out.set_len(n);
            std::ptr::copy_nonoverlapping(bytes.as_ptr(), out.as_mut_ptr() as *mut u8, bytes.len());
        }
        return Ok(out);
    }

    let mut out = Vec::with_capacity(n);
    for c in bytes.chunks_exact(8) {
        out.push(f64::from_le_bytes(c.try_into().unwrap()));
    }
    Ok(out)
}

fn bytes_to_f32_exact(bytes: &[u8]) -> Result<Vec<f32>, String> {
    if bytes.len() % 4 != 0 {
        return Err("Invalid f32 byte length".to_string());
    }
    let n = bytes.len() / 4;

    if cfg!(target_endian = "little") {
        let mut out: Vec<f32> = Vec::with_capacity(n);
        unsafe {
            out.set_len(n);
            std::ptr::copy_nonoverlapping(bytes.as_ptr(), out.as_mut_ptr() as *mut u8, bytes.len());
        }
        return Ok(out);
    }

    let mut out = Vec::with_capacity(n);
    for c in bytes.chunks_exact(4) {
        out.push(f32::from_le_bytes(c.try_into().unwrap()));
    }
    Ok(out)
}

#[inline]
fn cv_table_name(key: &str) -> Option<String> {
    cv_table::get(key)
        .and_then(|v| {
            v.as_str()
                .or_else(|| v.get("name").and_then(|n| n.as_str()))
        })
        .map(|s| s.to_string())
}

/// <cvParam>
#[inline]
fn ms_cv_param(accession_tail: u32) -> CvParam {
    let key = make_accession(Some("MS"), accession_tail)
        .unwrap_or_else(|| format!("MS:{:07}", accession_tail));
    let name = cv_table_name(&key).unwrap_or_default();
    CvParam {
        cv_ref: Some("MS".to_string()),
        accession: Some(key),
        name,
        value: Some(String::new()),
        ..Default::default()
    }
}

/// <cvParam>
fn decode_meta_block(
    bytes: &[u8],
    item_count: u32,
    meta_count: u32,
    num_count: u32,
    str_count: u32,
) -> Result<Vec<Vec<CvParam>>, String> {
    let mut offset = 0usize;

    let item_count = item_count as usize;
    let meta_count = meta_count as usize;
    let num_count = num_count as usize;
    let str_count = str_count as usize;

    let item_indices = read_u32_vec(
        read_slice(bytes, offset, (item_count + 1) * 4)?,
        item_count + 1,
    )?;
    offset += (item_count + 1) * 4;

    let meta_ref_codes = read_slice(bytes, offset, meta_count)?;
    offset += meta_count;

    let meta_accessions = read_u32_vec(read_slice(bytes, offset, meta_count * 4)?, meta_count)?;
    offset += meta_count * 4;

    let meta_unit_refs = read_slice(bytes, offset, meta_count)?;
    offset += meta_count;

    let meta_unit_accessions =
        read_u32_vec(read_slice(bytes, offset, meta_count * 4)?, meta_count)?;
    offset += meta_count * 4;

    let value_kinds = read_slice(bytes, offset, meta_count)?;
    offset += meta_count;

    let value_indices = read_u32_vec(read_slice(bytes, offset, meta_count * 4)?, meta_count)?;
    offset += meta_count * 4;

    let numeric_values = read_f64_vec(read_slice(bytes, offset, num_count * 8)?, num_count)?;
    offset += num_count * 8;

    let string_offsets = read_u32_vec(read_slice(bytes, offset, str_count * 4)?, str_count)?;
    offset += str_count * 4;

    let string_lengths = read_u32_vec(read_slice(bytes, offset, str_count * 4)?, str_count)?;
    offset += str_count * 4;

    let strings_data = bytes.get(offset..).ok_or_else(|| "EOF".to_string())?;

    let last = *item_indices.last().unwrap_or(&0) as usize;
    if last > meta_count {
        return Err("Invalid meta block indices".to_string());
    }

    let mut result = Vec::with_capacity(item_count);
    for i in 0..item_count {
        let start = item_indices[i] as usize;
        let end = item_indices[i + 1] as usize;
        if end > meta_count {
            return Err("Invalid meta block indices".to_string());
        }

        let mut item_params = Vec::with_capacity(end.saturating_sub(start));
        for m in start..end {
            let kind = value_kinds[m];
            let idx = value_indices[m] as usize;

            let value = if kind == 0 {
                numeric_values
                    .get(idx)
                    .map(|n| n.to_string())
                    .unwrap_or_default()
            } else if kind == 1 {
                let s_off = string_offsets.get(idx).copied().unwrap_or(0) as usize;
                let s_len = string_lengths.get(idx).copied().unwrap_or(0) as usize;
                if s_off + s_len <= strings_data.len() {
                    str::from_utf8(&strings_data[s_off..s_off + s_len])
                        .unwrap_or_default()
                        .to_string()
                } else {
                    String::new()
                }
            } else {
                String::new()
            };

            let cv_ref = cv_ref_from_code(meta_ref_codes[m]);
            let accession = make_accession(cv_ref, meta_accessions[m]);
            let name = accession
                .as_deref()
                .and_then(|k| cv_table_name(k))
                .unwrap_or_default();

            let unit_acc = meta_unit_accessions[m];
            let mut unit_ref = cv_ref_from_code(meta_unit_refs[m]);
            if unit_ref.is_none() && unit_acc != 0 {
                let uo_accession = make_accession(Some("UO"), unit_acc);
                unit_ref = if uo_accession
                    .as_deref()
                    .and_then(|k| cv_table::get(k))
                    .is_some()
                {
                    Some("UO")
                } else {
                    Some("MS")
                };
            }

            let unit_accession = make_accession(unit_ref, unit_acc);
            let unit_name = unit_accession.as_deref().and_then(|k| cv_table_name(k));

            item_params.push(CvParam {
                cv_ref: cv_ref.map(|s| s.to_string()),
                accession,
                name,
                value: Some(value),
                unit_cv_ref: unit_ref.map(|s| s.to_string()),
                unit_accession,
                unit_name,
            });
        }

        result.push(item_params);
    }

    Ok(result)
}

#[inline]
fn split_attr_value(params: &mut Vec<CvParam>, attr_tail: u32) -> Option<String> {
    let mut v: Option<String> = None;
    params.retain(|p| {
        if is_attr_param(p) && parse_acc_tail(p.accession.as_deref()) == attr_tail {
            v = attr_string_value(p);
            false
        } else {
            true
        }
    });
    v
}

#[inline]
fn split_id_attr(params: &mut Vec<CvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_ID)
}

#[inline]
fn split_name_attr(params: &mut Vec<CvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_NAME)
}

#[inline]
fn split_version_attr(params: &mut Vec<CvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_VERSION)
}

#[inline]
fn split_instrument_configuration_ref_attr(params: &mut Vec<CvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_INSTRUMENT_CONFIGURATION_REF)
}

#[inline]
fn collect_ref_attrs(params: &mut Vec<CvParam>) -> Vec<String> {
    let mut refs = Vec::<String>::new();
    params.retain(|p| {
        if is_attr_param(p) && parse_acc_tail(p.accession.as_deref()) == ACC_ATTR_REF {
            if let Some(v) = attr_string_value(p) {
                if !v.is_empty() {
                    refs.push(v);
                }
            }
            false
        } else {
            true
        }
    });
    refs
}

#[inline]
fn build_scan_settings_from_cv_params(mut params: Vec<CvParam>) -> ScanSettings {
    let id = split_id_attr(&mut params);
    let instrument_configuration_ref = split_instrument_configuration_ref_attr(&mut params);

    let source_file_refs = collect_ref_attrs(&mut params);
    let source_file_ref_list = (!source_file_refs.is_empty()).then(|| SourceFileRefList {
        count: Some(source_file_refs.len()),
        source_file_refs: source_file_refs
            .into_iter()
            .map(|r| SourceFileRef { r#ref: r })
            .collect(),
    });

    let mut target_cv = Vec::<CvParam>::new();
    let mut rest = Vec::with_capacity(params.len());
    for p in params.drain(..) {
        let t = parse_acc_tail(p.accession.as_deref());
        if matches!(
            t,
            ACC_ISO_TARGET_MZ | ACC_PRODUCT_ION_MZ | ACC_DWELL_TIME | ACC_COMPLETION_TIME
        ) {
            target_cv.push(p);
        } else {
            rest.push(p);
        }
    }

    let target_list = if target_cv.is_empty() {
        None
    } else {
        let mut targets: Vec<Target> = Vec::new();
        let mut cur: Vec<CvParam> = Vec::new();

        for p in target_cv {
            let t = parse_acc_tail(p.accession.as_deref());
            if t == ACC_ISO_TARGET_MZ && !cur.is_empty() {
                targets.push(Target {
                    cv_params: cur,
                    ..Default::default()
                });
                cur = Vec::new();
            }
            cur.push(p);
        }
        if !cur.is_empty() {
            targets.push(Target {
                cv_params: cur,
                ..Default::default()
            });
        }

        Some(TargetList {
            count: Some(targets.len()),
            targets,
        })
    };

    ScanSettings {
        id,
        instrument_configuration_ref,
        source_file_ref_list,
        target_list,
        cv_params: rest,
        ..Default::default()
    }
}

/// <cvList> <fileDescription> <referenceableParamGroupList> <sampleList>
/// <instrumentList> <softwareList> <dataProcessingList> <scanSettingsList>
fn decode_global_meta_structs(
    bytes: &[u8],
    m_cnt: u32,
    n_cnt: u32,
    s_cnt: u32,
) -> Result<
    (
        Option<CvList>,
        FileDescription,
        Option<ReferenceableParamGroupList>,
        Option<SampleList>,
        Option<InstrumentList>,
        Option<SoftwareList>,
        Option<DataProcessingList>,
        Option<ScanSettingsList>,
    ),
    String,
> {
    if bytes.len() < 32 {
        return Ok((
            None,
            FileDescription::default(),
            None,
            None,
            None,
            None,
            None,
            None,
        ));
    }

    let n_fd = read_u32_at(bytes, 0)?;
    let n_rpg = read_u32_at(bytes, 4)?;
    let n_samp = read_u32_at(bytes, 8)?;
    let n_inst = read_u32_at(bytes, 12)?;
    let n_soft = read_u32_at(bytes, 16)?;
    let n_dp = read_u32_at(bytes, 20)?;
    let n_acq = read_u32_at(bytes, 24)?;
    let n_cvs = read_u32_at(bytes, 28)?;

    let total = n_fd + n_rpg + n_samp + n_inst + n_soft + n_dp + n_acq + n_cvs;
    let items = decode_meta_block(&bytes[32..], total, m_cnt, n_cnt, s_cnt)?;
    let mut it = items.into_iter();

    let fd = if n_fd > 0 {
        split_file_description_from_cv_params(it.next().unwrap_or_default())
    } else {
        FileDescription::default()
    };

    let rpgs = if n_rpg > 0 {
        let mut groups = Vec::with_capacity(n_rpg as usize);
        for _ in 0..n_rpg {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();
            groups.push(ReferenceableParamGroup {
                id,
                cv_params: p,
                ..Default::default()
            });
        }
        Some(ReferenceableParamGroupList {
            count: Some(groups.len()),
            referenceable_param_groups: groups,
        })
    } else {
        None
    };

    let samps = if n_samp > 0 {
        let mut samples = Vec::with_capacity(n_samp as usize);
        for _ in 0..n_samp {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();
            let name = split_name_attr(&mut p).unwrap_or_default();
            samples.push(Sample {
                id,
                name,
                cv_params: p,
                ..Default::default()
            });
        }
        Some(SampleList {
            count: Some(samples.len().try_into().unwrap()),
            samples,
        })
    } else {
        None
    };

    let insts = if n_inst > 0 {
        let mut instruments = Vec::with_capacity(n_inst as usize);
        for _ in 0..n_inst {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();

            let mut component_cv = Vec::<CvParam>::new();
            let mut inst_cv = Vec::<CvParam>::new();
            for cv in p {
                let t = parse_acc_tail(cv.accession.as_deref());
                if t == 1_000_398 || t == 1_000_082 || t == 1_000_253 {
                    component_cv.push(cv);
                } else {
                    inst_cv.push(cv);
                }
            }

            let component_list = if component_cv.is_empty() {
                None
            } else {
                let mut src = Vec::<Source>::new();
                let mut an = Vec::<Analyzer>::new();
                let mut det = Vec::<Detector>::new();

                for cv in component_cv {
                    let t = parse_acc_tail(cv.accession.as_deref());
                    if t == 1_000_398 {
                        src.push(Source {
                            order: Some(1),
                            cv_param: vec![cv],
                            ..Default::default()
                        });
                    } else if t == 1_000_082 {
                        an.push(Analyzer {
                            order: Some(2),
                            cv_param: vec![cv],
                            ..Default::default()
                        });
                    } else if t == 1_000_253 {
                        det.push(Detector {
                            order: Some(3),
                            cv_param: vec![cv],
                            ..Default::default()
                        });
                    }
                }

                Some(ComponentList {
                    source: src,
                    analyzer: an,
                    detector: det,
                    ..Default::default()
                })
            };

            instruments.push(Instrument {
                id,
                cv_param: inst_cv,
                component_list,
                ..Default::default()
            });
        }
        Some(InstrumentList {
            count: Some(instruments.len()),
            instrument: instruments,
        })
    } else {
        None
    };

    let softs = if n_soft > 0 {
        let mut software = Vec::with_capacity(n_soft as usize);
        for _ in 0..n_soft {
            let mut p = it.next().unwrap_or_default();

            let id = split_id_attr(&mut p).unwrap_or_default();
            let version = split_version_attr(&mut p);

            let mut software_param = Vec::<SoftwareParam>::new();
            for cv in p {
                if is_attr_param(&cv) {
                    continue;
                }
                software_param.push(SoftwareParam {
                    cv_ref: cv.cv_ref,
                    accession: cv.accession.unwrap_or_default(),
                    name: cv.name,
                    version: version.clone(),
                    ..Default::default()
                });
            }

            software.push(Software {
                id,
                software_param,
                ..Default::default()
            });
        }
        Some(SoftwareList {
            count: Some(software.len()),
            software,
        })
    } else {
        None
    };

    let dps = if n_dp > 0 {
        let mut data_processing = Vec::with_capacity(n_dp as usize);
        for _ in 0..n_dp {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();
            data_processing.push(DataProcessing {
                id,
                processing_method: vec![ProcessingMethod {
                    cv_param: p,
                    ..Default::default()
                }],
                ..Default::default()
            });
        }
        Some(DataProcessingList {
            count: Some(data_processing.len()),
            data_processing,
        })
    } else {
        None
    };

    let acqs = if n_acq > 0 {
        let mut scan_settings = Vec::with_capacity(n_acq as usize);
        for _ in 0..n_acq {
            let p = it.next().unwrap_or_default();
            scan_settings.push(build_scan_settings_from_cv_params(p));
        }
        Some(ScanSettingsList {
            count: Some(scan_settings.len()),
            scan_settings,
        })
    } else {
        None
    };

    let cvs = if n_cvs > 0 {
        let mut cv = Vec::<Cv>::with_capacity(n_cvs as usize);
        for _ in 0..n_cvs {
            let p = it.next().unwrap_or_default();
            let mut c = Cv::default();
            for param in p {
                let tail = parse_acc_tail(param.accession.as_deref());
                if tail == 9_900_001 {
                    c.id = attr_string_value(&param).unwrap_or_default();
                } else if tail == 9_900_002 {
                    c.full_name = Some(attr_string_value(&param).unwrap_or_default());
                } else if tail == 9_900_003 {
                    c.version = Some(attr_string_value(&param).unwrap_or_default());
                } else if tail == 9_900_004 {
                    c.uri = Some(attr_string_value(&param).unwrap_or_default());
                }
            }
            cv.push(c);
        }
        Some(CvList {
            count: Some(cv.len()),
            cv,
        })
    } else {
        None
    };

    Ok((cvs, fd, rpgs, samps, insts, softs, dps, acqs))
}

#[inline]
fn read_u8_at(b: &[u8], o: usize) -> Result<u8, String> {
    b.get(o).copied().ok_or_else(|| "EOF".to_string())
}

#[inline]
fn read_u32_at(b: &[u8], o: usize) -> Result<u32, String> {
    let s = b.get(o..o + 4).ok_or_else(|| "EOF".to_string())?;
    Ok(u32::from_le_bytes(s.try_into().unwrap()))
}

#[inline]
fn read_u64_at(b: &[u8], o: usize) -> Result<u64, String> {
    let s = b.get(o..o + 8).ok_or_else(|| "EOF".to_string())?;
    Ok(u64::from_le_bytes(s.try_into().unwrap()))
}

#[inline]
fn read_slice(b: &[u8], o: usize, l: usize) -> Result<&[u8], String> {
    let end = o.checked_add(l).ok_or_else(|| "EOF".to_string())?;
    b.get(o..end).ok_or_else(|| "EOF".to_string())
}

fn read_u32_vec(b: &[u8], c: usize) -> Result<Vec<u32>, String> {
    let need = c.checked_mul(4).ok_or_else(|| "EOF".to_string())?;
    if b.len() < need {
        return Err("EOF".to_string());
    }
    let b = &b[..need];

    if cfg!(target_endian = "little") {
        let mut out: Vec<u32> = Vec::with_capacity(c);
        unsafe {
            out.set_len(c);
            std::ptr::copy_nonoverlapping(b.as_ptr(), out.as_mut_ptr() as *mut u8, need);
        }
        return Ok(out);
    }

    let mut out = Vec::with_capacity(c);
    for s in b.chunks_exact(4) {
        out.push(u32::from_le_bytes(s.try_into().unwrap()));
    }
    Ok(out)
}

fn read_f64_vec(b: &[u8], c: usize) -> Result<Vec<f64>, String> {
    let need = c.checked_mul(8).ok_or_else(|| "EOF".to_string())?;
    if b.len() < need {
        return Err("EOF".to_string());
    }
    let b = &b[..need];

    if cfg!(target_endian = "little") {
        let mut out: Vec<f64> = Vec::with_capacity(c);
        unsafe {
            out.set_len(c);
            std::ptr::copy_nonoverlapping(b.as_ptr(), out.as_mut_ptr() as *mut u8, need);
        }
        return Ok(out);
    }

    let mut out = Vec::with_capacity(c);
    for s in b.chunks_exact(8) {
        out.push(f64::from_le_bytes(s.try_into().unwrap()));
    }
    Ok(out)
}

fn cv_ref_from_code(c: u8) -> Option<&'static str> {
    match c {
        0 => Some("MS"),
        1 => Some("UO"),
        2 => Some("NCIT"),
        3 => Some("PEFF"),
        4 => Some(CV_REF_ATTR),
        _ => None,
    }
}

fn make_accession(r: Option<&str>, a: u32) -> Option<String> {
    if a == 0 {
        return None;
    }
    match r {
        Some("MS") | Some("UO") | Some("PEFF") => {
            let p = r.unwrap();
            Some(format!("{}:{:07}", p, a))
        }
        Some("NCIT") => Some(format!("NCIT:C{:05}", a)),
        Some(CV_REF_ATTR) => Some(format!("{}:{:07}", CV_REF_ATTR, a)),
        Some(cv) => Some(format!("{}:{}", cv, a)),
        None => Some(a.to_string()),
    }
}

fn parse_acc_tail(accession: Option<&str>) -> u32 {
    let s = accession.unwrap_or("");
    let tail = s.rsplit_once(':').map(|(_, t)| t).unwrap_or(s);

    let mut v: u32 = 0;
    let mut saw_digit = false;

    for b in tail.bytes() {
        if (b'0'..=b'9').contains(&b) {
            saw_digit = true;
            let d = (b - b'0') as u32;
            match v.checked_mul(10).and_then(|x| x.checked_add(d)) {
                Some(n) => v = n,
                None => return 0,
            }
        }
    }

    if saw_digit { v } else { 0 }
}
