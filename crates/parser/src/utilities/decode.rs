use std::io::Cursor;
use std::str;

use miniz_oxide::inflate::decompress_to_vec_zlib;

use zstd::bulk::decompress as zstd_decompress;
use zstd::stream::decode_all as zstd_decode_all;

use crate::utilities::{cv_table, mzml::*};

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

const ACC_SELECTED_ION_MZ: u32 = 1_000_744;
const ACC_CHARGE_STATE: u32 = 1_000_041;

const ACC_IN_SOURCE_CID: u32 = 1_001_880;
const ACC_COLLISION_ENERGY: u32 = 1_000_045;

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
        for i in 0..block_count {
            let base = i * BLOCK_DIR_ENTRY_SIZE;
            let comp_off = u64::from_le_bytes(container_bytes[base..base + 8].try_into().unwrap());
            let comp_size =
                u64::from_le_bytes(container_bytes[base + 8..base + 16].try_into().unwrap());
            let uncomp_bytes =
                u64::from_le_bytes(container_bytes[base + 16..base + 24].try_into().unwrap());
            dir.push(BlockDirEntry {
                comp_off,
                comp_size,
                uncomp_bytes,
            });
        }

        let compressed_region = &container_bytes[dir_bytes..];

        let mut block_start_elems = Vec::with_capacity(block_count + 1);
        block_start_elems.push(0);

        let elem_size_u64 = elem_size as u64;
        let mut acc = 0u64;
        for e in &dir {
            let elems = e.uncomp_bytes / elem_size_u64;
            acc = acc.saturating_add(elems);
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

        let comp = self
            .compressed_region
            .get(comp_off..end)
            .ok_or_else(|| "EOF".to_string())?;

        let needs_owned = self.compression_level != 0
            || (self.array_filter == ARRAY_FILTER_BYTE_SHUFFLE && self.elem_size > 1);

        if !needs_owned {
            return Ok(comp);
        }

        if self.cache[id].is_none() {
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

        let elem_size = self.elem_size;
        let local_elems = (global_elem_off - block_start) as usize;

        let byte_off = local_elems
            .checked_mul(elem_size)
            .ok_or_else(|| "Byte offset overflow".to_string())?;
        let byte_len = (elem_len as usize)
            .checked_mul(elem_size)
            .ok_or_else(|| "Byte length overflow".to_string())?;
        let end = byte_off
            .checked_add(byte_len)
            .ok_or_else(|| "Slice range overflow".to_string())?;

        let block = self.block_bytes(block_id)?;
        block.get(byte_off..end).ok_or_else(|| "EOF".to_string())
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
            BytesMaybeOwned::Borrowed(b) => b,
            BytesMaybeOwned::Owned(v) => v.as_slice(),
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
    matches!(t, ACC_SELECTED_ION_MZ | ACC_CHARGE_STATE)
}

#[inline]
fn is_activation_tail(t: u32) -> bool {
    matches!(t, ACC_IN_SOURCE_CID | ACC_COLLISION_ENERGY)
}

/// <precursorList>
fn infer_precursor_list_from_spectrum_cv(params: &mut Vec<CvParam>) -> Option<PrecursorList> {
    let mut iso = Vec::<CvParam>::new();
    let mut sel = Vec::<CvParam>::new();
    let mut act = Vec::<CvParam>::new();

    let mut rest = Vec::<CvParam>::with_capacity(params.len());
    for p in params.drain(..) {
        let tail = parse_acc_tail(p.accession.as_deref());
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

    if iso.is_empty() && sel.is_empty() && act.is_empty() {
        return None;
    }

    let isolation_window = if iso.is_empty() {
        None
    } else {
        Some(IsolationWindow {
            cv_params: iso,
            ..Default::default()
        })
    };

    let selected_ion_list = if sel.is_empty() {
        None
    } else {
        Some(SelectedIonList {
            count: Some(1),
            selected_ions: vec![SelectedIon {
                cv_params: sel,
                ..Default::default()
            }],
        })
    };

    let activation = if act.is_empty() {
        None
    } else {
        Some(Activation {
            cv_params: act,
            ..Default::default()
        })
    };

    Some(PrecursorList {
        count: Some(1),
        precursors: vec![Precursor {
            isolation_window,
            selected_ion_list,
            activation,
            ..Default::default()
        }],
    })
}

/// <mzML>
pub fn decode(bytes: &[u8]) -> Result<MzML, String> {
    if bytes.len() < HEADER_SIZE {
        return Err("Buffer too small for header".to_string());
    }

    let header = &bytes[..HEADER_SIZE];
    if &header[0..4] != b"B000" {
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

    let spectrum_index_bytes = read_slice(
        bytes,
        off_spec_index,
        spectrum_count as usize * INDEX_ENTRY_SIZE,
    )?;
    let chromatogram_index_bytes = read_slice(
        bytes,
        off_chrom_index,
        chrom_count as usize * INDEX_ENTRY_SIZE,
    )?;

    if off_chrom_meta < off_spec_meta || off_global_meta < off_chrom_meta {
        return Err("Invalid meta offsets".to_string());
    }

    let spec_meta_bytes = read_slice(bytes, off_spec_meta, off_chrom_meta - off_spec_meta)?;
    let chrom_meta_bytes = read_slice(bytes, off_chrom_meta, off_global_meta - off_chrom_meta)?;

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

    let global_meta_bytes = read_slice(
        bytes,
        off_global_meta,
        first_container_off - off_global_meta,
    )?;

    let spec_meta_bytes = decompress_meta_if_needed(
        codec,
        (codec_flags & HDR_FLAG_SPEC_META_COMP) != 0,
        spec_meta_bytes,
    )?;
    let chrom_meta_bytes = decompress_meta_if_needed(
        codec,
        (codec_flags & HDR_FLAG_CHROM_META_COMP) != 0,
        chrom_meta_bytes,
    )?;
    let global_meta_bytes = decompress_meta_if_needed(
        codec,
        (codec_flags & HDR_FLAG_GLOBAL_META_COMP) != 0,
        global_meta_bytes,
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

        let mut mz_ba = BinaryDataArray::default();
        mz_ba.array_length = Some(x_len as usize);
        mz_ba.is_f32 = Some(spect_x_fmt == 1);
        mz_ba.is_f64 = Some(spect_x_fmt == 2);
        mz_ba.cv_params.push(ms_cv_param(ACC_MZ_ARRAY));
        mz_ba.decoded_binary_f32 = mz_f32;
        mz_ba.decoded_binary_f64 = mz_f64;

        let mut inten_ba = BinaryDataArray::default();
        inten_ba.array_length = Some(y_len as usize);
        inten_ba.is_f32 = Some(spect_y_fmt == 1);
        inten_ba.is_f64 = Some(spect_y_fmt == 2);
        inten_ba.cv_params.push(ms_cv_param(ACC_INTENSITY_ARRAY));
        inten_ba.decoded_binary_f32 = in_f32;
        inten_ba.decoded_binary_f64 = in_f64;

        let mut spectrum_params = item_params;
        strip_binary_array_cv_params(&mut spectrum_params);

        let precursor_list = infer_precursor_list_from_spectrum_cv(&mut spectrum_params);

        spectra.push(Spectrum {
            id: format!("spectrum_{}", i),
            index: Some(i as u32),
            default_array_length: Some(x_len as usize),
            cv_params: spectrum_params,
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

        let mut time_ba = BinaryDataArray::default();
        time_ba.array_length = Some(x_len as usize);
        time_ba.is_f32 = Some(chrom_x_fmt == 1);
        time_ba.is_f64 = Some(chrom_x_fmt == 2);
        time_ba.cv_params.push(ms_cv_param(ACC_TIME_ARRAY));
        time_ba.decoded_binary_f32 = t_f32;
        time_ba.decoded_binary_f64 = t_f64;

        let mut inten_ba = BinaryDataArray::default();
        inten_ba.array_length = Some(y_len as usize);
        inten_ba.is_f32 = Some(chrom_y_fmt == 1);
        inten_ba.is_f64 = Some(chrom_y_fmt == 2);
        inten_ba.cv_params.push(ms_cv_param(ACC_INTENSITY_ARRAY));
        inten_ba.decoded_binary_f32 = in_f32;
        inten_ba.decoded_binary_f64 = in_f64;

        let mut chrom_params = item_params;
        strip_binary_array_cv_params(&mut chrom_params);

        chromatograms.push(Chromatogram {
            id: format!("chromatogram_{}", j),
            index: Some(j as u32),
            default_array_length: Some(x_len as usize),
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
        m = Some(match m {
            Some(cur) => cur.min(x),
            None => x,
        });
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

/// <cvParam>
#[inline]
fn ms_cv_param(accession_tail: u32) -> CvParam {
    let key = format!("MS:{:07}", accession_tail);
    let name = cv_table::get(&key)
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .to_string();
    CvParam {
        cv_ref: Some("MS".to_string()),
        accession: Some(key),
        name,
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

    let numeric_values = read_f64_vec(
        read_slice(bytes, offset, num_count as usize * 8)?,
        num_count as usize,
    )?;
    offset += num_count as usize * 8;

    let string_offsets = read_u32_vec(
        read_slice(bytes, offset, str_count as usize * 4)?,
        str_count as usize,
    )?;
    offset += str_count as usize * 4;

    let string_lengths = read_u32_vec(
        read_slice(bytes, offset, str_count as usize * 4)?,
        str_count as usize,
    )?;
    offset += str_count as usize * 4;

    let strings_data = bytes.get(offset..).ok_or_else(|| "EOF".to_string())?;

    let mut result = Vec::with_capacity(item_count);
    for i in 0..item_count {
        let start = item_indices[i] as usize;
        let end = item_indices[i + 1] as usize;

        let mut item_params = Vec::with_capacity(end.saturating_sub(start));
        for m in start..end {
            let kind = value_kinds[m];
            let idx = value_indices[m] as usize;

            let value = if kind == 0 && idx < numeric_values.len() {
                Some(numeric_values[idx].to_string())
            } else if kind == 1 && idx < string_offsets.len() {
                let s_off = string_offsets[idx] as usize;
                let s_len = string_lengths[idx] as usize;
                if s_off + s_len <= strings_data.len() {
                    Some(
                        str::from_utf8(&strings_data[s_off..s_off + s_len])
                            .unwrap_or_default()
                            .to_string(),
                    )
                } else {
                    Some(String::new())
                }
            } else {
                None
            };

            let cv_ref = cv_ref_from_code(meta_ref_codes[m]);
            let unit_ref = cv_ref_from_code(meta_unit_refs[m]);

            item_params.push(CvParam {
                cv_ref: cv_ref.map(|s| s.to_string()),
                accession: make_accession(cv_ref, meta_accessions[m]),
                name: cv_name_from_code(cv_ref, meta_accessions[m]).unwrap_or_default(),
                value,
                unit_cv_ref: unit_ref.map(|s| s.to_string()),
                unit_accession: make_accession(unit_ref, meta_unit_accessions[m]),
                unit_name: cv_name_from_code(unit_ref, meta_unit_accessions[m]),
            });
        }

        result.push(item_params);
    }

    Ok(result)
}

/// <cvList>
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

    let mut fd = FileDescription::default();
    if n_fd > 0 {
        fd.file_content.cv_params = it.next().unwrap_or_default();
    }

    let rpgs = if n_rpg > 0 {
        Some(ReferenceableParamGroupList {
            count: Some(n_rpg as usize),
            referenceable_param_groups: (0..n_rpg)
                .map(|_| ReferenceableParamGroup {
                    cv_params: it.next().unwrap_or_default(),
                    ..Default::default()
                })
                .collect(),
        })
    } else {
        None
    };

    let samps = if n_samp > 0 {
        Some(SampleList {
            count: Some(n_samp),
            samples: (0..n_samp)
                .map(|_| Sample {
                    cv_params: it.next().unwrap_or_default(),
                    ..Default::default()
                })
                .collect(),
        })
    } else {
        None
    };

    let insts = if n_inst > 0 {
        Some(InstrumentList {
            count: Some(n_inst as usize),
            instrument: (0..n_inst)
                .map(|_| Instrument {
                    cv_param: it.next().unwrap_or_default(),
                    ..Default::default()
                })
                .collect(),
        })
    } else {
        None
    };

    let softs = if n_soft > 0 {
        Some(SoftwareList {
            count: Some(n_soft as usize),
            software: (0..n_soft)
                .map(|_| Software {
                    cv_param: it.next().unwrap_or_default(),
                    ..Default::default()
                })
                .collect(),
        })
    } else {
        None
    };

    let dps = if n_dp > 0 {
        Some(DataProcessingList {
            count: Some(n_dp as usize),
            data_processing: (0..n_dp)
                .map(|_| DataProcessing {
                    processing_method: vec![ProcessingMethod {
                        cv_param: it.next().unwrap_or_default(),
                        ..Default::default()
                    }],
                    ..Default::default()
                })
                .collect(),
        })
    } else {
        None
    };

    let acqs = if n_acq > 0 {
        Some(ScanSettingsList {
            count: Some(n_acq as usize),
            scan_settings: (0..n_acq)
                .map(|_| ScanSettings {
                    cv_params: it.next().unwrap_or_default(),
                    ..Default::default()
                })
                .collect(),
        })
    } else {
        None
    };

    let cvs = if n_cvs > 0 {
        let cv: Vec<Cv> = (0..n_cvs)
            .map(|_| {
                let p = it.next().unwrap_or_default();
                let mut c = Cv::default();
                for param in p {
                    let tail = parse_acc_tail(param.accession.as_deref());
                    if tail == 9_900_001 {
                        c.id = param.value.unwrap_or_default();
                    } else if tail == 9_900_002 {
                        c.full_name = param.value;
                    } else if tail == 9_900_003 {
                        c.version = param.value;
                    } else if tail == 9_900_004 {
                        c.uri = param.value;
                    }
                }
                c
            })
            .collect();
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
    if b.len() < c * 4 {
        return Err("EOF".to_string());
    }
    let mut out = Vec::with_capacity(c);
    for s in b.chunks_exact(4).take(c) {
        out.push(u32::from_le_bytes(s.try_into().unwrap()));
    }
    Ok(out)
}

fn read_f64_vec(b: &[u8], c: usize) -> Result<Vec<f64>, String> {
    if b.len() < c * 8 {
        return Err("EOF".to_string());
    }
    let mut out = Vec::with_capacity(c);
    for s in b.chunks_exact(8).take(c) {
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
        _ => None,
    }
}

fn make_accession(r: Option<&str>, a: u32) -> Option<String> {
    if a == 0 {
        return None;
    }

    match r {
        Some("MS") | Some("UO") | Some("PEFF") => {
            let prefix = r.unwrap();
            Some(format!("{}:{:07}", prefix, a))
        }
        Some("NCIT") => Some(format!("NCIT:C{:05}", a)),
        Some(cv) => Some(format!("{}:{}", cv, a)),
        None => Some(a.to_string()),
    }
}

fn cv_name_from_code(r: Option<&str>, a: u32) -> Option<String> {
    if a == 0 || r.is_none() {
        return None;
    }
    cv_table::get(&format!("{}:{:07}", r.unwrap(), a))
        .and_then(|v| v.as_str())
        .map(|s| s.to_string())
}

fn parse_acc_tail(acc: Option<&str>) -> u32 {
    acc.and_then(|s| s.rsplit(':').next())
        .and_then(|s| s.parse().ok())
        .unwrap_or(0)
}
