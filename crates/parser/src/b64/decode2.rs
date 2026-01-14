use std::collections::{HashMap, HashSet};

use crate::{
    b64::utilities::{
        Header, common::*, parse_binary_data_array_list, parse_chromatogram_list,
        parse_cv_and_user_params, parse_header, parse_metadata, parse_precursor_list,
        parse_product_list, parse_scan_list, parse_spectrum_list,
    },
    mzml::{
        attr_meta::*,
        cv_table,
        schema::{SchemaTree as Schema, TagId, schema},
        structs::*,
    },
};

const INDEX_ENTRY_SIZE: usize = 32;
const BLOCK_DIR_ENTRY_SIZE: usize = 32;

const HDR_FLAG_SPEC_META_COMP: u8 = 1 << 4;
const HDR_FLAG_CHROM_META_COMP: u8 = 1 << 5;

const ARRAY_FILTER_NONE: u8 = 0;
const ARRAY_FILTER_BYTE_SHUFFLE: u8 = 1;

const ACC_MZ_ARRAY: u32 = 1_000_514;
const ACC_INTENSITY_ARRAY: u32 = 1_000_515;
const ACC_TIME_ARRAY: u32 = 1_000_595;

const ACC_32BIT_FLOAT: u32 = 1_000_521;
const ACC_64BIT_FLOAT: u32 = 1_000_523;

pub fn decode2(bytes: &[u8]) -> Result<MzML, String> {
    let schema = schema();
    let header = parse_header(bytes)?;
    Ok(MzML {
        cv_list: parse_cv_list(schema, bytes),
        file_description: parse_file_description(schema, bytes),
        referenceable_param_group_list: parse_referenceable_param_group_list(schema, bytes),
        sample_list: parse_sample_list(schema, bytes),
        instrument_list: parse_instrument_list(schema, bytes),
        software_list: parse_software_list(schema, bytes),
        data_processing_list: parse_data_processing_list(schema, bytes),
        scan_settings_list: parse_scan_settings_list(schema, bytes),
        run: parse_run(schema, bytes, &header)?,
    })
}

/// <run>
#[inline]
fn parse_run(schema: &Schema, bytes: &[u8], header: &Header) -> Result<Run, String> {
    let metadata = parse_metadata_section(
        bytes,
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        2,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
        42,
        "spectra",
    );

    let spec_child_index = ChildIndex::new(&metadata);

    Ok(Run {
        spectrum_list: parse_spectrum_list(schema, &metadata, &spec_child_index),
        chromatogram_list: parse_chromatogram_list(schema, &metadata, &spec_child_index),
        ..Default::default()
    })
}

#[derive(Clone, Copy, Debug)]
struct SpectrumIndexEntry {
    mz_element_off: u64,
    inten_element_off: u64,
    mz_element_len: u32,
    inten_element_len: u32,
    mz_block_id: u32,
    inten_block_id: u32,
}

#[derive(Clone, Copy, Debug)]
struct ChromIndexEntry {
    time_element_off: u64,
    inten_element_off: u64,
    time_element_len: u32,
    inten_element_len: u32,
    time_block_id: u32,
    inten_block_id: u32,
}

#[inline]
fn slice_at<'a>(
    bytes: &'a [u8],
    off: u64,
    len: u64,
    field: &'static str,
) -> Result<&'a [u8], String> {
    let off = usize::try_from(off).map_err(|_| format!("{field}: offset overflow"))?;
    let len = usize::try_from(len).map_err(|_| format!("{field}: len overflow"))?;
    let end = off
        .checked_add(len)
        .ok_or_else(|| format!("{field}: range overflow"))?;
    if end > bytes.len() {
        return Err(format!(
            "{field}: out of range (off={off}, len={len}, file_len={})",
            bytes.len()
        ));
    }
    Ok(&bytes[off..end])
}

#[inline]
fn read_u32_le_at(bytes: &[u8], pos: &mut usize, field: &'static str) -> Result<u32, String> {
    let s = take(bytes, pos, 4, field)?;
    Ok(u32::from_le_bytes(s.try_into().unwrap()))
}

#[inline]
fn read_u64_le_at(bytes: &[u8], pos: &mut usize, field: &'static str) -> Result<u64, String> {
    let s = take(bytes, pos, 8, field)?;
    Ok(u64::from_le_bytes(s.try_into().unwrap()))
}

#[inline]
fn parse_chrom_index(bytes: &[u8], header: &Header) -> Result<Vec<ChromIndexEntry>, String> {
    let count = header.chrom_count as usize;
    let off = header.off_chrom_index;
    let need = (count as u64)
        .checked_mul(INDEX_ENTRY_SIZE as u64)
        .ok_or_else(|| "chrom index size overflow".to_string())?;
    let raw = slice_at(bytes, off, need, "chrom index")?;

    let mut pos = 0usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let time_element_off = read_u64_le_at(raw, &mut pos, "time_element_off")?;
        let inten_element_off = read_u64_le_at(raw, &mut pos, "inten_element_off")?;
        let time_element_len = read_u32_le_at(raw, &mut pos, "time_element_len")?;
        let inten_element_len = read_u32_le_at(raw, &mut pos, "inten_element_len")?;
        let time_block_id = read_u32_le_at(raw, &mut pos, "time_block_id")?;
        let inten_block_id = read_u32_le_at(raw, &mut pos, "inten_block_id")?;
        out.push(ChromIndexEntry {
            time_element_off,
            inten_element_off,
            time_element_len,
            inten_element_len,
            time_block_id,
            inten_block_id,
        });
    }
    Ok(out)
}

#[derive(Clone, Copy)]
struct BlockDirEntry {
    comp_off: u64,
    comp_size: u64,
    uncomp_bytes: u64,
}

struct ContainerReader<'a> {
    bytes: &'a [u8],
    elem_size: usize,
    compression_level: u8,
    array_filter: u8,
    dir: Vec<BlockDirEntry>,
    comp_buf_start: usize,
    cache: Vec<Option<Vec<u8>>>,
    scratch: Vec<u8>,
}

impl<'a> ContainerReader<'a> {
    #[inline]
    fn new(
        bytes: &'a [u8],
        block_count: u32,
        elem_size: usize,
        compression_level: u8,
        array_filter: u8,
    ) -> Result<Self, String> {
        let bc = block_count as usize;
        let dir_bytes = bc
            .checked_mul(BLOCK_DIR_ENTRY_SIZE)
            .ok_or_else(|| "container dir size overflow".to_string())?;

        if bytes.len() < dir_bytes {
            return Err("container too small for directory".to_string());
        }

        let mut pos = 0usize;
        let dir_raw = &bytes[..dir_bytes];
        let mut dir = Vec::with_capacity(bc);
        for _ in 0..bc {
            let comp_off = read_u64_le_at(dir_raw, &mut pos, "comp_off")?;
            let comp_size = read_u64_le_at(dir_raw, &mut pos, "comp_size")?;
            let uncomp_bytes = read_u64_le_at(dir_raw, &mut pos, "uncomp_bytes")?;
            let _ = take(dir_raw, &mut pos, 8, "reserved")?;
            dir.push(BlockDirEntry {
                comp_off,
                comp_size,
                uncomp_bytes,
            });
        }

        Ok(Self {
            bytes,
            elem_size,
            compression_level,
            array_filter,
            dir,
            comp_buf_start: dir_bytes,
            cache: vec![None; bc],
            scratch: Vec::new(),
        })
    }

    #[inline]
    fn ensure_block(&mut self, block_id: u32) -> Result<(), String> {
        let i = block_id as usize;
        if i >= self.cache.len() {
            return Err(format!("block_id out of range: {block_id}"));
        }
        if self.cache[i].is_some() {
            return Ok(());
        }

        let e = self.dir[i];
        let start = self
            .comp_buf_start
            .checked_add(usize::try_from(e.comp_off).map_err(|_| "comp_off overflow".to_string())?)
            .ok_or_else(|| "comp start overflow".to_string())?;
        let size = usize::try_from(e.comp_size).map_err(|_| "comp_size overflow".to_string())?;
        let end = start
            .checked_add(size)
            .ok_or_else(|| "comp end overflow".to_string())?;

        if end > self.bytes.len() {
            return Err("container: block range out of bounds".to_string());
        }

        let comp = &self.bytes[start..end];
        let mut out = if self.compression_level == 0 {
            comp.to_vec()
        } else {
            decompress_zstd(comp)?
        };

        let expected =
            usize::try_from(e.uncomp_bytes).map_err(|_| "uncomp_bytes overflow".to_string())?;
        if out.len() != expected {
            return Err(format!(
                "container: bad block size (block_id={block_id}, got={}, expected={})",
                out.len(),
                expected
            ));
        }

        if self.array_filter == ARRAY_FILTER_BYTE_SHUFFLE && self.elem_size > 1 {
            self.scratch.resize(out.len(), 0);
            byte_unshuffle_into(&out, &mut self.scratch, self.elem_size);
            out.clear();
            out.extend_from_slice(&self.scratch);
        }

        self.cache[i] = Some(out);
        Ok(())
    }

    #[inline]
    fn block_bytes(&mut self, block_id: u32) -> Result<&[u8], String> {
        self.ensure_block(block_id)?;
        Ok(self.cache[block_id as usize].as_ref().unwrap().as_slice())
    }
}

#[inline]
fn byte_unshuffle_into(input: &[u8], output: &mut [u8], elem_size: usize) {
    let count = input.len() / elem_size;
    for b in 0..elem_size {
        let in_base = b * count;
        for e in 0..count {
            output[b + e * elem_size] = input[in_base + e];
        }
    }
}

#[derive(Clone, Debug)]
enum ArrayData {
    F32(Vec<f32>),
    F64(Vec<f64>),
}

#[inline]
fn bytes_to_f32_vec(raw: &[u8]) -> Vec<f32> {
    let mut out = Vec::with_capacity(raw.len() / 4);
    for c in raw.chunks_exact(4) {
        out.push(f32::from_le_bytes(c.try_into().unwrap()));
    }
    out
}

#[inline]
fn bytes_to_f64_vec(raw: &[u8]) -> Vec<f64> {
    let mut out = Vec::with_capacity(raw.len() / 8);
    for c in raw.chunks_exact(8) {
        out.push(f64::from_le_bytes(c.try_into().unwrap()));
    }
    out
}

#[inline]
fn compute_block_starts_for_x(
    index: &[SpectrumIndexEntry],
    block_count: u32,
) -> Result<Vec<u64>, String> {
    let mut starts = vec![u64::MAX; block_count as usize];
    for e in index {
        let bi = e.mz_block_id as usize;
        if bi >= starts.len() {
            return Err("mz_block_id out of range".to_string());
        }
        starts[bi] = starts[bi].min(e.mz_element_off);
    }
    Ok(starts)
}

#[inline]
fn compute_block_starts_for_y(
    index: &[SpectrumIndexEntry],
    block_count: u32,
) -> Result<Vec<u64>, String> {
    let mut starts = vec![u64::MAX; block_count as usize];
    for e in index {
        let bi = e.inten_block_id as usize;
        if bi >= starts.len() {
            return Err("inten_block_id out of range".to_string());
        }
        starts[bi] = starts[bi].min(e.inten_element_off);
    }
    Ok(starts)
}

#[inline]
fn compute_block_starts_for_cx(
    index: &[ChromIndexEntry],
    block_count: u32,
) -> Result<Vec<u64>, String> {
    let mut starts = vec![u64::MAX; block_count as usize];
    for e in index {
        let bi = e.time_block_id as usize;
        if bi >= starts.len() {
            return Err("time_block_id out of range".to_string());
        }
        starts[bi] = starts[bi].min(e.time_element_off);
    }
    Ok(starts)
}

#[inline]
fn compute_block_starts_for_cy(
    index: &[ChromIndexEntry],
    block_count: u32,
) -> Result<Vec<u64>, String> {
    let mut starts = vec![u64::MAX; block_count as usize];
    for e in index {
        let bi = e.inten_block_id as usize;
        if bi >= starts.len() {
            return Err("inten_block_id out of range".to_string());
        }
        starts[bi] = starts[bi].min(e.inten_element_off);
    }
    Ok(starts)
}

#[inline]
fn decode_item_array(
    reader: &mut ContainerReader<'_>,
    block_starts: &[u64],
    block_id: u32,
    global_off_elems: u64,
    len_elems: u32,
) -> Result<ArrayData, String> {
    let bi = block_id as usize;
    if bi >= block_starts.len() {
        return Err("block_id out of range for starts".to_string());
    }

    let start = block_starts[bi];
    if start == u64::MAX {
        return Err("block start unknown".to_string());
    }

    let local_off_elems = global_off_elems
        .checked_sub(start)
        .ok_or_else(|| "negative local offset".to_string())?;

    let elem_size = reader.elem_size;
    if elem_size != 4 && elem_size != 8 {
        return Err("unsupported elem_size".to_string());
    }

    let local_off_elems =
        usize::try_from(local_off_elems).map_err(|_| "local offset overflow".to_string())?;
    let len_elems = usize::try_from(len_elems).map_err(|_| "len overflow".to_string())?;

    let off_bytes = local_off_elems
        .checked_mul(elem_size)
        .ok_or_else(|| "offset bytes overflow".to_string())?;
    let len_bytes = len_elems
        .checked_mul(elem_size)
        .ok_or_else(|| "len bytes overflow".to_string())?;

    let raw = reader.block_bytes(block_id)?;
    let end = off_bytes
        .checked_add(len_bytes)
        .ok_or_else(|| "slice end overflow".to_string())?;
    if end > raw.len() {
        return Err("array slice out of bounds".to_string());
    }

    let slice = &raw[off_bytes..end];
    Ok(match elem_size {
        4 => ArrayData::F32(bytes_to_f32_vec(slice)),
        8 => ArrayData::F64(bytes_to_f64_vec(slice)),
        _ => unreachable!(),
    })
}

#[inline]
fn bda_has_array_kind(bda: &BinaryDataArray, tail: u32) -> bool {
    bda.cv_params.iter().any(|p| {
        p.accession
            .as_deref()
            .map(|a| parse_accession_tail_str(a) == tail)
            .unwrap_or(false)
    })
}

#[inline]
fn ensure_float_flag(bda: &mut BinaryDataArray, is_f32: bool) {
    if is_f32 {
        if !bda
            .cv_params
            .iter()
            .any(|p| p.accession.as_deref() == Some("MS:1000521"))
        {
            bda.cv_params.push(ms_float_param(ACC_32BIT_FLOAT));
        }
        bda.is_f32 = Some(true);
        bda.is_f64 = None;
    } else {
        if !bda
            .cv_params
            .iter()
            .any(|p| p.accession.as_deref() == Some("MS:1000523"))
        {
            bda.cv_params.push(ms_float_param(ACC_64BIT_FLOAT));
        }
        bda.is_f64 = Some(true);
        bda.is_f32 = None;
    }
}

#[inline]
fn attach_xy_arrays_to_bdal(
    list: &mut BinaryDataArrayList,
    x: &ArrayData,
    y: &ArrayData,
    x_kind: u32,
    y_kind: u32,
) {
    let mut x_i = None;
    let mut y_i = None;

    for (i, bda) in list.binary_data_arrays.iter().enumerate() {
        if x_i.is_none() && bda_has_array_kind(bda, x_kind) {
            x_i = Some(i);
        }
        if y_i.is_none() && bda_has_array_kind(bda, y_kind) {
            y_i = Some(i);
        }
    }

    let x_idx = x_i.unwrap_or(0);
    let y_idx = y_i.unwrap_or_else(|| {
        if list.binary_data_arrays.len() > 1 {
            1
        } else {
            0
        }
    });

    if list.binary_data_arrays.is_empty() {
        return;
    }

    if x_idx < list.binary_data_arrays.len() {
        let bda = &mut list.binary_data_arrays[x_idx];
        match x {
            ArrayData::F32(v) => {
                bda.decoded_binary_f32 = v.clone();
                bda.decoded_binary_f64.clear();
                ensure_float_flag(bda, true);
            }
            ArrayData::F64(v) => {
                bda.decoded_binary_f64 = v.clone();
                bda.decoded_binary_f32.clear();
                ensure_float_flag(bda, false);
            }
        }
    }

    if y_idx < list.binary_data_arrays.len() {
        let bda = &mut list.binary_data_arrays[y_idx];
        match y {
            ArrayData::F32(v) => {
                bda.decoded_binary_f32 = v.clone();
                bda.decoded_binary_f64.clear();
                ensure_float_flag(bda, true);
            }
            ArrayData::F64(v) => {
                bda.decoded_binary_f64 = v.clone();
                bda.decoded_binary_f32.clear();
                ensure_float_flag(bda, false);
            }
        }
    }

    list.count = Some(list.binary_data_arrays.len());
}

#[derive(Debug, Clone, PartialEq)]
pub enum MetadatumValue {
    Number(f64),
    Text(String),
    Empty,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Metadatum {
    pub item_index: u32,
    pub owner_id: u32,
    pub parent_index: u32,
    pub tag_id: TagId,
    pub accession: Option<String>,
    pub unit_accession: Option<String>,
    pub value: MetadatumValue,
}

/// <cvList>
#[inline]
fn parse_cv_list(_schema: &Schema, _bytes: &[u8]) -> Option<CvList> {
    None
}

/// <fileDescription>
#[inline]
fn parse_file_description(_schema: &Schema, _bytes: &[u8]) -> FileDescription {
    FileDescription::default()
}

/// <referenceableParamGroupList>
#[inline]
fn parse_referenceable_param_group_list(
    _schema: &Schema,
    _bytes: &[u8],
) -> Option<ReferenceableParamGroupList> {
    None
}

/// <sampleList>
#[inline]
fn parse_sample_list(_schema: &Schema, _bytes: &[u8]) -> Option<SampleList> {
    None
}

/// <instrumentList>
#[inline]
fn parse_instrument_list(_schema: &Schema, _bytes: &[u8]) -> Option<InstrumentList> {
    None
}

/// <softwareList>
#[inline]
fn parse_software_list(_schema: &Schema, _bytes: &[u8]) -> Option<SoftwareList> {
    None
}

/// <dataProcessingList>
#[inline]
fn parse_data_processing_list(_schema: &Schema, _bytes: &[u8]) -> Option<DataProcessingList> {
    None
}

/// <scanSettingsList>
#[inline]
fn parse_scan_settings_list(_schema: &Schema, _bytes: &[u8]) -> Option<ScanSettingsList> {
    None
}

#[inline]
fn ms_float_param(accession_tail: u32) -> CvParam {
    let name = if accession_tail == ACC_32BIT_FLOAT {
        "32-bit float"
    } else {
        "64-bit float"
    };
    CvParam {
        cv_ref: Some("MS".to_string()),
        accession: Some(format!("MS:{:07}", accession_tail)),
        name: name.to_string(),
        value: None,
        unit_cv_ref: None,
        unit_name: None,
        unit_accession: None,
    }
}

fn parse_metadata_section(
    bytes: &[u8],
    start_off: u64,
    end_off: u64,
    item_count: u32,
    expected_item_count: u32,
    meta_count: u32,
    num_count: u32,
    str_count: u32,
    compression_flag_bit: u8,
    expected_total_meta_len: usize,
    section_name: &str,
) -> Vec<Metadatum> {
    let header = parse_header(&bytes).expect("parse_header failed");

    let c0 = start_off as usize;
    let c1 = end_off as usize;

    assert!(
        c0 < c1,
        "invalid metadata offsets for {section_name}: start >= end"
    );
    assert!(
        c1 <= bytes.len(),
        "invalid metadata offsets for {section_name}: end out of bounds"
    );
    assert_eq!(
        item_count, expected_item_count,
        "test.b64 should contain {expected_item_count} {section_name} items"
    );

    let compressed = (header.reserved_flags & (1u8 << compression_flag_bit)) != 0;
    let slice = &bytes[c0..c1];

    let meta = parse_metadata(
        slice,
        item_count,
        meta_count,
        num_count,
        str_count,
        compressed,
        header.reserved_flags,
    )
    .expect("parse_metadata failed");

    assert_eq!(
        meta.len(),
        expected_total_meta_len,
        "unexpected {section_name} metadata count (expected {expected_total_meta_len} total items)"
    );

    meta
}
