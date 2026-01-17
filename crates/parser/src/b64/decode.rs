use std::collections::{HashMap, HashSet};

use crate::{
    b64::utilities::{
        Header,
        common::*,
        parse_chromatogram_list, parse_cv_and_user_params, parse_cv_list,
        parse_data_processing_list,
        parse_file_description::{
            allowed_from_rows, b000_attr_text, child_params_for_parent, is_child_of,
            parse_file_description, rows_for_owner, unique_ids,
        },
        parse_global_metadata::parse_global_metadata,
        parse_header, parse_instrument_list, parse_metadata, parse_referenceable_param_group_list,
        parse_sample_list, parse_scan_settings_list, parse_software_list, parse_spectrum_list,
    },
    mzml::{
        attr_meta::*,
        schema::{SchemaTree as Schema, TagId, schema},
        structs::*,
    },
};

const INDEX_ENTRY_SIZE: usize = 32;
const BLOCK_DIR_ENTRY_SIZE: usize = 32;

const ARRAY_FILTER_NONE: u8 = 0;
const ARRAY_FILTER_BYTE_SHUFFLE: u8 = 1;

const ACC_MZ_ARRAY: u32 = 1_000_514;
const ACC_INTENSITY_ARRAY: u32 = 1_000_515;
const ACC_TIME_ARRAY: u32 = 1_000_595;

const ACC_32BIT_FLOAT: u32 = 1_000_521;
const ACC_64BIT_FLOAT: u32 = 1_000_523;

const HDR_FLAG_GLOBAL_META_COMP: u8 = 1 << 6;

pub fn decode(bytes: &[u8]) -> Result<MzML, String> {
    let schema = schema();
    let header = parse_header(bytes)?;
    let global_meta = parse_global_metadata_section(bytes, &header)?;
    let global_child_index = ChildIndex::new(&global_meta);

    let cv_list = parse_cv_list(&schema, &global_meta, &global_child_index);

    Ok(MzML {
        cv_list,
        file_description: parse_file_description(&schema, &global_meta, &global_child_index)
            .expect("missing <fileDescription> in global metadata"),
        referenceable_param_group_list: parse_referenceable_param_group_list(
            &schema,
            &global_meta,
            &global_child_index,
        ),
        sample_list: parse_sample_list(&schema, &global_meta, &global_child_index),
        instrument_list: parse_instrument_list(&schema, &global_meta, &global_child_index),
        software_list: parse_software_list(&schema, &global_meta, &global_child_index),
        data_processing_list: parse_data_processing_list(
            &schema,
            &global_meta,
            &global_child_index,
        ),
        scan_settings_list: parse_scan_settings_list(&schema, &global_meta, &global_child_index),
        run: parse_run(&schema, bytes, &header, &global_meta)?,
    })
}

/// <run>
#[inline]
fn parse_run(
    schema: &Schema,
    bytes: &[u8],
    header: &Header,
    global_meta: &[Metadatum],
) -> Result<Run, String> {
    let spec_meta = parse_metadata_section(
        bytes,
        header,
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
    );

    let chrom_meta = parse_metadata_section(
        bytes,
        header,
        header.off_chrom_meta,
        header.off_global_meta,
        header.chrom_count,
        header.chrom_meta_count,
        header.chrom_num_count,
        header.chrom_str_count,
        5,
    );

    let run_child_index = ChildIndex::new(global_meta);

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> =
        HashMap::with_capacity(global_meta.len() / 2 + 1);
    for m in global_meta {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let run_id = global_meta
        .iter()
        .find(|m| m.tag_id == TagId::Run)
        .map(|m| m.owner_id)
        .unwrap_or(0);

    let run_rows = rows_for_owner(&owner_rows, run_id);

    let id = b000_attr_text(run_rows, ACC_ATTR_ID).unwrap_or_default();
    let start_time_stamp =
        b000_attr_text(run_rows, ACC_ATTR_START_TIME_STAMP).filter(|s| !s.is_empty());

    let default_instrument_configuration_ref =
        b000_attr_text(run_rows, ACC_ATTR_DEFAULT_INSTRUMENT_CONFIGURATION_REF)
            .or_else(|| b000_attr_text(run_rows, ACC_ATTR_INSTRUMENT_CONFIGURATION_REF))
            .filter(|s| !s.is_empty());

    let sample_ref = b000_attr_text(run_rows, ACC_ATTR_SAMPLE_REF).filter(|s| !s.is_empty());

    let allowed_run_schema: HashSet<&str> = find_node_by_tag(schema, TagId::Run)
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut params_meta = Vec::with_capacity(
        run_rows.len() + child_params_for_parent(&owner_rows, &run_child_index, run_id).len(),
    );
    params_meta.extend(run_rows.iter().copied());
    params_meta.extend(child_params_for_parent(
        &owner_rows,
        &run_child_index,
        run_id,
    ));

    let (cv_params, user_params) = if allowed_run_schema.is_empty() {
        let allowed_meta = allowed_from_rows(&params_meta);
        parse_cv_and_user_params(&allowed_meta, &params_meta)
    } else {
        parse_cv_and_user_params(&allowed_run_schema, &params_meta)
    };

    let source_file_ref_list =
        parse_source_file_ref_list(&owner_rows, &run_child_index, global_meta, run_id);

    let spec_child_index = ChildIndex::new(&spec_meta);
    let chrom_child_index = ChildIndex::new(&chrom_meta);

    Ok(Run {
        id,
        start_time_stamp,
        default_instrument_configuration_ref,
        sample_ref,
        cv_params,
        user_params,
        source_file_ref_list,
        spectrum_list: parse_spectrum_list(schema, &spec_meta, &spec_child_index),
        chromatogram_list: parse_chromatogram_list(schema, &chrom_meta, &chrom_child_index),
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
    let need = (count as u64)
        .checked_mul(INDEX_ENTRY_SIZE as u64)
        .ok_or_else(|| "chrom index size overflow".to_string())?;
    let raw = slice_at(bytes, header.off_chrom_index, need, "chrom index")?;

    let mut pos = 0usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        out.push(ChromIndexEntry {
            time_element_off: read_u64_le_at(raw, &mut pos, "time_element_off")?,
            inten_element_off: read_u64_le_at(raw, &mut pos, "inten_element_off")?,
            time_element_len: read_u32_le_at(raw, &mut pos, "time_element_len")?,
            inten_element_len: read_u32_le_at(raw, &mut pos, "inten_element_len")?,
            time_block_id: read_u32_le_at(raw, &mut pos, "time_block_id")?,
            inten_block_id: read_u32_le_at(raw, &mut pos, "inten_block_id")?,
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

        let comp_off = usize::try_from(e.comp_off).map_err(|_| "comp_off overflow".to_string())?;
        let comp_size =
            usize::try_from(e.comp_size).map_err(|_| "comp_size overflow".to_string())?;
        let expected =
            usize::try_from(e.uncomp_bytes).map_err(|_| "uncomp_bytes overflow".to_string())?;

        let start = self
            .comp_buf_start
            .checked_add(comp_off)
            .ok_or_else(|| "comp start overflow".to_string())?;
        let end = start
            .checked_add(comp_size)
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
            std::mem::swap(&mut out, &mut self.scratch);
            self.scratch.clear();
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
    let n = raw.len() / 4;
    let mut out = vec![0.0f32; n];
    for (i, c) in raw.chunks_exact(4).enumerate() {
        out[i] = f32::from_le_bytes(c.try_into().unwrap());
    }
    out
}

#[inline]
fn bytes_to_f64_vec(raw: &[u8]) -> Vec<f64> {
    let n = raw.len() / 8;
    let mut out = vec![0.0f64; n];
    for (i, c) in raw.chunks_exact(8).enumerate() {
        out[i] = f64::from_le_bytes(c.try_into().unwrap());
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
        let cur = starts[bi];
        let v = e.mz_element_off;
        starts[bi] = if cur < v { cur } else { v };
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
        let cur = starts[bi];
        let v = e.inten_element_off;
        starts[bi] = if cur < v { cur } else { v };
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
        let cur = starts[bi];
        let v = e.time_element_off;
        starts[bi] = if cur < v { cur } else { v };
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
        let cur = starts[bi];
        let v = e.inten_element_off;
        starts[bi] = if cur < v { cur } else { v };
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
    Ok(if elem_size == 4 {
        ArrayData::F32(bytes_to_f32_vec(slice))
    } else {
        ArrayData::F64(bytes_to_f64_vec(slice))
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
    if list.binary_data_arrays.is_empty() {
        return;
    }

    let mut x_i = None;
    let mut y_i = None;

    for (i, bda) in list.binary_data_arrays.iter().enumerate() {
        if x_i.is_none() && bda_has_array_kind(bda, x_kind) {
            x_i = Some(i);
        }
        if y_i.is_none() && bda_has_array_kind(bda, y_kind) {
            y_i = Some(i);
        }
        if x_i.is_some() && y_i.is_some() {
            break;
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

    if let Some(bda) = list.binary_data_arrays.get_mut(x_idx) {
        match x {
            ArrayData::F32(v) => {
                bda.decoded_binary_f32.clone_from(v);
                bda.decoded_binary_f64.clear();
                ensure_float_flag(bda, true);
            }
            ArrayData::F64(v) => {
                bda.decoded_binary_f64.clone_from(v);
                bda.decoded_binary_f32.clear();
                ensure_float_flag(bda, false);
            }
        }
    }

    if let Some(bda) = list.binary_data_arrays.get_mut(y_idx) {
        match y {
            ArrayData::F32(v) => {
                bda.decoded_binary_f32.clone_from(v);
                bda.decoded_binary_f64.clear();
                ensure_float_flag(bda, true);
            }
            ArrayData::F64(v) => {
                bda.decoded_binary_f64.clone_from(v);
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
    header: &Header,
    start_off: u64,
    end_off: u64,
    item_count: u32,
    meta_count: u32,
    num_count: u32,
    str_count: u32,
    compression_flag_bit: u8,
) -> Vec<Metadatum> {
    let c0 = start_off as usize;
    let c1 = end_off as usize;

    assert!(c0 < c1, "invalid metadata offsets: start >= end");
    assert!(
        c1 <= bytes.len(),
        "invalid metadata offsets: end out of bounds"
    );

    let compressed = (header.reserved_flags & (1u8 << compression_flag_bit)) != 0;
    let slice = &bytes[c0..c1];

    parse_metadata(
        slice,
        item_count,
        meta_count,
        num_count,
        str_count,
        compressed,
        header.reserved_flags,
    )
    .expect("parse_metadata failed")
}

fn parse_global_metadata_section(bytes: &[u8], header: &Header) -> Result<Vec<Metadatum>, String> {
    let start = header.off_global_meta as usize;
    let end = header.off_container_spect_x as usize;

    if start >= end {
        return Err("invalid global metadata offsets: start >= end".to_string());
    }
    if end > bytes.len() {
        return Err("invalid global metadata offsets: end out of bounds".to_string());
    }

    let slice = &bytes[start..end];
    let compressed = (header.reserved_flags & HDR_FLAG_GLOBAL_META_COMP) != 0;

    parse_global_metadata(
        slice,
        0,
        header.global_meta_count,
        header.global_num_count,
        header.global_str_count,
        compressed,
        header.reserved_flags,
    )
}

#[inline]
fn parse_source_file_ref_list(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    run_id: u32,
) -> Option<SourceFileRefList> {
    let mut list_ids = unique_ids(child_index.ids(run_id, TagId::SourceFileRefList));

    if list_ids.is_empty() {
        list_ids = ordered_unique_owner_ids(metadata, TagId::SourceFileRefList);

        if run_id != 0 && list_ids.len() > 1 {
            let filtered: Vec<u32> = list_ids
                .iter()
                .copied()
                .filter(|&id| is_child_of(owner_rows, id, run_id))
                .collect();
            if !filtered.is_empty() {
                list_ids = filtered;
            }
        }
    }

    let list_id = list_ids.first().copied()?;
    let list_rows = rows_for_owner(owner_rows, list_id);

    let mut count = b000_attr_text(list_rows, ACC_ATTR_COUNT).and_then(|s| s.parse::<usize>().ok());

    let mut ref_ids = unique_ids(child_index.ids(list_id, TagId::SourceFileRef));

    if ref_ids.is_empty() {
        ref_ids = ordered_unique_owner_ids(metadata, TagId::SourceFileRef);

        if ref_ids.len() > 1 {
            let filtered: Vec<u32> = ref_ids
                .iter()
                .copied()
                .filter(|&id| is_child_of(owner_rows, id, list_id))
                .collect();
            if !filtered.is_empty() {
                ref_ids = filtered;
            }
        }
    }

    let mut source_file_refs = Vec::with_capacity(ref_ids.len());
    for rid in ref_ids {
        let rows = rows_for_owner(owner_rows, rid);
        if let Some(r) = b000_attr_text(rows, ACC_ATTR_REF) {
            if !r.is_empty() {
                source_file_refs.push(SourceFileRef { r#ref: r });
            }
        }
    }

    if source_file_refs.is_empty() {
        return None;
    }

    if count.is_none() {
        count = Some(source_file_refs.len());
    }

    Some(SourceFileRefList {
        count,
        source_file_refs,
    })
}
