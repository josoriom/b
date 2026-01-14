use std::{io::Cursor, str};

use miniz_oxide::inflate::decompress_to_vec_zlib;
use zstd::{bulk::decompress as zstd_decompress, stream::decode_all as zstd_decode_all};

use crate::mzml::{
    attr_meta::*,
    cv_table,
    schema::{SchemaNode, SchemaTree, SchemaTree as Schema, TagId, schema},
    structs::*,
};

const HEADER_SIZE: usize = 192;
const INDEX_ENTRY_SIZE: usize = 32;
const BLOCK_DIR_ENTRY_SIZE: usize = 32;

const HDR_CODEC_MASK: u8 = 0x0F;
const HDR_CODEC_ZLIB: u8 = 0;
const HDR_CODEC_ZSTD: u8 = 1;

const HDR_FLAG_SPEC_META_COMP: u8 = 1 << 4;
const HDR_FLAG_CHROM_META_COMP: u8 = 1 << 5;
const HDR_FLAG_GLOBAL_META_COMP: u8 = 1 << 6;

const HDR_ARRAY_FILTER_OFF: usize = 178;
const ARRAY_FILTER_NONE: u8 = 0;
const ARRAY_FILTER_BYTE_SHUFFLE: u8 = 1;

const TAG_UNKNOWN: u8 = TagId::Unknown as u8;
const TAG_BINARY_DATA_ARRAY: u8 = TagId::BinaryDataArray as u8;
const TAG_FILE_DESCRIPTION: u8 = TagId::FileDescription as u8;
const TAG_SOURCE_FILE: u8 = TagId::SourceFile as u8;
const TAG_CONTACT: u8 = TagId::Contact as u8;
const TAG_SCAN_SETTINGS: u8 = TagId::ScanSettings as u8;
const TAG_TARGET: u8 = TagId::Target as u8;
const TAG_INSTRUMENT: u8 = TagId::Instrument as u8;
const TAG_COMPONENT_SOURCE: u8 = TagId::ComponentSource as u8;
const TAG_COMPONENT_ANALYZER: u8 = TagId::ComponentAnalyzer as u8;
const TAG_COMPONENT_DETECTOR: u8 = TagId::ComponentDetector as u8;
const TAG_SPECTRUM_DESCRIPTION: u8 = TagId::SpectrumDescription as u8;
const TAG_SCAN: u8 = TagId::Scan as u8;
const TAG_SCAN_WINDOW: u8 = TagId::ScanWindow as u8;
const TAG_ISOLATION_WINDOW: u8 = TagId::IsolationWindow as u8;
const TAG_SELECTED_ION: u8 = TagId::SelectedIon as u8;
const TAG_ACTIVATION: u8 = TagId::Activation as u8;
const TAG_BINARY_DATA_ARRAY_LIST: u8 = TagId::BinaryDataArrayList as u8;

#[inline]
fn tag_from_u8(b: u8) -> TagId {
    match b {
        TAG_FILE_DESCRIPTION => TagId::FileDescription,
        TAG_SOURCE_FILE => TagId::SourceFile,
        TAG_CONTACT => TagId::Contact,
        TAG_SCAN_SETTINGS => TagId::ScanSettings,
        TAG_TARGET => TagId::Target,
        TAG_INSTRUMENT => TagId::Instrument,
        TAG_COMPONENT_SOURCE => TagId::ComponentSource,
        TAG_COMPONENT_ANALYZER => TagId::ComponentAnalyzer,
        TAG_COMPONENT_DETECTOR => TagId::ComponentDetector,
        TAG_SPECTRUM_DESCRIPTION => TagId::SpectrumDescription,
        TAG_SCAN => TagId::Scan,
        TAG_SCAN_WINDOW => TagId::ScanWindow,
        TAG_ISOLATION_WINDOW => TagId::IsolationWindow,
        TAG_SELECTED_ION => TagId::SelectedIon,
        TAG_ACTIVATION => TagId::Activation,
        TAG_BINARY_DATA_ARRAY_LIST => TagId::BinaryDataArrayList,
        TAG_BINARY_DATA_ARRAY => TagId::BinaryDataArray,
        TAG_UNKNOWN => TagId::Unknown,
        _ => TagId::Unknown,
    }
}

#[derive(Clone)]
struct MetaCvParam {
    tag: TagId,
    cv_ref_code: u8,
    accession_tail: u32,
    cv: CvParam,
}

#[inline]
fn effective_tag(schema: Option<&Schema>, m: &MetaCvParam, hint: Option<TagId>) -> TagId {
    let _ = schema;
    let _ = hint;
    if m.tag != TagId::Unknown {
        m.tag
    } else {
        TagId::Unknown
    }
}

#[derive(Clone, Copy)]
struct SchemaNodes<'a> {
    file_description: Option<&'a SchemaNode>,
    source_file: Option<&'a SchemaNode>,
    contact: Option<&'a SchemaNode>,

    scan_settings: Option<&'a SchemaNode>,
    target: Option<&'a SchemaNode>,

    instrument: Option<&'a SchemaNode>,
    component_source: Option<&'a SchemaNode>,
    component_analyzer: Option<&'a SchemaNode>,
    component_detector: Option<&'a SchemaNode>,

    spectrum_description: Option<&'a SchemaNode>,
    scan: Option<&'a SchemaNode>,
    scan_window: Option<&'a SchemaNode>,

    isolation_window: Option<&'a SchemaNode>,
    selected_ion: Option<&'a SchemaNode>,
    activation: Option<&'a SchemaNode>,
}

#[inline]
fn child_node<'a>(parent: Option<&'a SchemaNode>, tag: TagId) -> Option<&'a SchemaNode> {
    let p = parent?;
    let key = p.child_key_for_tag(tag)?;
    p.children.get(key)
}

fn find_descendant<'a>(node: Option<&'a SchemaNode>, tag: TagId) -> Option<&'a SchemaNode> {
    let n = node?;
    if n.self_tags.iter().any(|&t| t == tag) {
        return Some(n);
    }
    for child in n.children.values() {
        if let Some(found) = find_descendant(Some(child), tag) {
            return Some(found);
        }
    }
    None
}

impl<'a> SchemaNodes<'a> {
    #[inline]
    fn new(tree: Option<&'a SchemaTree>) -> Self {
        let file_description = tree.and_then(|t| t.root_by_tag(TagId::FileDescription));
        let source_file = child_node(file_description, TagId::SourceFile)
            .or_else(|| tree.and_then(|t| t.root_by_tag(TagId::SourceFile)));
        let contact = child_node(file_description, TagId::Contact)
            .or_else(|| tree.and_then(|t| t.root_by_tag(TagId::Contact)));

        let scan_settings = tree.and_then(|t| t.root_by_tag(TagId::ScanSettings));
        let target = child_node(scan_settings, TagId::Target)
            .or_else(|| tree.and_then(|t| t.root_by_tag(TagId::Target)));

        let instrument = tree.and_then(|t| t.root_by_tag(TagId::Instrument));
        let component_source = tree
            .and_then(|t| t.root_by_tag(TagId::ComponentSource))
            .or_else(|| find_descendant(instrument, TagId::ComponentSource));
        let component_analyzer = tree
            .and_then(|t| t.root_by_tag(TagId::ComponentAnalyzer))
            .or_else(|| find_descendant(instrument, TagId::ComponentAnalyzer));
        let component_detector = tree
            .and_then(|t| t.root_by_tag(TagId::ComponentDetector))
            .or_else(|| find_descendant(instrument, TagId::ComponentDetector));

        let spectrum_description = tree.and_then(|t| t.root_by_tag(TagId::SpectrumDescription));
        let scan = tree
            .and_then(|t| t.root_by_tag(TagId::Scan))
            .or_else(|| find_descendant(spectrum_description, TagId::Scan));
        let scan_window = tree
            .and_then(|t| t.root_by_tag(TagId::ScanWindow))
            .or_else(|| find_descendant(scan, TagId::ScanWindow));

        let isolation_window = tree.and_then(|t| t.root_by_tag(TagId::IsolationWindow));
        let selected_ion = tree.and_then(|t| t.root_by_tag(TagId::SelectedIon));
        let activation = tree.and_then(|t| t.root_by_tag(TagId::Activation));

        Self {
            file_description,
            source_file,
            contact,
            scan_settings,
            target,
            instrument,
            component_source,
            component_analyzer,
            component_detector,
            spectrum_description,
            scan,
            scan_window,
            isolation_window,
            selected_ion,
            activation,
        }
    }
}

#[inline]
fn node_has_accession(node: &SchemaNode, acc: &str) -> bool {
    node.accessions.iter().any(|a| a == acc)
}

#[inline]
fn schema_tag_in_nodes(
    schema: Option<&Schema>,
    m: &MetaCvParam,
    hint: Option<TagId>,
    candidates: &[(TagId, Option<&SchemaNode>)],
) -> TagId {
    if m.tag != TagId::Unknown {
        return m.tag;
    }

    if let Some(acc) = m.cv.accession.as_deref() {
        let mut found = TagId::Unknown;
        for (tag, node) in candidates {
            let Some(node) = node else { continue };
            if node_has_accession(node, acc) {
                if found != TagId::Unknown && found != *tag {
                    found = TagId::Unknown;
                    break;
                }
                found = *tag;
            }
        }
        if found != TagId::Unknown {
            return found;
        }
    }

    effective_tag(schema, m, hint)
}

#[inline]
fn split_attr_value(params: &mut Vec<MetaCvParam>, attr_tail: u32) -> Option<String> {
    let mut v: Option<String> = None;
    params.retain(|m| {
        if is_attr_param(&m.cv) && parse_acc_tail(m.cv.accession.as_deref()) == attr_tail {
            v = attr_string_value(&m.cv);
            false
        } else {
            true
        }
    });
    v
}

#[inline]
fn split_id_attr(params: &mut Vec<MetaCvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_ID)
}

#[inline]
fn split_name_attr(params: &mut Vec<MetaCvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_NAME)
}

#[inline]
fn split_version_attr(params: &mut Vec<MetaCvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_VERSION)
}

#[inline]
fn split_instrument_configuration_ref_attr(params: &mut Vec<MetaCvParam>) -> Option<String> {
    split_attr_value(params, ACC_ATTR_INSTRUMENT_CONFIGURATION_REF)
}

#[inline]
fn collect_ref_attrs(params: &mut Vec<MetaCvParam>) -> Vec<String> {
    let mut refs = Vec::<String>::new();
    params.retain(|m| {
        if is_attr_param(&m.cv) && parse_acc_tail(m.cv.accession.as_deref()) == ACC_ATTR_REF {
            if let Some(v) = attr_string_value(&m.cv) {
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
fn base64_encoded_len(byte_len: usize) -> usize {
    ((byte_len + 2) / 3) * 4
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
fn cv_table_name(key: &str) -> Option<String> {
    cv_table::get(key)
        .and_then(|v| {
            v.as_str()
                .or_else(|| v.get("name").and_then(|n| n.as_str()))
        })
        .map(|s| s.to_string())
}

#[inline]
fn ms_tail(schema: Option<&Schema>, name: &str) -> Option<u32> {
    let _ = schema;
    match name {
        "m/z array" => Some(1_000_514),
        "intensity array" => Some(1_000_515),
        "time array" => Some(1_000_595),
        "32-bit float" => Some(1_000_521),
        "64-bit float" => Some(1_000_523),
        "no compression" => Some(1_000_576),
        "zlib compression" => Some(1_000_574),
        _ => None,
    }
}

#[derive(Clone, Copy)]
struct MsBinaryTails {
    mz_array: Option<u32>,
    intensity_array: Option<u32>,
    time_array: Option<u32>,
    f32: Option<u32>,
    f64: Option<u32>,
    no_compression: Option<u32>,
    zlib_compression: Option<u32>,
}

impl MsBinaryTails {
    #[inline]
    fn float_tail(self, fmt: u8) -> Option<u32> {
        match fmt {
            1 => self.f32,
            2 => self.f64,
            _ => None,
        }
    }
}

#[inline]
fn ms_cv_param(schema: Option<&Schema>, accession_tail: Option<u32>) -> CvParam {
    let _ = schema;
    let accession = accession_tail.map(|t| format!("MS:{:07}", t));
    let mut p = CvParam {
        cv_ref: Some("MS".to_string()),
        accession,
        name: String::new(),
        value: Some(String::new()),
        unit_cv_ref: None,
        unit_name: None,
        unit_accession: None,
    };

    if let Some(k) = p.accession.as_deref() {
        p.name = cv_table_name(k).unwrap_or_default();
    }

    p
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

#[inline]
fn filter_spectrum_top_level_cv_params(params: &mut Vec<CvParam>) {
    params.retain(|p| !(p.cv_ref.as_deref() == Some("MS") && p.name.is_empty()));
}

#[inline]
fn infer_term_cv_ref_code(mut code: u8, tail: u32) -> u8 {
    if tail == 0 || cv_ref_prefix_from_code(code).is_some() {
        return code;
    }
    let ms_key = format!("MS:{:07}", tail);
    let uo_key = format!("UO:{:07}", tail);
    if cv_table::get(&ms_key).is_some() {
        code = CV_CODE_MS;
    } else if cv_table::get(&uo_key).is_some() {
        code = CV_CODE_UO;
    } else {
        code = CV_CODE_B000;
    }
    code
}

#[inline]
fn infer_unit_cv_ref_code(mut code: u8, tail: u32) -> u8 {
    if tail == 0 || cv_ref_prefix_from_code(code).is_some() {
        return code;
    }
    let uo_key = format!("UO:{:07}", tail);
    let ms_key = format!("MS:{:07}", tail);
    if cv_table::get(&uo_key).is_some() {
        code = CV_CODE_UO;
    } else if cv_table::get(&ms_key).is_some() {
        code = CV_CODE_MS;
    } else {
        code = CV_CODE_B000;
    }
    code
}

#[inline]
fn attr_accession(tail: u32) -> Option<String> {
    (tail >= 9_900_000).then(|| format!("{}:{:07}", CV_REF_ATTR, tail))
}

#[inline]
fn cv_code_score(s: &[u8]) -> usize {
    let mut n = 0usize;
    for &b in s {
        if cv_ref_prefix_from_code(b).is_some() || b == CV_CODE_B000 {
            n += 1;
        }
    }
    n
}

#[inline]
fn is_target_boundary(p: &CvParam) -> bool {
    const ACC_ISO_TARGET_MZ: u32 = 1_000_827; // MS:1000827
    parse_acc_tail(p.accession.as_deref()) == ACC_ISO_TARGET_MZ
}

/// <scanSettings>
#[inline]
fn build_scan_settings_from_cv_params(
    schema: Option<&Schema>,
    schema_nodes: SchemaNodes<'_>,
    mut params: Vec<MetaCvParam>,
) -> ScanSettings {
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

    for m in params.drain(..) {
        let tag = schema_tag_in_nodes(
            schema,
            &m,
            Some(TagId::Target),
            &[(TagId::Target, schema_nodes.target)],
        );
        if tag == TagId::Target {
            target_cv.push(m.cv);
        } else {
            rest.push(m);
        }
    }
    let rest = rest.into_iter().map(|m| m.cv).collect::<Vec<_>>();

    let target_list = if target_cv.is_empty() {
        None
    } else {
        let mut targets: Vec<Target> = Vec::new();
        let mut cur: Vec<CvParam> = Vec::new();

        for p in target_cv {
            if is_target_boundary(&p) && !cur.is_empty() {
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

/// <fileDescription>
fn split_file_description_from_cv_params(
    schema: Option<&Schema>,
    schema_nodes: SchemaNodes<'_>,
    mut params: Vec<MetaCvParam>,
) -> FileDescription {
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

    for m in params.drain(..) {
        let tag = schema_tag_in_nodes(
            schema,
            &m,
            None,
            &[
                (TagId::SourceFile, schema_nodes.source_file),
                (TagId::Contact, schema_nodes.contact),
            ],
        );

        if tag == TagId::SourceFile {
            if is_attr_param(&m.cv) {
                let tail = parse_acc_tail(m.cv.accession.as_deref());
                let v = attr_string_value(&m.cv).unwrap_or_default();
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
                    continue;
                } else if tail == ACC_ATTR_NAME {
                    cur_name = v;
                    continue;
                } else if tail == ACC_ATTR_LOCATION {
                    cur_location = v;
                    continue;
                }
            }
            cur_cv.push(m.cv);
            continue;
        }

        if tag == TagId::Contact {
            contact_cv.push(m.cv);
            continue;
        }

        file_content_cv.push(m.cv);
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

/// <precursorList>
fn infer_precursor_from_cv(
    schema: Option<&Schema>,
    schema_nodes: SchemaNodes<'_>,
    params: &mut Vec<MetaCvParam>,
) -> Option<Precursor> {
    let mut spectrum_ref: Option<String> = None;
    let mut iso = Vec::<CvParam>::new();
    let mut sel = Vec::<CvParam>::new();
    let mut act = Vec::<CvParam>::new();

    let mut rest = Vec::<MetaCvParam>::with_capacity(params.len());
    for m in params.drain(..) {
        if is_attr_param(&m.cv)
            && parse_acc_tail(m.cv.accession.as_deref()) == ACC_ATTR_SPECTRUM_REF
        {
            if spectrum_ref.is_none() {
                spectrum_ref = attr_string_value(&m.cv);
            }
            continue;
        }

        let tag0 = schema_tag_in_nodes(
            schema,
            &m,
            None,
            &[
                (TagId::SelectedIon, schema_nodes.selected_ion),
                (TagId::IsolationWindow, schema_nodes.isolation_window),
                (TagId::Activation, schema_nodes.activation),
            ],
        );

        let tag = if tag0 == TagId::Unknown && m.cv_ref_code == CV_CODE_MS {
            match m.accession_tail {
                1_000_040 | 1_000_041 | 1_000_042 => TagId::SelectedIon,
                1_000_827 | 1_000_828 | 1_000_829 => TagId::IsolationWindow,
                1_000_133 | 1_000_045 => TagId::Activation,
                _ => TagId::Unknown,
            }
        } else {
            tag0
        };

        if tag == TagId::IsolationWindow {
            iso.push(m.cv);
        } else if tag == TagId::SelectedIon {
            sel.push(m.cv);
        } else if tag == TagId::Activation {
            act.push(m.cv);
        } else {
            rest.push(m);
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

    Some(Precursor {
        spectrum_ref,
        isolation_window,
        selected_ion_list,
        activation,
        ..Default::default()
    })
}

/// <precursorList>
fn infer_precursor_list_from_spectrum_cv(
    schema: Option<&Schema>,
    schema_nodes: SchemaNodes<'_>,
    params: &mut Vec<MetaCvParam>,
) -> Option<PrecursorList> {
    infer_precursor_from_cv(schema, schema_nodes, params).map(|p| PrecursorList {
        count: Some(1),
        precursors: vec![p],
        ..Default::default()
    })
}

/// <spectrumDescription>
fn infer_spectrum_description_from_spectrum_cv(
    schema: Option<&Schema>,
    schema_nodes: SchemaNodes<'_>,
    params: &mut Vec<MetaCvParam>,
) -> Option<SpectrumDescription> {
    let mut desc = Vec::<CvParam>::new();
    let mut scan = Vec::<CvParam>::new();
    let mut scan_window = Vec::<CvParam>::new();

    let mut rest = Vec::<MetaCvParam>::with_capacity(params.len());
    for m in params.drain(..) {
        let tag = schema_tag_in_nodes(
            schema,
            &m,
            None,
            &[
                (
                    TagId::SpectrumDescription,
                    schema_nodes.spectrum_description,
                ),
                (TagId::Scan, schema_nodes.scan),
                (TagId::ScanWindow, schema_nodes.scan_window),
            ],
        );

        if tag == TagId::SpectrumDescription {
            desc.push(m.cv);
        } else if tag == TagId::Scan {
            scan.push(m.cv);
        } else if tag == TagId::ScanWindow {
            scan_window.push(m.cv);
        } else {
            rest.push(m);
        }
    }
    *params = rest;

    if desc.is_empty() && scan.is_empty() && scan_window.is_empty() {
        return None;
    }

    let scan_window_list = (!scan_window.is_empty()).then(|| {
        let mut windows: Vec<ScanWindow> = Vec::new();
        let mut cur: Vec<CvParam> = Vec::new();

        for p in scan_window {
            if p.name == "scan window lower limit" && !cur.is_empty() {
                windows.push(ScanWindow {
                    cv_params: cur,
                    ..Default::default()
                });
                cur = Vec::new();
            }
            cur.push(p);
        }

        if !cur.is_empty() {
            windows.push(ScanWindow {
                cv_params: cur,
                ..Default::default()
            });
        }

        ScanWindowList {
            count: Some(windows.len()),
            scan_windows: windows,
        }
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
            ..Default::default()
        }),
        ..Default::default()
    })
}

/// <spectrum>
fn split_spectrum_attrs(
    item_idx: usize,
    x_len: u32,
    params: Vec<MetaCvParam>,
) -> (String, Option<u32>, Option<usize>, Vec<MetaCvParam>) {
    split_item_attrs(item_idx, x_len, params, "spectrum")
}

/// <chromatogram>
fn split_chromatogram_attrs(
    item_idx: usize,
    x_len: u32,
    params: Vec<MetaCvParam>,
) -> (String, Option<u32>, Option<usize>, Vec<MetaCvParam>) {
    split_item_attrs(item_idx, x_len, params, "chromatogram")
}

#[inline]
fn split_item_attrs(
    item_idx: usize,
    x_len: u32,
    params: Vec<MetaCvParam>,
    prefix: &str,
) -> (String, Option<u32>, Option<usize>, Vec<MetaCvParam>) {
    let mut id: Option<String> = None;
    let mut index: Option<u32> = None;
    let mut default_array_length: Option<usize> = None;

    let mut out = Vec::with_capacity(params.len());
    for m in params {
        if is_attr_param(&m.cv) {
            let tail = parse_acc_tail(m.cv.accession.as_deref());
            if tail == ACC_ATTR_ID {
                id = attr_string_value(&m.cv);
                continue;
            }
            if tail == ACC_ATTR_INDEX {
                index = Some(
                    attr_string_value(&m.cv)
                        .and_then(|s| s.parse::<u32>().ok())
                        .unwrap_or(item_idx as u32),
                );
                continue;
            }
            if tail == ACC_ATTR_DEFAULT_ARRAY_LENGTH {
                default_array_length = Some(
                    attr_string_value(&m.cv)
                        .and_then(|s| s.parse::<usize>().ok())
                        .unwrap_or(x_len as usize),
                );
                continue;
            }
        }
        out.push(m);
    }

    (
        id.unwrap_or_else(|| format!("{}_{}", prefix, item_idx)),
        Some(index.unwrap_or(item_idx as u32)),
        Some(default_array_length.unwrap_or(x_len as usize)),
        out,
    )
}

/// <mzML>
pub fn decode(bytes: &[u8]) -> Result<MzML, String> {
    let schema = Some(schema());
    let schema_nodes = SchemaNodes::new(schema);

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
        schema,
        spec_meta_bytes.as_slice(),
        spectrum_count,
        spec_meta_count,
        spec_num_count,
        spec_str_count,
    )?;
    let chrom_meta_by_item = decode_meta_block(
        schema,
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
        run,
    ) = decode_global_meta_structs(
        schema,
        schema_nodes,
        global_meta_bytes.as_slice(),
        global_meta_count,
        global_num_count,
        global_str_count,
    )?;

    let tails = MsBinaryTails {
        mz_array: ms_tail(schema, "m/z array"),
        intensity_array: ms_tail(schema, "intensity array"),
        time_array: ms_tail(schema, "time array"),
        f32: ms_tail(schema, "32-bit float"),
        f64: ms_tail(schema, "64-bit float"),
        no_compression: ms_tail(schema, "no compression"),
        zlib_compression: ms_tail(schema, "zlib compression"),
    };

    let mut spectra = Vec::with_capacity(spectrum_count as usize);
    for (i, mut item_params) in spec_meta_by_item.into_iter().enumerate() {
        let (x_off, y_off, x_len, y_len, x_block, y_block) =
            read_index_entry_with_blocks(spectrum_index_bytes, i)?;

        let mz_bytes = spect_x_container.slice_elems(x_block, x_off, x_len)?;
        let in_bytes = spect_y_container.slice_elems(y_block, y_off, y_len)?;

        let (mz_f32, mz_f64) = decode_array_by_fmt_from_bytes(mz_bytes, spect_x_fmt)?;
        let (in_f32, in_f64) = decode_array_by_fmt_from_bytes(in_bytes, spect_y_fmt)?;

        let arrays = [
            (BinaryArrayRole::Mz, mz_f32, mz_f64),
            (BinaryArrayRole::Intensity, in_f32, in_f64),
        ];

        let binary_data_array_list =
            parse_binary_data_array_list(schema, tails, &mut item_params, &arrays)?;

        let (id, index, default_array_length, mut item_params) =
            split_spectrum_attrs(i, x_len, item_params);

        let precursor_list =
            infer_precursor_list_from_spectrum_cv(schema, schema_nodes, &mut item_params);
        let spectrum_description =
            infer_spectrum_description_from_spectrum_cv(schema, schema_nodes, &mut item_params);

        let mut spectrum_params = item_params.into_iter().map(|m| m.cv).collect::<Vec<_>>();
        filter_spectrum_top_level_cv_params(&mut spectrum_params);

        spectra.push(Spectrum {
            id,
            index,
            default_array_length,
            cv_params: spectrum_params,
            spectrum_description,
            precursor_list,
            binary_data_array_list: binary_data_array_list,
            ..Default::default()
        });
    }

    let mut chromatograms = Vec::with_capacity(chrom_count as usize);
    for (j, mut item_params) in chrom_meta_by_item.into_iter().enumerate() {
        let (x_off, y_off, x_len, y_len, x_block, y_block) =
            read_index_entry_with_blocks(chromatogram_index_bytes, j)?;

        let t_bytes = chrom_x_container.slice_elems(x_block, x_off, x_len)?;
        let in_bytes = chrom_y_container.slice_elems(y_block, y_off, y_len)?;

        let (t_f32, t_f64) = decode_array_by_fmt_from_bytes(t_bytes, chrom_x_fmt)?;
        let (in_f32, in_f64) = decode_array_by_fmt_from_bytes(in_bytes, chrom_y_fmt)?;

        let arrays = [
            (BinaryArrayRole::Time, t_f32, t_f64),
            (BinaryArrayRole::Intensity, in_f32, in_f64),
        ];

        let binary_data_array_list =
            parse_binary_data_array_list(schema, tails, &mut item_params, &arrays)?;

        let (id, index, default_array_length, mut item_params) =
            split_chromatogram_attrs(j, x_len, item_params);

        let precursor = infer_precursor_from_cv(schema, schema_nodes, &mut item_params);

        let chrom_params = item_params.into_iter().map(|m| m.cv).collect::<Vec<_>>();

        chromatograms.push(Chromatogram {
            id,
            index,
            default_array_length,
            cv_params: chrom_params,
            precursor,
            binary_data_array_list: binary_data_array_list,
            ..Default::default()
        });
    }

    let source_file_ids: Vec<String> = file_description
        .source_file_list
        .source_file
        .iter()
        .map(|sf| sf.id.clone())
        .filter(|s| !s.is_empty())
        .collect();

    let run_source_file_ref_list = (!source_file_ids.is_empty()).then(|| SourceFileRefList {
        count: Some(source_file_ids.len()),
        source_file_refs: source_file_ids
            .iter()
            .cloned()
            .map(|id| SourceFileRef { r#ref: id })
            .collect(),
    });

    let default_source_file_ref = if file_description.source_file_list.source_file.len() == 1 {
        let id = file_description.source_file_list.source_file[0].id.clone();
        (!id.is_empty()).then_some(id)
    } else {
        None
    };

    let default_instrument_configuration_ref = instrument_list.as_ref().and_then(|il| {
        if il.instrument.len() == 1 {
            let id = il.instrument[0].id.clone();
            (!id.is_empty()).then_some(id)
        } else {
            None
        }
    });

    let sample_ref = sample_list.as_ref().and_then(|sl| {
        if sl.samples.len() == 1 {
            let id = sl.samples[0].id.clone();
            (!id.is_empty()).then_some(id)
        } else {
            None
        }
    });

    let run_result = run.unwrap();

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
            id: run_result.id,
            start_time_stamp: run_result.start_time_stamp,
            default_instrument_configuration_ref,
            default_source_file_ref,
            sample_ref,
            source_file_ref_list: run_source_file_ref_list,
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
fn decode_meta_block(
    schema: Option<&Schema>,
    bytes: &[u8],
    item_count: u32,
    meta_count: u32,
    num_count: u32,
    str_count: u32,
) -> Result<Vec<Vec<MetaCvParam>>, String> {
    let _ = schema;
    let item_count = item_count as usize;
    let meta_count = meta_count as usize;
    let num_count = num_count as usize;
    let str_count = str_count as usize;

    if item_count == 0 || meta_count == 0 {
        return Ok(vec![Vec::new(); item_count]);
    }

    let mut offset = 0usize;

    let item_indices = read_u32_vec(
        read_slice(bytes, offset, (item_count + 1) * 4)?,
        item_count + 1,
    )?;
    offset += (item_count + 1) * 4;

    let first = read_slice(bytes, offset, meta_count)?;
    let probe_n = meta_count.min(64);

    let second_off = offset
        .checked_add(meta_count)
        .ok_or_else(|| "Meta size overflow".to_string())?;
    let second_end = second_off
        .checked_add(meta_count)
        .ok_or_else(|| "Meta size overflow".to_string())?;

    let (tag_ids, meta_ref_codes) = if let Some(second) = bytes.get(second_off..second_end) {
        if cv_code_score(&second[..probe_n]) > cv_code_score(&first[..probe_n]) {
            offset += meta_count;
            (Some(first), read_slice(bytes, offset, meta_count)?)
        } else {
            (None, first)
        }
    } else {
        (None, first)
    };
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

            let accession_tail = meta_accessions[m];
            let unit_tail = meta_unit_accessions[m];

            let cv_ref_code = infer_term_cv_ref_code(meta_ref_codes[m], accession_tail);
            let unit_cv_ref_code = infer_unit_cv_ref_code(meta_unit_refs[m], unit_tail);

            let cv_ref = cv_ref_prefix_from_code(cv_ref_code)
                .map(|s| s.to_string())
                .or_else(|| (accession_tail >= 9_900_000).then(|| CV_REF_ATTR.to_string()));
            let accession = format_accession(cv_ref_code, accession_tail)
                .or_else(|| attr_accession(accession_tail));

            let unit_cv_ref = cv_ref_prefix_from_code(unit_cv_ref_code)
                .map(|s| s.to_string())
                .or_else(|| (unit_tail >= 9_900_000).then(|| CV_REF_ATTR.to_string()));
            let unit_accession =
                format_accession(unit_cv_ref_code, unit_tail).or_else(|| attr_accession(unit_tail));

            let mut cv = CvParam {
                cv_ref,
                accession,
                name: String::new(),
                value: Some(value),
                unit_cv_ref,
                unit_accession,
                unit_name: None,
            };

            if let Some(k) = cv.accession.as_deref() {
                cv.name = cv_table_name(k).unwrap_or_default();
            }
            if let Some(uk) = cv.unit_accession.as_deref() {
                cv.unit_name = cv_table_name(uk);
            }

            let tag = tag_ids
                .and_then(|t| t.get(m).copied())
                .map(tag_from_u8)
                .unwrap_or(TagId::Unknown);

            item_params.push(MetaCvParam {
                tag,
                cv_ref_code,
                accession_tail,
                cv,
            });
        }

        result.push(item_params);
    }

    Ok(result)
}

/// <cvList> <fileDescription> <referenceableParamGroupList> <sampleList>
/// <instrumentList> <softwareList> <dataProcessingList> <scanSettingsList>
fn decode_global_meta_structs(
    schema: Option<&Schema>,
    schema_nodes: SchemaNodes<'_>,
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
        Option<Run>,
    ),
    String,
> {
    let header_len = 9 * core::mem::size_of::<u32>();
    if bytes.len() < header_len {
        return Ok((
            None,
            FileDescription::default(),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ));
    }

    let n_fd = read_u32_at(bytes, 0)?;
    let n_run = read_u32_at(bytes, 4)?;
    let n_rpg = read_u32_at(bytes, 8)?;
    let n_samp = read_u32_at(bytes, 12)?;
    let n_inst = read_u32_at(bytes, 16)?;
    let n_soft = read_u32_at(bytes, 20)?;
    let n_dp = read_u32_at(bytes, 24)?;
    let n_acq = read_u32_at(bytes, 28)?;
    let n_cvs = read_u32_at(bytes, 32)?;

    let total = n_fd + n_run + n_rpg + n_samp + n_inst + n_soft + n_dp + n_acq + n_cvs;
    let items = decode_meta_block(schema, &bytes[header_len..], total, m_cnt, n_cnt, s_cnt)?;
    let mut it = items.into_iter();

    // Tag: FileDescription
    let fd = if n_fd > 0 {
        split_file_description_from_cv_params(schema, schema_nodes, it.next().unwrap_or_default())
    } else {
        FileDescription::default()
    };

    // Tag: Run
    let run = if n_run > 0 {
        let mut p = it.next().unwrap_or_default();

        let id = split_id_attr(&mut p).unwrap_or_default();
        let start_time_stamp = split_attr_value(&mut p, ACC_ATTR_START_TIME_STAMP);
        let default_instrument_configuration_ref =
            split_attr_value(&mut p, ACC_ATTR_DEFAULT_INSTRUMENT_CONFIGURATION_REF);
        let default_source_file_ref = split_attr_value(&mut p, ACC_ATTR_DEFAULT_SOURCE_FILE_REF);
        let sample_ref = split_attr_value(&mut p, ACC_ATTR_SAMPLE_REF);

        let source_file_refs = collect_ref_attrs(&mut p);
        let source_file_ref_list = (!source_file_refs.is_empty()).then(|| SourceFileRefList {
            count: Some(source_file_refs.len()),
            source_file_refs: source_file_refs
                .into_iter()
                .map(|r| SourceFileRef { r#ref: r })
                .collect(),
        });

        let cv_params = p.into_iter().map(|m| m.cv).collect::<Vec<_>>();

        for _ in 1..n_run {
            let _ = it.next();
        }

        Some(Run {
            id,
            start_time_stamp,
            default_instrument_configuration_ref,
            default_source_file_ref,
            sample_ref,
            source_file_ref_list,
            cv_params,
            ..Default::default()
        })
    } else {
        None
    };

    // Tag: ReferenceableParamGroupList
    let rpgs = if n_rpg > 0 {
        let mut groups = Vec::with_capacity(n_rpg as usize);
        for _ in 0..n_rpg {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();
            let cv_params = p.into_iter().map(|m| m.cv).collect::<Vec<_>>();
            groups.push(ReferenceableParamGroup {
                id,
                cv_params,
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

    // Tag: SampleList
    let samps = if n_samp > 0 {
        let mut samples = Vec::with_capacity(n_samp as usize);
        for _ in 0..n_samp {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();
            let name = split_name_attr(&mut p).unwrap_or_default();
            let cv_params = p.into_iter().map(|m| m.cv).collect::<Vec<_>>();
            samples.push(Sample {
                id,
                name,
                cv_params,
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

    // Tag: InstrumentList
    let insts = if n_inst > 0 {
        let mut instruments = Vec::with_capacity(n_inst as usize);
        for _ in 0..n_inst {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();

            let mut src_metas = Vec::<MetaCvParam>::new();
            let mut an_metas = Vec::<MetaCvParam>::new();
            let mut det_metas = Vec::<MetaCvParam>::new();
            let mut inst_cv = Vec::<CvParam>::new();

            for m in p {
                let tag = schema_tag_in_nodes(
                    schema,
                    &m,
                    None,
                    &[
                        (TagId::ComponentSource, schema_nodes.component_source),
                        (TagId::ComponentAnalyzer, schema_nodes.component_analyzer),
                        (TagId::ComponentDetector, schema_nodes.component_detector),
                    ],
                );

                if tag == TagId::ComponentSource {
                    src_metas.push(m);
                } else if tag == TagId::ComponentAnalyzer {
                    an_metas.push(m);
                } else if tag == TagId::ComponentDetector {
                    det_metas.push(m);
                } else {
                    inst_cv.push(m.cv);
                }
            }

            let source = split_components_by_order(src_metas, 1)
                .into_iter()
                .map(|(order, cv_param)| Source {
                    order, // adjust type if your field is Option<usize>
                    cv_param,
                    ..Default::default()
                })
                .collect::<Vec<_>>();

            let analyzer = split_components_by_order(an_metas, 2)
                .into_iter()
                .map(|(order, cv_param)| Analyzer {
                    order,
                    cv_param,
                    ..Default::default()
                })
                .collect::<Vec<_>>();

            let detector = split_components_by_order(det_metas, 3)
                .into_iter()
                .map(|(order, cv_param)| Detector {
                    order,
                    cv_param,
                    ..Default::default()
                })
                .collect::<Vec<_>>();

            let component_list = if source.is_empty() && analyzer.is_empty() && detector.is_empty()
            {
                None
            } else {
                Some(ComponentList {
                    source,
                    analyzer,
                    detector,
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

    // Tag: SoftwareList
    let softs = if n_soft > 0 {
        let mut software = Vec::with_capacity(n_soft as usize);

        for _ in 0..n_soft {
            let mut p = it.next().unwrap_or_default();

            let id = split_id_attr(&mut p).unwrap_or_default();
            let version = split_version_attr(&mut p);

            let cv_param = p
                .into_iter()
                .filter(|m| !is_attr_param(&m.cv))
                .map(|m| m.cv)
                .collect::<Vec<_>>();

            software.push(Software {
                id,
                version,
                cv_param,
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

    // Tag: DataProcessingList
    let dps = if n_dp > 0 {
        let mut data_processing = Vec::with_capacity(n_dp as usize);
        for _ in 0..n_dp {
            let mut p = it.next().unwrap_or_default();
            let id = split_id_attr(&mut p).unwrap_or_default();
            let cv_param = p.into_iter().map(|m| m.cv).collect::<Vec<_>>();
            data_processing.push(DataProcessing {
                id,
                processing_method: vec![ProcessingMethod {
                    cv_param,
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

    // Tag: ScanSettingsList
    let acqs = if n_acq > 0 {
        let mut scan_settings = Vec::with_capacity(n_acq as usize);
        for _ in 0..n_acq {
            let p = it.next().unwrap_or_default();
            scan_settings.push(build_scan_settings_from_cv_params(schema, schema_nodes, p));
        }
        Some(ScanSettingsList {
            count: Some(scan_settings.len()),
            scan_settings,
        })
    } else {
        None
    };

    // Tag: CvList
    let cvs = if n_cvs > 0 {
        let mut cv = Vec::<Cv>::with_capacity(n_cvs as usize);
        for _ in 0..n_cvs {
            let p = it.next().unwrap_or_default();
            let mut c = Cv::default();
            for m in p {
                let tail = parse_acc_tail(m.cv.accession.as_deref());
                if tail == 9_900_001 {
                    c.id = attr_string_value(&m.cv).unwrap_or_default();
                } else if tail == 9_900_002 {
                    c.full_name = Some(attr_string_value(&m.cv).unwrap_or_default());
                } else if tail == 9_900_003 {
                    c.version = Some(attr_string_value(&m.cv).unwrap_or_default());
                } else if tail == 9_900_004 {
                    c.uri = Some(attr_string_value(&m.cv).unwrap_or_default());
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

    Ok((cvs, fd, rpgs, samps, insts, softs, dps, acqs, run))
}

#[inline]
fn split_components_by_order(
    mut metas: Vec<MetaCvParam>,
    default_order: u32,
) -> Vec<(Option<u32>, Vec<CvParam>)> {
    let mut out: Vec<(Option<u32>, Vec<CvParam>)> = Vec::new();

    let mut cur_order: Option<u32> = None;
    let mut cur: Vec<CvParam> = Vec::new();

    for m in metas.drain(..) {
        if is_attr_param(&m.cv) {
            let tail = parse_acc_tail(m.cv.accession.as_deref());
            if tail == ACC_ATTR_ORDER {
                if !cur.is_empty() {
                    out.push((cur_order.or(Some(default_order)), std::mem::take(&mut cur)));
                }
                cur_order = attr_string_value(&m.cv).and_then(|s| s.parse::<u32>().ok());
            }
            continue;
        }

        cur.push(m.cv);
    }

    if !cur.is_empty() {
        out.push((cur_order.or(Some(default_order)), cur));
    }

    out
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum BinaryArrayRole {
    Mz,
    Time,
    Intensity,
}

#[inline]
fn split_count_attr(params: &mut Vec<MetaCvParam>) -> Option<usize> {
    split_attr_value(params, ACC_ATTR_COUNT).and_then(|s| s.parse::<usize>().ok())
}

#[inline]
fn split_array_length_attr(params: &mut Vec<MetaCvParam>) -> Option<usize> {
    split_attr_value(params, ACC_ATTR_DEFAULT_ARRAY_LENGTH).and_then(|s| s.parse::<usize>().ok())
}

#[inline]
fn split_encoded_length_attr(params: &mut Vec<MetaCvParam>) -> Option<usize> {
    split_attr_value(params, ACC_ATTR_ENCODED_LENGTH).and_then(|s| s.parse::<usize>().ok())
}

#[inline]
fn ensure_ms_cv(schema: Option<&Schema>, out: &mut Vec<CvParam>, tail: Option<u32>) {
    let Some(t) = tail else { return };
    if out
        .iter()
        .any(|p| parse_acc_tail(p.accession.as_deref()) == t)
    {
        return;
    }
    out.push(ms_cv_param(schema, Some(t)));
}

#[inline]
fn role_tail(tails: MsBinaryTails, role: BinaryArrayRole) -> Option<u32> {
    match role {
        BinaryArrayRole::Mz => tails.mz_array,
        BinaryArrayRole::Time => tails.time_array,
        BinaryArrayRole::Intensity => tails.intensity_array,
    }
}

/// <binaryDataArray>
fn parse_binary_data_array(
    schema: Option<&Schema>,
    tails: MsBinaryTails,
    role: BinaryArrayRole,
    mut params: Vec<MetaCvParam>,
    decoded_f32: Vec<f32>,
    decoded_f64: Vec<f64>,
) -> Result<BinaryDataArray, String> {
    let array_length = split_array_length_attr(&mut params);
    let _encoded_length = split_encoded_length_attr(&mut params);

    let fmt = if params
        .iter()
        .any(|m| m.accession_tail == tails.f32.unwrap_or(u32::MAX))
    {
        1
    } else if params
        .iter()
        .any(|m| m.accession_tail == tails.f64.unwrap_or(u32::MAX))
    {
        2
    } else if !decoded_f32.is_empty() {
        1
    } else {
        2
    };

    println!();

    let elem_size = fmt_elem_size(fmt)?;
    let len = if fmt == 1 {
        decoded_f32.len()
    } else {
        decoded_f64.len()
    };

    if let Some(al) = array_length {
        if al != len {
            return Err("binaryDataArray: arrayLength mismatch".to_string());
        }
    }

    let byte_len = len
        .checked_mul(elem_size)
        .ok_or_else(|| "binaryDataArray: byte length overflow".to_string())?;
    let enc_len = base64_encoded_len(byte_len);

    let mut cv_params = Vec::with_capacity(params.len() + 3);
    for m in params {
        if !is_attr_param(&m.cv) {
            cv_params.push(m.cv);
        }
    }

    ensure_ms_cv(schema, &mut cv_params, tails.float_tail(fmt));
    ensure_ms_cv(schema, &mut cv_params, tails.no_compression);
    ensure_ms_cv(schema, &mut cv_params, role_tail(tails, role));

    let mut ba = BinaryDataArray::default();
    ba.array_length = Some(len);
    ba.encoded_length = Some(enc_len);
    ba.is_f32 = Some(fmt == 1);
    ba.is_f64 = Some(fmt == 2);
    ba.cv_params = cv_params;
    ba.decoded_binary_f32 = decoded_f32;
    ba.decoded_binary_f64 = decoded_f64;
    Ok(ba)
}

/// <binaryDataArrayList>
fn parse_binary_data_array_list(
    schema: Option<&Schema>,
    tails: MsBinaryTails,
    params: &mut Vec<MetaCvParam>,
    arrays: &[(BinaryArrayRole, Vec<f32>, Vec<f64>)],
) -> Result<Option<BinaryDataArrayList>, String> {
    if arrays.is_empty() {
        return Ok(None);
    }

    let mut list_params = Vec::new();
    let mut bda_params = Vec::new();
    let mut rest = Vec::with_capacity(params.len());

    for m in params.drain(..) {
        match m.tag {
            TagId::BinaryDataArrayList => list_params.push(m),
            TagId::BinaryDataArray => bda_params.push(m),
            _ => rest.push(m),
        }
    }
    *params = rest;

    let count = split_count_attr(&mut list_params).unwrap_or(arrays.len());
    if count != arrays.len() {
        return Err("binaryDataArrayList: count mismatch".to_string());
    }

    let is_kind = |m: &MetaCvParam| {
        m.accession_tail == tails.mz_array.unwrap_or(u32::MAX)
            || m.accession_tail == tails.time_array.unwrap_or(u32::MAX)
            || m.accession_tail == tails.intensity_array.unwrap_or(u32::MAX)
    };

    let mut groups: Vec<Vec<MetaCvParam>> = Vec::new();
    let mut cur: Vec<MetaCvParam> = Vec::new();

    for m in bda_params {
        cur.push(m);
        if is_kind(cur.last().unwrap()) && groups.len() + 1 < count {
            groups.push(cur);
            cur = Vec::new();
        }
    }
    if !cur.is_empty() {
        groups.push(cur);
    }

    while groups.len() > count {
        let extra = groups.pop().unwrap();
        groups.last_mut().unwrap().extend(extra);
    }
    while groups.len() < count {
        groups.push(Vec::new());
    }

    let mut out = Vec::with_capacity(count);

    for (i, g) in groups.into_iter().enumerate() {
        let (role, f32v, f64v) = arrays[i].clone();
        let ba = parse_binary_data_array(schema, tails, role, g, f32v, f64v)?;
        out.push(ba);
    }

    Ok(Some(BinaryDataArrayList {
        count: Some(out.len()),
        binary_data_arrays: out,
    }))
}
