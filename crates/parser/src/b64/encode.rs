use serde::Serialize;

// use miniz_oxide::deflate::compress_to_vec_zlib;
use zstd::bulk::compress as zstd_compress;

use std::collections::HashMap;

use crate::{
    b64::utilities::assign_attributes,
    decode::MetadatumValue,
    mzml::{
        attr_meta::*,
        schema::TagId,
        structs::{
            BinaryDataArray, BinaryDataArrayList, Chromatogram, CvParam, MzML, Precursor,
            PrecursorList, Product, ProductList, ReferenceableParamGroup,
            ReferenceableParamGroupRef, ScanList, Spectrum,
        },
    },
};

#[derive(Debug)]
pub struct PackedMeta {
    pub index_offsets: Vec<u32>,          // CI
    pub owner_ids: Vec<u32>,              // MOI
    pub parent_indices: Vec<u32>,         // MPI
    pub tag_ids: Vec<u8>,                 // MTI
    pub ref_codes: Vec<u8>,               // MRI
    pub accession_numbers: Vec<u32>,      // MAN
    pub unit_ref_codes: Vec<u8>,          // MURI
    pub unit_accession_numbers: Vec<u32>, // MUAN
    pub value_kinds: Vec<u8>,             // VK
    pub value_indices: Vec<u32>,          // VI
    pub numeric_values: Vec<f64>,         // VN
    pub string_offsets: Vec<u32>,         // VOFF
    pub string_lengths: Vec<u32>,         // VLEN
    pub string_bytes: Vec<u8>,            // VS
}

#[derive(Debug)]
struct GlobalCounts {
    n_file_description: u32,
    n_run: u32,
    n_ref_param_groups: u32,
    n_samples: u32,
    n_instrument_configs: u32,
    n_software: u32,
    n_data_processing: u32,
    n_acquisition_settings: u32,
    n_cvs: u32,
}

#[derive(Debug)]
struct GlobalMetaItem {
    cvs: Vec<CvParam>,
    tags: Vec<u8>,
    owners: Vec<u32>,
    parents: Vec<u32>,
}

const HEADER_SIZE: usize = 192;
const INDEX_ENTRY_SIZE: usize = 32;
const BLOCK_DIR_ENTRY_SIZE: usize = 32;

const TARGET_BLOCK_UNCOMP_BYTES: usize = 512 * 1024 * 1024;

const ACC_MZ_ARRAY: u32 = 1_000_514;
const ACC_INTENSITY_ARRAY: u32 = 1_000_515;
const ACC_TIME_ARRAY: u32 = 1_000_595;
// const ACC_MS_LEVEL: u32 = 1_000_511;

const ACC_32BIT_FLOAT: u32 = 1_000_521;
const ACC_64BIT_FLOAT: u32 = 1_000_523;

// const HDR_CODEC_ZLIB: u8 = 0;
pub const HDR_CODEC_ZSTD: u8 = 1;
pub const HDR_CODEC_MASK: u8 = 0x0F;

const HDR_FLAG_SPEC_META_COMP: u8 = 1 << 4;
const HDR_FLAG_CHROM_META_COMP: u8 = 1 << 5;
const HDR_FLAG_GLOBAL_META_COMP: u8 = 1 << 6;

const HDR_ARRAY_FILTER_OFF: usize = 178;
const ARRAY_FILTER_NONE: u8 = 0;
const ARRAY_FILTER_BYTE_SHUFFLE: u8 = 1;

#[inline]
fn compress_bytes(input: &[u8], compression_level: u8) -> Vec<u8> {
    if compression_level == 0 {
        return input.to_vec();
    }
    zstd_compress(input, compression_level as i32).unwrap_or_else(|_| input.to_vec())
}

#[inline]
fn header_codec_and_flags(
    codec_id: u8,
    spec_meta_compressed: bool,
    chrom_meta_compressed: bool,
    global_meta_compressed: bool,
) -> u8 {
    let mut v = codec_id & 0x0F;
    if spec_meta_compressed {
        v |= HDR_FLAG_SPEC_META_COMP;
    }
    if chrom_meta_compressed {
        v |= HDR_FLAG_CHROM_META_COMP;
    }
    if global_meta_compressed {
        v |= HDR_FLAG_GLOBAL_META_COMP;
    }
    v
}

#[derive(Copy, Clone)]
enum ArrayRef<'a> {
    F32(&'a [f32]),
    F64(&'a [f64]),
}

impl<'a> ArrayRef<'a> {
    #[inline]
    fn len(self) -> usize {
        match self {
            ArrayRef::F32(s) => s.len(),
            ArrayRef::F64(s) => s.len(),
        }
    }
}

#[inline]
fn elem_size(store_f64: bool) -> usize {
    if store_f64 { 8 } else { 4 }
}

#[inline]
fn write_f32_slice_le(buf: &mut Vec<u8>, xs: &[f32]) {
    if cfg!(target_endian = "little") {
        unsafe {
            let p = xs.as_ptr() as *const u8;
            let b = std::slice::from_raw_parts(p, xs.len() * 4);
            buf.extend_from_slice(b);
        }
    } else {
        for &v in xs {
            write_f32_le(buf, v);
        }
    }
}

#[inline]
fn write_f64_slice_le(buf: &mut Vec<u8>, xs: &[f64]) {
    if cfg!(target_endian = "little") {
        unsafe {
            let p = xs.as_ptr() as *const u8;
            let b = std::slice::from_raw_parts(p, xs.len() * 8);
            buf.extend_from_slice(b);
        }
    } else {
        for &v in xs {
            write_f64_le(buf, v);
        }
    }
}

#[inline]
fn write_array_as_f64(buf: &mut Vec<u8>, arr: ArrayRef<'_>) {
    match arr {
        ArrayRef::F32(xs) => {
            for &v in xs {
                write_f64_le(buf, v as f64);
            }
        }
        ArrayRef::F64(xs) => write_f64_slice_le(buf, xs),
    }
}

#[inline]
fn write_array_as_f32(buf: &mut Vec<u8>, arr: ArrayRef<'_>) {
    match arr {
        ArrayRef::F32(xs) => write_f32_slice_le(buf, xs),
        ArrayRef::F64(xs) => {
            for &v in xs {
                write_f64_as_f32(buf, v);
            }
        }
    }
}

#[inline]
fn write_array(buf: &mut Vec<u8>, arr: ArrayRef<'_>, store_f64: bool) {
    if store_f64 {
        write_array_as_f64(buf, arr);
    } else {
        write_array_as_f32(buf, arr);
    }
}

#[inline]
fn byte_shuffle_into(input: &[u8], output: &mut [u8], elem_size: usize) {
    let count = input.len() / elem_size;
    for b in 0..elem_size {
        let out_base = b * count;
        let mut in_i = b;
        for e in 0..count {
            output[out_base + e] = input[in_i];
            in_i += elem_size;
        }
    }
}

#[inline]
fn merge_declared_width(cur: &mut Option<bool>, next_is_f64: bool, axis: &'static str) {
    match *cur {
        None => *cur = Some(next_is_f64),
        Some(prev) if prev == next_is_f64 => {}
        Some(prev) => {
            panic!(
                "Mixed float widths for {axis}: saw {} then {}. \
                 Your container format requires a single width per axis. \
                 Either set f32_compress=true (normalize to f32) or normalize the input mzML.",
                if prev { "f64" } else { "f32" },
                if next_is_f64 { "f64" } else { "f32" },
            );
        }
    }
}

#[derive(Clone, Copy)]
struct BlockDirEntry {
    comp_off: u64,
    comp_size: u64,
    uncomp_bytes: u64,
}

struct ContainerBuilder {
    target_uncomp_bytes: usize,
    compression_level: u8,
    elem_size: usize,
    do_shuffle: bool,
    current: Vec<u8>,
    entries: Vec<BlockDirEntry>,
    compressed: Vec<u8>,
    scratch: Vec<u8>,
}

impl ContainerBuilder {
    #[inline]
    fn new(
        target_uncomp_bytes: usize,
        compression_level: u8,
        elem_size: usize,
        do_shuffle: bool,
    ) -> Self {
        Self {
            target_uncomp_bytes,
            compression_level,
            elem_size,
            do_shuffle,
            current: Vec::new(),
            entries: Vec::new(),
            compressed: Vec::new(),
            scratch: Vec::new(),
        }
    }

    #[inline]
    fn current_block_id(&self) -> u32 {
        self.entries.len() as u32
    }

    #[inline]
    fn flush_current(&mut self) {
        if self.current.is_empty() {
            return;
        }

        let uncomp_bytes = self.current.len() as u64;
        let comp_off = self.compressed.len() as u64;

        if self.compression_level == 0 {
            self.entries.push(BlockDirEntry {
                comp_off,
                comp_size: uncomp_bytes,
                uncomp_bytes,
            });
            self.compressed.extend_from_slice(&self.current);
            self.current.clear();
            return;
        }

        let to_compress: &[u8] = if self.do_shuffle && self.elem_size > 1 {
            self.scratch.resize(self.current.len(), 0);
            byte_shuffle_into(
                self.current.as_slice(),
                self.scratch.as_mut_slice(),
                self.elem_size,
            );
            self.scratch.as_slice()
        } else {
            self.current.as_slice()
        };

        let comp = compress_bytes(to_compress, self.compression_level);
        let comp_size = comp.len() as u64;

        self.entries.push(BlockDirEntry {
            comp_off,
            comp_size,
            uncomp_bytes,
        });

        self.compressed.extend_from_slice(&comp);
        self.current.clear();
    }

    #[inline]
    fn ensure_room_for_item(&mut self, item_bytes: usize) {
        if !self.current.is_empty() && self.current.len() + item_bytes > self.target_uncomp_bytes {
            self.flush_current();
        }
    }

    #[inline]
    fn write_item<F>(&mut self, item_bytes: usize, write_fn: F) -> u32
    where
        F: FnOnce(&mut Vec<u8>),
    {
        if item_bytes > self.target_uncomp_bytes {
            if !self.current.is_empty() {
                self.flush_current();
            }
            let block_id = self.current_block_id();
            self.current.reserve(item_bytes);
            write_fn(&mut self.current);
            self.flush_current();
            return block_id;
        }

        self.ensure_room_for_item(item_bytes);
        let block_id = self.current_block_id();
        self.current.reserve(item_bytes);
        write_fn(&mut self.current);
        block_id
    }

    #[inline]
    fn finalize(mut self) -> (Vec<u8>, u32) {
        self.flush_current();

        let block_count = self.entries.len() as u32;
        let dir_bytes = self.entries.len() * BLOCK_DIR_ENTRY_SIZE;

        let mut container = Vec::with_capacity(dir_bytes + self.compressed.len());
        for e in &self.entries {
            write_u64_le(&mut container, e.comp_off);
            write_u64_le(&mut container, e.comp_size);
            write_u64_le(&mut container, e.uncomp_bytes);
            container.extend_from_slice(&[0u8; 8]);
        }
        container.extend_from_slice(&self.compressed);

        (container, block_count)
    }
}

#[derive(Debug, Default)]
struct NodeIdGen {
    next: u32,
}

impl NodeIdGen {
    #[inline]
    fn new() -> Self {
        Self { next: 1 }
    }

    #[inline]
    fn alloc(&mut self) -> u32 {
        let id = self.next;
        self.next += 1;
        id
    }
}

#[inline]
fn push_tagged_raw(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag_id: u8,
    owner_id: u32,
    parent_owner_id: u32,
    cv: CvParam,
) {
    out.push(cv);
    tags.push(tag_id);
    owners.push(owner_id);
    parents.push(parent_owner_id);
}

#[inline]
fn extend_tagged_raw(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag_id: u8,
    owner_id: u32,
    parent_owner_id: u32,
    cvs: &[CvParam],
) {
    out.extend_from_slice(cvs);
    tags.resize(tags.len() + cvs.len(), tag_id);
    owners.resize(owners.len() + cvs.len(), owner_id);
    parents.resize(parents.len() + cvs.len(), parent_owner_id);
}

#[inline]
fn push_tagged_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    cv: CvParam,
) {
    push_tagged_raw(
        out,
        tags,
        owners,
        parents,
        tag as u8,
        owner_id,
        parent_owner_id,
        cv,
    );
}

#[inline]
fn extend_tagged_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    cvs: &[CvParam],
) {
    extend_tagged_raw(
        out,
        tags,
        owners,
        parents,
        tag as u8,
        owner_id,
        parent_owner_id,
        cvs,
    );
}

#[inline]
fn push_attr_string_tagged_raw(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag_id: u8,
    owner_id: u32,
    parent_owner_id: u32,
    accession_tail: u32,
    value: &str,
) {
    if !value.is_empty() {
        push_tagged_raw(
            out,
            tags,
            owners,
            parents,
            tag_id,
            owner_id,
            parent_owner_id,
            attr_cv_param(accession_tail, value),
        );
    }
}

#[inline]
fn push_attr_u32_tagged_raw(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag_id: u8,
    owner_id: u32,
    parent_owner_id: u32,
    accession_tail: u32,
    value: Option<u32>,
) {
    if let Some(v) = value {
        push_tagged_raw(
            out,
            tags,
            owners,
            parents,
            tag_id,
            owner_id,
            parent_owner_id,
            attr_cv_param(accession_tail, &v.to_string()),
        );
    }
}

#[inline]
fn push_attr_string_tagged_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    accession_tail: u32,
    value: &str,
) {
    push_attr_string_tagged_raw(
        out,
        tags,
        owners,
        parents,
        tag as u8,
        owner_id,
        parent_owner_id,
        accession_tail,
        value,
    );
}

#[inline]
fn push_attr_usize_tagged_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    accession_tail: u32,
    value: Option<u32>,
) {
    push_attr_u32_tagged_raw(
        out,
        tags,
        owners,
        parents,
        tag as u8,
        owner_id,
        parent_owner_id,
        accession_tail,
        value,
    );
}

/// <referenceableParamGroupList>
fn build_ref_group_map<'a>(mzml: &'a MzML) -> HashMap<&'a str, &'a ReferenceableParamGroup> {
    let mut map = HashMap::new();
    if let Some(list) = &mzml.referenceable_param_group_list {
        for g in &list.referenceable_param_groups {
            map.insert(g.id.as_str(), g);
        }
    }
    map
}

/// <referenceableParamGroupRef>
#[inline]
fn extend_ref_group_cv_params_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    refs: &[ReferenceableParamGroupRef],
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for r in refs {
        if let Some(g) = ref_groups.get(r.r#ref.as_str()) {
            extend_tagged_ids(
                out,
                tags,
                owners,
                parents,
                tag,
                owner_id,
                parent_owner_id,
                &g.cv_params,
            );
        }
    }
}

/// <binaryDataArray>
fn extend_binary_data_array_cv_params_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    owner_id: u32,
    parent_owner_id: u32,
    ba: &BinaryDataArray,
    x_accession_tail: u32,
    y_accession_tail: u32,
    f32_compress: bool,
) {
    let mut is_x = false;
    let mut is_y = false;

    for cv in &ba.cv_params {
        let tail = parse_accession_tail(cv.accession.as_deref());
        if tail == x_accession_tail {
            is_x = true;
        } else if tail == y_accession_tail {
            is_y = true;
        }
    }

    let desired_float_tail = if f32_compress && (is_x || is_y) {
        Some(ACC_32BIT_FLOAT)
    } else {
        None
    };

    let mut wrote_float = false;

    for cv in &ba.cv_params {
        let tail = parse_accession_tail(cv.accession.as_deref());

        if tail == ACC_32BIT_FLOAT || tail == ACC_64BIT_FLOAT {
            if let Some(desired) = desired_float_tail {
                if !wrote_float {
                    push_tagged_ids(
                        out,
                        tags,
                        owners,
                        parents,
                        TagId::BinaryDataArray,
                        owner_id,
                        parent_owner_id,
                        ms_float_param(desired),
                    );
                    wrote_float = true;
                }
                continue;
            } else {
                push_tagged_ids(
                    out,
                    tags,
                    owners,
                    parents,
                    TagId::BinaryDataArray,
                    owner_id,
                    parent_owner_id,
                    cv.clone(),
                );
                wrote_float = true;
            }
        } else {
            push_tagged_ids(
                out,
                tags,
                owners,
                parents,
                TagId::BinaryDataArray,
                owner_id,
                parent_owner_id,
                cv.clone(),
            );
        }
    }

    if let Some(desired) = desired_float_tail {
        if !wrote_float {
            push_tagged_ids(
                out,
                tags,
                owners,
                parents,
                TagId::BinaryDataArray,
                owner_id,
                parent_owner_id,
                ms_float_param(desired),
            );
        }
    }
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

/// <scanList>
fn flatten_scan_list_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    scan_list: &ScanList,
    scan_list_owner_id: u32,
    id_gen: &mut NodeIdGen,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for scan in &scan_list.scans {
        let scan_id = id_gen.alloc();

        // <scan> attributes
        push_schema_attributes(
            out,
            tags,
            owners,
            parents,
            TagId::Scan,
            scan_id,
            scan_list_owner_id,
            scan,
        );

        // <scan> refGroup
        extend_ref_group_cv_params_ids(
            out,
            tags,
            owners,
            parents,
            TagId::Scan,
            scan_id,
            scan_list_owner_id,
            &scan.referenceable_param_group_refs,
            ref_groups,
        );

        // <scan> cvParams
        extend_tagged_ids(
            out,
            tags,
            owners,
            parents,
            TagId::Scan,
            scan_id,
            scan_list_owner_id,
            &scan.cv_params,
        );

        // <scanWindowList>/<scanWindow>
        if let Some(wl) = &scan.scan_window_list {
            for w in &wl.scan_windows {
                let win_id = id_gen.alloc();

                // <scanWindow> attributes
                push_schema_attributes(
                    out,
                    tags,
                    owners,
                    parents,
                    TagId::ScanWindow,
                    win_id,
                    scan_id,
                    w,
                );

                // <scanWindow> cvParams
                extend_tagged_ids(
                    out,
                    tags,
                    owners,
                    parents,
                    TagId::ScanWindow,
                    win_id,
                    scan_id,
                    &w.cv_params,
                );
            }
        }
    }
}

/// <precursorList>
fn flatten_precursor_list_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    precursor_list: &PrecursorList,
    parent_owner_id: u32,
    id_gen: &mut NodeIdGen,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for p in &precursor_list.precursors {
        flatten_precursor_ids(
            out,
            tags,
            owners,
            parents,
            p,
            parent_owner_id,
            id_gen,
            ref_groups,
        );
    }
}

/// <precursor>
fn flatten_precursor_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    precursor: &Precursor,
    parent_owner_id: u32,
    id_gen: &mut NodeIdGen,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    let precursor_id = id_gen.alloc();

    let base = out.len();

    push_schema_attributes(
        out,
        tags,
        owners,
        parents,
        TagId::Precursor,
        precursor_id,
        parent_owner_id,
        precursor,
    );

    if let Some(r) = precursor.spectrum_ref.as_deref() {
        let mut already = false;

        for i in base..out.len() {
            if owners[i] != precursor_id {
                continue;
            }
            if tags[i] != TagId::Precursor as u8 {
                continue;
            }

            let tail = parse_accession_tail(out[i].accession.as_deref());
            if tail == ACC_ATTR_SPECTRUM_REF {
                already = true;
                break;
            }
        }

        if !already {
            push_attr_string_tagged_ids(
                out,
                tags,
                owners,
                parents,
                TagId::Precursor,
                precursor_id,
                parent_owner_id,
                ACC_ATTR_SPECTRUM_REF,
                r,
            );
        }
    }

    if let Some(iw) = &precursor.isolation_window {
        let iw_id = id_gen.alloc();

        push_schema_attributes(
            out,
            tags,
            owners,
            parents,
            TagId::IsolationWindow,
            iw_id,
            precursor_id,
            iw,
        );

        extend_ref_group_cv_params_ids(
            out,
            tags,
            owners,
            parents,
            TagId::IsolationWindow,
            iw_id,
            precursor_id,
            &iw.referenceable_param_group_refs,
            ref_groups,
        );

        extend_tagged_ids(
            out,
            tags,
            owners,
            parents,
            TagId::IsolationWindow,
            iw_id,
            precursor_id,
            &iw.cv_params,
        );
    }

    if let Some(sil) = &precursor.selected_ion_list {
        let sil_id = id_gen.alloc();
        push_schema_attributes(
            out,
            tags,
            owners,
            parents,
            TagId::SelectedIonList,
            sil_id,
            precursor_id,
            sil,
        );
        for ion in &sil.selected_ions {
            let ion_id = id_gen.alloc();

            push_schema_attributes(
                out,
                tags,
                owners,
                parents,
                TagId::SelectedIon,
                ion_id,
                sil_id,
                ion,
            );

            extend_ref_group_cv_params_ids(
                out,
                tags,
                owners,
                parents,
                TagId::SelectedIon,
                ion_id,
                precursor_id,
                &ion.referenceable_param_group_refs,
                ref_groups,
            );

            extend_tagged_ids(
                out,
                tags,
                owners,
                parents,
                TagId::SelectedIon,
                ion_id,
                precursor_id,
                &ion.cv_params,
            );
        }
    }

    if let Some(act) = &precursor.activation {
        let act_id = id_gen.alloc();

        push_schema_attributes(
            out,
            tags,
            owners,
            parents,
            TagId::Activation,
            act_id,
            precursor_id,
            act,
        );

        extend_ref_group_cv_params_ids(
            out,
            tags,
            owners,
            parents,
            TagId::Activation,
            act_id,
            precursor_id,
            &act.referenceable_param_group_refs,
            ref_groups,
        );

        extend_tagged_ids(
            out,
            tags,
            owners,
            parents,
            TagId::Activation,
            act_id,
            precursor_id,
            &act.cv_params,
        );
    }
}

/// <productList>
fn flatten_product_list_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    product_list: &ProductList,
    parent_owner_id: u32,
    id_gen: &mut NodeIdGen,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for p in &product_list.products {
        flatten_product_ids(
            out,
            tags,
            owners,
            parents,
            p,
            parent_owner_id,
            id_gen,
            ref_groups,
        );
    }
}

/// <product>
fn flatten_product_ids(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    product: &Product,
    parent_owner_id: u32,
    id_gen: &mut NodeIdGen,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    let product_id = id_gen.alloc();
    push_schema_attributes(
        out,
        tags,
        owners,
        parents,
        TagId::Product,
        product_id,
        parent_owner_id,
        product,
    );

    if let Some(iw) = &product.isolation_window {
        let iw_id = id_gen.alloc();

        push_schema_attributes(
            out,
            tags,
            owners,
            parents,
            TagId::IsolationWindow,
            iw_id,
            product_id,
            iw,
        );

        // Tag: IsolationWindow
        extend_ref_group_cv_params_ids(
            out,
            tags,
            owners,
            parents,
            TagId::IsolationWindow,
            iw_id,
            product_id,
            &iw.referenceable_param_group_refs,
            ref_groups,
        );

        // Tag: IsolationWindow
        extend_tagged_ids(
            out,
            tags,
            owners,
            parents,
            TagId::IsolationWindow,
            iw_id,
            product_id,
            &iw.cv_params,
        );
    }
}

/// <spectrum>
fn flatten_spectrum_metadata_into_owned(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    spectrum: &Spectrum,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
    spectrum_id: u32,
    parent_owner_id: u32,
    id_gen: &mut NodeIdGen,
    x_accession_tail: u32,
    y_accession_tail: u32,
    x_store_f64: bool,
    y_store_f64: bool,
    f32_compress: bool,
) {
    let attrs = assign_attributes(spectrum, TagId::Spectrum, spectrum_id, parent_owner_id);
    push_assigned_attributes_as_cv_params(
        out,
        tags,
        owners,
        parents,
        TagId::Spectrum,
        spectrum_id,
        parent_owner_id,
        attrs,
    );

    extend_ref_group_cv_params_ids(
        out,
        tags,
        owners,
        parents,
        TagId::Spectrum,
        spectrum_id,
        parent_owner_id,
        &spectrum.referenceable_param_group_refs,
        ref_groups,
    );

    extend_tagged_ids(
        out,
        tags,
        owners,
        parents,
        TagId::Spectrum,
        spectrum_id,
        parent_owner_id,
        &spectrum.cv_params,
    );

    if let Some(sd) = &spectrum.spectrum_description {
        let sd_id = id_gen.alloc();
        extend_ref_group_cv_params_ids(
            out,
            tags,
            owners,
            parents,
            TagId::SpectrumDescription,
            sd_id,
            spectrum_id,
            &sd.referenceable_param_group_refs,
            ref_groups,
        );

        extend_tagged_ids(
            out,
            tags,
            owners,
            parents,
            TagId::SpectrumDescription,
            sd_id,
            spectrum_id,
            &sd.cv_params,
        );

        if let Some(sl) = &sd.scan_list {
            let sl_id = id_gen.alloc();

            let attrs = assign_attributes(sl, TagId::ScanList, sl_id, sd_id);
            push_assigned_attributes_as_cv_params(
                out,
                tags,
                owners,
                parents,
                TagId::ScanList,
                sl_id,
                sd_id,
                attrs,
            );

            extend_tagged_ids(
                out,
                tags,
                owners,
                parents,
                TagId::ScanList,
                sl_id,
                sd_id,
                &sl.cv_params,
            );

            flatten_scan_list_ids(out, tags, owners, parents, sl, sl_id, id_gen, ref_groups);
        }

        if let Some(pl) = &sd.precursor_list {
            let pl_id = id_gen.alloc();

            let attrs = assign_attributes(pl, TagId::PrecursorList, pl_id, sd_id);
            push_assigned_attributes_as_cv_params(
                out,
                tags,
                owners,
                parents,
                TagId::PrecursorList,
                pl_id,
                sd_id,
                attrs,
            );
            flatten_precursor_list_ids(out, tags, owners, parents, pl, pl_id, id_gen, ref_groups);
        }

        if let Some(pl) = &sd.product_list {
            let pl_id = id_gen.alloc();

            let attrs = assign_attributes(pl, TagId::ProductList, pl_id, sd_id);
            push_assigned_attributes_as_cv_params(
                out,
                tags,
                owners,
                parents,
                TagId::ProductList,
                pl_id,
                sd_id,
                attrs,
            );

            flatten_product_list_ids(out, tags, owners, parents, pl, pl_id, id_gen, ref_groups);
        }
    }

    if let Some(sl) = &spectrum.scan_list {
        let sl_id = id_gen.alloc();

        let attrs = assign_attributes(sl, TagId::ScanList, sl_id, spectrum_id);
        push_assigned_attributes_as_cv_params(
            out,
            tags,
            owners,
            parents,
            TagId::ScanList,
            sl_id,
            spectrum_id,
            attrs,
        );

        extend_tagged_ids(
            out,
            tags,
            owners,
            parents,
            TagId::ScanList,
            sl_id,
            spectrum_id,
            &sl.cv_params,
        );

        flatten_scan_list_ids(out, tags, owners, parents, sl, sl_id, id_gen, ref_groups);
    }

    if let Some(pl) = &spectrum.precursor_list {
        let pl_id = id_gen.alloc();

        let attrs = assign_attributes(pl, TagId::PrecursorList, pl_id, spectrum_id);
        push_assigned_attributes_as_cv_params(
            out,
            tags,
            owners,
            parents,
            TagId::PrecursorList,
            pl_id,
            spectrum_id,
            attrs,
        );

        flatten_precursor_list_ids(out, tags, owners, parents, pl, pl_id, id_gen, ref_groups);
    }

    if let Some(pl) = &spectrum.product_list {
        let pl_id = id_gen.alloc();

        let attrs = assign_attributes(pl, TagId::ProductList, pl_id, spectrum_id);
        push_assigned_attributes_as_cv_params(
            out,
            tags,
            owners,
            parents,
            TagId::ProductList,
            pl_id,
            spectrum_id,
            attrs,
        );

        flatten_product_list_ids(out, tags, owners, parents, pl, pl_id, id_gen, ref_groups);
    }

    if let Some(bal) = &spectrum.binary_data_array_list {
        let bal_id = id_gen.alloc();

        let attrs = assign_attributes(bal, TagId::BinaryDataArrayList, bal_id, spectrum_id);
        push_assigned_attributes_as_cv_params(
            out,
            tags,
            owners,
            parents,
            TagId::BinaryDataArrayList,
            bal_id,
            spectrum_id,
            attrs,
        );
        for ba in &bal.binary_data_arrays {
            let ba_id = id_gen.alloc();

            push_schema_attributes(
                out,
                tags,
                owners,
                parents,
                TagId::BinaryDataArray,
                ba_id,
                bal_id,
                ba,
            );
            extend_ref_group_cv_params_ids(
                out,
                tags,
                owners,
                parents,
                TagId::BinaryDataArray,
                ba_id,
                bal_id,
                &ba.referenceable_param_group_refs,
                ref_groups,
            );

            extend_binary_data_array_cv_params_ids(
                out,
                tags,
                owners,
                parents,
                ba_id,
                bal_id,
                ba,
                x_accession_tail,
                y_accession_tail,
                f32_compress,
            );

            let _ = (x_store_f64, y_store_f64);
        }
    }
}

/// <chromatogram>
fn flatten_chromatogram_metadata_into(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    chrom: &Chromatogram,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
    chrom_id: u32,
    parent_owner_id: u32,
    id_gen: &mut NodeIdGen,
    x_accession_tail: u32,
    y_accession_tail: u32,
    x_store_f64: bool,
    y_store_f64: bool,
    f32_compress: bool,
) {
    let attrs = assign_attributes(chrom, TagId::Chromatogram, chrom_id, parent_owner_id);
    push_assigned_attributes_as_cv_params(
        out,
        tags,
        owners,
        parents,
        TagId::Chromatogram,
        chrom_id,
        parent_owner_id,
        attrs,
    );

    extend_ref_group_cv_params_ids(
        out,
        tags,
        owners,
        parents,
        TagId::Chromatogram,
        chrom_id,
        parent_owner_id,
        &chrom.referenceable_param_group_refs,
        ref_groups,
    );

    extend_tagged_ids(
        out,
        tags,
        owners,
        parents,
        TagId::Chromatogram,
        chrom_id,
        parent_owner_id,
        &chrom.cv_params,
    );

    if let Some(p) = &chrom.precursor {
        flatten_precursor_ids(out, tags, owners, parents, p, chrom_id, id_gen, ref_groups);
    }

    if let Some(p) = &chrom.product {
        flatten_product_ids(out, tags, owners, parents, p, chrom_id, id_gen, ref_groups);
    }

    if let Some(bal) = &chrom.binary_data_array_list {
        let bal_id = id_gen.alloc();

        let attrs = assign_attributes(bal, TagId::BinaryDataArrayList, bal_id, chrom_id);
        push_assigned_attributes_as_cv_params(
            out,
            tags,
            owners,
            parents,
            TagId::BinaryDataArrayList,
            bal_id,
            chrom_id,
            attrs,
        );

        for ba in &bal.binary_data_arrays {
            let ba_id = id_gen.alloc();

            push_schema_attributes(
                out,
                tags,
                owners,
                parents,
                TagId::BinaryDataArray,
                ba_id,
                bal_id,
                ba,
            );
            extend_ref_group_cv_params_ids(
                out,
                tags,
                owners,
                parents,
                TagId::BinaryDataArray,
                ba_id,
                bal_id,
                &ba.referenceable_param_group_refs,
                ref_groups,
            );

            extend_binary_data_array_cv_params_ids(
                out,
                tags,
                owners,
                parents,
                ba_id,
                bal_id,
                ba,
                x_accession_tail,
                y_accession_tail,
                f32_compress,
            );

            let _ = (x_store_f64, y_store_f64);
        }
    }
}

/// <cvList>
fn build_global_meta_items(
    mzml: &MzML,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
    id_gen: &mut NodeIdGen,
) -> (Vec<GlobalMetaItem>, GlobalCounts) {
    fn cv_field(acc: u32, val: Option<&str>) -> CvParam {
        CvParam {
            cv_ref: Some(CV_REF_ATTR.to_string()),
            accession: Some(format!("{}:{:07}", CV_REF_ATTR, acc)),
            name: String::new(),
            value: val.map(|s| s.to_string()),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        }
    }

    let mut items: Vec<GlobalMetaItem> = Vec::new();

    {
        let fd = &mzml.file_description;

        let mut out = Vec::new();
        let mut tags = Vec::new();
        let mut owners = Vec::new();
        let mut parents = Vec::new();

        // <fileDescription>
        let file_desc_id = id_gen.alloc();

        // (optional but recommended) encode the <fileDescription> node itself as a tagged node.
        // If it has no attributes, you can just write a presence marker like count=0 (or skip).
        // Here we skip attributes because mzML fileDescription usually has none in 0.99.

        // <fileContent>
        let file_content_id = id_gen.alloc();

        extend_ref_group_cv_params_ids(
            &mut out,
            &mut tags,
            &mut owners,
            &mut parents,
            TagId::FileContent,
            file_content_id,
            file_desc_id,
            &fd.file_content.referenceable_param_group_refs,
            ref_groups,
        );

        extend_tagged_ids(
            &mut out,
            &mut tags,
            &mut owners,
            &mut parents,
            TagId::FileContent,
            file_content_id,
            file_desc_id,
            &fd.file_content.cv_params,
        );

        // <sourceFileList count="...">
        let sfl_id = id_gen.alloc();

        push_attr_usize_tagged_ids(
            &mut out,
            &mut tags,
            &mut owners,
            &mut parents,
            TagId::SourceFileList,
            sfl_id,
            file_desc_id,
            ACC_ATTR_COUNT,
            Some(fd.source_file_list.source_file.len() as u32),
        );

        // <sourceFile .../> children of sourceFileList
        for sf in &fd.source_file_list.source_file {
            let sf_id = id_gen.alloc();

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::SourceFile,
                sf_id,
                sfl_id,
                ACC_ATTR_ID,
                sf.id.as_str(),
            );

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::SourceFile,
                sf_id,
                sfl_id,
                ACC_ATTR_NAME,
                sf.name.as_str(),
            );

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::SourceFile,
                sf_id,
                sfl_id,
                ACC_ATTR_LOCATION,
                sf.location.as_str(),
            );

            extend_ref_group_cv_params_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::SourceFile,
                sf_id,
                sfl_id,
                &sf.referenceable_param_group_ref,
                ref_groups,
            );

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::SourceFile,
                sf_id,
                sfl_id,
                &sf.cv_param,
            );
        }

        for c in &fd.contacts {
            let contact_id = id_gen.alloc();

            extend_ref_group_cv_params_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Contact,
                contact_id,
                file_desc_id,
                &c.referenceable_param_group_refs,
                ref_groups,
            );

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Contact,
                contact_id,
                file_desc_id,
                &c.cv_params,
            );
        }

        items.push(GlobalMetaItem {
            cvs: out,
            tags,
            owners,
            parents,
        });
    }

    let n_file_description = 1u32;

    let n_run = 1u32;
    {
        let run = &mzml.run;

        let mut out: Vec<CvParam> = Vec::new();
        let mut tags: Vec<u8> = Vec::new();
        let mut owners: Vec<u32> = Vec::new();
        let mut parents: Vec<u32> = Vec::new();

        let run_id = id_gen.alloc();

        if !run.id.is_empty() {
            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Run,
                run_id,
                0,
                ACC_ATTR_ID,
                run.id.as_str(),
            );
        }

        if let Some(v) = run.start_time_stamp.as_deref() {
            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Run,
                run_id,
                0,
                ACC_ATTR_START_TIME_STAMP,
                v,
            );
        }

        if let Some(v) = run.default_instrument_configuration_ref.as_deref() {
            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Run,
                run_id,
                0,
                ACC_ATTR_DEFAULT_INSTRUMENT_CONFIGURATION_REF,
                v,
            );
        }

        if let Some(v) = run.default_source_file_ref.as_deref() {
            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Run,
                run_id,
                0,
                ACC_ATTR_DEFAULT_SOURCE_FILE_REF,
                v,
            );
        }

        if let Some(v) = run.sample_ref.as_deref() {
            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Run,
                run_id,
                0,
                ACC_ATTR_SAMPLE_REF,
                v,
            );
        }

        if let Some(sfrl) = &run.source_file_ref_list {
            let sfrl_id = id_gen.alloc();

            push_attr_usize_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::SourceFileRefList,
                sfrl_id,
                run_id,
                ACC_ATTR_COUNT,
                Some(sfrl.source_file_refs.len() as u32),
            );

            for sref in &sfrl.source_file_refs {
                let sref_id = id_gen.alloc();
                push_attr_string_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::SourceFileRef,
                    sref_id,
                    sfrl_id,
                    ACC_ATTR_REF,
                    sref.r#ref.as_str(),
                );
            }
        }

        for r in &run.referenceable_param_group_refs {
            let rgr_id = id_gen.alloc();
            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::ReferenceableParamGroupRef,
                rgr_id,
                run_id,
                ACC_ATTR_REF,
                r.r#ref.as_str(),
            );
        }

        if !run.cv_params.is_empty() {
            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Run,
                run_id,
                0,
                &run.cv_params,
            );
        }

        items.push(GlobalMetaItem {
            cvs: out,
            tags,
            owners,
            parents,
        });
    }

    let ref_start = items.len();
    if let Some(rpgl) = &mzml.referenceable_param_group_list {
        for g in &rpgl.referenceable_param_groups {
            let mut out = Vec::new();
            let mut tags = Vec::new();
            let mut owners = Vec::new();
            let mut parents = Vec::new();

            let g_id = id_gen.alloc();

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::ReferenceableParamGroup,
                g_id,
                0,
                ACC_ATTR_ID,
                g.id.as_str(),
            );

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::ReferenceableParamGroup,
                g_id,
                0,
                &g.cv_params,
            );

            items.push(GlobalMetaItem {
                cvs: out,
                tags,
                owners,
                parents,
            });
        }
    }
    let n_ref_param_groups = (items.len() - ref_start) as u32;

    let samples_start = items.len();
    if let Some(sl) = &mzml.sample_list {
        for s in &sl.samples {
            let mut out = Vec::new();
            let mut tags = Vec::new();
            let mut owners = Vec::new();
            let mut parents = Vec::new();

            let sample_id = id_gen.alloc();

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Sample,
                sample_id,
                0,
                ACC_ATTR_ID,
                s.id.as_str(),
            );

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Sample,
                sample_id,
                0,
                ACC_ATTR_NAME,
                s.name.as_str(),
            );

            if let Some(r) = &s.referenceable_param_group_ref {
                extend_ref_group_cv_params_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::Sample,
                    sample_id,
                    0,
                    std::slice::from_ref(r),
                    ref_groups,
                );
            }

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Sample,
                sample_id,
                0,
                &[],
            );

            items.push(GlobalMetaItem {
                cvs: out,
                tags,
                owners,
                parents,
            });
        }
    }
    let n_samples = (items.len() - samples_start) as u32;

    let instr_start = items.len();
    if let Some(il) = &mzml.instrument_list {
        for ic in &il.instrument {
            let mut out = Vec::new();
            let mut tags = Vec::new();
            let mut owners = Vec::new();
            let mut parents = Vec::new();

            let inst_id = id_gen.alloc();

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Instrument,
                inst_id,
                0,
                ACC_ATTR_ID,
                ic.id.as_str(),
            );

            extend_ref_group_cv_params_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Instrument,
                inst_id,
                0,
                &ic.referenceable_param_group_ref,
                ref_groups,
            );

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Instrument,
                inst_id,
                0,
                &ic.cv_param,
            );

            if let Some(cl) = &ic.component_list {
                for s in &cl.source {
                    let cid = id_gen.alloc();
                    push_attr_usize_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentSource,
                        cid,
                        inst_id,
                        ACC_ATTR_ORDER,
                        s.order,
                    );

                    extend_ref_group_cv_params_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentSource,
                        cid,
                        inst_id,
                        &s.referenceable_param_group_ref,
                        ref_groups,
                    );

                    extend_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentSource,
                        cid,
                        inst_id,
                        &s.cv_param,
                    );
                }

                for a in &cl.analyzer {
                    let cid = id_gen.alloc();

                    push_attr_usize_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentAnalyzer,
                        cid,
                        inst_id,
                        ACC_ATTR_ORDER,
                        a.order,
                    );

                    extend_ref_group_cv_params_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentAnalyzer,
                        cid,
                        inst_id,
                        &a.referenceable_param_group_ref,
                        ref_groups,
                    );

                    extend_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentAnalyzer,
                        cid,
                        inst_id,
                        &a.cv_param,
                    );
                }

                for d in &cl.detector {
                    let cid = id_gen.alloc();

                    push_attr_usize_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentDetector,
                        cid,
                        inst_id,
                        ACC_ATTR_ORDER,
                        d.order,
                    );

                    extend_ref_group_cv_params_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentDetector,
                        cid,
                        inst_id,
                        &d.referenceable_param_group_ref,
                        ref_groups,
                    );

                    extend_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::ComponentDetector,
                        cid,
                        inst_id,
                        &d.cv_param,
                    );
                }
            }

            items.push(GlobalMetaItem {
                cvs: out,
                tags,
                owners,
                parents,
            });
        }
    }
    let n_instrument_configs = (items.len() - instr_start) as u32;

    let sw_start = items.len();
    if let Some(sw) = &mzml.software_list {
        for s in &sw.software {
            let mut out = Vec::new();
            let mut tags = Vec::new();
            let mut owners = Vec::new();
            let mut parents = Vec::new();

            let sw_id = id_gen.alloc();

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Software,
                sw_id,
                0,
                ACC_ATTR_ID,
                s.id.as_str(),
            );

            let ver = s
                .version
                .as_deref()
                .or_else(|| s.software_param.first().and_then(|p| p.version.as_deref()));
            if let Some(ver) = ver {
                push_attr_string_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::Software,
                    sw_id,
                    0,
                    ACC_ATTR_VERSION,
                    ver,
                );
            }

            for p in &s.software_param {
                let sp_id = id_gen.alloc();

                push_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::SoftwareParam,
                    sp_id,
                    sw_id,
                    CvParam {
                        cv_ref: p.cv_ref.clone(),
                        accession: Some(p.accession.clone()),
                        name: p.name.clone(),
                        value: Some(String::new()),
                        unit_cv_ref: None,
                        unit_name: None,
                        unit_accession: None,
                    },
                );
            }

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::Software,
                sw_id,
                0,
                &s.cv_param,
            );

            items.push(GlobalMetaItem {
                cvs: out,
                tags,
                owners,
                parents,
            });
        }
    }
    let n_software = (items.len() - sw_start) as u32;

    let dp_start = items.len();
    if let Some(dpl) = &mzml.data_processing_list {
        for dp in &dpl.data_processing {
            let mut out = Vec::new();
            let mut tags = Vec::new();
            let mut owners = Vec::new();
            let mut parents = Vec::new();

            let dp_id = id_gen.alloc();

            push_attr_string_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::DataProcessing,
                dp_id,
                0,
                ACC_ATTR_ID,
                dp.id.as_str(),
            );

            for m in &dp.processing_method {
                let pm_id = id_gen.alloc();

                extend_ref_group_cv_params_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::ProcessingMethod,
                    pm_id,
                    dp_id,
                    &m.referenceable_param_group_ref,
                    ref_groups,
                );

                extend_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::ProcessingMethod,
                    pm_id,
                    dp_id,
                    &m.cv_param,
                );
            }

            items.push(GlobalMetaItem {
                cvs: out,
                tags,
                owners,
                parents,
            });
        }
    }
    let n_data_processing = (items.len() - dp_start) as u32;

    let acq_start = items.len();
    if let Some(ssl) = &mzml.scan_settings_list {
        for ss in &ssl.scan_settings {
            let mut out = Vec::new();
            let mut tags = Vec::new();
            let mut owners = Vec::new();
            let mut parents = Vec::new();

            let ss_id = id_gen.alloc();

            if let Some(id) = ss.id.as_deref() {
                push_attr_string_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::ScanSettings,
                    ss_id,
                    0,
                    ACC_ATTR_ID,
                    id,
                );
            }

            if let Some(icr) = ss.instrument_configuration_ref.as_deref() {
                push_attr_string_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::ScanSettings,
                    ss_id,
                    0,
                    ACC_ATTR_INSTRUMENT_CONFIGURATION_REF,
                    icr,
                );
            }

            if let Some(sfrl) = &ss.source_file_ref_list {
                let sfrl_id = id_gen.alloc();

                push_attr_usize_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::SourceFileRefList,
                    sfrl_id,
                    ss_id,
                    ACC_ATTR_COUNT,
                    Some(sfrl.source_file_refs.len() as u32),
                );

                for sref in &sfrl.source_file_refs {
                    let sfr_id = id_gen.alloc();
                    push_attr_string_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::SourceFileRef,
                        sfr_id,
                        sfrl_id,
                        ACC_ATTR_REF,
                        sref.r#ref.as_str(),
                    );
                }
            }

            extend_ref_group_cv_params_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::ScanSettings,
                ss_id,
                0,
                &ss.referenceable_param_group_refs,
                ref_groups,
            );

            extend_tagged_ids(
                &mut out,
                &mut tags,
                &mut owners,
                &mut parents,
                TagId::ScanSettings,
                ss_id,
                0,
                &ss.cv_params,
            );

            if let Some(tl) = &ss.target_list {
                for t in &tl.targets {
                    let tgt_id = id_gen.alloc();

                    extend_ref_group_cv_params_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::Target,
                        tgt_id,
                        ss_id,
                        &t.referenceable_param_group_refs,
                        ref_groups,
                    );

                    extend_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::Target,
                        tgt_id,
                        ss_id,
                        &t.cv_params,
                    );
                }
            }

            items.push(GlobalMetaItem {
                cvs: out,
                tags,
                owners,
                parents,
            });
        }
    }
    let n_acquisition_settings = (items.len() - acq_start) as u32;

    let cv_start = items.len();
    if let Some(cl) = &mzml.cv_list {
        let cv_count = cl.cv.len() as u32;

        if cv_count != 0 {
            let cv_list_id = id_gen.alloc();

            for (i, cv) in cl.cv.iter().enumerate() {
                let mut out = Vec::new();
                let mut tags = Vec::new();
                let mut owners = Vec::new();
                let mut parents = Vec::new();
                if i == 0 {
                    push_attr_usize_tagged_ids(
                        &mut out,
                        &mut tags,
                        &mut owners,
                        &mut parents,
                        TagId::CvList,
                        cv_list_id,
                        0,
                        ACC_ATTR_COUNT,
                        Some(cv_count),
                    );
                }

                let cv_id = id_gen.alloc();

                push_attr_string_tagged_ids(
                    &mut out,
                    &mut tags,
                    &mut owners,
                    &mut parents,
                    TagId::Cv,
                    cv_id,
                    cv_list_id,
                    ACC_ATTR_LABEL,
                    cv.id.as_str(),
                );

                if let Some(v) = cv.full_name.as_deref() {
                    if !v.is_empty() {
                        push_attr_string_tagged_ids(
                            &mut out,
                            &mut tags,
                            &mut owners,
                            &mut parents,
                            TagId::Cv,
                            cv_id,
                            cv_list_id,
                            ACC_ATTR_CV_FULL_NAME,
                            v,
                        );
                    }
                }

                if let Some(v) = cv.version.as_deref() {
                    if !v.is_empty() {
                        push_attr_string_tagged_ids(
                            &mut out,
                            &mut tags,
                            &mut owners,
                            &mut parents,
                            TagId::Cv,
                            cv_id,
                            cv_list_id,
                            ACC_ATTR_CV_VERSION,
                            v,
                        );
                    }
                }

                if let Some(v) = cv.uri.as_deref() {
                    if !v.is_empty() {
                        push_attr_string_tagged_ids(
                            &mut out,
                            &mut tags,
                            &mut owners,
                            &mut parents,
                            TagId::Cv,
                            cv_id,
                            cv_list_id,
                            ACC_ATTR_CV_URI,
                            v,
                        );
                    }
                }

                items.push(GlobalMetaItem {
                    cvs: out,
                    tags,
                    owners,
                    parents,
                });
            }
        }
    }

    let n_cvs = (items.len() - cv_start) as u32;

    (
        items,
        GlobalCounts {
            n_file_description,
            n_run,
            n_ref_param_groups,
            n_samples,
            n_instrument_configs,
            n_software,
            n_data_processing,
            n_acquisition_settings,
            n_cvs,
        },
    )
}

#[inline]
fn pack_cv_param(
    tag_id: u8,
    owner_id: u32,
    parent_owner_id: u32,
    cv: &CvParam,
    tag_ids: &mut Vec<u8>,
    owner_ids: &mut Vec<u32>,
    parent_indices: &mut Vec<u32>,
    ref_codes: &mut Vec<u8>,
    accession_numbers: &mut Vec<u32>,
    unit_ref_codes: &mut Vec<u8>,
    unit_accession_numbers: &mut Vec<u32>,
    value_kinds: &mut Vec<u8>,
    value_indices: &mut Vec<u32>,
    numeric_values: &mut Vec<f64>,
    string_offsets: &mut Vec<u32>,
    string_lengths: &mut Vec<u32>,
    string_bytes: &mut Vec<u8>,
    numeric_index: &mut u32,
    string_index: &mut u32,
) {
    tag_ids.push(tag_id);
    owner_ids.push(owner_id);
    parent_indices.push(parent_owner_id);

    let cv_ref = cv_ref_from_accession(cv.accession.as_deref()).or(cv.cv_ref.as_deref());
    ref_codes.push(cv_ref_code_from_str(cv_ref));

    accession_numbers.push(parse_accession_tail(cv.accession.as_deref()));

    let unit_ref =
        cv_ref_from_accession(cv.unit_accession.as_deref()).or(cv.unit_cv_ref.as_deref());
    unit_ref_codes.push(cv_ref_code_from_str(unit_ref));

    unit_accession_numbers.push(parse_accession_tail(cv.unit_accession.as_deref()));

    let (kind, idx) = match cv.value.as_deref() {
        None | Some("") => (2u8, 0u32),
        Some(val) => {
            if let Ok(num) = val.parse::<f64>() {
                let i = *numeric_index;
                numeric_values.push(num);
                *numeric_index += 1;
                (0u8, i)
            } else {
                let i = *string_index;
                let bytes = val.as_bytes();
                let off = string_bytes.len() as u32;
                let len = bytes.len() as u32;

                string_bytes.extend_from_slice(bytes);
                string_offsets.push(off);
                string_lengths.push(len);

                *string_index += 1;
                (1u8, i)
            }
        }
    };

    value_kinds.push(kind);
    value_indices.push(idx);
}

#[inline]
fn cv_ref_from_accession<'a>(acc: Option<&'a str>) -> Option<&'a str> {
    acc.and_then(|s| s.split_once(':').map(|(p, _)| p))
}

/// <cvParam>
fn pack_meta_streaming<T, F>(items: &[T], mut fill: F) -> PackedMeta
where
    F: FnMut(&mut Vec<CvParam>, &mut Vec<u8>, &mut Vec<u32>, &mut Vec<u32>, &T),
{
    let item_count = items.len();

    let mut index_offsets = Vec::with_capacity(item_count + 1);

    let mut owner_ids: Vec<u32> = Vec::new();
    let mut parent_indices: Vec<u32> = Vec::new();

    let mut tag_ids: Vec<u8> = Vec::new();
    let mut ref_codes: Vec<u8> = Vec::new();
    let mut accession_numbers: Vec<u32> = Vec::new();
    let mut unit_ref_codes: Vec<u8> = Vec::new();
    let mut unit_accession_numbers: Vec<u32> = Vec::new();
    let mut value_kinds: Vec<u8> = Vec::new();
    let mut value_indices: Vec<u32> = Vec::new();

    let mut numeric_values: Vec<f64> = Vec::new();
    let mut string_offsets: Vec<u32> = Vec::new();
    let mut string_lengths: Vec<u32> = Vec::new();
    let mut string_bytes: Vec<u8> = Vec::new();

    let mut scratch: Vec<CvParam> = Vec::new();
    let mut scratch_tags: Vec<u8> = Vec::new();
    let mut scratch_owners: Vec<u32> = Vec::new();
    let mut scratch_parents: Vec<u32> = Vec::new();

    let mut numeric_index: u32 = 0;
    let mut string_index: u32 = 0;
    let mut meta_index: u32 = 0;

    index_offsets.push(0);

    for item in items {
        scratch.clear();
        scratch_tags.clear();
        scratch_owners.clear();
        scratch_parents.clear();

        fill(
            &mut scratch,
            &mut scratch_tags,
            &mut scratch_owners,
            &mut scratch_parents,
            item,
        );

        debug_assert_eq!(scratch.len(), scratch_tags.len());
        debug_assert_eq!(scratch.len(), scratch_owners.len());
        debug_assert_eq!(scratch.len(), scratch_parents.len());

        let n = scratch.len();

        for i in 0..n {
            pack_cv_param(
                scratch_tags[i],
                scratch_owners[i],
                scratch_parents[i],
                &scratch[i],
                &mut tag_ids,
                &mut owner_ids,
                &mut parent_indices,
                &mut ref_codes,
                &mut accession_numbers,
                &mut unit_ref_codes,
                &mut unit_accession_numbers,
                &mut value_kinds,
                &mut value_indices,
                &mut numeric_values,
                &mut string_offsets,
                &mut string_lengths,
                &mut string_bytes,
                &mut numeric_index,
                &mut string_index,
            );
            meta_index += 1;
        }

        index_offsets.push(meta_index);
    }

    PackedMeta {
        index_offsets,
        owner_ids,
        parent_indices,
        tag_ids,
        ref_codes,
        accession_numbers,
        unit_ref_codes,
        unit_accession_numbers,
        value_kinds,
        value_indices,
        numeric_values,
        string_offsets,
        string_lengths,
        string_bytes,
    }
}

fn pack_meta_slices<T, F>(items: &[T], meta_of: F) -> PackedMeta
where
    F: Fn(&T) -> (&[CvParam], &[u8], &[u32], &[u32]),
{
    let item_count = items.len();

    let mut total_meta_count = 0usize;
    for item in items {
        total_meta_count += meta_of(item).0.len();
    }

    let mut index_offsets = Vec::with_capacity(item_count + 1);

    let mut owner_ids = Vec::with_capacity(total_meta_count);
    let mut parent_indices = Vec::with_capacity(total_meta_count);

    let mut tag_ids = Vec::with_capacity(total_meta_count);
    let mut ref_codes = Vec::with_capacity(total_meta_count);
    let mut accession_numbers = Vec::with_capacity(total_meta_count);
    let mut unit_ref_codes = Vec::with_capacity(total_meta_count);
    let mut unit_accession_numbers = Vec::with_capacity(total_meta_count);
    let mut value_kinds = Vec::with_capacity(total_meta_count);
    let mut value_indices = Vec::with_capacity(total_meta_count);

    let mut numeric_values = Vec::with_capacity(total_meta_count);
    let mut string_offsets = Vec::with_capacity(total_meta_count);
    let mut string_lengths = Vec::with_capacity(total_meta_count);
    let mut string_bytes = Vec::new();

    let mut numeric_index: u32 = 0;
    let mut string_index: u32 = 0;
    let mut meta_index: u32 = 0;

    index_offsets.push(0);

    for item in items {
        let (xs, ts, os, ps) = meta_of(item);
        debug_assert_eq!(xs.len(), ts.len());
        debug_assert_eq!(xs.len(), os.len());
        debug_assert_eq!(xs.len(), ps.len());

        for i in 0..xs.len() {
            pack_cv_param(
                ts[i],
                os[i],
                ps[i],
                &xs[i],
                &mut tag_ids,
                &mut owner_ids,
                &mut parent_indices,
                &mut ref_codes,
                &mut accession_numbers,
                &mut unit_ref_codes,
                &mut unit_accession_numbers,
                &mut value_kinds,
                &mut value_indices,
                &mut numeric_values,
                &mut string_offsets,
                &mut string_lengths,
                &mut string_bytes,
                &mut numeric_index,
                &mut string_index,
            );
            meta_index += 1;
        }

        index_offsets.push(meta_index);
    }

    PackedMeta {
        index_offsets,
        owner_ids,
        parent_indices,
        tag_ids,
        ref_codes,
        accession_numbers,
        unit_ref_codes,
        unit_accession_numbers,
        value_kinds,
        value_indices,
        numeric_values,
        string_offsets,
        string_lengths,
        string_bytes,
    }
}

fn packed_meta_byte_len(meta: &PackedMeta) -> usize {
    meta.index_offsets.len() * 4
        + meta.owner_ids.len() * 4
        + meta.parent_indices.len() * 4
        + meta.tag_ids.len()
        + meta.ref_codes.len()
        + meta.accession_numbers.len() * 4
        + meta.unit_ref_codes.len()
        + meta.unit_accession_numbers.len() * 4
        + meta.value_kinds.len()
        + meta.value_indices.len() * 4
        + meta.numeric_values.len() * 8
        + meta.string_offsets.len() * 4
        + meta.string_lengths.len() * 4
        + meta.string_bytes.len()
}

fn write_packed_meta_into(buf: &mut Vec<u8>, meta: &PackedMeta) {
    // CI
    for &v in &meta.index_offsets {
        write_u32_le(buf, v);
    }
    // MOI
    for &v in &meta.owner_ids {
        write_u32_le(buf, v);
    }
    // MPI
    for &v in &meta.parent_indices {
        write_u32_le(buf, v);
    }

    // MTI, MRI
    buf.extend_from_slice(&meta.tag_ids);
    buf.extend_from_slice(&meta.ref_codes);

    // MAN
    for &v in &meta.accession_numbers {
        write_u32_le(buf, v);
    }

    // MURI
    buf.extend_from_slice(&meta.unit_ref_codes);

    // MUAN
    for &v in &meta.unit_accession_numbers {
        write_u32_le(buf, v);
    }

    // VK
    buf.extend_from_slice(&meta.value_kinds);

    // VI
    for &v in &meta.value_indices {
        write_u32_le(buf, v);
    }

    // VN
    for &v in &meta.numeric_values {
        write_f64_le(buf, v);
    }

    // VOFF / VLEN / VS
    for &v in &meta.string_offsets {
        write_u32_le(buf, v);
    }
    for &v in &meta.string_lengths {
        write_u32_le(buf, v);
    }
    buf.extend_from_slice(&meta.string_bytes);
}

fn write_packed_meta_bytes(meta: &PackedMeta) -> Vec<u8> {
    let mut buf = Vec::with_capacity(packed_meta_byte_len(meta));
    write_packed_meta_into(&mut buf, meta);
    buf
}

fn write_global_meta_bytes(counts: &GlobalCounts, meta: &PackedMeta) -> Vec<u8> {
    let mut buf = Vec::with_capacity(9 * 4 + packed_meta_byte_len(meta));

    write_u32_le(&mut buf, counts.n_file_description);
    write_u32_le(&mut buf, counts.n_run);
    write_u32_le(&mut buf, counts.n_ref_param_groups);
    write_u32_le(&mut buf, counts.n_samples);
    write_u32_le(&mut buf, counts.n_instrument_configs);
    write_u32_le(&mut buf, counts.n_software);
    write_u32_le(&mut buf, counts.n_data_processing);
    write_u32_le(&mut buf, counts.n_acquisition_settings);
    write_u32_le(&mut buf, counts.n_cvs);

    write_packed_meta_into(&mut buf, meta);
    buf
}

fn append_aligned_8(output: &mut Vec<u8>, bytes: &[u8]) -> u64 {
    let aligned = align_to_8(output.len());
    if aligned > output.len() {
        output.resize(aligned, 0);
    }
    let offset = output.len() as u64;
    output.extend_from_slice(bytes);
    offset
}

/// <binaryDataArray>
#[inline]
fn array_ref<'a>(ba: &'a BinaryDataArray) -> Option<ArrayRef<'a>> {
    match bda_declared_is_f64(ba) {
        Some(false) => {
            if !ba.decoded_binary_f32.is_empty() {
                return Some(ArrayRef::F32(&ba.decoded_binary_f32));
            }
            if !ba.decoded_binary_f64.is_empty() {
                return Some(ArrayRef::F64(&ba.decoded_binary_f64));
            }
        }
        Some(true) => {
            if !ba.decoded_binary_f64.is_empty() {
                return Some(ArrayRef::F64(&ba.decoded_binary_f64));
            }
            if !ba.decoded_binary_f32.is_empty() {
                return Some(ArrayRef::F32(&ba.decoded_binary_f32));
            }
        }
        None => {
            if !ba.decoded_binary_f64.is_empty() {
                return Some(ArrayRef::F64(&ba.decoded_binary_f64));
            }
            if !ba.decoded_binary_f32.is_empty() {
                return Some(ArrayRef::F32(&ba.decoded_binary_f32));
            }
        }
    }
    None
}

#[inline]
fn parse_accession_tail(accession: Option<&str>) -> u32 {
    let s = accession.unwrap_or("");
    let tail = match s.rsplit_once(':') {
        Some((_, t)) => t,
        None => s,
    };

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

#[inline]
fn align_to_8(x: usize) -> usize {
    (x + 7) & !7
}

#[inline]
fn write_u32_le(buf: &mut Vec<u8>, value: u32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

#[inline]
fn write_u64_le(buf: &mut Vec<u8>, value: u64) {
    buf.extend_from_slice(&value.to_le_bytes());
}

#[inline]
fn write_f64_le(buf: &mut Vec<u8>, value: f64) {
    buf.extend_from_slice(&value.to_le_bytes());
}

#[inline]
fn write_f32_le(buf: &mut Vec<u8>, value: f32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

#[inline]
fn write_f64_as_f32(buf: &mut Vec<u8>, value: f64) {
    buf.extend_from_slice(&(value as f32).to_le_bytes());
}

#[inline]
fn set_u8_at(buf: &mut [u8], offset: usize, value: u8) {
    buf[offset] = value;
}

#[inline]
fn set_u32_at(buf: &mut [u8], offset: usize, value: u32) {
    buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

#[inline]
fn set_u64_at(buf: &mut [u8], offset: usize, value: u64) {
    buf[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

#[inline]
fn bda_declared_is_f64(ba: &BinaryDataArray) -> Option<bool> {
    if ba.is_f64 == Some(true) {
        return Some(true);
    }
    if ba.is_f32 == Some(true) {
        return Some(false);
    }

    let mut saw32 = false;
    let mut saw64 = false;

    for cv in &ba.cv_params {
        match parse_accession_tail(cv.accession.as_deref()) {
            ACC_32BIT_FLOAT => saw32 = true,
            ACC_64BIT_FLOAT => saw64 = true,
            _ => {}
        }
    }

    match (saw32, saw64) {
        (true, false) => Some(false),
        (false, true) => Some(true),
        _ => None,
    }
}

#[inline]
fn find_xy_ba<'a>(
    list: Option<&'a BinaryDataArrayList>,
    x_accession_tail: u32,
    y_accession_tail: u32,
) -> (Option<&'a BinaryDataArray>, Option<&'a BinaryDataArray>) {
    let mut x: Option<&'a BinaryDataArray> = None;
    let mut y: Option<&'a BinaryDataArray> = None;

    if let Some(list) = list {
        for ba in &list.binary_data_arrays {
            let mut is_x = false;
            let mut is_y = false;

            for cv in &ba.cv_params {
                let tail = parse_accession_tail(cv.accession.as_deref());
                if tail == x_accession_tail {
                    is_x = true;
                } else if tail == y_accession_tail {
                    is_y = true;
                }
            }

            if x.is_none() && is_x {
                x = Some(ba);
            }
            if y.is_none() && is_y {
                y = Some(ba);
            }

            if x.is_some() && y.is_some() {
                break;
            }
        }
    }

    (x, y)
}

/// <mzML>
pub fn encode(mzml: &MzML, compression_level: u8, f32_compress: bool) -> Vec<u8> {
    assert!(compression_level <= 22);

    #[inline]
    fn fix_attr_values(out: &mut Vec<CvParam>) {
        for cv in out.iter_mut() {
            if cv.cv_ref.as_deref() == Some(CV_REF_ATTR) {
                let empty_val = cv.value.as_deref().map_or(true, |s| s.is_empty());
                if empty_val && !cv.name.is_empty() {
                    cv.value = Some(std::mem::take(&mut cv.name));
                }
            }
        }
    }

    let compress_meta = compression_level != 0;
    let do_shuffle = compress_meta;

    let array_filter_id = if do_shuffle {
        ARRAY_FILTER_BYTE_SHUFFLE
    } else {
        ARRAY_FILTER_NONE
    };

    let run = &mzml.run;

    let spectra: &[Spectrum] = run
        .spectrum_list
        .as_ref()
        .map(|sl| sl.spectra.as_slice())
        .unwrap_or(&[]);

    let chromatograms: &[Chromatogram] = run
        .chromatogram_list
        .as_ref()
        .map(|cl| cl.chromatograms.as_slice())
        .unwrap_or(&[]);

    let spectrum_count = spectra.len() as u32;
    let chrom_count = chromatograms.len() as u32;

    let mut spectrum_x_decl_f64: Option<bool> = None;
    let mut spectrum_y_decl_f64: Option<bool> = None;
    let mut spectrum_xy_cache: Vec<(ArrayRef<'_>, ArrayRef<'_>)> =
        Vec::with_capacity(spectra.len());

    for s in spectra {
        let (xba, yba) = find_xy_ba(
            s.binary_data_array_list.as_ref(),
            ACC_MZ_ARRAY,
            ACC_INTENSITY_ARRAY,
        );

        let x_arr = xba.and_then(array_ref).unwrap_or(ArrayRef::F32(&[]));
        let y_arr = yba.and_then(array_ref).unwrap_or(ArrayRef::F32(&[]));

        if let Some(xba) = xba {
            let x_decl_f64 = bda_declared_is_f64(xba).unwrap_or(matches!(x_arr, ArrayRef::F64(_)));
            merge_declared_width(&mut spectrum_x_decl_f64, x_decl_f64, "spectrum x");
        }

        if let Some(yba) = yba {
            let y_decl_f64 = bda_declared_is_f64(yba).unwrap_or(matches!(y_arr, ArrayRef::F64(_)));
            merge_declared_width(&mut spectrum_y_decl_f64, y_decl_f64, "spectrum y");
        }

        spectrum_xy_cache.push((x_arr, y_arr));
    }

    let mut chrom_x_decl_f64: Option<bool> = None;
    let mut chrom_y_decl_f64: Option<bool> = None;
    let mut chrom_xy_cache: Vec<(ArrayRef<'_>, ArrayRef<'_>)> =
        Vec::with_capacity(chromatograms.len());

    for c in chromatograms {
        let (xba, yba) = find_xy_ba(
            c.binary_data_array_list.as_ref(),
            ACC_TIME_ARRAY,
            ACC_INTENSITY_ARRAY,
        );

        let x_arr = xba.and_then(array_ref).unwrap_or(ArrayRef::F32(&[]));
        let y_arr = yba.and_then(array_ref).unwrap_or(ArrayRef::F32(&[]));

        if let Some(xba) = xba {
            let x_decl_f64 = bda_declared_is_f64(xba).unwrap_or(matches!(x_arr, ArrayRef::F64(_)));
            merge_declared_width(&mut chrom_x_decl_f64, x_decl_f64, "chromatogram x");
        }

        if let Some(yba) = yba {
            let y_decl_f64 = bda_declared_is_f64(yba).unwrap_or(matches!(y_arr, ArrayRef::F64(_)));
            merge_declared_width(&mut chrom_y_decl_f64, y_decl_f64, "chromatogram y");
        }

        chrom_xy_cache.push((x_arr, y_arr));
    }

    let spect_x_store_f64 = if f32_compress {
        false
    } else {
        spectrum_x_decl_f64.unwrap_or(false)
    };

    let spect_y_store_f64 = if f32_compress {
        false
    } else {
        spectrum_y_decl_f64.unwrap_or(false)
    };

    let chrom_x_store_f64 = if f32_compress {
        false
    } else {
        chrom_x_decl_f64.unwrap_or(false)
    };

    let chrom_y_store_f64 = if f32_compress {
        false
    } else {
        chrom_y_decl_f64.unwrap_or(false)
    };

    let spec_x_elem_size = elem_size(spect_x_store_f64);
    let spec_y_elem_size = elem_size(spect_y_store_f64);
    let chrom_x_elem_size = elem_size(chrom_x_store_f64);
    let chrom_y_elem_size = elem_size(chrom_y_store_f64);

    let ref_groups = build_ref_group_map(mzml);

    let mut id_gen = NodeIdGen::new();

    let spectrum_list_owner_id: u32 = if run.spectrum_list.is_some() {
        id_gen.alloc()
    } else {
        0
    };
    let chromatogram_list_owner_id: u32 = if run.chromatogram_list.is_some() {
        id_gen.alloc()
    } else {
        0
    };

    let (mut global_items, global_counts) = build_global_meta_items(mzml, &ref_groups, &mut id_gen);
    for item in &mut global_items {
        fix_attr_values(&mut item.cvs);
    }

    // Spectra
    let mut spec_i: usize = 0;

    let spectrum_meta = pack_meta_streaming(spectra, |out, tags, owners, parents, s| {
        let idx = spec_i;
        spec_i += 1;
        if idx == 0 {
            if let Some(sl) = run.spectrum_list.as_ref() {
                if spectrum_list_owner_id != 0 {
                    push_schema_attributes(
                        out,
                        tags,
                        owners,
                        parents,
                        TagId::SpectrumList,
                        spectrum_list_owner_id,
                        0,
                        sl,
                    );
                }
            }
        }

        let spectrum_id = id_gen.alloc();

        flatten_spectrum_metadata_into_owned(
            out,
            tags,
            owners,
            parents,
            s,
            &ref_groups,
            spectrum_id,
            0,
            &mut id_gen,
            ACC_MZ_ARRAY,
            ACC_INTENSITY_ARRAY,
            spect_x_store_f64,
            spect_y_store_f64,
            f32_compress,
        );

        fix_attr_values(out);
    });

    let mut chrom_i: usize = 0;

    let chromatogram_meta = pack_meta_streaming(chromatograms, |out, tags, owners, parents, c| {
        let idx = chrom_i;
        chrom_i += 1;
        if idx == 0 {
            if let Some(cl) = run.chromatogram_list.as_ref() {
                if chromatogram_list_owner_id != 0 {
                    push_schema_attributes(
                        out,
                        tags,
                        owners,
                        parents,
                        TagId::ChromatogramList,
                        chromatogram_list_owner_id,
                        0,
                        cl,
                    );
                }
            }
        }

        let chrom_id = id_gen.alloc();

        flatten_chromatogram_metadata_into(
            out,
            tags,
            owners,
            parents,
            c,
            &ref_groups,
            chrom_id,
            0,
            &mut id_gen,
            ACC_TIME_ARRAY,
            ACC_INTENSITY_ARRAY,
            chrom_x_store_f64,
            chrom_y_store_f64,
            f32_compress,
        );

        fix_attr_values(out);
    });

    let global_meta = pack_meta_slices(&global_items, |m| {
        (
            m.cvs.as_slice(),
            m.tags.as_slice(),
            m.owners.as_slice(),
            m.parents.as_slice(),
        )
    });

    let spec_meta_count = spectrum_meta.ref_codes.len() as u32;
    let spec_num_count = spectrum_meta.numeric_values.len() as u32;
    let spec_str_count = spectrum_meta.string_offsets.len() as u32;

    let chrom_meta_count = chromatogram_meta.ref_codes.len() as u32;
    let chrom_num_count = chromatogram_meta.numeric_values.len() as u32;
    let chrom_str_count = chromatogram_meta.string_offsets.len() as u32;

    let global_meta_count = global_meta.ref_codes.len() as u32;
    let global_num_count = global_meta.numeric_values.len() as u32;
    let global_str_count = global_meta.string_offsets.len() as u32;

    let mut spectrum_meta_bytes = write_packed_meta_bytes(&spectrum_meta);
    let mut chromatogram_meta_bytes = write_packed_meta_bytes(&chromatogram_meta);
    let mut global_meta_bytes = write_global_meta_bytes(&global_counts, &global_meta);

    if compress_meta {
        spectrum_meta_bytes = compress_bytes(&spectrum_meta_bytes, compression_level);
        chromatogram_meta_bytes = compress_bytes(&chromatogram_meta_bytes, compression_level);
        global_meta_bytes = compress_bytes(&global_meta_bytes, compression_level);
    }

    let mut spec_x_builder = ContainerBuilder::new(
        TARGET_BLOCK_UNCOMP_BYTES,
        compression_level,
        spec_x_elem_size,
        do_shuffle,
    );
    let mut spec_y_builder = ContainerBuilder::new(
        TARGET_BLOCK_UNCOMP_BYTES,
        compression_level,
        spec_y_elem_size,
        do_shuffle,
    );
    let mut chrom_x_builder = ContainerBuilder::new(
        TARGET_BLOCK_UNCOMP_BYTES,
        compression_level,
        chrom_x_elem_size,
        do_shuffle,
    );
    let mut chrom_y_builder = ContainerBuilder::new(
        TARGET_BLOCK_UNCOMP_BYTES,
        compression_level,
        chrom_y_elem_size,
        do_shuffle,
    );

    let mut spec_index_bytes = Vec::with_capacity(spectra.len() * INDEX_ENTRY_SIZE);
    let mut chrom_index_bytes = Vec::with_capacity(chromatograms.len() * INDEX_ENTRY_SIZE);

    let mut spec_x_off_elems: u64 = 0;
    let mut spec_y_off_elems: u64 = 0;

    for &(x, y) in &spectrum_xy_cache {
        let x_len = x.len() as u32;
        let y_len = y.len() as u32;

        let x_item_bytes = x.len() * spec_x_elem_size;
        let y_item_bytes = y.len() * spec_y_elem_size;

        let x_block_id = spec_x_builder.write_item(x_item_bytes, |buf| {
            write_array(buf, x, spect_x_store_f64);
        });
        let y_block_id = spec_y_builder.write_item(y_item_bytes, |buf| {
            write_array(buf, y, spect_y_store_f64);
        });

        write_u64_le(&mut spec_index_bytes, spec_x_off_elems);
        write_u64_le(&mut spec_index_bytes, spec_y_off_elems);
        write_u32_le(&mut spec_index_bytes, x_len);
        write_u32_le(&mut spec_index_bytes, y_len);
        write_u32_le(&mut spec_index_bytes, x_block_id);
        write_u32_le(&mut spec_index_bytes, y_block_id);

        spec_x_off_elems += x_len as u64;
        spec_y_off_elems += y_len as u64;
    }

    let mut chrom_x_off_elems: u64 = 0;
    let mut chrom_y_off_elems: u64 = 0;

    for &(x, y) in &chrom_xy_cache {
        let x_len = x.len() as u32;
        let y_len = y.len() as u32;

        let x_item_bytes = x.len() * chrom_x_elem_size;
        let y_item_bytes = y.len() * chrom_y_elem_size;

        let x_block_id = chrom_x_builder.write_item(x_item_bytes, |buf| {
            write_array(buf, x, chrom_x_store_f64);
        });
        let y_block_id = chrom_y_builder.write_item(y_item_bytes, |buf| {
            write_array(buf, y, chrom_y_store_f64);
        });

        write_u64_le(&mut chrom_index_bytes, chrom_x_off_elems);
        write_u64_le(&mut chrom_index_bytes, chrom_y_off_elems);
        write_u32_le(&mut chrom_index_bytes, x_len);
        write_u32_le(&mut chrom_index_bytes, y_len);
        write_u32_le(&mut chrom_index_bytes, x_block_id);
        write_u32_le(&mut chrom_index_bytes, y_block_id);

        chrom_x_off_elems += x_len as u64;
        chrom_y_off_elems += y_len as u64;
    }

    let (container_spect_x, block_count_spect_x) = spec_x_builder.finalize();
    let (container_spect_y, block_count_spect_y) = spec_y_builder.finalize();
    let (container_chrom_x, block_count_chrom_x) = chrom_x_builder.finalize();
    let (container_chrom_y, block_count_chrom_y) = chrom_y_builder.finalize();

    let mut output = Vec::with_capacity(
        HEADER_SIZE
            + spec_index_bytes.len()
            + chrom_index_bytes.len()
            + spectrum_meta_bytes.len()
            + chromatogram_meta_bytes.len()
            + global_meta_bytes.len()
            + container_spect_x.len()
            + container_spect_y.len()
            + container_chrom_x.len()
            + container_chrom_y.len()
            + 64,
    );

    output.resize(HEADER_SIZE, 0);

    let off_spec_index = HEADER_SIZE as u64;
    output.extend_from_slice(&spec_index_bytes);

    let off_chrom_index = output.len() as u64;
    output.extend_from_slice(&chrom_index_bytes);

    let off_spec_meta = append_aligned_8(&mut output, &spectrum_meta_bytes);
    let off_chrom_meta = append_aligned_8(&mut output, &chromatogram_meta_bytes);
    let off_global_meta = append_aligned_8(&mut output, &global_meta_bytes);

    let off_container_spect_x = append_aligned_8(&mut output, &container_spect_x);
    let size_container_spect_x = container_spect_x.len() as u64;

    let off_container_spect_y = append_aligned_8(&mut output, &container_spect_y);
    let size_container_spect_y = container_spect_y.len() as u64;

    let off_container_chrom_x = append_aligned_8(&mut output, &container_chrom_x);
    let size_container_chrom_x = container_chrom_x.len() as u64;

    let off_container_chrom_y = append_aligned_8(&mut output, &container_chrom_y);
    let size_container_chrom_y = container_chrom_y.len() as u64;

    {
        let header = &mut output[0..HEADER_SIZE];

        header[0..4].copy_from_slice(b"B000");
        set_u8_at(header, 4, 0);
        header[5] = 1;
        header[6] = 0;
        header[7] = 0;

        set_u64_at(header, 8, off_spec_index);
        set_u64_at(header, 16, off_chrom_index);
        set_u64_at(header, 24, off_spec_meta);
        set_u64_at(header, 32, off_chrom_meta);
        set_u64_at(header, 40, off_global_meta);

        set_u64_at(header, 48, size_container_spect_x);
        set_u64_at(header, 56, off_container_spect_x);

        set_u64_at(header, 64, size_container_spect_y);
        set_u64_at(header, 72, off_container_spect_y);

        set_u64_at(header, 80, size_container_chrom_x);
        set_u64_at(header, 88, off_container_chrom_x);

        set_u64_at(header, 96, size_container_chrom_y);
        set_u64_at(header, 104, off_container_chrom_y);

        set_u32_at(header, 112, spectrum_count);
        set_u32_at(header, 116, chrom_count);

        set_u32_at(header, 120, spec_meta_count);
        set_u32_at(header, 124, spec_num_count);
        set_u32_at(header, 128, spec_str_count);

        set_u32_at(header, 132, chrom_meta_count);
        set_u32_at(header, 136, chrom_num_count);
        set_u32_at(header, 140, chrom_str_count);

        set_u32_at(header, 144, global_meta_count);
        set_u32_at(header, 148, global_num_count);
        set_u32_at(header, 152, global_str_count);

        set_u32_at(header, 156, block_count_spect_x);
        set_u32_at(header, 160, block_count_spect_y);
        set_u32_at(header, 164, block_count_chrom_x);
        set_u32_at(header, 168, block_count_chrom_y);

        set_u8_at(
            header,
            172,
            header_codec_and_flags(HDR_CODEC_ZSTD, compress_meta, compress_meta, compress_meta),
        );

        set_u8_at(header, 173, if chrom_x_store_f64 { 2 } else { 1 });
        set_u8_at(header, 174, if chrom_y_store_f64 { 2 } else { 1 });
        set_u8_at(header, 175, if spect_x_store_f64 { 2 } else { 1 });
        set_u8_at(header, 176, if spect_y_store_f64 { 2 } else { 1 });

        set_u8_at(header, 177, compression_level);
        set_u8_at(header, HDR_ARRAY_FILTER_OFF, array_filter_id);
    }

    output
}

#[inline]
fn push_assigned_attributes_as_cv_params(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    attrs: Vec<crate::b64::decode::Metadatum>,
) {
    for m in attrs {
        let tail = parse_accession_tail(m.accession.as_deref());
        if tail == 0 {
            continue;
        }

        match m.value {
            MetadatumValue::Text(v) => {
                if v.is_empty() {
                    continue;
                }
                push_tagged_ids(
                    out,
                    tags,
                    owners,
                    parents,
                    tag,
                    owner_id,
                    parent_owner_id,
                    attr_cv_param(tail, &v),
                );
            }
            MetadatumValue::Number(n) => {
                let s = n.to_string();
                if s.is_empty() {
                    continue;
                }
                push_tagged_ids(
                    out,
                    tags,
                    owners,
                    parents,
                    tag,
                    owner_id,
                    parent_owner_id,
                    attr_cv_param(tail, &s),
                );
            }
            MetadatumValue::Empty => {
                push_tagged_ids(
                    out,
                    tags,
                    owners,
                    parents,
                    tag,
                    owner_id,
                    parent_owner_id,
                    attr_cv_param(tail, ""),
                );
            }
        }
    }
}

#[inline]
fn push_schema_attributes<T: Serialize>(
    out: &mut Vec<CvParam>,
    tags: &mut Vec<u8>,
    owners: &mut Vec<u32>,
    parents: &mut Vec<u32>,
    tag: TagId,
    owner_id: u32,
    parent_owner_id: u32,
    expected: &T,
) {
    let attrs = assign_attributes(expected, tag, owner_id, parent_owner_id);
    for m in attrs {
        let tail = parse_accession_tail(m.accession.as_deref());
        if tail == 0 {
            continue;
        }

        let s = match m.value {
            MetadatumValue::Text(v) => v,
            MetadatumValue::Number(n) => n.to_string(),
            MetadatumValue::Empty => continue,
        };

        if s.is_empty() {
            continue;
        }

        push_tagged_ids(
            out,
            tags,
            owners,
            parents,
            tag,
            owner_id,
            parent_owner_id,
            attr_cv_param(tail, &s),
        );
    }
}
