// use miniz_oxide::deflate::compress_to_vec_zlib;
use zstd::bulk::compress as zstd_compress;

use std::collections::HashMap;

use crate::utilities::attr_meta::*;
use crate::utilities::mzml::*;

#[derive(Debug)]
pub struct PackedMeta {
    pub index_offsets: Vec<u32>,
    pub ref_codes: Vec<u8>,
    pub accession_numbers: Vec<u32>,
    pub unit_ref_codes: Vec<u8>,
    pub unit_accession_numbers: Vec<u32>,
    pub value_kinds: Vec<u8>,
    pub value_indices: Vec<u32>,
    pub numeric_values: Vec<f64>,
    pub string_offsets: Vec<u32>,
    pub string_lengths: Vec<u32>,
    pub string_bytes: Vec<u8>,
}

#[derive(Debug)]
struct GlobalCounts {
    n_file_description: u32,
    n_ref_param_groups: u32,
    n_samples: u32,
    n_instrument_configs: u32,
    n_software: u32,
    n_data_processing: u32,
    n_acquisition_settings: u32,
    n_cvs: u32,
}

const HEADER_SIZE: usize = 192;
const INDEX_ENTRY_SIZE: usize = 32;
const BLOCK_DIR_ENTRY_SIZE: usize = 32;

const TARGET_BLOCK_UNCOMP_BYTES: usize = 512 * 1024 * 1024;

const ACC_MZ_ARRAY: u32 = 1_000_514;
const ACC_INTENSITY_ARRAY: u32 = 1_000_515;
const ACC_TIME_ARRAY: u32 = 1_000_595;

const ACC_32BIT_FLOAT: u32 = 1_000_521;
const ACC_64BIT_FLOAT: u32 = 1_000_523;

// const HDR_CODEC_ZLIB: u8 = 0;
const HDR_CODEC_ZSTD: u8 = 1;

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

type ArrayInfo<'a> = ArrayRef<'a>;

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

    let mut spectrum_x_has_f64 = false;
    let mut spectrum_y_has_f64 = false;
    let mut spectrum_xy_cache: Vec<(ArrayInfo<'_>, ArrayInfo<'_>)> =
        Vec::with_capacity(spectra.len());

    for s in spectra {
        let (x, y) = spectrum_xy(s);
        spectrum_x_has_f64 |= matches!(x, ArrayRef::F64(_));
        spectrum_y_has_f64 |= matches!(y, ArrayRef::F64(_));
        spectrum_xy_cache.push((x, y));
    }

    let mut chrom_x_has_f64 = false;
    let mut chrom_y_has_f64 = false;
    let mut chrom_xy_cache: Vec<(ArrayInfo<'_>, ArrayInfo<'_>)> =
        Vec::with_capacity(chromatograms.len());

    for c in chromatograms {
        let (x, y) = chromatogram_xy(c);
        chrom_x_has_f64 |= matches!(x, ArrayRef::F64(_));
        chrom_y_has_f64 |= matches!(y, ArrayRef::F64(_));
        chrom_xy_cache.push((x, y));
    }

    let spect_x_store_f64 = spectrum_x_has_f64 && !f32_compress;
    let spect_y_store_f64 = spectrum_y_has_f64 && !f32_compress;
    let chrom_x_store_f64 = chrom_x_has_f64 && !f32_compress;
    let chrom_y_store_f64 = chrom_y_has_f64 && !f32_compress;

    let spec_x_elem_size = elem_size(spect_x_store_f64);
    let spec_y_elem_size = elem_size(spect_y_store_f64);
    let chrom_x_elem_size = elem_size(chrom_x_store_f64);
    let chrom_y_elem_size = elem_size(chrom_y_store_f64);

    let ref_groups = build_ref_group_map(mzml);

    let (mut global_items, global_counts) = build_global_meta_items(mzml, &ref_groups);
    for v in &mut global_items {
        fix_attr_values(v);
    }

    let spectrum_meta = pack_meta_streaming(spectra, |out, s| {
        flatten_spectrum_metadata_into(
            out,
            s,
            &ref_groups,
            ACC_MZ_ARRAY,
            ACC_INTENSITY_ARRAY,
            spect_x_store_f64,
            spect_y_store_f64,
        );
        fix_attr_values(out);
    });

    let chromatogram_meta = pack_meta_streaming(chromatograms, |out, c| {
        flatten_chromatogram_metadata_into(
            out,
            c,
            &ref_groups,
            ACC_TIME_ARRAY,
            ACC_INTENSITY_ARRAY,
            chrom_x_store_f64,
            chrom_y_store_f64,
        );
        fix_attr_values(out);
    });

    let global_meta = pack_meta(&global_items, |m| m.as_slice());

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
        header[5] = 0;
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

/// <cvList>
fn build_global_meta_items(
    mzml: &MzML,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) -> (Vec<Vec<CvParam>>, GlobalCounts) {
    const CV_REF: &str = "NCIT";
    const ACC_CV_ID: u32 = 9_900_001;
    const ACC_CV_FULL_NAME: u32 = 9_900_002;
    const ACC_CV_VERSION: u32 = 9_900_003;
    const ACC_CV_URI: u32 = 9_900_004;

    fn cv_field(acc: u32, val: Option<&str>) -> CvParam {
        CvParam {
            cv_ref: Some(CV_REF.to_string()),
            accession: Some(format!("{}:{:07}", CV_REF, acc)),
            name: String::new(),
            value: val.map(|s| s.to_string()),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        }
    }

    let mut items: Vec<Vec<CvParam>> = Vec::new();

    {
        let fd = &mzml.file_description;
        let mut v = Vec::new();

        extend_ref_group_cv_params(
            &mut v,
            &fd.file_content.referenceable_param_group_refs,
            ref_groups,
        );
        v.extend_from_slice(&fd.file_content.cv_params);

        for sf in &fd.source_file_list.source_file {
            push_attr_string(&mut v, ACC_ATTR_ID, sf.id.as_str()); // ACC_ATTR_ID
            push_attr_string(&mut v, ACC_ATTR_NAME, sf.name.as_str()); // ACC_ATTR_NAME
            push_attr_string(&mut v, ACC_ATTR_LOCATION, sf.location.as_str()); // ACC_ATTR_LOCATION

            extend_ref_group_cv_params(&mut v, &sf.referenceable_param_group_ref, ref_groups);
            v.extend_from_slice(&sf.cv_param);
        }

        for c in &fd.contacts {
            extend_ref_group_cv_params(&mut v, &c.referenceable_param_group_refs, ref_groups);
            v.extend_from_slice(&c.cv_params);
        }

        items.push(v);
    }
    let n_file_description = 1u32;

    let ref_start = items.len();
    if let Some(rpgl) = &mzml.referenceable_param_group_list {
        for g in &rpgl.referenceable_param_groups {
            let mut v = Vec::new();
            push_attr_string(&mut v, ACC_ATTR_ID, g.id.as_str()); // ACC_ATTR_ID
            v.extend_from_slice(&g.cv_params);
            items.push(v);
        }
    }
    let n_ref_param_groups = (items.len() - ref_start) as u32;

    let samples_start = items.len();
    if let Some(sl) = &mzml.sample_list {
        for s in &sl.samples {
            let mut v = Vec::new();
            push_attr_string(&mut v, ACC_ATTR_ID, s.id.as_str()); // ACC_ATTR_ID
            push_attr_string(&mut v, ACC_ATTR_NAME, s.name.as_str()); // ACC_ATTR_NAME

            if let Some(r) = &s.referenceable_param_group_ref {
                extend_ref_group_cv_params(&mut v, std::slice::from_ref(r), ref_groups);
            }
            v.extend_from_slice(&s.cv_params);
            items.push(v);
        }
    }
    let n_samples = (items.len() - samples_start) as u32;

    let instr_start = items.len();
    if let Some(il) = &mzml.instrument_list {
        for ic in &il.instrument {
            let mut v = Vec::new();
            push_attr_string(&mut v, ACC_ATTR_ID, ic.id.as_str()); // ACC_ATTR_ID

            extend_ref_group_cv_params(&mut v, &ic.referenceable_param_group_ref, ref_groups);
            v.extend_from_slice(&ic.cv_param);

            if let Some(cl) = &ic.component_list {
                for s in &cl.source {
                    extend_ref_group_cv_params(
                        &mut v,
                        &s.referenceable_param_group_ref,
                        ref_groups,
                    );
                    v.extend_from_slice(&s.cv_param);
                }
                for a in &cl.analyzer {
                    extend_ref_group_cv_params(
                        &mut v,
                        &a.referenceable_param_group_ref,
                        ref_groups,
                    );
                    v.extend_from_slice(&a.cv_param);
                }
                for d in &cl.detector {
                    extend_ref_group_cv_params(
                        &mut v,
                        &d.referenceable_param_group_ref,
                        ref_groups,
                    );
                    v.extend_from_slice(&d.cv_param);
                }
            }

            items.push(v);
        }
    }
    let n_instrument_configs = (items.len() - instr_start) as u32;

    let sw_start = items.len();
    if let Some(sw) = &mzml.software_list {
        for s in &sw.software {
            let mut v = Vec::new();
            push_attr_string(&mut v, ACC_ATTR_ID, s.id.as_str()); // ACC_ATTR_ID
            let ver = s
                .version
                .as_deref()
                .or_else(|| s.software_param.first().and_then(|p| p.version.as_deref()));
            if let Some(ver) = ver {
                push_attr_string(&mut v, ACC_ATTR_VERSION, ver); // ACC_ATTR_VERSION
            }

            for p in &s.software_param {
                v.push(CvParam {
                    cv_ref: p.cv_ref.clone(),
                    accession: Some(p.accession.clone()),
                    name: p.name.clone(),
                    value: Some(String::new()),
                    unit_cv_ref: None,
                    unit_name: None,
                    unit_accession: None,
                });
            }
            v.extend_from_slice(&s.cv_param);

            items.push(v);
        }
    }
    let n_software = (items.len() - sw_start) as u32;

    let dp_start = items.len();
    if let Some(dpl) = &mzml.data_processing_list {
        for dp in &dpl.data_processing {
            let mut v = Vec::new();
            push_attr_string(&mut v, ACC_ATTR_ID, dp.id.as_str()); // ACC_ATTR_ID

            for m in &dp.processing_method {
                extend_ref_group_cv_params(&mut v, &m.referenceable_param_group_ref, ref_groups);
                v.extend_from_slice(&m.cv_param);
            }
            items.push(v);
        }
    }
    let n_data_processing = (items.len() - dp_start) as u32;

    let acq_start = items.len();
    if let Some(ssl) = &mzml.scan_settings_list {
        for ss in &ssl.scan_settings {
            let mut v = Vec::new();

            if let Some(id) = ss.id.as_deref() {
                push_attr_string(&mut v, ACC_ATTR_ID, id); // ACC_ATTR_ID
            }
            if let Some(icr) = ss.instrument_configuration_ref.as_deref() {
                push_attr_string(&mut v, ACC_ATTR_INSTRUMENT_CONFIGURATION_REF, icr); // ACC_ATTR_INSTRUMENT_CONFIGURATION_REF
            }

            if let Some(sfrl) = &ss.source_file_ref_list {
                for sref in &sfrl.source_file_refs {
                    push_attr_string(&mut v, ACC_ATTR_REF, sref.r#ref.as_str()); // ACC_ATTR_REF
                }
            }

            extend_ref_group_cv_params(&mut v, &ss.referenceable_param_group_refs, ref_groups);
            v.extend_from_slice(&ss.cv_params);

            if let Some(tl) = &ss.target_list {
                for t in &tl.targets {
                    extend_ref_group_cv_params(
                        &mut v,
                        &t.referenceable_param_group_refs,
                        ref_groups,
                    );
                    v.extend_from_slice(&t.cv_params);
                }
            }

            items.push(v);
        }
    }
    let n_acquisition_settings = (items.len() - acq_start) as u32;

    let cv_start = items.len();
    if let Some(cl) = &mzml.cv_list {
        for cv in &cl.cv {
            let mut v = Vec::new();
            v.push(cv_field(ACC_CV_ID, Some(cv.id.as_str())));
            v.push(cv_field(ACC_CV_FULL_NAME, cv.full_name.as_deref()));
            v.push(cv_field(ACC_CV_VERSION, cv.version.as_deref()));
            v.push(cv_field(ACC_CV_URI, cv.uri.as_deref()));
            items.push(v);
        }
    }
    let n_cvs = (items.len() - cv_start) as u32;

    (
        items,
        GlobalCounts {
            n_file_description,
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
fn extend_ref_group_cv_params(
    out: &mut Vec<CvParam>,
    refs: &[ReferenceableParamGroupRef],
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for r in refs {
        if let Some(g) = ref_groups.get(r.r#ref.as_str()) {
            out.extend_from_slice(&g.cv_params);
        }
    }
}

/// <spectrum>
fn flatten_spectrum_metadata_into(
    out: &mut Vec<CvParam>,
    spectrum: &Spectrum,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
    x_accession_tail: u32,
    y_accession_tail: u32,
    x_store_f64: bool,
    y_store_f64: bool,
) {
    push_attr_string(out, ACC_ATTR_ID, spectrum.id.as_str()); // ACC_ATTR_ID
    push_attr_u32(out, ACC_ATTR_INDEX, spectrum.index); // ACC_ATTR_INDEX
    push_attr_usize(
        out,
        ACC_ATTR_DEFAULT_ARRAY_LENGTH,
        spectrum.default_array_length,
    ); // ACC_ATTR_DEFAULT_ARRAY_LENGTH

    extend_ref_group_cv_params(out, &spectrum.referenceable_param_group_refs, ref_groups);
    out.extend_from_slice(&spectrum.cv_params);

    if let Some(sd) = &spectrum.spectrum_description {
        extend_ref_group_cv_params(out, &sd.referenceable_param_group_refs, ref_groups);
        out.extend_from_slice(&sd.cv_params);

        if let Some(sl) = &sd.scan_list {
            flatten_scan_list(out, sl, ref_groups);
        }
        if let Some(pl) = &sd.precursor_list {
            flatten_precursor_list(out, pl, ref_groups);
        }
        if let Some(pl) = &sd.product_list {
            flatten_product_list(out, pl, ref_groups);
        }
    }

    if let Some(sl) = &spectrum.scan_list {
        flatten_scan_list(out, sl, ref_groups);
    }
    if let Some(pl) = &spectrum.precursor_list {
        flatten_precursor_list(out, pl, ref_groups);
    }
    if let Some(pl) = &spectrum.product_list {
        flatten_product_list(out, pl, ref_groups);
    }

    if let Some(bal) = &spectrum.binary_data_array_list {
        for ba in &bal.binary_data_arrays {
            extend_ref_group_cv_params(out, &ba.referenceable_param_group_refs, ref_groups);
            extend_binary_data_array_cv_params(
                out,
                ba,
                x_accession_tail,
                y_accession_tail,
                x_store_f64,
                y_store_f64,
            );
        }
    }
}

/// <chromatogram>
fn flatten_chromatogram_metadata_into(
    out: &mut Vec<CvParam>,
    chrom: &Chromatogram,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
    x_accession_tail: u32,
    y_accession_tail: u32,
    x_store_f64: bool,
    y_store_f64: bool,
) {
    push_attr_string(out, ACC_ATTR_ID, chrom.id.as_str()); // ACC_ATTR_ID
    push_attr_u32(out, ACC_ATTR_INDEX, chrom.index); // ACC_ATTR_INDEX
    push_attr_usize(
        out,
        ACC_ATTR_DEFAULT_ARRAY_LENGTH,
        chrom.default_array_length,
    ); // ACC_ATTR_DEFAULT_ARRAY_LENGTH

    extend_ref_group_cv_params(out, &chrom.referenceable_param_group_refs, ref_groups);
    out.extend_from_slice(&chrom.cv_params);

    if let Some(p) = &chrom.precursor {
        flatten_precursor(out, p, ref_groups);
    }
    if let Some(p) = &chrom.product {
        flatten_product(out, p, ref_groups);
    }

    if let Some(bal) = &chrom.binary_data_array_list {
        for ba in &bal.binary_data_arrays {
            extend_ref_group_cv_params(out, &ba.referenceable_param_group_refs, ref_groups);
            extend_binary_data_array_cv_params(
                out,
                ba,
                x_accession_tail,
                y_accession_tail,
                x_store_f64,
                y_store_f64,
            );
        }
    }
}

/// <scanList>
fn flatten_scan_list(
    out: &mut Vec<CvParam>,
    scan_list: &ScanList,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for scan in &scan_list.scans {
        extend_ref_group_cv_params(out, &scan.referenceable_param_group_refs, ref_groups);
        out.extend_from_slice(&scan.cv_params);

        if let Some(wl) = &scan.scan_window_list {
            for w in &wl.scan_windows {
                out.extend_from_slice(&w.cv_params);
            }
        }
    }
}

/// <precursorList>
fn flatten_precursor_list(
    out: &mut Vec<CvParam>,
    precursor_list: &PrecursorList,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for p in &precursor_list.precursors {
        flatten_precursor(out, p, ref_groups);
    }
}

/// <precursor>
fn flatten_precursor(
    out: &mut Vec<CvParam>,
    precursor: &Precursor,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    if let Some(r) = precursor.spectrum_ref.as_deref() {
        push_attr_string(out, ACC_ATTR_SPECTRUM_REF, r); // ACC_ATTR_SPECTRUM_REF
    }

    if let Some(iw) = &precursor.isolation_window {
        extend_ref_group_cv_params(out, &iw.referenceable_param_group_refs, ref_groups);
        out.extend_from_slice(&iw.cv_params);
    }
    if let Some(sil) = &precursor.selected_ion_list {
        for ion in &sil.selected_ions {
            extend_ref_group_cv_params(out, &ion.referenceable_param_group_refs, ref_groups);
            out.extend_from_slice(&ion.cv_params);
        }
    }
    if let Some(act) = &precursor.activation {
        extend_ref_group_cv_params(out, &act.referenceable_param_group_refs, ref_groups);
        out.extend_from_slice(&act.cv_params);
    }
}

/// <productList>
fn flatten_product_list(
    out: &mut Vec<CvParam>,
    product_list: &ProductList,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    for p in &product_list.products {
        flatten_product(out, p, ref_groups);
    }
}

/// <product>
fn flatten_product(
    out: &mut Vec<CvParam>,
    product: &Product,
    ref_groups: &HashMap<&str, &ReferenceableParamGroup>,
) {
    if let Some(iw) = &product.isolation_window {
        extend_ref_group_cv_params(out, &iw.referenceable_param_group_refs, ref_groups);
        out.extend_from_slice(&iw.cv_params);
    }
}

/// <binaryDataArray>
fn extend_binary_data_array_cv_params(
    out: &mut Vec<CvParam>,
    ba: &BinaryDataArray,
    x_accession_tail: u32,
    y_accession_tail: u32,
    x_store_f64: bool,
    y_store_f64: bool,
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

    let desired_float_tail = if is_x {
        Some(if x_store_f64 {
            ACC_64BIT_FLOAT
        } else {
            ACC_32BIT_FLOAT
        })
    } else if is_y {
        Some(if y_store_f64 {
            ACC_64BIT_FLOAT
        } else {
            ACC_32BIT_FLOAT
        })
    } else {
        None
    };

    let mut wrote_float = false;

    for cv in &ba.cv_params {
        let tail = parse_accession_tail(cv.accession.as_deref());
        if tail == ACC_32BIT_FLOAT || tail == ACC_64BIT_FLOAT {
            if let Some(desired) = desired_float_tail {
                if !wrote_float {
                    out.push(ms_float_param(desired));
                    wrote_float = true;
                }
            } else {
                out.push(cv.clone());
            }
        } else {
            out.push(cv.clone());
        }
    }

    if let Some(desired) = desired_float_tail {
        if !wrote_float {
            out.push(ms_float_param(desired));
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

#[inline]
fn pack_cv_param(
    cv: &CvParam,
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
    ref_codes.push(opt_cv_ref_code(cv.cv_ref.as_deref()));
    accession_numbers.push(parse_accession_tail(cv.accession.as_deref()));

    unit_ref_codes.push(opt_cv_ref_code(cv.unit_cv_ref.as_deref()));
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

/// <cvParam>
fn pack_meta<T, F>(items: &[T], cv_params_of: F) -> PackedMeta
where
    F: Fn(&T) -> &[CvParam],
{
    let item_count = items.len();

    let mut total_meta_count = 0usize;
    for item in items {
        total_meta_count += cv_params_of(item).len();
    }

    let mut index_offsets = Vec::with_capacity(item_count + 1);
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
        for cv in cv_params_of(item) {
            pack_cv_param(
                cv,
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

fn pack_meta_streaming<T, F>(items: &[T], mut fill: F) -> PackedMeta
where
    F: FnMut(&mut Vec<CvParam>, &T),
{
    let item_count = items.len();

    let mut index_offsets = Vec::with_capacity(item_count + 1);
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

    let mut numeric_index: u32 = 0;
    let mut string_index: u32 = 0;
    let mut meta_index: u32 = 0;

    index_offsets.push(0);

    for item in items {
        scratch.clear();
        fill(&mut scratch, item);

        let n = scratch.len();
        ref_codes.reserve(n);
        accession_numbers.reserve(n);
        unit_ref_codes.reserve(n);
        unit_accession_numbers.reserve(n);
        value_kinds.reserve(n);
        value_indices.reserve(n);

        for cv in &scratch {
            pack_cv_param(
                cv,
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
    for &v in &meta.index_offsets {
        write_u32_le(buf, v);
    }
    buf.extend_from_slice(&meta.ref_codes);
    for &v in &meta.accession_numbers {
        write_u32_le(buf, v);
    }
    buf.extend_from_slice(&meta.unit_ref_codes);
    for &v in &meta.unit_accession_numbers {
        write_u32_le(buf, v);
    }
    buf.extend_from_slice(&meta.value_kinds);
    for &v in &meta.value_indices {
        write_u32_le(buf, v);
    }
    for &v in &meta.numeric_values {
        write_f64_le(buf, v);
    }
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
    let mut buf = Vec::with_capacity(8 * 4 + packed_meta_byte_len(meta));

    write_u32_le(&mut buf, counts.n_file_description);
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

/// <spectrum>
fn spectrum_xy<'a>(spectrum: &'a Spectrum) -> (ArrayInfo<'a>, ArrayInfo<'a>) {
    find_xy_in_bda_list(
        spectrum.binary_data_array_list.as_ref(),
        ACC_MZ_ARRAY,
        ACC_INTENSITY_ARRAY,
    )
}

/// <chromatogram>
fn chromatogram_xy<'a>(chrom: &'a Chromatogram) -> (ArrayInfo<'a>, ArrayInfo<'a>) {
    find_xy_in_bda_list(
        chrom.binary_data_array_list.as_ref(),
        ACC_TIME_ARRAY,
        ACC_INTENSITY_ARRAY,
    )
}

/// <binaryDataArrayList>
fn find_xy_in_bda_list<'a>(
    list: Option<&'a BinaryDataArrayList>,
    x_accession_tail: u32,
    y_accession_tail: u32,
) -> (ArrayInfo<'a>, ArrayInfo<'a>) {
    let mut x: Option<ArrayInfo<'a>> = None;
    let mut y: Option<ArrayInfo<'a>> = None;

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
                x = Some(array_ref(ba).unwrap_or(ArrayRef::F32(&[])));
            }
            if y.is_none() && is_y {
                y = Some(array_ref(ba).unwrap_or(ArrayRef::F32(&[])));
            }

            if x.is_some() && y.is_some() {
                break;
            }
        }
    }

    (
        x.unwrap_or(ArrayRef::F32(&[])),
        y.unwrap_or(ArrayRef::F32(&[])),
    )
}

/// <binaryDataArray>
#[inline]
fn array_ref<'a>(ba: &'a BinaryDataArray) -> Option<ArrayRef<'a>> {
    if !ba.decoded_binary_f64.is_empty() {
        Some(ArrayRef::F64(ba.decoded_binary_f64.as_slice()))
    } else if !ba.decoded_binary_f32.is_empty() {
        Some(ArrayRef::F32(ba.decoded_binary_f32.as_slice()))
    } else {
        None
    }
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

/// <cvParam>
#[inline]
fn opt_cv_ref_code(cv_ref: Option<&str>) -> u8 {
    match cv_ref {
        Some("MS") => 0,
        Some("UO") => 1,
        Some("NCIT") => 2,
        Some("PEFF") => 3,
        Some(CV_REF_ATTR) => 4,
        _ => 255,
    }
}
