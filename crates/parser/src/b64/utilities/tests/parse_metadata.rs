use std::{fs, path::PathBuf};

use crate::b64::{
    decode::{Metadatum, MetadatumValue},
    utilities::{parse_header, parse_metadata},
};
use crate::mzml::{attr_meta::*, schema::TagId};

const PATH: &str = "data/b64/test.b64";

fn read_bytes(path: &str) -> Vec<u8> {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(path);
    fs::read(&full).unwrap_or_else(|e| panic!("cannot read {:?}: {}", full, e))
}

fn parse_metadata_section_from_test_file(
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
    let bytes = read_bytes(PATH);
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

#[test]
fn check_first_spectrum() {
    let header_bytes = read_bytes(PATH);
    let header = parse_header(&header_bytes).expect("parse_header failed");
    let spec_meta = parse_metadata_section_from_test_file(
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        2,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
        54,
        "spectra",
    );

    assert_eq!(item_meta_count(&spec_meta, 0), 24);
    expect_text_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_B000,
        ACC_ATTR_ID,
        CV_CODE_UNKNOWN,
        0,
        "scan=1",
    );
    expect_number_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_B000,
        ACC_ATTR_INDEX,
        CV_CODE_UNKNOWN,
        0,
        0.0,
    );
    expect_number_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_B000,
        ACC_ATTR_DEFAULT_ARRAY_LENGTH,
        CV_CODE_UNKNOWN,
        0,
        340032.0,
    );

    expect_number_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_MS,
        1000511,
        CV_CODE_UNKNOWN,
        0,
        1.0,
    );
    expect_empty_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_MS,
        1000579,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_MS,
        1000130,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_number_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_MS,
        1000505,
        CV_CODE_UNKNOWN,
        0,
        24998.0,
    );
    expect_number_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_MS,
        1000285,
        CV_CODE_UNKNOWN,
        0,
        440132.0,
    );
    expect_empty_one(
        &spec_meta,
        0,
        TagId::Spectrum,
        CV_CODE_MS,
        1000128,
        CV_CODE_UNKNOWN,
        0,
    );

    // ScanList
    expect_number_one(
        &spec_meta,
        0,
        TagId::ScanList,
        CV_CODE_B000,
        ACC_ATTR_COUNT,
        CV_CODE_UNKNOWN,
        0,
        1.0,
    );

    // Scan
    expect_number_one(
        &spec_meta,
        0,
        TagId::Scan,
        CV_CODE_MS,
        1000016,
        CV_CODE_UO,
        10,
        0.191,
    );

    // ScanWindow
    expect_number_one(
        &spec_meta,
        0,
        TagId::ScanWindow,
        CV_CODE_MS,
        1000501,
        CV_CODE_MS,
        1000040,
        30.0,
    );
    expect_number_one(
        &spec_meta,
        0,
        TagId::ScanWindow,
        CV_CODE_MS,
        1000500,
        CV_CODE_MS,
        1000040,
        1000.0,
    );

    // BinaryDataArrayList
    expect_number_one(
        &spec_meta,
        0,
        TagId::BinaryDataArrayList,
        CV_CODE_B000,
        ACC_ATTR_COUNT,
        CV_CODE_UNKNOWN,
        0,
        2.0,
    );

    // BinaryDataArray
    expect_empty_one(
        &spec_meta,
        0,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000523,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_count(
        &spec_meta,
        0,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000576,
        CV_CODE_UNKNOWN,
        0,
        2,
    );
    expect_empty_one(
        &spec_meta,
        0,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000514,
        CV_CODE_MS,
        1000040,
    );
    expect_empty_one(
        &spec_meta,
        0,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000521,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_one(
        &spec_meta,
        0,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000515,
        CV_CODE_MS,
        1000131,
    );
}

#[test]
fn check_second_spectrum() {
    let header_bytes = read_bytes(PATH);
    let header = parse_header(&header_bytes).expect("parse_header failed");
    let spec_meta = parse_metadata_section_from_test_file(
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        2,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
        54,
        "spectra",
    );

    assert_eq!(item_meta_count(&spec_meta, 1), 30);

    expect_text_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_B000,
        ACC_ATTR_ID,
        CV_CODE_UNKNOWN,
        0,
        "scan=3476",
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_B000,
        ACC_ATTR_INDEX,
        CV_CODE_UNKNOWN,
        0,
        3475.0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_B000,
        ACC_ATTR_DEFAULT_ARRAY_LENGTH,
        CV_CODE_UNKNOWN,
        0,
        4340.0,
    );

    expect_number_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_MS,
        1000511,
        CV_CODE_UNKNOWN,
        0,
        2.0,
    );
    expect_empty_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_MS,
        1000580,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_MS,
        1000130,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_MS,
        1000505,
        CV_CODE_UNKNOWN,
        0,
        20032.0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_MS,
        1000285,
        CV_CODE_UNKNOWN,
        0,
        359026.0,
    );
    expect_empty_one(
        &spec_meta,
        1,
        TagId::Spectrum,
        CV_CODE_MS,
        1000127,
        CV_CODE_UNKNOWN,
        0,
    );

    expect_number_one(
        &spec_meta,
        1,
        TagId::Scan,
        CV_CODE_MS,
        1000016,
        CV_CODE_UO,
        10,
        452.262,
    );

    // ScanWindow
    expect_number_one(
        &spec_meta,
        1,
        TagId::ScanWindow,
        CV_CODE_MS,
        1000501,
        CV_CODE_MS,
        1000040,
        30.0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::ScanWindow,
        CV_CODE_MS,
        1000500,
        CV_CODE_MS,
        1000040,
        1000.0,
    );

    // PrecursorList
    expect_number_one(
        &spec_meta,
        1,
        TagId::PrecursorList,
        CV_CODE_B000,
        ACC_ATTR_COUNT,
        CV_CODE_UNKNOWN,
        0,
        1.0,
    );
    // IsolationWindow
    expect_number_one(
        &spec_meta,
        1,
        TagId::IsolationWindow,
        CV_CODE_MS,
        1000827,
        CV_CODE_MS,
        1000040,
        515.0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::IsolationWindow,
        CV_CODE_MS,
        1000828,
        CV_CODE_MS,
        1000040,
        485.0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::IsolationWindow,
        CV_CODE_MS,
        1000829,
        CV_CODE_MS,
        1000040,
        485.0,
    );
    // SelectedIonList
    expect_number_one(
        &spec_meta,
        1,
        TagId::SelectedIonList,
        CV_CODE_B000,
        ACC_ATTR_COUNT,
        CV_CODE_UNKNOWN,
        0,
        1.0,
    );
    // SelectedIon
    expect_number_one(
        &spec_meta,
        1,
        TagId::SelectedIon,
        CV_CODE_MS,
        1000744,
        CV_CODE_MS,
        1000040,
        515.0,
    );

    // Activation
    expect_empty_one(
        &spec_meta,
        1,
        TagId::Activation,
        CV_CODE_MS,
        1001880,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_number_one(
        &spec_meta,
        1,
        TagId::Activation,
        CV_CODE_MS,
        1000045,
        CV_CODE_UNKNOWN,
        0,
        20.0,
    );

    // BinaryDataArray
    expect_empty_one(
        &spec_meta,
        1,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000523,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_count(
        &spec_meta,
        1,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000576,
        CV_CODE_UNKNOWN,
        0,
        2,
    );
    expect_empty_one(
        &spec_meta,
        1,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000514,
        CV_CODE_MS,
        1000040,
    );
    expect_empty_one(
        &spec_meta,
        1,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000521,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_one(
        &spec_meta,
        1,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000515,
        CV_CODE_MS,
        1000131,
    );
}

fn assert_chromatogram_binary_data_array_list(meta: &[Metadatum], item_index: u32) {
    // BinaryDataArray
    expect_empty_one(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000523,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_count(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000576,
        CV_CODE_UNKNOWN,
        0,
        3,
    );
    expect_empty_one(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000595,
        CV_CODE_UO,
        10,
    );

    expect_empty_one(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000521,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_empty_one(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000515,
        CV_CODE_MS,
        1000131,
    );

    expect_empty_one(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000522,
        CV_CODE_UNKNOWN,
        0,
    );
    expect_text_one(
        meta,
        item_index,
        TagId::BinaryDataArray,
        CV_CODE_MS,
        1000786,
        CV_CODE_UO,
        186,
        "ms level",
    );
}

#[test]
fn check_first_chromatogram() {
    let header_bytes = read_bytes(PATH);
    let header = parse_header(&header_bytes).expect("parse_header failed");
    let chrom_meta = parse_metadata_section_from_test_file(
        header.off_chrom_meta,
        header.off_global_meta,
        header.chrom_count,
        2,
        header.chrom_meta_count,
        header.chrom_num_count,
        header.chrom_str_count,
        5,
        38,
        "chromatograms",
    );
    assert_eq!(chrom_meta.len(), 38);
    assert_eq!(item_meta_count(&chrom_meta, 0), 20);

    // Chromatogram
    expect_text_one(
        &chrom_meta,
        0,
        TagId::Chromatogram,
        CV_CODE_B000,
        ACC_ATTR_ID,
        CV_CODE_UNKNOWN,
        0,
        "TIC",
    );
    expect_number_one(
        &chrom_meta,
        0,
        TagId::Chromatogram,
        CV_CODE_B000,
        ACC_ATTR_INDEX,
        CV_CODE_UNKNOWN,
        0,
        0.0,
    );
    expect_number_one(
        &chrom_meta,
        0,
        TagId::Chromatogram,
        CV_CODE_B000,
        ACC_ATTR_DEFAULT_ARRAY_LENGTH,
        CV_CODE_UNKNOWN,
        0,
        3476.0,
    );

    expect_empty_one(
        &chrom_meta,
        0,
        TagId::Chromatogram,
        CV_CODE_MS,
        1000235,
        CV_CODE_UNKNOWN,
        0,
    );

    assert_chromatogram_binary_data_array_list(&chrom_meta, 0);
}

#[test]
fn check_second_chromatogram() {
    let header_bytes = read_bytes(PATH);
    let header = parse_header(&header_bytes).expect("parse_header failed");
    let chrom_meta = parse_metadata_section_from_test_file(
        header.off_chrom_meta,
        header.off_global_meta,
        header.chrom_count,
        2,
        header.chrom_meta_count,
        header.chrom_num_count,
        header.chrom_str_count,
        5,
        38,
        "chromatograms",
    );
    assert_eq!(chrom_meta.len(), 38);
    assert_eq!(item_meta_count(&chrom_meta, 1), 18);

    // Chromatogram
    expect_text_one(
        &chrom_meta,
        1,
        TagId::Chromatogram,
        CV_CODE_B000,
        ACC_ATTR_ID,
        CV_CODE_UNKNOWN,
        0,
        "BPC",
    );
    expect_number_one(
        &chrom_meta,
        1,
        TagId::Chromatogram,
        CV_CODE_B000,
        ACC_ATTR_INDEX,
        CV_CODE_UNKNOWN,
        0,
        1.0,
    );
    expect_number_one(
        &chrom_meta,
        1,
        TagId::Chromatogram,
        CV_CODE_B000,
        ACC_ATTR_DEFAULT_ARRAY_LENGTH,
        CV_CODE_UNKNOWN,
        0,
        3476.0,
    );

    expect_empty_one(
        &chrom_meta,
        1,
        TagId::Chromatogram,
        CV_CODE_MS,
        1000628,
        CV_CODE_UNKNOWN,
        0,
    );

    assert_chromatogram_binary_data_array_list(&chrom_meta, 1);
}

fn item_meta_count(meta: &[Metadatum], item_index: u32) -> usize {
    let mut n = 0;
    for m in meta {
        if m.item_index == item_index {
            n += 1;
        }
    }
    n
}

fn find_meta_all<'a>(
    meta: &'a [Metadatum],
    item_index: u32,
    tag_id: TagId,
    ref_id: u8,
    accession_tail: u32,
) -> Vec<&'a Metadatum> {
    let expected_accession = format_accession(ref_id, accession_tail);

    let mut out = Vec::new();
    for m in meta {
        if m.item_index == item_index
            && m.tag_id == tag_id
            && m.accession.as_deref() == expected_accession.as_deref()
        {
            out.push(m);
        }
    }
    out
}

fn find_meta_one<'a>(
    meta: &'a [Metadatum],
    item_index: u32,
    tag_id: TagId,
    ref_id: u8,
    accession_tail: u32,
) -> &'a Metadatum {
    let hits = find_meta_all(meta, item_index, tag_id, ref_id, accession_tail);
    if hits.len() != 1 {
        let expected_accession = format_accession(ref_id, accession_tail);
        panic!(
            "expected exactly 1 metadatum, found {}: item_index={}, tag_id={:?}, accession={:?}",
            hits.len(),
            item_index,
            tag_id,
            expected_accession
        );
    }
    hits[0]
}

fn expect_text_one(
    meta: &[Metadatum],
    item_index: u32,
    tag_id: TagId,
    ref_id: u8,
    accession_tail: u32,
    unit_ref_id: u8,
    unit_accession_tail: u32,
    expected: &str,
) {
    let m = find_meta_one(meta, item_index, tag_id, ref_id, accession_tail);

    let expected_unit_accession = format_accession(unit_ref_id, unit_accession_tail);
    assert_eq!(
        m.unit_accession.as_deref(),
        expected_unit_accession.as_deref()
    );

    match &m.value {
        MetadatumValue::Text(s) => assert_eq!(s.as_str(), expected),
        other => panic!("expected Text({expected:?}), got {other:?} for {m:?}"),
    }
}

fn expect_number_one(
    meta: &[Metadatum],
    item_index: u32,
    tag_id: TagId,
    ref_id: u8,
    accession_tail: u32,
    unit_ref_id: u8,
    unit_accession_tail: u32,
    expected: f64,
) {
    let m = find_meta_one(meta, item_index, tag_id, ref_id, accession_tail);

    let expected_unit_accession = format_accession(unit_ref_id, unit_accession_tail);
    assert_eq!(
        m.unit_accession.as_deref(),
        expected_unit_accession.as_deref()
    );

    match &m.value {
        MetadatumValue::Number(v) => {
            let tol = 1e-9_f64.max(expected.abs() * 1e-9);
            assert!(
                (v - expected).abs() <= tol,
                "expected Number({expected}), got Number({v}) for {m:?}"
            );
        }
        other => panic!("expected Number({expected}), got {other:?} for {m:?}"),
    }
}

fn expect_empty_one(
    meta: &[Metadatum],
    item_index: u32,
    tag_id: TagId,
    ref_id: u8,
    accession_tail: u32,
    unit_ref_id: u8,
    unit_accession_tail: u32,
) {
    let m = find_meta_one(meta, item_index, tag_id, ref_id, accession_tail);

    let expected_unit_accession = format_accession(unit_ref_id, unit_accession_tail);
    assert_eq!(
        m.unit_accession.as_deref(),
        expected_unit_accession.as_deref()
    );

    match &m.value {
        MetadatumValue::Empty => {}
        other => panic!("expected Empty, got {other:?} for {m:?}"),
    }
}

fn expect_empty_count(
    meta: &[Metadatum],
    item_index: u32,
    tag_id: TagId,
    ref_id: u8,
    accession_tail: u32,
    unit_ref_id: u8,
    unit_accession_tail: u32,
    expected_count: usize,
) {
    let hits = find_meta_all(meta, item_index, tag_id, ref_id, accession_tail);

    assert_eq!(
        hits.len(),
        expected_count,
        "unexpected metadatum count for item_index={}, tag_id={:?}, accession={:?}",
        item_index,
        tag_id,
        format_accession(ref_id, accession_tail)
    );

    let expected_unit_accession = format_accession(unit_ref_id, unit_accession_tail);
    for m in hits {
        assert_eq!(
            m.unit_accession.as_deref(),
            expected_unit_accession.as_deref()
        );
        match &m.value {
            MetadatumValue::Empty => {}
            other => panic!("expected Empty, got {other:?} for {m:?}"),
        }
    }
}
