use std::{fs, path::PathBuf};

use crate::b64::utilities::common::ChildIndex;
use crate::mzml::{schema::TagId, schema::schema};
use crate::{
    CvParam,
    b64::decode::Metadatum,
    b64::utilities::{parse_header, parse_metadata, parse_scan_list},
};

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

fn assert_cv_param(
    p: &CvParam,
    cv_ref: Option<&str>,
    accession: Option<&str>,
    name: &str,
    value: Option<&str>,
    unit_cv_ref: Option<&str>,
    unit_name: Option<&str>,
    unit_accession: Option<&str>,
) {
    assert_eq!(p.cv_ref.as_deref(), cv_ref);
    assert_eq!(p.accession.as_deref(), accession);
    assert_eq!(p.name.as_str(), name);
    assert_eq!(p.value.as_deref(), value);
    assert_eq!(p.unit_cv_ref.as_deref(), unit_cv_ref);
    assert_eq!(p.unit_name.as_deref(), unit_name);
    assert_eq!(p.unit_accession.as_deref(), unit_accession);
}

#[test]
fn first_spectrum_scan_list_cv_params_item_by_item() {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let meta = parse_metadata_section_from_test_file(
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        2,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
        header.spec_meta_count as usize,
        "spectra",
    );

    let scan_item_index = meta
        .iter()
        .find(|m| m.tag_id == TagId::Scan)
        .map(|m| m.item_index)
        .expect("no Scan entries found in spectra metadata");

    let scoped: Vec<Metadatum> = meta
        .into_iter()
        .filter(|m| m.item_index == scan_item_index)
        .collect();

    let child_index = ChildIndex::new(&scoped);

    let scan_list =
        parse_scan_list(schema(), &scoped, &child_index).expect("parse_scan_list returned None");
    assert_eq!(scan_list.count, Some(1));
    assert_eq!(scan_list.scans.len(), 1);

    let scan = &scan_list.scans[0];

    let p = scan
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000016"))
        .expect("missing MS:1000016 (scan start time)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000016"),
        "scan start time",
        Some("0.191"),
        Some("UO"),
        Some("second"),
        Some("UO:0000010"),
    );

    let swl = scan
        .scan_window_list
        .as_ref()
        .expect("missing scanWindowList");
    assert_eq!(swl.count, Some(1));
    assert_eq!(swl.scan_windows.len(), 1);

    let sw = &swl.scan_windows[0];

    let p = sw
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000501"))
        .expect("missing MS:1000501 (scan window lower limit)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000501"),
        "scan window lower limit",
        Some("30"),
        Some("MS"),
        Some("m/z"),
        Some("MS:1000040"),
    );

    let p = sw
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000500"))
        .expect("missing MS:1000500 (scan window upper limit)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000500"),
        "scan window upper limit",
        Some("1000"),
        Some("MS"),
        Some("m/z"),
        Some("MS:1000040"),
    );
}

#[test]
fn second_spectrum_scan_list_cv_params_item_by_item() {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let meta = parse_metadata_section_from_test_file(
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        2,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
        header.spec_meta_count as usize,
        "spectra",
    );

    let mut scan_item_indices: Vec<_> = meta
        .iter()
        .filter(|m| m.tag_id == TagId::Scan)
        .map(|m| m.item_index)
        .collect();

    scan_item_indices.sort_unstable();
    scan_item_indices.dedup();

    let scan_item_index = scan_item_indices
        .get(1)
        .copied()
        .expect("no second Scan item_index found in spectra metadata");

    let scoped: Vec<Metadatum> = meta
        .into_iter()
        .filter(|m| m.item_index == scan_item_index)
        .collect();

    let child_index = ChildIndex::new(&scoped);

    let scan_list =
        parse_scan_list(schema(), &scoped, &child_index).expect("parse_scan_list returned None");

    assert_eq!(scan_list.count, Some(1));
    assert_eq!(scan_list.scans.len(), 1);

    let scan = &scan_list.scans[0];

    let p = scan
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000016"))
        .expect("missing MS:1000016 (scan start time)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000016"),
        "scan start time",
        Some("452.262"),
        Some("UO"),
        Some("second"),
        Some("UO:0000010"),
    );

    let swl = scan
        .scan_window_list
        .as_ref()
        .expect("missing scanWindowList");
    assert_eq!(swl.count, Some(1));
    assert_eq!(swl.scan_windows.len(), 1);

    let sw = &swl.scan_windows[0];

    let p = sw
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000501"))
        .expect("missing MS:1000501 (scan window lower limit)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000501"),
        "scan window lower limit",
        Some("30"),
        Some("MS"),
        Some("m/z"),
        Some("MS:1000040"),
    );

    let p = sw
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000500"))
        .expect("missing MS:1000500 (scan window upper limit)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000500"),
        "scan window upper limit",
        Some("1000"),
        Some("MS"),
        Some("m/z"),
        Some("MS:1000040"),
    );
}
