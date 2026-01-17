use std::{fs, path::PathBuf};

use crate::b64::decode::Metadatum;
use crate::b64::utilities::common::ChildIndex;
use crate::b64::utilities::{parse_chromatogram_list, parse_header, parse_metadata};
use crate::mzml::schema::schema;
use crate::{ChromatogramList, CvParam};

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

fn parse_chromatogram_list_from_test_file() -> ChromatogramList {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let meta = parse_metadata_section_from_test_file(
        header.off_chrom_meta,
        header.off_global_meta,
        header.chrom_count,
        2,
        header.chrom_meta_count,
        header.chrom_num_count,
        header.chrom_str_count,
        5,
        header.chrom_meta_count as usize,
        "chromatograms",
    );

    let child_index = ChildIndex::new(&meta);

    parse_chromatogram_list(schema(), &meta, &child_index)
        .expect("parse_chromatogram_list returned None")
}

#[test]
fn chromatogram_list_emits_count_and_default_dp_ref() {
    let cl = parse_chromatogram_list_from_test_file();

    assert_eq!(cl.chromatograms.len(), 2);
    assert_eq!(cl.count, Some(2));
    assert_eq!(
        cl.default_data_processing_ref.as_deref(),
        Some("pwiz_Reader_Bruker_conversion")
    );
}

#[test]
fn chromatogram_0_tic_cv_params_and_bdal_shape() {
    let cl = parse_chromatogram_list_from_test_file();
    let c = &cl.chromatograms[0];

    assert_eq!(c.id.as_str(), "TIC");
    assert_eq!(c.index, Some(0));
    assert_eq!(c.default_array_length, Some(3476));

    // <cvParam accession="MS:1000235" name="total ion current chromatogram" value=""/>
    let p = c
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000235"))
        .expect("missing MS:1000235 (total ion current chromatogram)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000235"),
        "total ion current chromatogram",
        None,
        None,
        None,
        None,
    );

    let bal = c
        .binary_data_array_list
        .as_ref()
        .expect("missing <binaryDataArrayList> on TIC");
    assert_eq!(bal.count, Some(3));
    assert_eq!(bal.binary_data_arrays.len(), 3);

    let ba0 = &bal.binary_data_arrays[0];
    assert_eq!(ba0.encoded_length, Some(37080));
    let p = ba0
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000523"))
        .expect("missing MS:1000523 (64-bit float) on time array");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000523"),
        "64-bit float",
        None,
        None,
        None,
        None,
    );

    let p = ba0
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000576"))
        .expect("missing MS:1000576 (no compression) on time array");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000576"),
        "no compression",
        None,
        None,
        None,
        None,
    );

    let p = ba0
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000595"))
        .expect("missing MS:1000595 (time array)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000595"),
        "time array",
        None,
        Some("UO"),
        Some("second"),
        Some("UO:0000010"),
    );

    let ba1 = &bal.binary_data_arrays[1];
    assert_eq!(ba1.encoded_length, Some(18540));

    let p = ba1
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000521"))
        .expect("missing MS:1000521 (32-bit float) on intensity array");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000521"),
        "32-bit float",
        None,
        None,
        None,
        None,
    );

    let p = ba1
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000576"))
        .expect("missing MS:1000576 (no compression) on intensity array");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000576"),
        "no compression",
        None,
        None,
        None,
        None,
    );

    let p = ba1
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000515"))
        .expect("missing MS:1000515 (intensity array)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000515"),
        "intensity array",
        None,
        Some("MS"),
        Some("number of detector counts"),
        Some("MS:1000131"),
    );

    let ba2 = &bal.binary_data_arrays[2];
    assert_eq!(ba2.array_length, Some(3476));
    assert_eq!(ba2.encoded_length, Some(37080));

    let p = ba2
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000522"))
        .expect("missing MS:1000522 (64-bit integer) on non-standard array");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000522"),
        "64-bit integer",
        None,
        None,
        None,
        None,
    );

    let p = ba2
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000576"))
        .expect("missing MS:1000576 (no compression) on non-standard array");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000576"),
        "no compression",
        None,
        None,
        None,
        None,
    );

    let p = ba2
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000786"))
        .expect("missing MS:1000786 (non-standard data array)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000786"),
        "non-standard data array",
        Some("ms level"),
        Some("UO"),
        Some("dimensionless unit"),
        Some("UO:0000186"),
    );
}

#[test]
fn chromatogram_1_bpc_has_expected_cv_and_bdal_lengths_match() {
    let cl = parse_chromatogram_list_from_test_file();
    let c = &cl.chromatograms[1];

    assert_eq!(c.id.as_str(), "BPC");
    assert_eq!(c.index, Some(1));
    assert_eq!(c.default_array_length, Some(3476));

    // <cvParam accession="MS:1000628" name="basepeak chromatogram" value=""/>
    let p = c
        .cv_params
        .iter()
        .find(|p| p.accession.as_deref() == Some("MS:1000628"))
        .expect("missing MS:1000628 (basepeak chromatogram)");
    assert_cv_param(
        p,
        Some("MS"),
        Some("MS:1000628"),
        "basepeak chromatogram",
        None,
        None,
        None,
        None,
    );

    let bal = c
        .binary_data_array_list
        .as_ref()
        .expect("missing <binaryDataArrayList> on BPC");
    assert_eq!(bal.count, Some(3));
    assert_eq!(bal.binary_data_arrays.len(), 3);

    assert_eq!(bal.binary_data_arrays[0].encoded_length, Some(37080));
    assert_eq!(bal.binary_data_arrays[1].encoded_length, Some(18540));
    assert_eq!(bal.binary_data_arrays[2].array_length, Some(3476));
    assert_eq!(bal.binary_data_arrays[2].encoded_length, Some(37080));
}
