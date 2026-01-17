use std::{
    collections::{HashMap, HashSet},
    fs,
    path::PathBuf,
};

use crate::{
    CvParam,
    b64::decode::Metadatum,
    b64::utilities::{parse_binary_data_array_list, parse_header, parse_metadata},
    mzml::schema::{TagId, schema},
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
fn first_chrom_cv_params_item_by_item() {
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

    let mut by_parent: HashMap<u32, HashSet<u32>> = HashMap::new();
    for m in &meta {
        if m.tag_id == TagId::BinaryDataArray {
            by_parent
                .entry(m.parent_index)
                .or_default()
                .insert(m.owner_id);
        }
    }

    let mut parent_ids: Vec<u32> = by_parent
        .iter()
        .filter(|(_, owners)| owners.len() == 3)
        .map(|(pid, _)| *pid)
        .collect();
    parent_ids.sort_unstable();
    let pid = parent_ids[0];

    let scoped: Vec<Metadatum> = meta
        .into_iter()
        .filter(|m| m.tag_id == TagId::BinaryDataArray && m.parent_index == pid)
        .collect();

    let bdal = parse_binary_data_array_list(&scoped).unwrap();
    assert_eq!(bdal.count, Some(3));
    assert_eq!(bdal.binary_data_arrays.len(), 3);

    for bda in &bdal.binary_data_arrays {
        let accs: HashSet<&str> = bda
            .cv_params
            .iter()
            .filter_map(|p| p.accession.as_deref())
            .collect();

        if accs.contains("MS:1000595") {
            let exp: HashSet<&str> = ["MS:1000523", "MS:1000576", "MS:1000595"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000523"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000595"))
                .unwrap();
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
        } else if accs.contains("MS:1000786") {
            let exp: HashSet<&str> = ["MS:1000522", "MS:1000576", "MS:1000786"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000522"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000786"))
                .unwrap();
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
        } else {
            let exp: HashSet<&str> = ["MS:1000521", "MS:1000576", "MS:1000515"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000521"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000515"))
                .unwrap();
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
        }
    }
}

#[test]
fn second_chrom_cv_params_item_by_item() {
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

    let mut by_parent: HashMap<u32, HashSet<u32>> = HashMap::new();
    for m in &meta {
        if m.tag_id == TagId::BinaryDataArray {
            by_parent
                .entry(m.parent_index)
                .or_default()
                .insert(m.owner_id);
        }
    }

    let mut parent_ids: Vec<u32> = by_parent
        .iter()
        .filter(|(_, owners)| owners.len() == 3)
        .map(|(pid, _)| *pid)
        .collect();
    parent_ids.sort_unstable();
    let pid = parent_ids[parent_ids.len() - 1];

    let scoped: Vec<Metadatum> = meta
        .into_iter()
        .filter(|m| m.tag_id == TagId::BinaryDataArray && m.parent_index == pid)
        .collect();

    let bdal = parse_binary_data_array_list(&scoped).unwrap();
    assert_eq!(bdal.count, Some(3));
    assert_eq!(bdal.binary_data_arrays.len(), 3);

    for bda in &bdal.binary_data_arrays {
        let accs: HashSet<&str> = bda
            .cv_params
            .iter()
            .filter_map(|p| p.accession.as_deref())
            .collect();

        if accs.contains("MS:1000595") {
            let exp: HashSet<&str> = ["MS:1000523", "MS:1000576", "MS:1000595"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000523"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000595"))
                .unwrap();
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
        } else if accs.contains("MS:1000786") {
            let exp: HashSet<&str> = ["MS:1000522", "MS:1000576", "MS:1000786"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000522"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000786"))
                .unwrap();
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
        } else {
            let exp: HashSet<&str> = ["MS:1000521", "MS:1000576", "MS:1000515"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000521"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000515"))
                .unwrap();
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
        }
    }
}

#[test]
fn first_spectrum_cv_params_item_by_item() {
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

    let mut by_parent: HashMap<u32, HashSet<u32>> = HashMap::new();
    for m in &meta {
        if m.tag_id == TagId::BinaryDataArray {
            by_parent
                .entry(m.parent_index)
                .or_default()
                .insert(m.owner_id);
        }
    }

    let mut parent_ids: Vec<u32> = by_parent
        .iter()
        .filter(|(_, owners)| owners.len() == 2)
        .map(|(pid, _)| *pid)
        .collect();
    parent_ids.sort_unstable();
    let pid = parent_ids[0];

    let scoped: Vec<Metadatum> = meta
        .into_iter()
        .filter(|m| m.tag_id == TagId::BinaryDataArray && m.parent_index == pid)
        .collect();

    let bdal = parse_binary_data_array_list(&scoped).unwrap();
    assert_eq!(bdal.count, Some(2));
    assert_eq!(bdal.binary_data_arrays.len(), 2);

    for bda in &bdal.binary_data_arrays {
        let accs: HashSet<&str> = bda
            .cv_params
            .iter()
            .filter_map(|p| p.accession.as_deref())
            .collect();

        if accs.contains("MS:1000514") {
            let exp: HashSet<&str> = ["MS:1000523", "MS:1000576", "MS:1000514"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000523"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000514"))
                .unwrap();
            assert_cv_param(
                p,
                Some("MS"),
                Some("MS:1000514"),
                "m/z array",
                None,
                Some("MS"),
                Some("m/z"),
                Some("MS:1000040"),
            );
        } else {
            let exp: HashSet<&str> = ["MS:1000521", "MS:1000576", "MS:1000515"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000521"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000515"))
                .unwrap();
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
        }
    }
}

#[test]
fn second_spectrum_cv_params_item_by_item() {
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

    let mut by_parent: HashMap<u32, HashSet<u32>> = HashMap::new();
    for m in &meta {
        if m.tag_id == TagId::BinaryDataArray {
            by_parent
                .entry(m.parent_index)
                .or_default()
                .insert(m.owner_id);
        }
    }

    let mut parent_ids: Vec<u32> = by_parent
        .iter()
        .filter(|(_, owners)| owners.len() == 2)
        .map(|(pid, _)| *pid)
        .collect();
    parent_ids.sort_unstable();
    let pid = parent_ids[parent_ids.len() - 1];

    let scoped: Vec<Metadatum> = meta
        .into_iter()
        .filter(|m| m.tag_id == TagId::BinaryDataArray && m.parent_index == pid)
        .collect();

    let bdal = parse_binary_data_array_list(&scoped).unwrap();
    assert_eq!(bdal.count, Some(2));
    assert_eq!(bdal.binary_data_arrays.len(), 2);

    for bda in &bdal.binary_data_arrays {
        let accs: HashSet<&str> = bda
            .cv_params
            .iter()
            .filter_map(|p| p.accession.as_deref())
            .collect();

        if accs.contains("MS:1000514") {
            let exp: HashSet<&str> = ["MS:1000523", "MS:1000576", "MS:1000514"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000523"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000514"))
                .unwrap();
            assert_cv_param(
                p,
                Some("MS"),
                Some("MS:1000514"),
                "m/z array",
                None,
                Some("MS"),
                Some("m/z"),
                Some("MS:1000040"),
            );
        } else {
            let exp: HashSet<&str> = ["MS:1000521", "MS:1000576", "MS:1000515"]
                .into_iter()
                .collect();
            assert_eq!(accs, exp);

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000521"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000576"))
                .unwrap();
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

            let p = bda
                .cv_params
                .iter()
                .find(|p| p.accession.as_deref() == Some("MS:1000515"))
                .unwrap();
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
        }
    }
}
