// crates/parser/src/b64/utilities/tests/parse_file_description.rs

use std::{fs, path::PathBuf};

use serde_json::Value;

use crate::b64::{
    decode::Metadatum,
    utilities::{
        common::ChildIndex, parse_file_description::parse_file_description,
        parse_global_metadata::parse_global_metadata, parse_header,
    },
};
use crate::mzml::schema::schema;

const PATH: &str = "data/b64/test.b64";
const GLOBAL_META_COMPRESSION_BIT: u8 = 6;

fn read_bytes(path: &str) -> Vec<u8> {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(path);
    fs::read(&full).unwrap_or_else(|e| panic!("cannot read {:?}: {}", full, e))
}

fn read_u32_at(bytes: &[u8], off: usize, little: bool) -> u32 {
    let b = [bytes[off], bytes[off + 1], bytes[off + 2], bytes[off + 3]];
    if little {
        u32::from_le_bytes(b)
    } else {
        u32::from_be_bytes(b)
    }
}

fn parse_global_metadata_from_test_file() -> Vec<Metadatum> {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let start = header.off_global_meta as usize;
    let end = header.off_container_spect_x as usize;

    assert!(start < end, "invalid global metadata offsets: start >= end");
    assert!(
        end <= bytes.len(),
        "invalid global metadata offsets: end out of bounds"
    );

    let slice = &bytes[start..end];
    let compressed = (header.reserved_flags & (1u8 << GLOBAL_META_COMPRESSION_BIT)) != 0;

    let meta = parse_global_metadata(
        slice,
        0, // ignored/derived from the global header inside the section
        header.global_meta_count,
        header.global_num_count,
        header.global_str_count,
        compressed,
        header.reserved_flags,
    )
    .expect("parse_global_metadata failed");

    meta
}

fn must_obj<'a>(v: &'a Value, key: &str) -> &'a serde_json::Map<String, Value> {
    v.get(key)
        .unwrap_or_else(|| panic!("missing key {key:?}"))
        .as_object()
        .unwrap_or_else(|| panic!("key {key:?} is not an object"))
}

fn must_arr<'a>(v: &'a Value, key: &str) -> &'a Vec<Value> {
    v.get(key)
        .unwrap_or_else(|| panic!("missing key {key:?}"))
        .as_array()
        .unwrap_or_else(|| panic!("key {key:?} is not an array"))
}

fn must_str<'a>(v: &'a Value, key: &str) -> &'a str {
    v.get(key)
        .unwrap_or_else(|| panic!("missing key {key:?}"))
        .as_str()
        .unwrap_or_else(|| panic!("key {key:?} is not a string"))
}

fn find_cv_by_accession<'a>(arr: &'a [Value], accession: &str) -> &'a Value {
    for v in arr {
        if v.get("accession").and_then(|x| x.as_str()) == Some(accession) {
            return v;
        }
    }
    panic!("missing cvParam with accession {accession:?}");
}

#[test]
fn parse_file_description_matches_test_mzml() {
    let meta = parse_global_metadata_from_test_file();

    let child_index = ChildIndex::new(&meta);
    let schema = schema();

    let fd = parse_file_description(&schema, &meta, &child_index)
        .expect("parse_file_description returned None");

    assert_eq!(fd.source_file_list.count, Some(1));
    assert_eq!(fd.source_file_list.source_file.len(), 1);

    let sf = &fd.source_file_list.source_file[0];
    assert_eq!(sf.id, "anpc_file.d_x005c_Analysis.baf");
    assert_eq!(sf.name, "Analysis.baf");
    assert_eq!(sf.location, r"file://Z:\inputDirectory\anpc_file.d");

    let json = serde_json::to_value(&fd).expect("serialize FileDescription failed");

    let file_content_map = must_obj(&json, "file_content");
    let file_content_val = Value::Object(file_content_map.clone());
    let fc_cv = must_arr(&file_content_val, "cv_params");

    find_cv_by_accession(fc_cv, "MS:1000579");
    find_cv_by_accession(fc_cv, "MS:1000580");

    let sfl_map = must_obj(&json, "source_file_list");
    let sfl_val = Value::Object(sfl_map.clone());
    let sfl_files = must_arr(&sfl_val, "source_file");
    assert_eq!(sfl_files.len(), 1);

    let sf0 = sfl_files[0]
        .as_object()
        .unwrap_or_else(|| panic!("source_file[0] is not an object"));
    let sf_cv = sf0
        .get("cv_param")
        .and_then(|v| v.as_array())
        .unwrap_or_else(|| panic!("source_file[0].cv_param is not an array"));

    find_cv_by_accession(sf_cv, "MS:1000772");
    find_cv_by_accession(sf_cv, "MS:1000815");

    let sha = find_cv_by_accession(sf_cv, "MS:1000569");
    let sha_val = sha
        .get("value")
        .and_then(|v| v.as_str())
        .unwrap_or_default();
    assert_eq!(sha_val, "36a9346b9d32b3ef5b30e48d1a20cf1515232083");
}
