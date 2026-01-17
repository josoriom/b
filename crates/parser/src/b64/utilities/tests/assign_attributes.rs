use crate::b64::decode::{Metadatum, MetadatumValue};
use crate::b64::utilities::assign_attributes;
use crate::mzml::attr_meta::*;
use crate::mzml::schema::TagId;
use crate::mzml::structs::{Chromatogram, Spectrum};

fn b000_accession(tail: u32) -> String {
    format!("{CV_REF_ATTR}:{tail:07}")
}

fn find_by_tail<'a>(xs: &'a [Metadatum], tail: u32) -> Option<&'a Metadatum> {
    let acc = b000_accession(tail);
    xs.iter()
        .find(|m| m.accession.as_deref() == Some(acc.as_str()))
}

fn assert_has_b000_tail(xs: &[Metadatum], tail: u32) {
    let acc = b000_accession(tail);
    assert!(
        xs.iter()
            .any(|m| m.accession.as_deref() == Some(acc.as_str())),
        "expected generated meta to contain {acc}"
    );
}

fn assert_missing_b000_tail(xs: &[Metadatum], tail: u32) {
    let acc = b000_accession(tail);
    assert!(
        xs.iter()
            .all(|m| m.accession.as_deref() != Some(acc.as_str())),
        "expected generated meta to NOT contain {acc}"
    );
}

fn assert_text(xs: &[Metadatum], tail: u32, s: &str) {
    let m = find_by_tail(xs, tail).unwrap_or_else(|| panic!("missing tail={tail}"));
    assert_eq!(m.value, MetadatumValue::Text(s.to_string()));
}

fn assert_num(xs: &[Metadatum], tail: u32, n: f64) {
    let m = find_by_tail(xs, tail).unwrap_or_else(|| panic!("missing tail={tail}"));
    assert_eq!(m.value, MetadatumValue::Number(n));
}

#[test]
fn assign_attrs_spectrum_emits_schema_b000_only() {
    let mut s = Spectrum::default();

    s.id = "scan=1".to_string();
    s.index = Some(0);
    s.default_array_length = Some(340032);
    s.native_id = Some("controllerType=0 controllerNumber=1 scan=1".to_string());
    s.data_processing_ref = Some("dp1".to_string());
    s.source_file_ref = Some("sf1".to_string());
    s.spot_id = Some("spotA".to_string());
    s.ms_level = Some(2);

    let owner_id = 1u32;
    let parent_index = 0u32;

    let out = assign_attributes(&s, TagId::Spectrum, owner_id, parent_index);

    assert_has_b000_tail(&out, ACC_ATTR_ID);
    assert_has_b000_tail(&out, ACC_ATTR_INDEX);
    assert_has_b000_tail(&out, ACC_ATTR_DEFAULT_ARRAY_LENGTH);
    assert_has_b000_tail(&out, ACC_ATTR_NATIVE_ID);
    assert_has_b000_tail(&out, ACC_ATTR_DATA_PROCESSING_REF);
    assert_has_b000_tail(&out, ACC_ATTR_SOURCE_FILE_REF);
    assert_has_b000_tail(&out, ACC_ATTR_SPOT_ID);

    assert_text(&out, ACC_ATTR_ID, "scan=1");
    assert_num(&out, ACC_ATTR_INDEX, 0.0);
    assert_num(&out, ACC_ATTR_DEFAULT_ARRAY_LENGTH, 340032.0);
    assert_text(
        &out,
        ACC_ATTR_NATIVE_ID,
        "controllerType=0 controllerNumber=1 scan=1",
    );
    assert_text(&out, ACC_ATTR_DATA_PROCESSING_REF, "dp1");
    assert_text(&out, ACC_ATTR_SOURCE_FILE_REF, "sf1");
    assert_text(&out, ACC_ATTR_SPOT_ID, "spotA");

    assert_missing_b000_tail(&out, ACC_ATTR_MS_LEVEL);
    assert_missing_b000_tail(&out, ACC_ATTR_SCAN_NUMBER);
}

#[test]
fn assign_attrs_emits_only_b000_accessions() {
    let s = Spectrum::default();

    let owner_id = 1u32;
    let parent_index = 0u32;

    let out = assign_attributes(&s, TagId::Spectrum, owner_id, parent_index);
    assert_has_b000_tail(&out, ACC_ATTR_ID);
    assert_text(&out, ACC_ATTR_ID, "");
    assert_missing_b000_tail(&out, ACC_ATTR_INDEX);
    assert_missing_b000_tail(&out, ACC_ATTR_DEFAULT_ARRAY_LENGTH);
    assert_missing_b000_tail(&out, ACC_ATTR_NATIVE_ID);
    assert!(out.iter().all(|m| {
        m.accession
            .as_deref()
            .map(|a| a.starts_with(CV_REF_ATTR))
            .unwrap_or(false)
    }));
}

#[test]
fn assign_attrs_spectrum_id_is_preserved() {
    let mut s = Spectrum::default();
    s.id = "ok".to_string();

    let owner_id = 1u32;
    let parent_index = 0u32;

    let out = assign_attributes(&s, TagId::Spectrum, owner_id, parent_index);

    assert_has_b000_tail(&out, ACC_ATTR_ID);
    assert_text(&out, ACC_ATTR_ID, "ok");
}

#[test]
fn assign_attrs_chromatogram_emits_schema_b000_only() {
    let mut c = Chromatogram::default();

    c.id = "TIC".to_string();
    c.index = Some(0);
    c.default_array_length = Some(3476);
    c.native_id = Some("nativeX".to_string());
    c.data_processing_ref = Some("dpX".to_string());

    let owner_id = 7u32;
    let parent_index = 0u32;

    let out = assign_attributes(&c, TagId::Chromatogram, owner_id, parent_index);

    assert_has_b000_tail(&out, ACC_ATTR_ID);
    assert_has_b000_tail(&out, ACC_ATTR_INDEX);
    assert_has_b000_tail(&out, ACC_ATTR_DEFAULT_ARRAY_LENGTH);
    assert_has_b000_tail(&out, ACC_ATTR_NATIVE_ID);
    assert_has_b000_tail(&out, ACC_ATTR_DATA_PROCESSING_REF);

    assert_text(&out, ACC_ATTR_ID, "TIC");
    assert_num(&out, ACC_ATTR_INDEX, 0.0);
    assert_num(&out, ACC_ATTR_DEFAULT_ARRAY_LENGTH, 3476.0);
    assert_text(&out, ACC_ATTR_NATIVE_ID, "nativeX");
    assert_text(&out, ACC_ATTR_DATA_PROCESSING_REF, "dpX");
}
