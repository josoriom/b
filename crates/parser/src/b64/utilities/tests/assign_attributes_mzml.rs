use std::{
    collections::{HashMap, HashSet},
    fs,
    path::PathBuf,
    sync::OnceLock,
};

use crate::{
    b64::{
        decode::{Metadatum, MetadatumValue},
        utilities::common::find_node_by_tag,
        utilities::{assign_attributes, parse_header, parse_metadata},
    },
    mzml::{
        attr_meta::*,
        schema::{TagId, schema},
        structs::{MzML, Spectrum},
    },
    utilities::test::mzml as parse_mzml_cached,
};

static MZML_CACHE: OnceLock<MzML> = OnceLock::new();

const MZML_PATH: &str = "data/mzml/test.mzML";
const B64_PATH: &str = "data/b64/test.b64";

fn read_bytes(path: &str) -> Vec<u8> {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(path);
    fs::read(&full).unwrap_or_else(|e| panic!("cannot read {:?}: {}", full, e))
}

fn spectra_meta_from_test_b64() -> Vec<Metadatum> {
    let bytes = read_bytes(B64_PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let c0 = header.off_spec_meta as usize;
    let c1 = header.off_chrom_meta as usize;

    assert!(c0 < c1, "invalid spectra meta offsets (start >= end)");
    assert!(
        c1 <= bytes.len(),
        "invalid spectra meta offsets (end out of bounds)"
    );

    let compressed = (header.reserved_flags & (1u8 << 4)) != 0;

    parse_metadata(
        &bytes[c0..c1],
        header.spectrum_count,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        compressed,
        header.reserved_flags,
    )
    .expect("parse_metadata(spectra) failed")
}

fn schema_attrs_for_tag(tag: TagId) -> &'static HashMap<String, Vec<String>> {
    let s = schema();
    let node =
        find_node_by_tag(s, tag).unwrap_or_else(|| panic!("schema missing node for {tag:?}"));
    &node.attributes
}

fn parse_b000_tail(accession: &str) -> Option<u32> {
    let (cv, tail) = accession.split_once(':')?;
    if cv != CV_REF_ATTR {
        return None;
    }
    tail.parse::<u32>().ok()
}

fn tail_to_field_key(tail: u32) -> Option<&'static str> {
    match tail {
        ACC_ATTR_ID => Some("id"),
        ACC_ATTR_REF => Some("ref"),
        ACC_ATTR_LOCATION => Some("location"),
        ACC_ATTR_START_TIME_STAMP => Some("start_time_stamp"),
        ACC_ATTR_DEFAULT_INSTRUMENT_CONFIGURATION_REF => {
            Some("default_instrument_configuration_ref")
        }
        ACC_ATTR_DEFAULT_SOURCE_FILE_REF => Some("default_source_file_ref"),
        ACC_ATTR_SAMPLE_REF => Some("sample_ref"),
        ACC_ATTR_DEFAULT_DATA_PROCESSING_REF => Some("default_data_processing_ref"),
        ACC_ATTR_DATA_PROCESSING_REF => Some("data_processing_ref"),
        ACC_ATTR_SOURCE_FILE_REF => Some("source_file_ref"),
        ACC_ATTR_NATIVE_ID => Some("native_id"),
        ACC_ATTR_SPOT_ID => Some("spot_id"),
        ACC_ATTR_EXTERNAL_SPECTRUM_ID => Some("external_spectrum_id"),
        ACC_ATTR_SPECTRUM_REF => Some("spectrum_ref"),
        ACC_ATTR_SCAN_SETTINGS_REF => Some("scan_settings_ref"),
        ACC_ATTR_INSTRUMENT_CONFIGURATION_REF => Some("instrument_configuration_ref"),
        ACC_ATTR_SOFTWARE_REF => Some("software_ref"),
        ACC_ATTR_VERSION => Some("version"),
        ACC_ATTR_COUNT => Some("count"),
        ACC_ATTR_ORDER => Some("order"),
        ACC_ATTR_INDEX => Some("index"),
        ACC_ATTR_SCAN_NUMBER => Some("scan_number"),
        ACC_ATTR_DEFAULT_ARRAY_LENGTH => Some("default_array_length"),
        ACC_ATTR_ARRAY_LENGTH => Some("array_length"),
        ACC_ATTR_ENCODED_LENGTH => Some("encoded_length"),
        ACC_ATTR_MS_LEVEL => Some("ms_level"),
        _ => None,
    }
}

fn find_by_tail<'a>(meta: &'a [Metadatum], tail: u32) -> Option<&'a Metadatum> {
    meta.iter().find(|m| {
        m.accession
            .as_deref()
            .and_then(parse_b000_tail)
            .is_some_and(|t| t == tail)
    })
}

fn assert_meta_value_equiv(a: &MetadatumValue, b: &MetadatumValue, ctx: &str) {
    match (a, b) {
        (MetadatumValue::Text(x), MetadatumValue::Text(y)) => assert_eq!(x, y, "{ctx}"),
        (MetadatumValue::Number(x), MetadatumValue::Number(y)) => {
            let dx = (x - y).abs();
            assert!(dx <= 1e-9, "{ctx}: {x} != {y}");
        }
        (MetadatumValue::Text(x), MetadatumValue::Number(y)) => {
            let px = x.parse::<f64>().unwrap_or(f64::NAN);
            let dx = (px - y).abs();
            assert!(dx <= 1e-9, "{ctx}: Text({x}) != Number({y})");
        }
        (MetadatumValue::Number(x), MetadatumValue::Text(y)) => {
            let py = y.parse::<f64>().unwrap_or(f64::NAN);
            let dx = (x - py).abs();
            assert!(dx <= 1e-9, "{ctx}: Number({x}) != Text({y})");
        }
        _ => panic!("{ctx}: type mismatch: {a:?} vs {b:?}"),
    }
}

fn candidate_schema_attr_tails_for_owner(
    meta: &[Metadatum],
    owner_id: u32,
    schema_attrs: &HashMap<String, Vec<String>>,
) -> HashSet<u32> {
    let mut tails = HashSet::new();

    for m in meta.iter().filter(|m| m.owner_id == owner_id) {
        let Some(acc) = m.accession.as_deref() else {
            continue;
        };
        let Some(tail) = parse_b000_tail(acc) else {
            continue;
        };
        let Some(k) = tail_to_field_key(tail) else {
            continue;
        };

        if schema_attrs.contains_key(k) {
            tails.insert(tail);
        }
    }

    tails
}

fn find_owner_id_by_tail_text(meta: &[Metadatum], tag: TagId, tail: u32, text: &str) -> u32 {
    meta.iter()
        .find(|m| {
            m.tag_id == tag
                && m.accession.as_deref().and_then(parse_b000_tail) == Some(tail)
                && matches!(&m.value, MetadatumValue::Text(t) if t == text)
        })
        .unwrap_or_else(|| panic!("cannot find owner_id for {tag:?} tail={tail} text={text:?}"))
        .owner_id
}

#[test]
fn assign_attrs_spectrum_list_b64_equiv_schema_subset() {
    let mzml = parse_mzml_cached(&MZML_CACHE, MZML_PATH);
    let sl = mzml
        .run
        .spectrum_list
        .as_ref()
        .expect("spectrumList parsed");

    let b64_meta = spectra_meta_from_test_b64();
    let schema_attrs = schema_attrs_for_tag(TagId::SpectrumList);

    let expected_ddpr = sl
        .default_data_processing_ref
        .as_deref()
        .expect("mzML spectrumList default_data_processing_ref present");

    let owner_id = find_owner_id_by_tail_text(
        &b64_meta,
        TagId::SpectrumList,
        ACC_ATTR_DEFAULT_DATA_PROCESSING_REF,
        expected_ddpr,
    );

    let parent_index = 0u32;
    let generated = assign_attributes(sl, TagId::SpectrumList, owner_id, parent_index);

    let tails = candidate_schema_attr_tails_for_owner(&b64_meta, owner_id, schema_attrs);
    assert!(
        !tails.is_empty(),
        "no schema-listed B000 attrs found for SpectrumList owner_id={owner_id}"
    );

    for tail in tails {
        let b = b64_meta
            .iter()
            .find(|m| {
                m.owner_id == owner_id
                    && m.accession.as_deref().and_then(parse_b000_tail) == Some(tail)
            })
            .unwrap_or_else(|| panic!("b64 missing tail={tail} for owner_id={owner_id}"));

        let g = find_by_tail(&generated, tail)
            .unwrap_or_else(|| panic!("generated meta missing tail={tail}"));

        assert_eq!(
            g.tag_id,
            TagId::SpectrumList,
            "generated tag_id wrong for tail={tail}"
        );
        assert_meta_value_equiv(&b.value, &g.value, &format!("tail={tail}"));
    }

    assert!(
        generated.iter().all(|m| {
            m.accession
                .as_deref()
                .map(|a| a.starts_with(CV_REF_ATTR))
                .unwrap_or(true)
        }),
        "generated assign_attributes should only emit {CV_REF_ATTR}:* accessions: {generated:#?}"
    );
}

#[test]
fn assign_attrs_spectrum_first_b64_equiv_schema_subset() {
    let mzml = parse_mzml_cached(&MZML_CACHE, MZML_PATH);
    let sl = mzml
        .run
        .spectrum_list
        .as_ref()
        .expect("spectrumList parsed");
    let s0 = &sl.spectra[0];

    let b64_meta = spectra_meta_from_test_b64();
    let schema_attrs = schema_attrs_for_tag(TagId::Spectrum);

    let owner_id =
        find_owner_id_by_tail_text(&b64_meta, TagId::Spectrum, ACC_ATTR_ID, s0.id.as_str());

    let parent_index = 0u32;
    let generated = assign_attributes(s0, TagId::Spectrum, owner_id, parent_index);

    let tails = candidate_schema_attr_tails_for_owner(&b64_meta, owner_id, schema_attrs);
    assert!(
        !tails.is_empty(),
        "no schema-listed B000 attrs found for Spectrum owner_id={owner_id}"
    );

    for tail in tails {
        let b = b64_meta
            .iter()
            .find(|m| {
                m.owner_id == owner_id
                    && m.accession.as_deref().and_then(parse_b000_tail) == Some(tail)
            })
            .unwrap_or_else(|| panic!("b64 missing tail={tail} for owner_id={owner_id}"));

        let g = find_by_tail(&generated, tail)
            .unwrap_or_else(|| panic!("generated meta missing tail={tail}"));

        assert_eq!(
            g.tag_id,
            TagId::Spectrum,
            "generated tag_id wrong for tail={tail}"
        );
        assert_meta_value_equiv(&b.value, &g.value, &format!("tail={tail}"));
    }
}

#[test]
fn assign_attrs_spectrum_last_b64_equiv_schema_subset() {
    let mzml = parse_mzml_cached(&MZML_CACHE, MZML_PATH);
    let sl = mzml
        .run
        .spectrum_list
        .as_ref()
        .expect("spectrumList parsed");
    let s_last: &Spectrum = sl.spectra.last().expect("last spectrum");

    let b64_meta = spectra_meta_from_test_b64();
    let schema_attrs = schema_attrs_for_tag(TagId::Spectrum);

    let owner_id =
        find_owner_id_by_tail_text(&b64_meta, TagId::Spectrum, ACC_ATTR_ID, s_last.id.as_str());

    let parent_index = 0u32;
    let generated = assign_attributes(s_last, TagId::Spectrum, owner_id, parent_index);

    let tails = candidate_schema_attr_tails_for_owner(&b64_meta, owner_id, schema_attrs);
    assert!(
        !tails.is_empty(),
        "no schema-listed B000 attrs found for Spectrum owner_id={owner_id}"
    );

    for tail in tails {
        let b = b64_meta
            .iter()
            .find(|m| {
                m.owner_id == owner_id
                    && m.accession.as_deref().and_then(parse_b000_tail) == Some(tail)
            })
            .unwrap_or_else(|| panic!("b64 missing tail={tail} for owner_id={owner_id}"));

        let g = find_by_tail(&generated, tail)
            .unwrap_or_else(|| panic!("generated meta missing tail={tail}"));

        assert_eq!(
            g.tag_id,
            TagId::Spectrum,
            "generated tag_id wrong for tail={tail}"
        );
        assert_meta_value_equiv(&b.value, &g.value, &format!("tail={tail}"));
    }
}
