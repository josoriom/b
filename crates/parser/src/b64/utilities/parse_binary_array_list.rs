use std::collections::HashMap;

use crate::{
    BinaryDataArray, BinaryDataArrayList, CvParam, UserParam,
    b64::utilities::common::{is_cv_prefix, unit_cv_ref, value_to_opt_string},
    decode2::{Metadatum, MetadatumValue},
    mzml::{
        attr_meta::{
            ACC_ATTR_ARRAY_LENGTH, ACC_ATTR_COUNT, ACC_ATTR_DATA_PROCESSING_REF,
            ACC_ATTR_ENCODED_LENGTH, CV_REF_ATTR,
        },
        cv_table,
        schema::TagId,
    },
};

pub fn parse_binary_data_array_list(metadata: &[Metadatum]) -> Option<BinaryDataArrayList> {
    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::BinaryDataArrayList)
        .map(|m| m.owner_id)?;

    let count = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::BinaryDataArrayList && m.owner_id == list_id)
        .find_map(|m| match b000_tail(m.accession.as_deref()) {
            Some(tail) if tail == ACC_ATTR_COUNT => as_u32(&m.value).map(|v| v as usize),
            _ => None,
        });

    let mut groups: HashMap<u32, Vec<&Metadatum>> = HashMap::new();
    for m in metadata {
        if m.tag_id == TagId::BinaryDataArray && m.parent_index == list_id {
            groups.entry(m.owner_id).or_default().push(m);
        }
    }

    let mut entries: Vec<(u32, Vec<&Metadatum>)> = groups.into_iter().collect();
    entries.sort_unstable_by_key(|(k, _)| *k);

    let mut binary_data_arrays = Vec::with_capacity(entries.len());
    for (_, group) in entries {
        binary_data_arrays.push(parse_binary_data_array(&group));
    }

    Some(BinaryDataArrayList {
        count: count.or(Some(binary_data_arrays.len())),
        binary_data_arrays,
    })
}

fn parse_binary_data_array(metadata: &[&Metadatum]) -> BinaryDataArray {
    let mut out = BinaryDataArray {
        array_length: None,
        encoded_length: None,
        data_processing_ref: None,
        referenceable_param_group_refs: Vec::new(),
        cv_params: Vec::with_capacity(metadata.len()),
        user_params: Vec::new(),
        is_f32: None,
        is_f64: None,
        decoded_binary_f32: Vec::new(),
        decoded_binary_f64: Vec::new(),
    };

    for m in metadata {
        let Some(acc) = m.accession.as_deref() else {
            continue;
        };
        let Some((prefix, _)) = acc.split_once(':') else {
            continue;
        };

        if prefix == CV_REF_ATTR {
            if let Some(tail) = b000_tail(Some(acc)) {
                match tail {
                    ACC_ATTR_ARRAY_LENGTH => {
                        out.array_length = as_u32(&m.value).map(|v| v as usize);
                    }
                    ACC_ATTR_ENCODED_LENGTH => {
                        out.encoded_length = as_u32(&m.value).map(|v| v as usize);
                    }
                    ACC_ATTR_DATA_PROCESSING_REF => {
                        out.data_processing_ref = as_string(&m.value);
                    }
                    _ => {}
                }
            }
            continue;
        }

        let value = value_to_opt_string(&m.value);
        let unit_accession = m.unit_accession.clone();
        let unit_cv_ref = unit_cv_ref(&unit_accession);

        if is_cv_prefix(prefix) {
            if acc == "MS:1000521" {
                out.is_f32 = Some(true);
            } else if acc == "MS:1000523" {
                out.is_f64 = Some(true);
            }

            let unit_name = unit_accession
                .as_deref()
                .and_then(|ua| cv_table::get(ua).and_then(|v| v.as_str()))
                .map(|s| s.to_string());

            out.cv_params.push(CvParam {
                cv_ref: Some(prefix.to_string()),
                accession: Some(acc.to_string()),
                name: cv_table::get(acc)
                    .and_then(|v| v.as_str())
                    .unwrap_or(acc)
                    .to_string(),
                value,
                unit_cv_ref,
                unit_name,
                unit_accession,
            });
        } else {
            out.user_params.push(UserParam {
                name: acc.to_string(),
                r#type: None,
                unit_accession,
                unit_cv_ref,
                unit_name: None,
                value,
            });
        }
    }

    out
}

fn b000_tail(acc: Option<&str>) -> Option<u32> {
    let acc = acc?;
    let (prefix, tail) = acc.split_once(':')?;
    if prefix != CV_REF_ATTR {
        return None;
    }
    tail.parse::<u32>().ok()
}

fn as_u32(v: &MetadatumValue) -> Option<u32> {
    match v {
        MetadatumValue::Number(f) => {
            if f.is_finite() && f.fract() == 0.0 && *f >= 0.0 && *f <= (u32::MAX as f64) {
                Some(*f as u32)
            } else {
                None
            }
        }
        MetadatumValue::Text(s) => s.parse::<u32>().ok(),
        MetadatumValue::Empty => None,
    }
}

fn as_string(v: &MetadatumValue) -> Option<String> {
    match v {
        MetadatumValue::Text(s) => Some(s.clone()),
        MetadatumValue::Number(f) => Some(f.to_string()),
        MetadatumValue::Empty => None,
    }
}
