use std::collections::HashMap;

use crate::{
    BinaryDataArray, BinaryDataArrayList, CvParam, UserParam,
    b64::utilities::common::{is_cv_prefix, unit_cv_ref, value_to_opt_string},
    decode::{Metadatum, MetadatumValue},
    mzml::{
        attr_meta::{
            ACC_ATTR_ARRAY_LENGTH, ACC_ATTR_COUNT, ACC_ATTR_DATA_PROCESSING_REF,
            ACC_ATTR_DEFAULT_ARRAY_LENGTH, ACC_ATTR_ENCODED_LENGTH, CV_REF_ATTR,
        },
        cv_table,
        schema::TagId,
    },
};

/// <binaryDataArrayList>
#[inline]
pub fn parse_binary_data_array_list(metadata: &[Metadatum]) -> Option<BinaryDataArrayList> {
    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::BinaryDataArrayList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::BinaryDataArray)
                .map(|m| m.parent_index)
        })?;

    let mut count: Option<usize> = None;

    let mut groups: HashMap<u32, Vec<&Metadatum>> = HashMap::new();
    for m in metadata {
        if count.is_none()
            && m.tag_id == TagId::BinaryDataArrayList
            && m.owner_id == list_id
            && b000_tail(m.accession.as_deref()) == Some(ACC_ATTR_COUNT)
        {
            count = as_u32(&m.value).map(|v| v as usize);
        }

        if m.tag_id == TagId::BinaryDataArray && m.parent_index == list_id {
            groups.entry(m.owner_id).or_default().push(m);
        }
    }

    if groups.is_empty() {
        return Some(BinaryDataArrayList {
            count: count.or(Some(0)),
            binary_data_arrays: Vec::new(),
        });
    }

    let mut ids: Vec<u32> = groups.keys().copied().collect();
    ids.sort_unstable();

    let mut binary_data_arrays = Vec::with_capacity(ids.len());
    for id in ids {
        if let Some(group) = groups.get(&id) {
            binary_data_arrays.push(parse_binary_data_array(group));
        }
    }

    inherit_array_length_from_parent(metadata, &mut binary_data_arrays);

    Some(BinaryDataArrayList {
        count: count.or(Some(binary_data_arrays.len())),
        binary_data_arrays,
    })
}

/// <binaryDataArrayList>
#[inline]
fn inherit_array_length_from_parent(metadata: &[Metadatum], bdas: &mut [BinaryDataArray]) {
    let parent_default = metadata
        .iter()
        .find_map(|m| match b000_tail(m.accession.as_deref()) {
            Some(tail) if tail == ACC_ATTR_DEFAULT_ARRAY_LENGTH => {
                as_u32(&m.value).map(|v| v as usize)
            }
            _ => None,
        });

    if let Some(len) = parent_default {
        for bda in bdas {
            if bda.array_length.is_none() {
                bda.array_length = Some(len);
            }
        }
        return;
    }

    for bda in bdas {
        if bda.array_length.is_some() {
            continue;
        }
        let Some(enc_chars) = bda.encoded_length else {
            continue;
        };
        if enc_chars % 4 != 0 {
            continue;
        }

        let elem_bytes = if bda.is_f64 == Some(true) {
            8usize
        } else if bda.is_f32 == Some(true) {
            4usize
        } else {
            continue;
        };

        let decoded_bytes = (enc_chars / 4) * 3;
        if decoded_bytes % elem_bytes == 0 {
            bda.array_length = Some(decoded_bytes / elem_bytes);
        }
    }
}

/// <binaryDataArray>
#[inline]
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
        let Some((prefix, tail)) = acc.split_once(':') else {
            continue;
        };

        if prefix == CV_REF_ATTR {
            if let Ok(tail_u32) = tail.parse::<u32>() {
                match tail_u32 {
                    ACC_ATTR_ARRAY_LENGTH => {
                        out.array_length = as_u32(&m.value).map(|v| v as usize)
                    }
                    ACC_ATTR_ENCODED_LENGTH => {
                        out.encoded_length = as_u32(&m.value).map(|v| v as usize)
                    }
                    ACC_ATTR_DATA_PROCESSING_REF => out.data_processing_ref = as_string(&m.value),
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

#[inline]
fn b000_tail(acc: Option<&str>) -> Option<u32> {
    let acc = acc?;
    let (prefix, tail) = acc.split_once(':')?;
    if prefix != CV_REF_ATTR {
        return None;
    }
    tail.parse::<u32>().ok()
}

#[inline]
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

#[inline]
fn as_string(v: &MetadatumValue) -> Option<String> {
    match v {
        MetadatumValue::Text(s) => Some(s.clone()),
        MetadatumValue::Number(f) => Some(f.to_string()),
        MetadatumValue::Empty => None,
    }
}
