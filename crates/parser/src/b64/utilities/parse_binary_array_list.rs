use std::collections::HashMap;

use crate::{
    BinaryDataArray, BinaryDataArrayList, CvParam, NumericType, UserParam,
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

    let mut index: HashMap<u32, usize> = HashMap::with_capacity(metadata.len() / 8 + 1);
    let mut tmp: Vec<(u32, BinaryDataArray)> = Vec::new();

    for m in metadata {
        if count.is_none()
            && m.tag_id == TagId::BinaryDataArrayList
            && m.owner_id == list_id
            && b000_tail(m.accession.as_deref()) == Some(ACC_ATTR_COUNT)
        {
            count = as_u32(&m.value).map(|v| v as usize);
            continue;
        }

        if m.tag_id != TagId::BinaryDataArray || m.parent_index != list_id {
            continue;
        }

        let id = m.owner_id;
        let at = match index.get(&id).copied() {
            Some(i) => i,
            None => {
                let i = tmp.len();
                tmp.push((id, new_binary_data_array()));
                index.insert(id, i);
                i
            }
        };

        apply_binary_data_array_metadatum(&mut tmp[at].1, m);
    }

    if tmp.is_empty() {
        return Some(BinaryDataArrayList {
            count: count.or(Some(0)),
            binary_data_arrays: Vec::new(),
        });
    }

    tmp.sort_unstable_by_key(|(id, _)| *id);

    let mut binary_data_arrays: Vec<BinaryDataArray> =
        tmp.into_iter().map(|(_, bda)| bda).collect();

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

    let Some(len) = parent_default else {
        return;
    };

    for bda in bdas {
        if bda.array_length.is_none() {
            bda.array_length = Some(len);
        }
    }
}

/// <binaryDataArray>
#[inline]
fn new_binary_data_array() -> BinaryDataArray {
    BinaryDataArray {
        array_length: None,
        encoded_length: None,
        data_processing_ref: None,
        referenceable_param_group_refs: Vec::new(),
        cv_params: Vec::with_capacity(8),
        user_params: Vec::with_capacity(2),
        numeric_type: None,
        binary: None,
    }
}

/// <binaryDataArray>
#[inline]
fn parse_binary_data_array(metadata: &[&Metadatum]) -> BinaryDataArray {
    let mut out = new_binary_data_array();
    out.cv_params
        .reserve(metadata.len().saturating_sub(out.cv_params.capacity()));
    for m in metadata {
        apply_binary_data_array_metadatum(&mut out, m);
    }
    out
}

#[inline]
fn apply_binary_data_array_metadatum(out: &mut BinaryDataArray, m: &Metadatum) {
    let Some(acc) = m.accession.as_deref() else {
        return;
    };
    let Some((prefix, tail)) = acc.split_once(':') else {
        return;
    };

    if prefix == CV_REF_ATTR {
        let Ok(tail_u32) = tail.parse::<u32>() else {
            return;
        };
        match tail_u32 {
            ACC_ATTR_ARRAY_LENGTH => out.array_length = as_u32(&m.value).map(|v| v as usize),
            ACC_ATTR_ENCODED_LENGTH => out.encoded_length = as_u32(&m.value).map(|v| v as usize),
            ACC_ATTR_DATA_PROCESSING_REF => out.data_processing_ref = as_string(&m.value),
            _ => {}
        }
        return;
    }

    let value = value_to_opt_string(&m.value);
    let unit_accession = m.unit_accession.clone();
    let unit_cv_ref = unit_cv_ref(&unit_accession);

    if is_cv_prefix(prefix) {
        if prefix == "MS" {
            if let Ok(tail_u32) = tail.parse::<u32>() {
                let new_ty = match tail_u32 {
                    1000519 => Some(NumericType::Int32),
                    1000521 => Some(NumericType::Float32),
                    1000522 => Some(NumericType::Int64),
                    1000523 => Some(NumericType::Float64),
                    _ => None,
                };

                if let Some(nt) = new_ty {
                    out.numeric_type = match out.numeric_type {
                        None => Some(nt),
                        Some(cur) if cur == nt => Some(cur),
                        Some(_) => None,
                    };
                }
            }
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
