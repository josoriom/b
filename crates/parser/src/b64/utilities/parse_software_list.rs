use std::collections::{HashMap, HashSet};

use crate::{
    b64::utilities::{
        common::{ChildIndex, ordered_unique_owner_ids},
        parse_cv_and_user_params,
        parse_file_description::{
            b000_attr_text, child_params_for_parent, is_child_of, rows_for_owner, unique_ids,
        },
    },
    decode::Metadatum,
    mzml::{
        attr_meta::{ACC_ATTR_ID, ACC_ATTR_NAME, ACC_ATTR_REF, ACC_ATTR_VERSION},
        cv_table,
        schema::{SchemaTree as Schema, TagId},
        structs::{Software, SoftwareList, SoftwareParam},
    },
};

/// <softwareList>
#[inline]
pub fn parse_software_list(
    _schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<SoftwareList> {
    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::SoftwareList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::Software)
                .map(|m| m.parent_index)
        })?;

    let mut software_ids = unique_ids(child_index.ids(list_id, TagId::Software));
    if software_ids.is_empty() {
        software_ids = ordered_unique_owner_ids(metadata, TagId::Software);
        software_ids.retain(|&id| is_child_of(&owner_rows, id, list_id));
    }

    if software_ids.is_empty() {
        return None;
    }

    let empty_allowed: HashSet<&str> = HashSet::new();

    let mut software = Vec::with_capacity(software_ids.len());
    for id in software_ids {
        software.push(parse_software(
            &owner_rows,
            child_index,
            metadata,
            id,
            &empty_allowed,
        ));
    }

    Some(SoftwareList {
        count: Some(software.len()),
        software,
    })
}

/// <software>
#[inline]
fn parse_software(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    software_id: u32,
    empty_allowed: &HashSet<&str>,
) -> Software {
    let rows = rows_for_owner(owner_rows, software_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).unwrap_or_default();
    let version = b000_attr_text(rows, ACC_ATTR_VERSION).filter(|s| !s.is_empty());

    let software_param = parse_software_params(
        owner_rows,
        child_index,
        metadata,
        software_id,
        version.as_deref(),
        empty_allowed,
    );

    let child_meta = child_params_for_parent(owner_rows, child_index, software_id);
    let mut params_meta = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, _) = parse_cv_and_user_params(empty_allowed, &params_meta);

    Software {
        id,
        version,
        software_param,
        cv_param,
    }
}

/// <softwareParam>
#[inline]
fn parse_software_params(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    software_id: u32,
    parent_version: Option<&str>,
    empty_allowed: &HashSet<&str>,
) -> Vec<SoftwareParam> {
    let mut ids = unique_ids(child_index.ids(software_id, TagId::SoftwareParam));
    if ids.is_empty() {
        ids = ordered_unique_owner_ids(metadata, TagId::SoftwareParam);
        ids.retain(|&id| is_child_of(owner_rows, id, software_id));
    }

    if ids.is_empty() {
        return Vec::new();
    }

    let mut out = Vec::with_capacity(ids.len());
    for id in ids {
        out.push(parse_software_param(
            owner_rows,
            id,
            parent_version,
            empty_allowed,
        ));
    }
    out
}

/// <softwareParam>
#[inline]
fn parse_software_param(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    software_param_id: u32,
    parent_version: Option<&str>,
    empty_allowed: &HashSet<&str>,
) -> SoftwareParam {
    let rows = rows_for_owner(owner_rows, software_param_id);

    let version = b000_attr_text(rows, ACC_ATTR_VERSION)
        .filter(|s| !s.is_empty())
        .or_else(|| parent_version.map(|s| s.to_string()));

    let (cv_params, _user_params) = parse_cv_and_user_params(empty_allowed, rows);

    if let Some(cv) = cv_params.into_iter().next() {
        let accession = cv.accession.unwrap_or_default();

        let cv_ref = b000_attr_text(rows, ACC_ATTR_REF)
            .filter(|s| !s.is_empty())
            .or(cv.cv_ref);

        return SoftwareParam {
            cv_ref,
            accession,
            name: cv.name,
            version,
        };
    }

    let accession = software_param_accession(rows).unwrap_or_default();
    let cv_ref = b000_attr_text(rows, ACC_ATTR_REF).filter(|s| !s.is_empty());

    let name = b000_attr_text(rows, ACC_ATTR_NAME)
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| {
            cv_table::get(accession.as_str())
                .and_then(|v| v.as_str())
                .unwrap_or(accession.as_str())
                .to_string()
        });

    SoftwareParam {
        cv_ref,
        accession,
        name,
        version,
    }
}

#[inline]
fn software_param_accession(rows: &[&Metadatum]) -> Option<String> {
    for m in rows {
        let acc = m.accession.as_deref()?;
        if acc.starts_with("B000:") {
            continue;
        }
        if acc.contains(':') {
            return Some(acc.to_string());
        }
    }
    None
}
