use std::collections::HashMap;

use crate::{
    b64::utilities::{
        common::{ChildIndex, child_node, find_node_by_tag, ordered_unique_owner_ids},
        parse_file_description::{b000_attr_text, is_child_of, rows_for_owner, unique_ids},
    },
    decode::Metadatum,
    mzml::{
        attr_meta::{
            ACC_ATTR_CV_FULL_NAME, ACC_ATTR_CV_URI, ACC_ATTR_CV_VERSION, ACC_ATTR_ID,
            ACC_ATTR_LABEL,
        },
        schema::{SchemaTree as Schema, TagId},
        structs::{Cv, CvList},
    },
};

/// <cvList>
#[inline]
pub fn parse_cv_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<CvList> {
    if let Some(list_node) = find_node_by_tag(schema, TagId::CvList) {
        let _ = child_node(Some(list_node), TagId::Cv);
    }

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::CvList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::Cv)
                .map(|m| m.parent_index)
        })?;

    let mut cv_ids = unique_ids(child_index.ids(list_id, TagId::Cv));
    if cv_ids.is_empty() {
        cv_ids = ordered_unique_owner_ids(metadata, TagId::Cv);
        cv_ids.retain(|&id| is_child_of(&owner_rows, id, list_id));
    }

    if cv_ids.is_empty() {
        return None;
    }

    let mut cv = Vec::with_capacity(cv_ids.len());
    for id in cv_ids {
        cv.push(parse_cv(&owner_rows, id));
    }

    Some(CvList {
        count: Some(cv.len()),
        cv,
    })
}

/// <cv>
#[inline]
fn parse_cv(owner_rows: &HashMap<u32, Vec<&Metadatum>>, cv_id: u32) -> Cv {
    let rows = rows_for_owner(owner_rows, cv_id);

    let id = b000_attr_text(rows, ACC_ATTR_LABEL)
        .or_else(|| b000_attr_text(rows, ACC_ATTR_ID))
        .unwrap_or_default();

    let full_name = b000_attr_text(rows, ACC_ATTR_CV_FULL_NAME).filter(|s| !s.is_empty());
    let version = b000_attr_text(rows, ACC_ATTR_CV_VERSION).filter(|s| !s.is_empty());
    let uri = b000_attr_text(rows, ACC_ATTR_CV_URI).filter(|s| !s.is_empty());

    Cv {
        id,
        full_name,
        version,
        uri,
    }
}
