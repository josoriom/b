use std::collections::{HashMap, HashSet};

use crate::{
    b64::utilities::{
        common::{ChildIndex, child_node, find_node_by_tag, ordered_unique_owner_ids},
        parse_cv_and_user_params,
        parse_file_description::{
            allowed_from_rows, b000_attr_text, child_params_for_parent, is_child_of,
            rows_for_owner, unique_ids,
        },
    },
    decode::Metadatum,
    mzml::{
        attr_meta::{ACC_ATTR_COUNT, ACC_ATTR_ID},
        schema::{SchemaTree as Schema, TagId},
        structs::{ReferenceableParamGroup, ReferenceableParamGroupList},
    },
};

#[inline]
fn parse_params(
    allowed_schema: &HashSet<&str>,
    rows: &[&Metadatum],
) -> (
    Vec<crate::mzml::structs::CvParam>,
    Vec<crate::mzml::structs::UserParam>,
) {
    if allowed_schema.is_empty() {
        let allowed_meta = allowed_from_rows(rows);
        parse_cv_and_user_params(&allowed_meta, rows)
    } else {
        parse_cv_and_user_params(allowed_schema, rows)
    }
}

/// <referenceableParamGroupList>
#[inline]
pub fn parse_referenceable_param_group_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<ReferenceableParamGroupList> {
    let list_node = find_node_by_tag(schema, TagId::ReferenceableParamGroupList)?;
    let group_node = child_node(Some(list_node), TagId::ReferenceableParamGroup)?;

    let allowed_schema: HashSet<&str> = child_node(Some(group_node), TagId::CvParam)
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::ReferenceableParamGroupList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::ReferenceableParamGroup)
                .map(|m| m.parent_index)
        })?;

    let mut group_ids = unique_ids(child_index.ids(list_id, TagId::ReferenceableParamGroup));
    if group_ids.is_empty() {
        group_ids = ordered_unique_owner_ids(metadata, TagId::ReferenceableParamGroup);
        group_ids.retain(|&id| is_child_of(&owner_rows, id, list_id));
    }

    if group_ids.is_empty() {
        return None;
    }

    let mut referenceable_param_groups = Vec::with_capacity(group_ids.len());
    for group_id in group_ids {
        referenceable_param_groups.push(parse_referenceable_param_group(
            &allowed_schema,
            &owner_rows,
            child_index,
            group_id,
        ));
    }

    let count = b000_attr_text(rows_for_owner(&owner_rows, list_id), ACC_ATTR_COUNT)
        .and_then(|s| s.parse::<usize>().ok())
        .or(Some(referenceable_param_groups.len()));

    Some(ReferenceableParamGroupList {
        count,
        referenceable_param_groups,
    })
}

/// <referenceableParamGroup>
#[inline]
fn parse_referenceable_param_group(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    group_id: u32,
) -> ReferenceableParamGroup {
    let rows = rows_for_owner(owner_rows, group_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).unwrap_or_default();

    let child_meta = child_params_for_parent(owner_rows, child_index, group_id);
    let mut params_meta = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_params, user_params) = parse_params(allowed_schema, &params_meta);

    ReferenceableParamGroup {
        id,
        cv_params,
        user_params,
    }
}
