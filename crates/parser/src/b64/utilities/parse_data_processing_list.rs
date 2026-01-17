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
        attr_meta::{ACC_ATTR_ID, ACC_ATTR_ORDER, ACC_ATTR_REF, ACC_ATTR_SOFTWARE_REF},
        schema::{SchemaTree as Schema, TagId},
        structs::{
            DataProcessing, DataProcessingList, ProcessingMethod, ReferenceableParamGroupRef,
        },
    },
};

/// <dataProcessingList>
#[inline]
pub fn parse_data_processing_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<DataProcessingList> {
    let list_node = find_node_by_tag(schema, TagId::DataProcessingList)?;
    let dp_node = child_node(Some(list_node), TagId::DataProcessing)?;

    let pm_node = child_node(Some(dp_node), TagId::ProcessingMethod)?;
    let allowed_processing_method_schema: HashSet<&str> = child_node(Some(pm_node), TagId::CvParam)
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::DataProcessingList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::DataProcessing)
                .map(|m| m.parent_index)
        })?;

    let mut data_processing_ids = unique_ids(child_index.ids(list_id, TagId::DataProcessing));
    if data_processing_ids.is_empty() {
        data_processing_ids = ordered_unique_owner_ids(metadata, TagId::DataProcessing);
        data_processing_ids.retain(|&id| is_child_of(&owner_rows, id, list_id));
    }

    if data_processing_ids.is_empty() {
        return None;
    }

    let mut data_processing = Vec::with_capacity(data_processing_ids.len());
    for id in data_processing_ids {
        data_processing.push(parse_data_processing(
            &allowed_processing_method_schema,
            &owner_rows,
            child_index,
            metadata,
            id,
        ));
    }

    Some(DataProcessingList {
        count: Some(data_processing.len()),
        data_processing,
    })
}

/// <dataProcessing>
#[inline]
fn parse_data_processing(
    allowed_processing_method_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    data_processing_id: u32,
) -> DataProcessing {
    let rows = rows_for_owner(owner_rows, data_processing_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).unwrap_or_default();
    let software_ref = b000_attr_text(rows, ACC_ATTR_SOFTWARE_REF).filter(|s| !s.is_empty());

    let processing_method = parse_processing_methods(
        allowed_processing_method_schema,
        owner_rows,
        child_index,
        metadata,
        data_processing_id,
    );

    DataProcessing {
        id,
        software_ref,
        processing_method,
    }
}

/// <processingMethod>
#[inline]
fn parse_processing_methods(
    allowed_processing_method_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    data_processing_id: u32,
) -> Vec<ProcessingMethod> {
    let mut ids = unique_ids(child_index.ids(data_processing_id, TagId::ProcessingMethod));
    if ids.is_empty() {
        ids = ordered_unique_owner_ids(metadata, TagId::ProcessingMethod);
        ids.retain(|&id| is_child_of(owner_rows, id, data_processing_id));
    }

    let mut out = Vec::with_capacity(ids.len());
    for id in ids {
        out.push(parse_processing_method(
            allowed_processing_method_schema,
            owner_rows,
            child_index,
            id,
        ));
    }
    out
}

/// <processingMethod>
#[inline]
fn parse_processing_method(
    allowed_processing_method_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    processing_method_id: u32,
) -> ProcessingMethod {
    let rows = rows_for_owner(owner_rows, processing_method_id);

    let order = b000_attr_text(rows, ACC_ATTR_ORDER).and_then(|s| s.parse::<u32>().ok());
    let software_ref = b000_attr_text(rows, ACC_ATTR_SOFTWARE_REF).filter(|s| !s.is_empty());

    let referenceable_param_group_ref =
        parse_referenceable_param_group_refs(owner_rows, child_index, processing_method_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, processing_method_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, user_param) = if allowed_processing_method_schema.is_empty() {
        let allowed_meta = allowed_from_rows(&params_meta);
        parse_cv_and_user_params(&allowed_meta, &params_meta)
    } else {
        parse_cv_and_user_params(allowed_processing_method_schema, &params_meta)
    };

    ProcessingMethod {
        order,
        software_ref,
        referenceable_param_group_ref,
        cv_param,
        user_param,
    }
}

/// <referenceableParamGroupRef>
#[inline]
fn parse_referenceable_param_group_refs(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    parent_id: u32,
) -> Vec<ReferenceableParamGroupRef> {
    let mut ref_ids = unique_ids(child_index.ids(parent_id, TagId::ReferenceableParamGroupRef));
    if ref_ids.is_empty() {
        return Vec::new();
    }

    let mut out: Vec<ReferenceableParamGroupRef> = Vec::with_capacity(ref_ids.len());
    for ref_id in ref_ids.drain(..) {
        let ref_rows = rows_for_owner(owner_rows, ref_id);
        if let Some(r) = b000_attr_text(ref_rows, ACC_ATTR_REF) {
            if !r.is_empty() {
                out.push(ReferenceableParamGroupRef { r#ref: r });
            }
        }
    }

    out
}
