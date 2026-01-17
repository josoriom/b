use std::collections::HashMap;

use crate::{
    b64::utilities::{
        common::{ChildIndex, child_node, find_node_by_tag, ordered_unique_owner_ids},
        parse_file_description::{b000_attr_text, is_child_of, rows_for_owner, unique_ids},
    },
    decode::Metadatum,
    mzml::{
        attr_meta::{ACC_ATTR_ID, ACC_ATTR_NAME, ACC_ATTR_REF},
        schema::{SchemaTree as Schema, TagId},
        structs::{ReferenceableParamGroupRef, Sample, SampleList},
    },
};

/// <sampleList>
#[inline]
pub fn parse_sample_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<SampleList> {
    let list_node = find_node_by_tag(schema, TagId::SampleList)?;
    let _ = child_node(Some(list_node), TagId::Sample)?;

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let sample_list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::SampleList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::Sample)
                .map(|m| m.parent_index)
        })?;

    let mut sample_ids = unique_ids(child_index.ids(sample_list_id, TagId::Sample));
    if sample_ids.is_empty() {
        sample_ids = ordered_unique_owner_ids(metadata, TagId::Sample);
        sample_ids.retain(|&id| is_child_of(&owner_rows, id, sample_list_id));
    }

    if sample_ids.is_empty() {
        return None;
    }

    let mut samples = Vec::with_capacity(sample_ids.len());
    for sample_id in sample_ids {
        samples.push(parse_sample(&owner_rows, child_index, sample_id));
    }

    Some(SampleList {
        count: Some(samples.len() as u32),
        samples,
    })
}

/// <sample>
#[inline]
fn parse_sample(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    sample_id: u32,
) -> Sample {
    let rows = rows_for_owner(owner_rows, sample_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).unwrap_or_default();
    let name = b000_attr_text(rows, ACC_ATTR_NAME).expect("name not found");

    let mut referenceable_param_group_ref = None;

    let ref_ids = unique_ids(child_index.ids(sample_id, TagId::ReferenceableParamGroupRef));
    if let Some(&ref_id) = ref_ids.first() {
        let ref_rows = rows_for_owner(owner_rows, ref_id);
        if let Some(r) = b000_attr_text(ref_rows, ACC_ATTR_REF) {
            if !r.is_empty() {
                referenceable_param_group_ref = Some(ReferenceableParamGroupRef { r#ref: r });
            }
        }
    }

    Sample {
        id,
        name,
        referenceable_param_group_ref,
    }
}
