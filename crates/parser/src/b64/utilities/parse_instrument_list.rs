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
        attr_meta::{ACC_ATTR_ID, ACC_ATTR_ORDER, ACC_ATTR_REF, ACC_ATTR_SCAN_SETTINGS_REF},
        schema::{SchemaTree as Schema, TagId},
        structs::{
            Analyzer, ComponentList, Detector, Instrument, InstrumentList, InstrumentSoftwareRef,
            ReferenceableParamGroupRef, ScanSettingsRef, Source,
        },
    },
};

#[inline]
fn build_owner_rows<'a>(metadata: &'a [Metadatum]) -> HashMap<u32, Vec<&'a Metadatum>> {
    let mut owner_rows: HashMap<u32, Vec<&'a Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }
    owner_rows
}

#[inline]
fn ids_for_parent(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    parent_id: u32,
    tag: TagId,
) -> Vec<u32> {
    let mut ids = unique_ids(child_index.ids(parent_id, tag));
    if ids.is_empty() {
        ids = ordered_unique_owner_ids(metadata, tag);
        ids.retain(|&id| is_child_of(owner_rows, id, parent_id));
    }
    ids
}

#[inline]
fn parse_params(
    allowed_schema: &HashSet<&str>,
    params_meta: &[&Metadatum],
) -> (
    Vec<crate::mzml::structs::CvParam>,
    Vec<crate::mzml::structs::UserParam>,
) {
    if allowed_schema.is_empty() {
        let allowed_meta = allowed_from_rows(params_meta);
        parse_cv_and_user_params(&allowed_meta, params_meta)
    } else {
        parse_cv_and_user_params(allowed_schema, params_meta)
    }
}

/// <instrumentConfigurationList>
#[inline]
pub fn parse_instrument_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<InstrumentList> {
    let list_node = find_node_by_tag(schema, TagId::InstrumentConfigurationList)?;
    let instrument_node = child_node(Some(list_node), TagId::Instrument)?;

    let allowed_instrument_schema: HashSet<&str> =
        child_node(Some(instrument_node), TagId::CvParam)
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let component_list_node = child_node(Some(instrument_node), TagId::ComponentList)?;

    let allowed_source_schema: HashSet<&str> =
        child_node(Some(component_list_node), TagId::ComponentSource)
            .and_then(|n| child_node(Some(n), TagId::CvParam))
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let allowed_analyzer_schema: HashSet<&str> =
        child_node(Some(component_list_node), TagId::ComponentAnalyzer)
            .and_then(|n| child_node(Some(n), TagId::CvParam))
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let allowed_detector_schema: HashSet<&str> =
        child_node(Some(component_list_node), TagId::ComponentDetector)
            .and_then(|n| child_node(Some(n), TagId::CvParam))
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let owner_rows = build_owner_rows(metadata);

    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::InstrumentConfigurationList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::Instrument)
                .map(|m| m.parent_index)
        })?;

    let instrument_ids = ids_for_parent(
        &owner_rows,
        child_index,
        metadata,
        list_id,
        TagId::Instrument,
    );
    if instrument_ids.is_empty() {
        return None;
    }

    let mut instrument = Vec::with_capacity(instrument_ids.len());
    for id in instrument_ids {
        instrument.push(parse_instrument(
            &allowed_instrument_schema,
            &allowed_source_schema,
            &allowed_analyzer_schema,
            &allowed_detector_schema,
            &owner_rows,
            child_index,
            metadata,
            id,
        ));
    }

    Some(InstrumentList {
        count: Some(instrument.len()),
        instrument,
    })
}

/// <instrumentConfiguration>
#[inline]
fn parse_instrument(
    allowed_instrument_schema: &HashSet<&str>,
    allowed_source_schema: &HashSet<&str>,
    allowed_analyzer_schema: &HashSet<&str>,
    allowed_detector_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    instrument_id: u32,
) -> Instrument {
    let rows = rows_for_owner(owner_rows, instrument_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).unwrap_or_default();

    let scan_settings_ref = b000_attr_text(rows, ACC_ATTR_SCAN_SETTINGS_REF)
        .filter(|s| !s.is_empty())
        .map(|s| ScanSettingsRef { r#ref: s });

    let referenceable_param_group_ref =
        parse_referenceable_param_group_refs(owner_rows, child_index, instrument_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, instrument_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, user_param) = parse_params(allowed_instrument_schema, &params_meta);

    let component_list = parse_component_list(
        allowed_source_schema,
        allowed_analyzer_schema,
        allowed_detector_schema,
        owner_rows,
        child_index,
        metadata,
        instrument_id,
    );

    let software_ref =
        parse_instrument_software_ref(owner_rows, child_index, metadata, instrument_id);

    Instrument {
        id,
        scan_settings_ref,
        cv_param,
        user_param,
        referenceable_param_group_ref,
        component_list,
        software_ref,
    }
}

/// <componentList>
#[inline]
fn parse_component_list(
    allowed_source_schema: &HashSet<&str>,
    allowed_analyzer_schema: &HashSet<&str>,
    allowed_detector_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    instrument_id: u32,
) -> Option<ComponentList> {
    let component_list_ids = ids_for_parent(
        owner_rows,
        child_index,
        metadata,
        instrument_id,
        TagId::ComponentList,
    );
    let component_list_id = component_list_ids.first().copied().unwrap_or(0);

    let primary_parent = if component_list_id != 0 {
        component_list_id
    } else {
        instrument_id
    };

    let mut source_ids = ids_for_parent(
        owner_rows,
        child_index,
        metadata,
        primary_parent,
        TagId::ComponentSource,
    );
    let mut analyzer_ids = ids_for_parent(
        owner_rows,
        child_index,
        metadata,
        primary_parent,
        TagId::ComponentAnalyzer,
    );
    let mut detector_ids = ids_for_parent(
        owner_rows,
        child_index,
        metadata,
        primary_parent,
        TagId::ComponentDetector,
    );

    if source_ids.is_empty()
        && analyzer_ids.is_empty()
        && detector_ids.is_empty()
        && primary_parent != instrument_id
    {
        source_ids = ids_for_parent(
            owner_rows,
            child_index,
            metadata,
            instrument_id,
            TagId::ComponentSource,
        );
        analyzer_ids = ids_for_parent(
            owner_rows,
            child_index,
            metadata,
            instrument_id,
            TagId::ComponentAnalyzer,
        );
        detector_ids = ids_for_parent(
            owner_rows,
            child_index,
            metadata,
            instrument_id,
            TagId::ComponentDetector,
        );
    }

    if source_ids.is_empty() && analyzer_ids.is_empty() && detector_ids.is_empty() {
        return None;
    }

    let mut source = Vec::with_capacity(source_ids.len());
    for id in source_ids {
        source.push(parse_source(
            allowed_source_schema,
            owner_rows,
            child_index,
            id,
        ));
    }

    let mut analyzer = Vec::with_capacity(analyzer_ids.len());
    for id in analyzer_ids {
        analyzer.push(parse_analyzer(
            allowed_analyzer_schema,
            owner_rows,
            child_index,
            id,
        ));
    }

    let mut detector = Vec::with_capacity(detector_ids.len());
    for id in detector_ids {
        detector.push(parse_detector(
            allowed_detector_schema,
            owner_rows,
            child_index,
            id,
        ));
    }

    Some(ComponentList {
        count: Some(source.len() + analyzer.len() + detector.len()),
        source,
        analyzer,
        detector,
    })
}

/// <source>
#[inline]
fn parse_source(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    source_id: u32,
) -> Source {
    let rows = rows_for_owner(owner_rows, source_id);

    let order = b000_attr_text(rows, ACC_ATTR_ORDER).and_then(|s| s.parse::<u32>().ok());

    let referenceable_param_group_ref =
        parse_referenceable_param_group_refs(owner_rows, child_index, source_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, source_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, user_param) = parse_params(allowed_schema, &params_meta);

    Source {
        order,
        referenceable_param_group_ref,
        cv_param,
        user_param,
    }
}

/// <analyzer>
#[inline]
fn parse_analyzer(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    analyzer_id: u32,
) -> Analyzer {
    let rows = rows_for_owner(owner_rows, analyzer_id);

    let order = b000_attr_text(rows, ACC_ATTR_ORDER).and_then(|s| s.parse::<u32>().ok());

    let referenceable_param_group_ref =
        parse_referenceable_param_group_refs(owner_rows, child_index, analyzer_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, analyzer_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, user_param) = parse_params(allowed_schema, &params_meta);

    Analyzer {
        order,
        referenceable_param_group_ref,
        cv_param,
        user_param,
    }
}

/// <detector>
#[inline]
fn parse_detector(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    detector_id: u32,
) -> Detector {
    let rows = rows_for_owner(owner_rows, detector_id);

    let order = b000_attr_text(rows, ACC_ATTR_ORDER).and_then(|s| s.parse::<u32>().ok());

    let referenceable_param_group_ref =
        parse_referenceable_param_group_refs(owner_rows, child_index, detector_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, detector_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, user_param) = parse_params(allowed_schema, &params_meta);

    Detector {
        order,
        referenceable_param_group_ref,
        cv_param,
        user_param,
    }
}

/// <instrumentSoftwareRef> / <softwareRef>
#[inline]
fn parse_instrument_software_ref(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    instrument_id: u32,
) -> Option<InstrumentSoftwareRef> {
    let software_ref_ids = ids_for_parent(
        owner_rows,
        child_index,
        metadata,
        instrument_id,
        TagId::SoftwareRef,
    );
    let software_ref_id = software_ref_ids.first().copied().unwrap_or(0);
    if software_ref_id == 0 {
        return None;
    }

    let rows = rows_for_owner(owner_rows, software_ref_id);
    let r#ref = b000_attr_text(rows, ACC_ATTR_REF)?;
    (!r#ref.is_empty()).then(|| InstrumentSoftwareRef { r#ref })
}

/// <referenceableParamGroupRef>
#[inline]
fn parse_referenceable_param_group_refs(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    parent_id: u32,
) -> Vec<ReferenceableParamGroupRef> {
    let ref_ids = unique_ids(child_index.ids(parent_id, TagId::ReferenceableParamGroupRef));
    if ref_ids.is_empty() {
        return Vec::new();
    }

    let mut out: Vec<ReferenceableParamGroupRef> = Vec::with_capacity(ref_ids.len());
    for ref_id in ref_ids {
        let ref_rows = rows_for_owner(owner_rows, ref_id);
        if let Some(r) = b000_attr_text(ref_rows, ACC_ATTR_REF) {
            if !r.is_empty() {
                out.push(ReferenceableParamGroupRef { r#ref: r });
            }
        }
    }

    out
}
