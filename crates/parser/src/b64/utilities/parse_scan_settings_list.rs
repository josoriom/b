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
        attr_meta::{ACC_ATTR_ID, ACC_ATTR_INSTRUMENT_CONFIGURATION_REF, ACC_ATTR_REF},
        schema::{SchemaTree as Schema, TagId},
        structs::{
            ReferenceableParamGroupRef, ScanSettings, ScanSettingsList, SourceFileRef,
            SourceFileRefList, Target, TargetList,
        },
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

/// <scanSettingsList> / <acquisitionSettingsList>
#[inline]
pub fn parse_scan_settings_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<ScanSettingsList> {
    let list_node = find_node_by_tag(schema, TagId::ScanSettingsList)
        .or_else(|| find_node_by_tag(schema, TagId::AcquisitionSettingsList))?;

    let settings_node = child_node(Some(list_node), TagId::ScanSettings)
        .or_else(|| child_node(Some(list_node), TagId::AcquisitionSettings))?;

    let allowed_scan_settings_schema: HashSet<&str> =
        child_node(Some(settings_node), TagId::CvParam)
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let allowed_target_schema: HashSet<&str> = child_node(Some(settings_node), TagId::TargetList)
        .and_then(|n| child_node(Some(n), TagId::Target))
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::ScanSettingsList || m.tag_id == TagId::AcquisitionSettingsList)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::ScanSettings || m.tag_id == TagId::AcquisitionSettings)
                .map(|m| m.parent_index)
        })?;

    let mut scan_settings_ids: Vec<u32> = Vec::new();
    scan_settings_ids.extend_from_slice(child_index.ids(list_id, TagId::ScanSettings));
    scan_settings_ids.extend_from_slice(child_index.ids(list_id, TagId::AcquisitionSettings));
    scan_settings_ids = unique_ids(&scan_settings_ids);

    if scan_settings_ids.is_empty() {
        scan_settings_ids = ordered_unique_owner_ids(metadata, TagId::ScanSettings);
        scan_settings_ids.extend(ordered_unique_owner_ids(
            metadata,
            TagId::AcquisitionSettings,
        ));
        scan_settings_ids.retain(|&id| is_child_of(&owner_rows, id, list_id));
        scan_settings_ids.sort_unstable();
        scan_settings_ids.dedup();
    }

    if scan_settings_ids.is_empty() {
        return None;
    }

    let mut scan_settings = Vec::with_capacity(scan_settings_ids.len());
    for id in scan_settings_ids {
        scan_settings.push(parse_scan_settings(
            &allowed_scan_settings_schema,
            &allowed_target_schema,
            &owner_rows,
            child_index,
            metadata,
            id,
        ));
    }

    Some(ScanSettingsList {
        count: Some(scan_settings.len()),
        scan_settings,
    })
}

/// <scanSettings> / <acquisitionSettings>
#[inline]
fn parse_scan_settings(
    allowed_scan_settings_schema: &HashSet<&str>,
    allowed_target_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    scan_settings_id: u32,
) -> ScanSettings {
    let rows = rows_for_owner(owner_rows, scan_settings_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).filter(|s| !s.is_empty());
    let instrument_configuration_ref =
        b000_attr_text(rows, ACC_ATTR_INSTRUMENT_CONFIGURATION_REF).filter(|s| !s.is_empty());

    let referenceable_param_group_refs =
        parse_referenceable_param_group_refs(owner_rows, child_index, scan_settings_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, scan_settings_id);
    let mut params_meta = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_params, user_params) = parse_params(allowed_scan_settings_schema, &params_meta);

    let source_file_ref_list =
        parse_source_file_ref_list(owner_rows, child_index, metadata, scan_settings_id);

    let target_list = parse_target_list(
        allowed_target_schema,
        owner_rows,
        child_index,
        metadata,
        scan_settings_id,
    );

    ScanSettings {
        id,
        instrument_configuration_ref,
        referenceable_param_group_refs,
        cv_params,
        user_params,
        source_file_ref_list,
        target_list,
    }
}

/// <sourceFileRefList> OR (mzML 0.99.x) <sourceFileList>
#[inline]
fn parse_source_file_ref_list(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    scan_settings_id: u32,
) -> Option<SourceFileRefList> {
    let mut list_ids = unique_ids(child_index.ids(scan_settings_id, TagId::SourceFileRefList));
    if list_ids.is_empty() {
        list_ids = ordered_unique_owner_ids(metadata, TagId::SourceFileRefList);
        list_ids.retain(|&id| is_child_of(owner_rows, id, scan_settings_id));
    }

    if let Some(list_id) = list_ids.first().copied() {
        let mut ref_ids = unique_ids(child_index.ids(list_id, TagId::SourceFileRef));
        if ref_ids.is_empty() {
            ref_ids = ordered_unique_owner_ids(metadata, TagId::SourceFileRef);
            ref_ids.retain(|&id| is_child_of(owner_rows, id, list_id));
        }

        if !ref_ids.is_empty() {
            let mut source_file_refs = Vec::with_capacity(ref_ids.len());
            for id in ref_ids {
                let rows = rows_for_owner(owner_rows, id);
                if let Some(r) = b000_attr_text(rows, ACC_ATTR_REF) {
                    if !r.is_empty() {
                        source_file_refs.push(SourceFileRef { r#ref: r });
                    }
                }
            }

            if !source_file_refs.is_empty() {
                return Some(SourceFileRefList {
                    count: Some(source_file_refs.len()),
                    source_file_refs,
                });
            }
        }
    }

    let mut sf_list_ids = unique_ids(child_index.ids(scan_settings_id, TagId::SourceFileList));
    if sf_list_ids.is_empty() {
        sf_list_ids = ordered_unique_owner_ids(metadata, TagId::SourceFileList);
        sf_list_ids.retain(|&id| is_child_of(owner_rows, id, scan_settings_id));
    }

    let list_id = sf_list_ids.first().copied()?;

    let mut sf_ids = unique_ids(child_index.ids(list_id, TagId::SourceFile));
    if sf_ids.is_empty() {
        sf_ids = ordered_unique_owner_ids(metadata, TagId::SourceFile);
        sf_ids.retain(|&id| is_child_of(owner_rows, id, list_id));
    }

    if sf_ids.is_empty() {
        return None;
    }

    let mut source_file_refs = Vec::with_capacity(sf_ids.len());
    for id in sf_ids {
        let rows = rows_for_owner(owner_rows, id);
        if let Some(sf_id) = b000_attr_text(rows, ACC_ATTR_ID) {
            if !sf_id.is_empty() {
                source_file_refs.push(SourceFileRef { r#ref: sf_id });
            }
        }
    }

    (!source_file_refs.is_empty()).then(|| SourceFileRefList {
        count: Some(source_file_refs.len()),
        source_file_refs,
    })
}

/// <targetList>
#[inline]
fn parse_target_list(
    allowed_target_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    scan_settings_id: u32,
) -> Option<TargetList> {
    let mut list_ids = unique_ids(child_index.ids(scan_settings_id, TagId::TargetList));
    if list_ids.is_empty() {
        list_ids = ordered_unique_owner_ids(metadata, TagId::TargetList);
        list_ids.retain(|&id| is_child_of(owner_rows, id, scan_settings_id));
    }

    let mut target_ids: Vec<u32> = Vec::new();

    if let Some(list_id) = list_ids.first().copied() {
        target_ids = unique_ids(child_index.ids(list_id, TagId::Target));
        if target_ids.is_empty() {
            target_ids = ordered_unique_owner_ids(metadata, TagId::Target);
            target_ids.retain(|&id| is_child_of(owner_rows, id, list_id));
        }
    }

    if target_ids.is_empty() {
        target_ids = unique_ids(child_index.ids(scan_settings_id, TagId::Target));
        if target_ids.is_empty() {
            target_ids = ordered_unique_owner_ids(metadata, TagId::Target);
            target_ids.retain(|&id| is_child_of(owner_rows, id, scan_settings_id));
        }
    }

    if !target_ids.is_empty() {
        let mut targets = Vec::with_capacity(target_ids.len());
        for id in target_ids {
            targets.push(parse_target(
                allowed_target_schema,
                owner_rows,
                child_index,
                id,
            ));
        }
        return Some(TargetList {
            count: Some(targets.len()),
            targets,
        });
    }

    let mut cv_ids = ordered_unique_owner_ids(metadata, TagId::CvParam);
    cv_ids.retain(|&id| is_child_of(owner_rows, id, scan_settings_id));
    if cv_ids.is_empty() {
        return None;
    }

    let mut target_cv: Vec<crate::mzml::structs::CvParam> = Vec::new();
    for cv_id in cv_ids {
        let rows = rows_for_owner(owner_rows, cv_id);
        let (mut cv_params, _user_params) = parse_cv_and_user_params(allowed_target_schema, rows);
        target_cv.append(&mut cv_params);
    }

    target_cv.retain(|p| {
        p.accession
            .as_deref()
            .map(|a| {
                a.ends_with("1000827")
                    || a.ends_with("1001225")
                    || a.ends_with("1000502")
                    || a.ends_with("1000747")
            })
            .unwrap_or(false)
    });

    if target_cv.is_empty() {
        return None;
    }

    let is_boundary = |p: &crate::mzml::structs::CvParam| {
        p.accession
            .as_deref()
            .map(|a| a.ends_with("1000827"))
            .unwrap_or(false)
    };

    let mut targets: Vec<Target> = Vec::new();
    let mut cur: Vec<crate::mzml::structs::CvParam> = Vec::new();

    for p in target_cv {
        if is_boundary(&p) && !cur.is_empty() {
            targets.push(Target {
                referenceable_param_group_refs: Vec::new(),
                cv_params: cur,
                user_params: Vec::new(),
            });
            cur = Vec::new();
        }
        cur.push(p);
    }

    if !cur.is_empty() {
        targets.push(Target {
            referenceable_param_group_refs: Vec::new(),
            cv_params: cur,
            user_params: Vec::new(),
        });
    }

    (!targets.is_empty()).then(|| TargetList {
        count: Some(targets.len()),
        targets,
    })
}

/// <target>
#[inline]
fn parse_target(
    allowed_target_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    target_id: u32,
) -> Target {
    let rows = rows_for_owner(owner_rows, target_id);

    let referenceable_param_group_refs =
        parse_referenceable_param_group_refs(owner_rows, child_index, target_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, target_id);
    let mut params_meta = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_params, user_params) = parse_params(allowed_target_schema, &params_meta);

    Target {
        referenceable_param_group_refs,
        cv_params,
        user_params,
    }
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
