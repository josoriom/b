use std::collections::{HashMap, HashSet};

use crate::{
    Activation, IsolationWindow, Precursor, PrecursorList, SelectedIon, SelectedIonList,
    b64::utilities::{
        common::{
            ChildIndex, child_node, find_node_by_tag, get_attr_text, ordered_unique_owner_ids,
        },
        parse_cv_and_user_params,
    },
    decode::Metadatum,
    mzml::{
        attr_meta::{
            ACC_ATTR_EXTERNAL_SPECTRUM_ID, ACC_ATTR_SOURCE_FILE_REF, ACC_ATTR_SPECTRUM_REF,
        },
        schema::{SchemaTree as Schema, TagId},
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
    parse_cv_and_user_params(allowed_schema, rows)
}

#[inline]
fn params_for_tagged_owner<'a>(
    metas_by_owner: &HashMap<u32, Vec<&'a Metadatum>>,
    child_index: &ChildIndex,
    owner_id: u32,
    tag: TagId,
) -> Vec<&'a Metadatum> {
    let mut out: Vec<&'a Metadatum> = Vec::new();

    if let Some(xs) = metas_by_owner.get(&owner_id) {
        for &m in xs {
            if m.tag_id == tag {
                out.push(m);
            }
        }
    }

    let cv_ids = child_index.ids(owner_id, TagId::CvParam);
    let up_ids = child_index.ids(owner_id, TagId::UserParam);
    let rpg_ids = child_index.ids(owner_id, TagId::ReferenceableParamGroupRef);

    out.reserve(cv_ids.len() + up_ids.len() + rpg_ids.len());

    for &id in cv_ids {
        if let Some(xs) = metas_by_owner.get(&id) {
            out.extend(xs.iter().copied());
        }
    }
    for &id in up_ids {
        if let Some(xs) = metas_by_owner.get(&id) {
            out.extend(xs.iter().copied());
        }
    }
    for &id in rpg_ids {
        if let Some(xs) = metas_by_owner.get(&id) {
            out.extend(xs.iter().copied());
        }
    }

    out
}

/// <precursorList>
#[inline]
pub fn parse_precursor_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<PrecursorList> {
    let list_node = find_node_by_tag(schema, TagId::PrecursorList)?;
    let precursor_node = child_node(Some(list_node), TagId::Precursor)?;

    let allowed_isolation_window: HashSet<&str> =
        child_node(Some(precursor_node), TagId::IsolationWindow)
            .and_then(|n| child_node(Some(n), TagId::CvParam))
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let allowed_selected_ion: HashSet<&str> =
        child_node(Some(precursor_node), TagId::SelectedIonList)
            .and_then(|n| child_node(Some(n), TagId::SelectedIon))
            .and_then(|n| child_node(Some(n), TagId::CvParam))
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let allowed_activation: HashSet<&str> = child_node(Some(precursor_node), TagId::Activation)
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut metas_by_owner: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    let mut parent_by_owner: HashMap<u32, u32> = HashMap::with_capacity(metadata.len());

    for m in metadata {
        metas_by_owner.entry(m.owner_id).or_default().push(m);
        parent_by_owner.entry(m.owner_id).or_insert(m.parent_index);
    }

    let mut precursor_ids: Vec<u32> = Vec::new();
    if let Some(list_id) = metadata
        .iter()
        .find(|m| m.tag_id == TagId::PrecursorList)
        .map(|m| m.owner_id)
    {
        precursor_ids.extend_from_slice(child_index.ids(list_id, TagId::Precursor));
    }

    if precursor_ids.is_empty() {
        precursor_ids = ordered_unique_owner_ids(metadata, TagId::Precursor);
    }

    if precursor_ids.is_empty() {
        let mut set: HashSet<u32> = HashSet::new();
        for m in metadata {
            match m.tag_id {
                TagId::IsolationWindow | TagId::SelectedIon | TagId::Activation => {
                    if m.parent_index != 0 {
                        set.insert(m.parent_index);
                    }
                }
                _ => {}
            }
        }
        if !set.is_empty() {
            precursor_ids = set.into_iter().collect();
            precursor_ids.sort_unstable();
        }
    }

    if precursor_ids.is_empty() {
        return None;
    }

    let mut precursors = Vec::with_capacity(precursor_ids.len());
    for precursor_id in precursor_ids {
        precursors.push(parse_precursor(
            precursor_id,
            &allowed_isolation_window,
            &allowed_selected_ion,
            &allowed_activation,
            &metas_by_owner,
            &parent_by_owner,
            child_index,
        ));
    }

    Some(PrecursorList {
        count: Some(precursors.len()),
        cv_params: Vec::new(),
        user_params: Vec::new(),
        precursors,
    })
}

/// <precursor>
#[inline]
fn parse_precursor(
    precursor_id: u32,
    allowed_isolation_window: &HashSet<&str>,
    allowed_selected_ion: &HashSet<&str>,
    allowed_activation: &HashSet<&str>,
    metas_by_owner: &HashMap<u32, Vec<&Metadatum>>,
    parent_by_owner: &HashMap<u32, u32>,
    child_index: &ChildIndex,
) -> Precursor {
    let precursor_parent = parent_by_owner.get(&precursor_id).copied().unwrap_or(0);

    let rows: &[&Metadatum] = metas_by_owner
        .get(&precursor_id)
        .map(|v| v.as_slice())
        .unwrap_or(&[]);

    let spectrum_ref = get_attr_text(rows, ACC_ATTR_SPECTRUM_REF);
    let source_file_ref = get_attr_text(rows, ACC_ATTR_SOURCE_FILE_REF);
    let external_spectrum_id = get_attr_text(rows, ACC_ATTR_EXTERNAL_SPECTRUM_ID);

    let isolation_window = parse_isolation_window(
        allowed_isolation_window,
        metas_by_owner,
        child_index,
        precursor_id,
        precursor_parent,
    );

    let selected_ion_list = parse_selected_ion_list(
        allowed_selected_ion,
        metas_by_owner,
        child_index,
        precursor_id,
        precursor_parent,
    );

    let activation = parse_activation(
        allowed_activation,
        metas_by_owner,
        child_index,
        precursor_id,
        precursor_parent,
    );

    Precursor {
        spectrum_ref,
        source_file_ref,
        external_spectrum_id,
        isolation_window,
        selected_ion_list,
        activation,
    }
}

/// <isolationWindow>
#[inline]
fn parse_isolation_window(
    allowed_isolation_window: &HashSet<&str>,
    metas_by_owner: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    precursor_id: u32,
    precursor_parent: u32,
) -> Option<IsolationWindow> {
    let isolation_id = child_index
        .first_id(precursor_id, TagId::IsolationWindow)
        .or_else(|| child_index.first_id(precursor_parent, TagId::IsolationWindow));

    let mut meta: Vec<&Metadatum> = Vec::new();

    if let Some(iw_id) = isolation_id {
        meta = params_for_tagged_owner(metas_by_owner, child_index, iw_id, TagId::IsolationWindow);
    }
    if meta.is_empty() {
        meta = params_for_tagged_owner(
            metas_by_owner,
            child_index,
            precursor_id,
            TagId::IsolationWindow,
        );
    }
    if meta.is_empty() && precursor_parent != 0 && precursor_parent != precursor_id {
        meta = params_for_tagged_owner(
            metas_by_owner,
            child_index,
            precursor_parent,
            TagId::IsolationWindow,
        );
    }

    if meta.is_empty() && isolation_id.is_none() {
        return None;
    }

    let (cv_params, user_params) = parse_params(allowed_isolation_window, &meta);

    Some(IsolationWindow {
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
    })
}

/// <selectedIonList>
#[inline]
fn parse_selected_ion_list(
    allowed_selected_ion: &HashSet<&str>,
    metas_by_owner: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    precursor_id: u32,
    precursor_parent: u32,
) -> Option<SelectedIonList> {
    let mut selected_ion_ids: Vec<u32> = Vec::new();

    if let Some(list_id) = child_index.first_id(precursor_id, TagId::SelectedIonList) {
        selected_ion_ids.extend_from_slice(child_index.ids(list_id, TagId::SelectedIon));
    }
    if selected_ion_ids.is_empty() {
        selected_ion_ids.extend_from_slice(child_index.ids(precursor_id, TagId::SelectedIon));
    }
    if selected_ion_ids.is_empty() {
        if let Some(list_id) = child_index.first_id(precursor_parent, TagId::SelectedIonList) {
            selected_ion_ids.extend_from_slice(child_index.ids(list_id, TagId::SelectedIon));
        }
    }
    if selected_ion_ids.is_empty() {
        selected_ion_ids.extend_from_slice(child_index.ids(precursor_parent, TagId::SelectedIon));
    }

    if selected_ion_ids.is_empty() {
        return None;
    }

    selected_ion_ids.sort_unstable();
    selected_ion_ids.dedup();

    let mut selected_ions = Vec::with_capacity(selected_ion_ids.len());
    for sid in selected_ion_ids {
        let meta = params_for_tagged_owner(metas_by_owner, child_index, sid, TagId::SelectedIon);
        let (cv_params, user_params) = parse_params(allowed_selected_ion, &meta);

        selected_ions.push(SelectedIon {
            referenceable_param_group_refs: Vec::new(),
            cv_params,
            user_params,
        });
    }

    Some(SelectedIonList {
        count: Some(selected_ions.len()),
        selected_ions,
    })
}

/// <activation>
#[inline]
fn parse_activation(
    allowed_activation: &HashSet<&str>,
    metas_by_owner: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    precursor_id: u32,
    precursor_parent: u32,
) -> Option<Activation> {
    let activation_id = child_index
        .first_id(precursor_id, TagId::Activation)
        .or_else(|| child_index.first_id(precursor_parent, TagId::Activation));

    let mut meta: Vec<&Metadatum> = Vec::new();

    if let Some(act_id) = activation_id {
        meta = params_for_tagged_owner(metas_by_owner, child_index, act_id, TagId::Activation);
    }
    if meta.is_empty() {
        meta =
            params_for_tagged_owner(metas_by_owner, child_index, precursor_id, TagId::Activation);
    }
    if meta.is_empty() && precursor_parent != 0 && precursor_parent != precursor_id {
        meta = params_for_tagged_owner(
            metas_by_owner,
            child_index,
            precursor_parent,
            TagId::Activation,
        );
    }

    if meta.is_empty() && activation_id.is_none() {
        return None;
    }

    let (cv_params, user_params) = parse_params(allowed_activation, &meta);

    Some(Activation {
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
    })
}
