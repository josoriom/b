use std::collections::{HashMap, HashSet};

use crate::{
    Chromatogram, ChromatogramList,
    b64::utilities::{
        common::{
            ChildIndex, child_node, find_node_by_tag, get_attr_text, get_attr_u32,
            ordered_unique_owner_ids, xy_lengths_from_bdal,
        },
        parse_binary_data_array_list, parse_cv_and_user_params,
    },
    decode::Metadatum,
    mzml::{
        schema::{SchemaTree as Schema, TagId},
        structs::{Activation, IsolationWindow, Precursor, Product, SelectedIon, SelectedIonList},
    },
};

use crate::mzml::attr_meta::{
    ACC_ATTR_COUNT, ACC_ATTR_DATA_PROCESSING_REF, ACC_ATTR_DEFAULT_ARRAY_LENGTH,
    ACC_ATTR_DEFAULT_DATA_PROCESSING_REF, ACC_ATTR_EXTERNAL_SPECTRUM_ID, ACC_ATTR_ID,
    ACC_ATTR_INDEX, ACC_ATTR_NATIVE_ID, ACC_ATTR_SOURCE_FILE_REF, ACC_ATTR_SPECTRUM_REF,
};

/// <chromatogramList>
#[inline]
pub fn parse_chromatogram_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<ChromatogramList> {
    let list_node = find_node_by_tag(schema, TagId::ChromatogramList)?;
    let chromatogram_node = child_node(Some(list_node), TagId::Chromatogram)?;

    let allowed_chromatogram: HashSet<&str> = child_node(Some(chromatogram_node), TagId::CvParam)
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let precursor_node = child_node(Some(chromatogram_node), TagId::Precursor);
    let product_node = child_node(Some(chromatogram_node), TagId::Product);

    let allowed_iso_precursor = allowed_isolation_window_cv(precursor_node);
    let allowed_sel_ion = allowed_selected_ion_cv(precursor_node);
    let allowed_activation = allowed_activation_cv(precursor_node);
    let allowed_iso_product = allowed_isolation_window_cv(product_node);

    let chromatogram_list_rows: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::ChromatogramList)
        .collect();

    if chromatogram_list_rows.is_empty() {
        return None;
    }

    let default_data_processing_ref = get_attr_text(
        &chromatogram_list_rows,
        ACC_ATTR_DEFAULT_DATA_PROCESSING_REF,
    );

    let count_attr = get_attr_u32(&chromatogram_list_rows, ACC_ATTR_COUNT).map(|v| v as usize);

    let chromatogram_ids = ordered_unique_owner_ids(metadata, TagId::Chromatogram);
    if chromatogram_ids.is_empty() {
        return None;
    }

    let mut chromatogram_owner_to_item: HashMap<u32, u32> =
        HashMap::with_capacity(chromatogram_ids.len());
    for m in metadata {
        if m.tag_id == TagId::Chromatogram {
            chromatogram_owner_to_item
                .entry(m.owner_id)
                .or_insert(m.item_index);
        }
    }

    let mut want_items: HashSet<u32> = HashSet::with_capacity(chromatogram_ids.len());
    for &id in &chromatogram_ids {
        if let Some(&item_idx) = chromatogram_owner_to_item.get(&id) {
            want_items.insert(item_idx);
        }
    }

    let mut by_item_index: HashMap<u32, Vec<Metadatum>> = HashMap::with_capacity(want_items.len());
    for m in metadata {
        if want_items.contains(&m.item_index) {
            by_item_index
                .entry(m.item_index)
                .or_default()
                .push(m.clone());
        }
    }

    let mut chromatograms = Vec::with_capacity(chromatogram_ids.len());

    for (fallback_index, chromatogram_id) in chromatogram_ids.into_iter().enumerate() {
        let item_idx = match chromatogram_owner_to_item.get(&chromatogram_id).copied() {
            Some(v) => v,
            None => continue,
        };

        let item_meta = match by_item_index.get(&item_idx) {
            Some(v) => v,
            None => continue,
        };

        let mut has_other_root = false;
        for m in item_meta {
            if m.tag_id == TagId::Chromatogram && m.owner_id != chromatogram_id {
                has_other_root = true;
                break;
            }
        }

        if has_other_root {
            let ids = subtree_owner_ids(chromatogram_id, child_index);
            let keep: HashSet<u32> = ids.into_iter().collect();

            let mut scoped: Vec<Metadatum> = Vec::with_capacity(item_meta.len());
            for m in item_meta {
                if keep.contains(&m.owner_id) {
                    scoped.push(m.clone());
                }
            }

            chromatograms.push(parse_chromatogram(
                &scoped,
                chromatogram_id,
                child_index,
                fallback_index as u32,
                &allowed_chromatogram,
                default_data_processing_ref.as_deref(),
                &allowed_iso_precursor,
                &allowed_sel_ion,
                &allowed_activation,
                &allowed_iso_product,
            ));
        } else {
            chromatograms.push(parse_chromatogram(
                item_meta,
                chromatogram_id,
                child_index,
                fallback_index as u32,
                &allowed_chromatogram,
                default_data_processing_ref.as_deref(),
                &allowed_iso_precursor,
                &allowed_sel_ion,
                &allowed_activation,
                &allowed_iso_product,
            ));
        }
    }

    Some(ChromatogramList {
        count: count_attr.or(Some(chromatograms.len())),
        default_data_processing_ref,
        chromatograms,
    })
}

/// <chromatogram>
#[inline]
fn parse_chromatogram(
    metadata: &[Metadatum],
    chromatogram_id: u32,
    child_index: &ChildIndex,
    fallback_index: u32,
    allowed_chromatogram: &HashSet<&str>,
    default_data_processing_ref: Option<&str>,
    allowed_iso_precursor: &HashSet<&str>,
    allowed_sel_ion: &HashSet<&str>,
    allowed_activation: &HashSet<&str>,
    allowed_iso_product: &HashSet<&str>,
) -> Chromatogram {
    let mut chromatogram_rows: Vec<&Metadatum> = Vec::new();
    let mut chromatogram_params_meta: Vec<&Metadatum> = Vec::new();

    for m in metadata {
        if m.owner_id == chromatogram_id && m.tag_id == TagId::Chromatogram {
            chromatogram_rows.push(m);
            chromatogram_params_meta.push(m);
            continue;
        }

        if m.parent_index == chromatogram_id
            && (m.tag_id == TagId::CvParam || m.tag_id == TagId::UserParam)
        {
            chromatogram_params_meta.push(m);
        }
    }

    let id = get_attr_text(&chromatogram_rows, ACC_ATTR_ID).unwrap_or_default();
    let index = get_attr_u32(&chromatogram_rows, ACC_ATTR_INDEX).or(Some(fallback_index));
    let native_id = get_attr_text(&chromatogram_rows, ACC_ATTR_NATIVE_ID);

    let data_processing_ref = get_attr_text(&chromatogram_rows, ACC_ATTR_DATA_PROCESSING_REF)
        .or_else(|| default_data_processing_ref.map(|s| s.to_string()));

    let default_array_length_attr =
        get_attr_u32(&chromatogram_rows, ACC_ATTR_DEFAULT_ARRAY_LENGTH).map(|v| v as usize);

    let (cv_params, user_params) =
        parse_cv_and_user_params(allowed_chromatogram, &chromatogram_params_meta);

    let binary_data_array_list = parse_binary_data_array_list(metadata);

    let (x_len, y_len) = xy_lengths_from_bdal(binary_data_array_list.as_ref());
    let default_array_length: Option<usize> =
        default_array_length_attr.or(x_len).or(y_len).or(Some(0));

    let precursor = parse_precursor_for_chromatogram(
        metadata,
        chromatogram_id,
        child_index,
        allowed_iso_precursor,
        allowed_sel_ion,
        allowed_activation,
    );

    let product =
        parse_product_for_chromatogram(metadata, chromatogram_id, child_index, allowed_iso_product);

    Chromatogram {
        id,
        native_id,
        index,
        default_array_length,
        data_processing_ref,
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
        precursor,
        product,
        binary_data_array_list,
    }
}

/// <precursor>
#[inline]
fn parse_precursor_for_chromatogram(
    metadata: &[Metadatum],
    chromatogram_id: u32,
    child_index: &ChildIndex,
    allowed_iso: &HashSet<&str>,
    allowed_sel_ion: &HashSet<&str>,
    allowed_activation: &HashSet<&str>,
) -> Option<Precursor> {
    let precursor_id = child_index.first_id(chromatogram_id, TagId::Precursor)?;

    let precursor_rows: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::Precursor && m.owner_id == precursor_id)
        .collect();

    let spectrum_ref = get_attr_text(&precursor_rows, ACC_ATTR_SPECTRUM_REF);
    let source_file_ref = get_attr_text(&precursor_rows, ACC_ATTR_SOURCE_FILE_REF);
    let external_spectrum_id = get_attr_text(&precursor_rows, ACC_ATTR_EXTERNAL_SPECTRUM_ID);

    Some(Precursor {
        spectrum_ref,
        source_file_ref,
        external_spectrum_id,
        isolation_window: parse_isolation_window(metadata, precursor_id, child_index, allowed_iso),
        selected_ion_list: parse_selected_ion_list(
            metadata,
            precursor_id,
            child_index,
            allowed_sel_ion,
        ),
        activation: parse_activation(metadata, precursor_id, child_index, allowed_activation),
    })
}

/// <product>
#[inline]
fn parse_product_for_chromatogram(
    metadata: &[Metadatum],
    chromatogram_id: u32,
    child_index: &ChildIndex,
    allowed_iso: &HashSet<&str>,
) -> Option<Product> {
    let product_id = child_index.first_id(chromatogram_id, TagId::Product)?;

    Some(Product {
        spectrum_ref: None,
        source_file_ref: None,
        external_spectrum_id: None,
        isolation_window: parse_isolation_window(metadata, product_id, child_index, allowed_iso),
    })
}

/// <isolationWindow>
#[inline]
fn parse_isolation_window(
    metadata: &[Metadatum],
    parent_id: u32,
    child_index: &ChildIndex,
    allowed_iso: &HashSet<&str>,
) -> Option<IsolationWindow> {
    let isolation_id = child_index.first_id(parent_id, TagId::IsolationWindow)?;

    let iso_params_meta: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| {
            m.parent_index == isolation_id
                && (m.tag_id == TagId::CvParam || m.tag_id == TagId::UserParam)
        })
        .collect();

    let (cv_params, user_params) = parse_cv_and_user_params(allowed_iso, &iso_params_meta);

    Some(IsolationWindow {
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
    })
}

/// <selectedIonList>
#[inline]
fn parse_selected_ion_list(
    metadata: &[Metadatum],
    precursor_id: u32,
    child_index: &ChildIndex,
    allowed_sel_ion: &HashSet<&str>,
) -> Option<SelectedIonList> {
    let list_id = child_index.first_id(precursor_id, TagId::SelectedIonList)?;

    let list_rows: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::SelectedIonList && m.owner_id == list_id)
        .collect();

    let count_attr = get_attr_u32(&list_rows, ACC_ATTR_COUNT).map(|v| v as usize);

    let mut selected_ions = Vec::new();
    let mut seen = HashSet::new();

    for &ion_id in child_index.ids(list_id, TagId::SelectedIon) {
        if !seen.insert(ion_id) {
            continue;
        }

        let ion_params_meta: Vec<&Metadatum> = metadata
            .iter()
            .filter(|m| {
                m.parent_index == ion_id
                    && (m.tag_id == TagId::CvParam || m.tag_id == TagId::UserParam)
            })
            .collect();

        let (cv_params, user_params) = parse_cv_and_user_params(allowed_sel_ion, &ion_params_meta);

        selected_ions.push(SelectedIon {
            referenceable_param_group_refs: Vec::new(),
            cv_params,
            user_params,
        });
    }

    Some(SelectedIonList {
        count: count_attr.or(Some(selected_ions.len())),
        selected_ions,
    })
}

/// <activation>
#[inline]
fn parse_activation(
    metadata: &[Metadatum],
    precursor_id: u32,
    child_index: &ChildIndex,
    allowed_activation: &HashSet<&str>,
) -> Option<Activation> {
    let activation_id = child_index.first_id(precursor_id, TagId::Activation)?;

    let activation_params_meta: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| {
            m.parent_index == activation_id
                && (m.tag_id == TagId::CvParam || m.tag_id == TagId::UserParam)
        })
        .collect();

    let (cv_params, user_params) =
        parse_cv_and_user_params(allowed_activation, &activation_params_meta);

    Some(Activation {
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
    })
}

#[inline]
fn allowed_isolation_window_cv(parent: Option<&crate::mzml::schema::SchemaNode>) -> HashSet<&str> {
    parent
        .and_then(|p| child_node(Some(p), TagId::IsolationWindow))
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default()
}

#[inline]
fn allowed_selected_ion_cv(parent: Option<&crate::mzml::schema::SchemaNode>) -> HashSet<&str> {
    parent
        .and_then(|p| child_node(Some(p), TagId::SelectedIonList))
        .and_then(|n| child_node(Some(n), TagId::SelectedIon))
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default()
}

#[inline]
fn allowed_activation_cv(parent: Option<&crate::mzml::schema::SchemaNode>) -> HashSet<&str> {
    parent
        .and_then(|p| child_node(Some(p), TagId::Activation))
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default()
}

#[inline]
fn subtree_owner_ids(root_id: u32, child_index: &ChildIndex) -> Vec<u32> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    let mut stack = vec![root_id];

    while let Some(id) = stack.pop() {
        if !seen.insert(id) {
            continue;
        }

        out.push(id);

        for &child_id in child_index.children(id) {
            stack.push(child_id);
        }
    }

    out
}
