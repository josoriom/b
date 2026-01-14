use std::collections::{HashMap, HashSet};

use crate::{
    Spectrum, SpectrumList,
    b64::decode2::Metadatum,
    b64::utilities::{
        common::{
            ChildIndex, child_node, find_node_by_tag, get_attr_text, get_attr_u32,
            xy_lengths_from_bdal,
        },
        parse_binary_data_array_list, parse_cv_and_user_params, parse_precursor_list,
        parse_product_list, parse_scan_list,
    },
    mzml::{
        schema::{SchemaTree as Schema, TagId},
        structs::BinaryDataArrayList,
    },
};

use crate::mzml::attr_meta::{
    ACC_ATTR_COUNT, ACC_ATTR_DATA_PROCESSING_REF, ACC_ATTR_DEFAULT_ARRAY_LENGTH,
    ACC_ATTR_DEFAULT_DATA_PROCESSING_REF, ACC_ATTR_ID, ACC_ATTR_INDEX, ACC_ATTR_MS_LEVEL,
    ACC_ATTR_NATIVE_ID, ACC_ATTR_SCAN_NUMBER, ACC_ATTR_SOURCE_FILE_REF, ACC_ATTR_SPOT_ID,
};

#[inline]
pub fn parse_spectrum_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<SpectrumList> {
    // <spectrumList>/<spectrum>/<cvParam>
    let allowed_spectrum: HashSet<&str> = find_node_by_tag(schema, TagId::SpectrumList)
        .and_then(|n| child_node(Some(n), TagId::Spectrum))
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let spectrum_list_rows: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::SpectrumList)
        .collect();

    let all_rows: Vec<&Metadatum> = metadata.iter().collect();

    let default_data_processing_ref =
        get_attr_text(&spectrum_list_rows, ACC_ATTR_DEFAULT_DATA_PROCESSING_REF)
            .or_else(|| get_attr_text(&spectrum_list_rows, ACC_ATTR_DATA_PROCESSING_REF))
            .or_else(|| get_attr_text(&all_rows, ACC_ATTR_DEFAULT_DATA_PROCESSING_REF))
            .or_else(|| get_attr_text(&all_rows, ACC_ATTR_DATA_PROCESSING_REF));

    let count_attr = get_attr_u32(&spectrum_list_rows, ACC_ATTR_COUNT)
        .or_else(|| get_attr_u32(&all_rows, ACC_ATTR_COUNT))
        .map(|v| v as usize);

    let spectrum_rows: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::Spectrum)
        .collect();

    if spectrum_rows.is_empty() {
        return None;
    }

    let mut spectrum_owner_to_item: HashMap<u32, u32> = HashMap::with_capacity(spectrum_rows.len());
    let mut spectrum_item_in_order: Vec<u32> = Vec::with_capacity(spectrum_rows.len());
    let mut seen_item: HashSet<u32> = HashSet::with_capacity(spectrum_rows.len());

    for m in &spectrum_rows {
        spectrum_owner_to_item.insert(m.owner_id, m.item_index);
        if seen_item.insert(m.item_index) {
            spectrum_item_in_order.push(m.item_index);
        }
    }

    let spectrum_item_indices: Vec<u32> = if spectrum_list_rows.is_empty() {
        spectrum_item_in_order.clone()
    } else {
        let spectrum_list_id = spectrum_list_rows[0].owner_id;
        let direct = child_index.ids(spectrum_list_id, TagId::Spectrum);

        if direct.is_empty() {
            spectrum_item_in_order.clone()
        } else {
            let mut out = Vec::with_capacity(direct.len());
            let mut seen = HashSet::with_capacity(direct.len());

            for &spectrum_id in direct {
                if let Some(&item_index) = spectrum_owner_to_item.get(&spectrum_id) {
                    if seen.insert(item_index) {
                        out.push(item_index);
                    }
                }
            }

            if out.is_empty() {
                spectrum_item_in_order.clone()
            } else {
                out
            }
        }
    };

    let want_items: HashSet<u32> = spectrum_item_indices.iter().copied().collect();
    let mut by_item_index: HashMap<u32, Vec<Metadatum>> = HashMap::with_capacity(want_items.len());

    for m in metadata {
        if want_items.contains(&m.item_index) {
            by_item_index
                .entry(m.item_index)
                .or_default()
                .push(m.clone());
        }
    }

    let mut spectra = Vec::with_capacity(spectrum_item_indices.len());

    for (fallback_index, item_index) in spectrum_item_indices.into_iter().enumerate() {
        let spectrum_meta = match by_item_index.get(&item_index) {
            Some(v) => v,
            None => continue,
        };

        let local_child_index = ChildIndex::new(spectrum_meta);

        spectra.push(parse_spectrum(
            schema,
            spectrum_meta,
            fallback_index as u32,
            &allowed_spectrum,
            default_data_processing_ref.as_deref(),
            &local_child_index,
        ));
    }

    if spectra.is_empty() {
        return None;
    }

    Some(SpectrumList {
        count: count_attr.or(Some(spectra.len())),
        default_data_processing_ref,
        spectra,
    })
}

#[inline]
fn parse_spectrum(
    schema: &Schema,
    metadata: &[Metadatum],
    fallback_index: u32,
    allowed_spectrum: &HashSet<&str>,
    default_data_processing_ref: Option<&str>,
    child_index: &ChildIndex,
) -> Spectrum {
    // <spectrum>
    let spectrum_rows: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::Spectrum)
        .collect();

    let spectrum_id = spectrum_rows.first().map(|m| m.owner_id).unwrap_or(0);

    let id = get_attr_text(&spectrum_rows, ACC_ATTR_ID).unwrap_or_default();
    let index = get_attr_u32(&spectrum_rows, ACC_ATTR_INDEX).or(Some(fallback_index));

    let scan_number = get_attr_u32(&spectrum_rows, ACC_ATTR_SCAN_NUMBER);
    let ms_level = get_attr_u32(&spectrum_rows, ACC_ATTR_MS_LEVEL);

    let native_id = get_attr_text(&spectrum_rows, ACC_ATTR_NATIVE_ID);
    let source_file_ref = get_attr_text(&spectrum_rows, ACC_ATTR_SOURCE_FILE_REF);
    let spot_id = get_attr_text(&spectrum_rows, ACC_ATTR_SPOT_ID);

    let data_processing_ref = get_attr_text(&spectrum_rows, ACC_ATTR_DATA_PROCESSING_REF)
        .or_else(|| default_data_processing_ref.map(|s| s.to_string()));

    let spectrum_params_meta: Vec<&Metadatum> = metadata
        .iter()
        .filter(|m| m.tag_id == TagId::Spectrum)
        .filter(|m| m.owner_id == spectrum_id)
        .filter(|m| {
            !m.accession
                .as_deref()
                .is_some_and(|a| a.starts_with("B000:"))
        })
        .collect();

    let (cv_params, user_params) = if allowed_spectrum.is_empty() {
        let mut allow_all: HashSet<&str> = HashSet::new();
        allow_all.insert("");
        parse_cv_and_user_params(&allow_all, &spectrum_params_meta)
    } else {
        parse_cv_and_user_params(allowed_spectrum, &spectrum_params_meta)
    };

    let scan_list = parse_scan_list(schema, metadata, child_index);
    let product_list = parse_product_list(schema, metadata, child_index);
    let precursor_list = parse_precursor_list(schema, metadata, child_index);

    let binary_data_array_list: Option<BinaryDataArrayList> =
        parse_binary_data_array_list(metadata);

    let (x_len, y_len) = xy_lengths_from_bdal(binary_data_array_list.as_ref());

    let default_array_length_attr =
        get_attr_u32(&spectrum_rows, ACC_ATTR_DEFAULT_ARRAY_LENGTH).map(|v| v as usize);

    let default_array_length = default_array_length_attr.or(x_len).or(y_len).or(Some(0));

    Spectrum {
        id,
        index,
        scan_number,
        default_array_length,
        native_id,
        data_processing_ref,
        source_file_ref,
        spot_id,
        ms_level,
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
        spectrum_description: None,
        scan_list,
        precursor_list,
        product_list,
        binary_data_array_list,
    }
}
