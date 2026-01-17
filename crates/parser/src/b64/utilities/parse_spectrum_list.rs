use std::collections::{HashMap, HashSet};

use crate::{
    Spectrum, SpectrumDescription, SpectrumList,
    b64::{
        decode::Metadatum,
        utilities::{
            common::{
                ChildIndex, child_node, find_node_by_tag, get_attr_text, get_attr_u32,
                xy_lengths_from_bdal,
            },
            parse_binary_data_array_list, parse_cv_and_user_params, parse_precursor_list,
            parse_product_list, parse_scan_list,
        },
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
fn parse_params(
    allowed_schema: &HashSet<&str>,
    rows: &[&Metadatum],
) -> (
    Vec<crate::mzml::structs::CvParam>,
    Vec<crate::mzml::structs::UserParam>,
) {
    if allowed_schema.is_empty() {
        let allowed_meta = crate::b64::utilities::parse_file_description::allowed_from_rows(rows);
        parse_cv_and_user_params(&allowed_meta, rows)
    } else {
        parse_cv_and_user_params(allowed_schema, rows)
    }
}

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

    let mut all_rows: Vec<&Metadatum> = Vec::with_capacity(metadata.len());
    let mut spectrum_list_rows: Vec<&Metadatum> = Vec::new();
    let mut spectrum_rows: Vec<&Metadatum> = Vec::new();

    for m in metadata {
        all_rows.push(m);
        match m.tag_id {
            TagId::SpectrumList => spectrum_list_rows.push(m),
            TagId::Spectrum => spectrum_rows.push(m),
            _ => {}
        }
    }

    let default_data_processing_ref =
        get_attr_text(&spectrum_list_rows, ACC_ATTR_DEFAULT_DATA_PROCESSING_REF)
            .or_else(|| get_attr_text(&spectrum_list_rows, ACC_ATTR_DATA_PROCESSING_REF))
            .or_else(|| get_attr_text(&all_rows, ACC_ATTR_DEFAULT_DATA_PROCESSING_REF))
            .or_else(|| get_attr_text(&all_rows, ACC_ATTR_DATA_PROCESSING_REF));

    let count_attr = get_attr_u32(&spectrum_list_rows, ACC_ATTR_COUNT)
        .or_else(|| get_attr_u32(&all_rows, ACC_ATTR_COUNT))
        .map(|v| v as usize);

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

    let n = spectrum_item_indices.len();
    if n == 0 {
        return None;
    }

    let mut pos_by_item: HashMap<u32, usize> = HashMap::with_capacity(n);
    for (pos, item_index) in spectrum_item_indices.iter().copied().enumerate() {
        pos_by_item.insert(item_index, pos);
    }

    let mut buckets: Vec<Vec<Metadatum>> = vec![Vec::new(); n];
    for m in metadata {
        if let Some(&pos) = pos_by_item.get(&m.item_index) {
            buckets[pos].push(m.clone());
        }
    }

    let mut spectra = Vec::with_capacity(n);

    for (pos, spectrum_meta) in buckets.iter().enumerate() {
        if spectrum_meta.is_empty() {
            continue;
        }
        let local_child_index = ChildIndex::new(spectrum_meta);
        spectra.push(parse_spectrum(
            schema,
            spectrum_meta,
            pos as u32,
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
    let mut spectrum_rows: Vec<&Metadatum> = Vec::new();
    let mut spectrum_id: u32 = 0;

    for m in metadata {
        if m.tag_id == TagId::Spectrum {
            if spectrum_id == 0 {
                spectrum_id = m.owner_id;
            }
            spectrum_rows.push(m);
        }
    }

    let id = get_attr_text(&spectrum_rows, ACC_ATTR_ID).unwrap_or_default();
    let index = get_attr_u32(&spectrum_rows, ACC_ATTR_INDEX).or(Some(fallback_index));

    let scan_number = get_attr_u32(&spectrum_rows, ACC_ATTR_SCAN_NUMBER);
    let ms_level = get_attr_u32(&spectrum_rows, ACC_ATTR_MS_LEVEL);

    let native_id = get_attr_text(&spectrum_rows, ACC_ATTR_NATIVE_ID);
    let source_file_ref = get_attr_text(&spectrum_rows, ACC_ATTR_SOURCE_FILE_REF);
    let spot_id = get_attr_text(&spectrum_rows, ACC_ATTR_SPOT_ID);

    let data_processing_ref = get_attr_text(&spectrum_rows, ACC_ATTR_DATA_PROCESSING_REF)
        .or_else(|| default_data_processing_ref.map(|s| s.to_string()));

    let mut spectrum_params_meta: Vec<&Metadatum> = Vec::new();
    if spectrum_id != 0 {
        for m in metadata {
            if m.tag_id != TagId::Spectrum || m.owner_id != spectrum_id {
                continue;
            }
            if m.accession
                .as_deref()
                .is_some_and(|a| a.starts_with("B000:"))
            {
                continue;
            }
            spectrum_params_meta.push(m);
        }
    }

    let (cv_params, user_params) = parse_params(allowed_spectrum, &spectrum_params_meta);

    let scan_list = parse_scan_list(schema, metadata, child_index);
    let product_list = parse_product_list(schema, metadata, child_index);
    let precursor_list = parse_precursor_list(schema, metadata, child_index);

    let binary_data_array_list: Option<BinaryDataArrayList> =
        parse_binary_data_array_list(metadata);

    let (x_len, y_len) = xy_lengths_from_bdal(binary_data_array_list.as_ref());

    let default_array_length_attr =
        get_attr_u32(&spectrum_rows, ACC_ATTR_DEFAULT_ARRAY_LENGTH).map(|v| v as usize);

    let default_array_length = default_array_length_attr.or(x_len).or(y_len).or(Some(0));
    let spectrum_description = parse_spectrum_description(schema, metadata, child_index);

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
        spectrum_description,
        scan_list,
        precursor_list,
        product_list,
        binary_data_array_list,
    }
}

#[inline]
fn parse_spectrum_description(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<SpectrumDescription> {
    // <spectrumDescription>
    let mut spectrum_id: u32 = 0;
    for m in metadata {
        if m.tag_id == TagId::Spectrum {
            spectrum_id = m.owner_id;
            break;
        }
    }
    if spectrum_id == 0 {
        return None;
    }

    let sd_id = child_index
        .first_id(spectrum_id, TagId::SpectrumDescription)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::SpectrumDescription && m.parent_index == spectrum_id)
                .map(|m| m.owner_id)
        })?;

    if sd_id == 0 {
        return None;
    }

    let allowed_sd: HashSet<&str> = find_node_by_tag(schema, TagId::SpectrumDescription)
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let sd_meta = collect_subtree_metadata(metadata, child_index, sd_id);
    let sd_child_index = ChildIndex::new(&sd_meta);

    let mut sd_rows: Vec<&Metadatum> = Vec::new();
    for m in &sd_meta {
        if m.tag_id != TagId::SpectrumDescription || m.owner_id != sd_id {
            continue;
        }
        if m.accession
            .as_deref()
            .is_some_and(|a| a.starts_with("B000:"))
        {
            continue;
        }
        sd_rows.push(m);
    }

    let (cv_params, user_params) = parse_params(&allowed_sd, &sd_rows);

    let scan_list = parse_scan_list(schema, &sd_meta, &sd_child_index);
    let product_list = parse_product_list(schema, &sd_meta, &sd_child_index);
    let precursor_list = parse_precursor_list(schema, &sd_meta, &sd_child_index);

    Some(SpectrumDescription {
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
        scan_list,
        precursor_list,
        product_list,
    })
}

#[inline]
fn collect_subtree_metadata(
    metadata: &[Metadatum],
    child_index: &ChildIndex,
    root_owner_id: u32,
) -> Vec<Metadatum> {
    let mut keep: HashSet<u32> = HashSet::new();
    let mut stack = vec![root_owner_id];

    while let Some(id) = stack.pop() {
        if !keep.insert(id) {
            continue;
        }
        for &ch in child_index.children(id) {
            stack.push(ch);
        }
    }

    let mut out = Vec::new();
    out.reserve(metadata.len().min(256));

    for m in metadata {
        if keep.contains(&m.owner_id) {
            out.push(m.clone());
        }
    }

    out
}
