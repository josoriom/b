use std::collections::{HashMap, HashSet};

use crate::{
    ScanList, ScanWindow, ScanWindowList,
    b64::utilities::{
        common::{ChildIndex, child_node, find_node_by_tag, ordered_unique_owner_ids},
        parse_cv_and_user_params,
    },
    decode::Metadatum,
    mzml::{
        schema::{SchemaTree as Schema, TagId},
        structs::Scan,
    },
};

#[inline]
pub fn parse_scan_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<ScanList> {
    let list_node = find_node_by_tag(schema, TagId::ScanList)?;
    let scan_node = child_node(Some(list_node), TagId::Scan)?;

    let mut allowed_scan_list: HashSet<&str> = child_node(Some(list_node), TagId::CvParam)
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let allowed_scan: HashSet<&str> = child_node(Some(scan_node), TagId::CvParam)
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let scan_window_node = child_node(Some(scan_node), TagId::ScanWindowList)
        .and_then(|n| child_node(Some(n), TagId::ScanWindow));

    let allowed_scan_window: HashSet<&str> = scan_window_node
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::new();
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let scan_list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::ScanList)
        .map(|m| m.owner_id)?;

    let scan_list_rows = rows_for_owner(&owner_rows, scan_list_id);

    if allowed_scan_list.is_empty() {
        for m in scan_list_rows {
            if let Some(acc) = m.accession.as_deref() {
                if !acc.starts_with("B000:") {
                    allowed_scan_list.insert(acc);
                }
            }
        }
    }

    let mut scan_list_params_meta: Vec<&Metadatum> = Vec::new();
    scan_list_params_meta.extend(scan_list_rows.iter().copied());
    scan_list_params_meta.extend(child_params_for_parent(
        &owner_rows,
        child_index,
        scan_list_id,
    ));

    let (cv_params, user_params) =
        parse_cv_and_user_params(&allowed_scan_list, &scan_list_params_meta);

    let mut scan_ids: Vec<u32> = unique_ids(child_index.ids(scan_list_id, TagId::Scan));
    if scan_ids.is_empty() {
        scan_ids = ordered_unique_owner_ids(metadata, TagId::Scan);
    }
    if scan_ids.is_empty() {
        return None;
    }

    let single_scan_id = if scan_ids.len() == 1 {
        Some(scan_ids[0])
    } else {
        None
    };

    let all_scan_window_ids = ordered_unique_owner_ids(metadata, TagId::ScanWindow);

    let mut scans = Vec::with_capacity(scan_ids.len());
    for scan_id in scan_ids {
        let scan_rows = rows_for_owner(&owner_rows, scan_id);
        let scan_parent = scan_rows.first().map(|m| m.parent_index).unwrap_or(0);

        let scan_window_ids = scan_window_ids_for_scan(
            child_index,
            scan_id,
            scan_parent,
            single_scan_id,
            &all_scan_window_ids,
        );

        let scan_window_list =
            parse_scan_window_list(&allowed_scan_window, &owner_rows, &scan_window_ids);

        scans.push(parse_scan(&allowed_scan, scan_rows, scan_window_list));
    }

    Some(ScanList {
        count: Some(scans.len()),
        cv_params,
        user_params,
        scans,
    })
}

/// <scan>
#[inline]
fn parse_scan(
    allowed_scan: &HashSet<&str>,
    scan_rows: &[&Metadatum],
    scan_window_list: Option<ScanWindowList>,
) -> Scan {
    let (cv_params, user_params) = parse_cv_and_user_params(allowed_scan, scan_rows);

    Scan {
        instrument_configuration_ref: None,
        external_spectrum_id: None,
        source_file_ref: None,
        spectrum_ref: None,
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
        scan_window_list,
    }
}

/// <scanWindowList>
#[inline]
fn parse_scan_window_list(
    allowed_scan_window: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    scan_window_ids: &[u32],
) -> Option<ScanWindowList> {
    if scan_window_ids.is_empty() {
        return None;
    }

    let mut scan_windows = Vec::with_capacity(scan_window_ids.len());
    for &sw_id in scan_window_ids {
        let sw_rows = rows_for_owner(owner_rows, sw_id);
        scan_windows.push(parse_scan_window(allowed_scan_window, sw_rows));
    }

    Some(ScanWindowList {
        count: Some(scan_windows.len()),
        scan_windows,
    })
}

/// <scanWindow>
#[inline]
fn parse_scan_window(
    allowed_scan_window: &HashSet<&str>,
    scan_window_rows: &[&Metadatum],
) -> ScanWindow {
    let (cv_params, user_params) = parse_cv_and_user_params(allowed_scan_window, scan_window_rows);

    ScanWindow {
        cv_params,
        user_params,
        ..Default::default()
    }
}

#[inline]
fn rows_for_owner<'a>(
    owner_rows: &'a HashMap<u32, Vec<&'a Metadatum>>,
    owner_id: u32,
) -> &'a [&'a Metadatum] {
    owner_rows
        .get(&owner_id)
        .map(|v| v.as_slice())
        .unwrap_or(&[])
}

#[inline]
fn child_params_for_parent<'a>(
    owner_rows: &HashMap<u32, Vec<&'a Metadatum>>,
    child_index: &ChildIndex,
    parent_id: u32,
) -> Vec<&'a Metadatum> {
    let cv_ids = child_index.ids(parent_id, TagId::CvParam);
    let up_ids = child_index.ids(parent_id, TagId::UserParam);

    let mut out = Vec::with_capacity(cv_ids.len() + up_ids.len());

    for &id in cv_ids {
        if let Some(rows) = owner_rows.get(&id) {
            out.extend(rows.iter().copied());
        }
    }
    for &id in up_ids {
        if let Some(rows) = owner_rows.get(&id) {
            out.extend(rows.iter().copied());
        }
    }

    out
}

#[inline]
fn scan_window_ids_for_scan(
    child_index: &ChildIndex,
    scan_id: u32,
    scan_parent: u32,
    single_scan_id: Option<u32>,
    all_scan_window_ids: &[u32],
) -> Vec<u32> {
    let mut out: Vec<u32> = Vec::new();

    out.extend_from_slice(child_index.ids(scan_id, TagId::ScanWindow));
    if out.is_empty() {
        out.extend_from_slice(child_index.ids(scan_parent, TagId::ScanWindow));
    }

    if out.is_empty() {
        for parent in [scan_id, scan_parent] {
            for &swl_id in child_index.ids(parent, TagId::ScanWindowList) {
                out.extend_from_slice(child_index.ids(swl_id, TagId::ScanWindow));
            }
        }
    }

    if out.is_empty() && single_scan_id == Some(scan_id) {
        out.extend_from_slice(all_scan_window_ids);
    }

    let mut seen = HashSet::with_capacity(out.len());
    out.retain(|id| seen.insert(*id));
    out
}

#[inline]
fn unique_ids(ids: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(ids.len());
    let mut seen = HashSet::with_capacity(ids.len());
    for &id in ids {
        if seen.insert(id) {
            out.push(id);
        }
    }
    out
}
