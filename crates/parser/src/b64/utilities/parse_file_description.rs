use std::collections::{HashMap, HashSet};

use crate::{
    b64::utilities::{
        common::{ChildIndex, child_node, find_node_by_tag, ordered_unique_owner_ids},
        parse_cv_and_user_params,
    },
    decode::{Metadatum, MetadatumValue},
    mzml::{
        attr_meta::{ACC_ATTR_COUNT, ACC_ATTR_ID, ACC_ATTR_LOCATION, ACC_ATTR_NAME},
        schema::{SchemaTree as Schema, TagId},
        structs::{
            Contact, FileContent, FileDescription, ReferenceableParamGroupRef, SourceFile,
            SourceFileList,
        },
    },
};

/// <fileDescription>
#[inline]
pub fn parse_file_description(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<FileDescription> {
    let fd_node = find_node_by_tag(schema, TagId::FileDescription)?;

    let file_content_node = child_node(Some(fd_node), TagId::FileContent)?;
    let source_file_node = child_node(
        child_node(Some(fd_node), TagId::SourceFileList),
        TagId::SourceFile,
    );
    let contact_node = child_node(Some(fd_node), TagId::Contact);

    let allowed_file_content_schema: HashSet<&str> =
        child_node(Some(file_content_node), TagId::CvParam)
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let allowed_source_file_schema: HashSet<&str> = source_file_node
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let allowed_contact_schema: HashSet<&str> = contact_node
        .and_then(|n| child_node(Some(n), TagId::CvParam))
        .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
        .unwrap_or_default();

    let mut owner_rows: HashMap<u32, Vec<&Metadatum>> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        owner_rows.entry(m.owner_id).or_default().push(m);
    }

    let file_desc_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::FileDescription)
        .map(|m| m.owner_id)
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::FileContent)
                .map(|m| m.parent_index)
        })
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::SourceFileList)
                .map(|m| m.parent_index)
        })
        .or_else(|| {
            let sfl_id = metadata
                .iter()
                .find(|m| m.tag_id == TagId::SourceFile)
                .map(|m| m.parent_index)?;
            owner_rows
                .get(&sfl_id)
                .and_then(|rows| rows.first())
                .map(|m| m.parent_index)
        })
        .or_else(|| {
            metadata
                .iter()
                .find(|m| m.tag_id == TagId::Contact)
                .map(|m| m.parent_index)
        })?;

    let file_content = parse_file_content(
        &allowed_file_content_schema,
        &owner_rows,
        child_index,
        metadata,
        file_desc_id,
    );

    let source_file_list = parse_source_file_list(
        &allowed_source_file_schema,
        &owner_rows,
        child_index,
        metadata,
        file_desc_id,
    );

    let contacts = parse_contacts(
        &allowed_contact_schema,
        &owner_rows,
        child_index,
        metadata,
        file_desc_id,
    );

    Some(FileDescription {
        file_content,
        source_file_list,
        contacts,
    })
}

/// <fileContent>
#[inline]
fn parse_file_content(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    file_desc_id: u32,
) -> FileContent {
    let mut ids = unique_ids(child_index.ids(file_desc_id, TagId::FileContent));
    if ids.is_empty() {
        ids = ordered_unique_owner_ids(metadata, TagId::FileContent);
        ids.retain(|&id| is_child_of(owner_rows, id, file_desc_id));
    }

    let file_content_id = ids.first().copied().unwrap_or(0);
    if file_content_id == 0 {
        return FileContent::default();
    }

    let rows = rows_for_owner(owner_rows, file_content_id);
    let child_meta = child_params_for_parent(owner_rows, child_index, file_content_id);

    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_params, user_params) = if allowed_schema.is_empty() {
        let allowed_meta = allowed_from_rows(&params_meta);
        parse_cv_and_user_params(&allowed_meta, &params_meta)
    } else {
        parse_cv_and_user_params(allowed_schema, &params_meta)
    };

    FileContent {
        referenceable_param_group_refs: Vec::<ReferenceableParamGroupRef>::new(),
        cv_params,
        user_params,
    }
}

/// <sourceFileList>
#[inline]
fn parse_source_file_list(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    file_desc_id: u32,
) -> SourceFileList {
    let mut sfl_ids = unique_ids(child_index.ids(file_desc_id, TagId::SourceFileList));
    if sfl_ids.is_empty() {
        sfl_ids = ordered_unique_owner_ids(metadata, TagId::SourceFileList);
        sfl_ids.retain(|&id| is_child_of(owner_rows, id, file_desc_id));
    }
    let sfl_id = sfl_ids.first().copied().unwrap_or(0);
    if sfl_id == 0 {
        return SourceFileList {
            count: Some(0),
            source_file: Vec::new(),
        };
    }

    let sfl_rows = rows_for_owner(owner_rows, sfl_id);
    let count_attr = b000_attr_text(sfl_rows, ACC_ATTR_COUNT).and_then(|s| s.parse::<usize>().ok());

    let mut source_file_ids = unique_ids(child_index.ids(sfl_id, TagId::SourceFile));
    if source_file_ids.is_empty() {
        source_file_ids = ordered_unique_owner_ids(metadata, TagId::SourceFile);
        source_file_ids.retain(|&id| is_child_of(owner_rows, id, sfl_id));
    }

    let mut source_file = Vec::with_capacity(source_file_ids.len());
    for id in source_file_ids {
        source_file.push(parse_source_file(
            allowed_schema,
            owner_rows,
            child_index,
            id,
        ));
    }

    SourceFileList {
        count: count_attr.or(Some(source_file.len())),
        source_file,
    }
}

/// <sourceFile>
#[inline]
fn parse_source_file(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    source_file_id: u32,
) -> SourceFile {
    let rows = rows_for_owner(owner_rows, source_file_id);

    let id = b000_attr_text(rows, ACC_ATTR_ID).unwrap_or_default();
    let name = b000_attr_text(rows, ACC_ATTR_NAME).unwrap_or_default();
    let location = b000_attr_text(rows, ACC_ATTR_LOCATION).unwrap_or_default();

    let child_meta = child_params_for_parent(owner_rows, child_index, source_file_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_param, user_param) = if allowed_schema.is_empty() {
        let allowed_meta = allowed_from_rows(&params_meta);
        parse_cv_and_user_params(&allowed_meta, &params_meta)
    } else {
        parse_cv_and_user_params(allowed_schema, &params_meta)
    };

    SourceFile {
        id,
        name,
        location,
        referenceable_param_group_ref: Vec::<ReferenceableParamGroupRef>::new(),
        cv_param,
        user_param,
    }
}

/// <contact>
#[inline]
fn parse_contacts(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    metadata: &[Metadatum],
    file_desc_id: u32,
) -> Vec<Contact> {
    let mut contact_ids = unique_ids(child_index.ids(file_desc_id, TagId::Contact));
    if contact_ids.is_empty() {
        contact_ids = ordered_unique_owner_ids(metadata, TagId::Contact);
        contact_ids.retain(|&id| is_child_of(owner_rows, id, file_desc_id));
    }

    let mut contacts = Vec::with_capacity(contact_ids.len());
    for id in contact_ids {
        contacts.push(parse_contact(allowed_schema, owner_rows, child_index, id));
    }
    contacts
}

/// <contact>
#[inline]
fn parse_contact(
    allowed_schema: &HashSet<&str>,
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_index: &ChildIndex,
    contact_id: u32,
) -> Contact {
    let rows = rows_for_owner(owner_rows, contact_id);

    let child_meta = child_params_for_parent(owner_rows, child_index, contact_id);
    let mut params_meta: Vec<&Metadatum> = Vec::with_capacity(rows.len() + child_meta.len());
    params_meta.extend_from_slice(rows);
    params_meta.extend(child_meta);

    let (cv_params, user_params) = if allowed_schema.is_empty() {
        let allowed_meta = allowed_from_rows(&params_meta);
        parse_cv_and_user_params(&allowed_meta, &params_meta)
    } else {
        parse_cv_and_user_params(allowed_schema, &params_meta)
    };

    Contact {
        referenceable_param_group_refs: Vec::<ReferenceableParamGroupRef>::new(),
        cv_params,
        user_params,
    }
}

#[inline]
pub fn rows_for_owner<'a>(
    owner_rows: &'a HashMap<u32, Vec<&'a Metadatum>>,
    owner_id: u32,
) -> &'a [&'a Metadatum] {
    owner_rows
        .get(&owner_id)
        .map(|v| v.as_slice())
        .unwrap_or(&[])
}

#[inline]
pub fn child_params_for_parent<'a>(
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
pub fn allowed_from_rows<'a>(rows: &[&'a Metadatum]) -> HashSet<&'a str> {
    let mut allowed = HashSet::new();
    for m in rows {
        if let Some(acc) = m.accession.as_deref() {
            if !acc.starts_with("B000:") {
                allowed.insert(acc);
            }
        }
    }
    allowed
}

#[inline]
pub fn is_child_of(
    owner_rows: &HashMap<u32, Vec<&Metadatum>>,
    child_id: u32,
    parent_id: u32,
) -> bool {
    rows_for_owner(owner_rows, child_id)
        .first()
        .map(|m| m.parent_index == parent_id)
        .unwrap_or(false)
}

#[inline]
pub fn b000_attr_text(rows: &[&Metadatum], accession_tail: u32) -> Option<String> {
    for m in rows {
        let acc = m.accession.as_deref()?;
        if !acc.starts_with("B000:") {
            continue;
        }
        if parse_accession_tail(Some(acc)) != accession_tail {
            continue;
        }
        return match &m.value {
            MetadatumValue::Text(s) => Some(s.clone()),
            MetadatumValue::Number(n) => Some(n.to_string()),
            _ => None,
        };
    }
    None
}

#[inline]
pub fn parse_accession_tail(accession: Option<&str>) -> u32 {
    let s = accession.unwrap_or("");
    let tail = match s.rsplit_once(':') {
        Some((_, t)) => t,
        None => s,
    };

    let mut v: u32 = 0;
    let mut saw_digit = false;

    for b in tail.bytes() {
        if (b'0'..=b'9').contains(&b) {
            saw_digit = true;
            let d = (b - b'0') as u32;
            match v.checked_mul(10).and_then(|x| x.checked_add(d)) {
                Some(n) => v = n,
                None => return 0,
            }
        }
    }

    if saw_digit { v } else { 0 }
}

#[inline]
pub fn unique_ids(ids: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(ids.len());
    let mut seen = HashSet::with_capacity(ids.len());
    for &id in ids {
        if seen.insert(id) {
            out.push(id);
        }
    }
    out
}
