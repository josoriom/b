use std::{
    collections::{HashMap, HashSet},
    io::Read,
};

use crate::{
    BinaryDataArray, BinaryDataArrayList,
    decode2::{Metadatum, MetadatumValue},
    mzml::{
        attr_meta::CV_REF_ATTR,
        schema::{SchemaNode, SchemaTree as Schema, TagId},
    },
};

pub const ACC_Y_INTENSITY: &str = "MS:1000515";
pub const ACC_Y_SNR: &str = "MS:1000786";

#[inline]
pub fn take<'a>(
    bytes: &'a [u8],
    pos: &mut usize,
    n: usize,
    field: &'static str,
) -> Result<&'a [u8], String> {
    let end = pos
        .checked_add(n)
        .ok_or_else(|| format!("overflow while reading {field}"))?;
    if end > bytes.len() {
        return Err(format!(
            "unexpected EOF while reading {field}: need {n} bytes at pos {pos}, len {}",
            bytes.len()
        ));
    }
    let out = &bytes[*pos..end];
    *pos = end;
    Ok(out)
}

#[inline]
pub fn read_u32_vec(bytes: &[u8], pos: &mut usize, n: usize) -> Result<Vec<u32>, String> {
    let raw = take(bytes, pos, n * 4, "u32 vector")?;
    let mut out = Vec::with_capacity(n);
    for chunk in raw.chunks_exact(4) {
        out.push(u32::from_le_bytes(chunk.try_into().unwrap()));
    }
    Ok(out)
}

#[inline]
pub fn read_f64_vec(bytes: &[u8], pos: &mut usize, n: usize) -> Result<Vec<f64>, String> {
    let raw = take(bytes, pos, n * 8, "f64 vector")?;
    let mut out = Vec::with_capacity(n);
    for chunk in raw.chunks_exact(8) {
        out.push(f64::from_le_bytes(chunk.try_into().unwrap()));
    }
    Ok(out)
}

#[inline]
pub fn vs_len_bytes(vk: &[u8], vi: &[u32], voff: &[u32], vlen: &[u32]) -> Result<usize, String> {
    let mut max_end = 0usize;

    for j in 0..vk.len() {
        if vk[j] != 1 {
            continue;
        }
        let idx = vi[j] as usize;
        if idx >= voff.len() || idx >= vlen.len() {
            return Err("string VI out of range".to_string());
        }
        let end = (voff[idx] as usize)
            .checked_add(vlen[idx] as usize)
            .ok_or_else(|| "VOFF+VLEN overflow".to_string())?;
        if end > max_end {
            max_end = end;
        }
    }

    Ok(max_end)
}

#[inline]
pub fn decompress_zstd_allow_aligned_padding(input: &[u8]) -> Result<Vec<u8>, String> {
    match decompress_zstd(input) {
        Ok(v) => Ok(v),
        Err(first_err) => {
            let mut trimmed = input;
            for _ in 0..7 {
                if trimmed.is_empty() || *trimmed.last().unwrap() != 0 {
                    break;
                }
                trimmed = &trimmed[..trimmed.len() - 1];
                if let Ok(v) = decompress_zstd(trimmed) {
                    return Ok(v);
                }
            }
            Err(first_err)
        }
    }
}

#[inline]
pub fn decompress_zstd(mut input: &[u8]) -> Result<Vec<u8>, String> {
    let mut dec = zstd::Decoder::new(&mut input).map_err(|e| format!("zstd decoder init: {e}"))?;
    let mut out = Vec::new();
    dec.read_to_end(&mut out)
        .map_err(|e| format!("zstd decode: {e}"))?;
    Ok(out)
}

#[inline]
pub fn find_node_by_tag<'a>(schema: &'a Schema, tag: TagId) -> Option<&'a SchemaNode> {
    if let Some(n) = schema.root_by_tag(tag) {
        return Some(n);
    }
    for root in schema.roots.values() {
        if let Some(n) = find_node_by_tag_rec(root, tag) {
            return Some(n);
        }
    }
    None
}

#[inline]
pub fn find_node_by_tag_rec<'a>(node: &'a SchemaNode, tag: TagId) -> Option<&'a SchemaNode> {
    if node.self_tags.iter().any(|&t| t == tag) {
        return Some(node);
    }
    for child in node.children.values() {
        if let Some(n) = find_node_by_tag_rec(child, tag) {
            return Some(n);
        }
    }
    None
}

#[inline]
pub fn child_node<'a>(parent: Option<&'a SchemaNode>, tag: TagId) -> Option<&'a SchemaNode> {
    let p = parent?;
    let key = p.child_key_for_tag(tag)?;
    p.children.get(key)
}

#[inline]
pub fn value_to_opt_string(v: &MetadatumValue) -> Option<String> {
    match v {
        MetadatumValue::Empty => None,
        MetadatumValue::Text(s) => Some(s.clone()),
        MetadatumValue::Number(x) => Some(x.to_string()),
    }
}

#[inline]
pub fn is_cv_prefix(p: &str) -> bool {
    matches!(p, "MS" | "UO" | "NCIT" | "PEFF")
}

#[inline]
pub fn unit_cv_ref(unit_accession: &Option<String>) -> Option<String> {
    unit_accession
        .as_deref()
        .and_then(|u| u.split_once(':'))
        .map(|(prefix, _)| prefix.to_string())
}

#[inline]
pub fn get_attr_u32(rows: &[&Metadatum], tail: u32) -> Option<u32> {
    let s = get_attr_text(rows, tail)?;
    s.parse::<u32>().ok()
}

#[inline]
pub fn get_attr_text(rows: &[&Metadatum], tail: u32) -> Option<String> {
    for m in rows {
        let acc = m.accession.as_deref()?;
        let (p, _) = split_prefix(acc)?;
        if p != CV_REF_ATTR {
            continue;
        }
        if parse_accession_tail_str(acc) != tail {
            continue;
        }
        return value_to_opt_string(&m.value);
    }
    None
}

#[inline]
pub fn split_prefix(acc: &str) -> Option<(&str, &str)> {
    acc.split_once(':')
}

#[inline]
pub fn parse_accession_tail_str(acc: &str) -> u32 {
    let tail = match acc.rsplit_once(':') {
        Some((_, t)) => t,
        None => acc,
    };

    let mut v: u32 = 0;
    let mut saw = false;
    for b in tail.bytes() {
        if (b'0'..=b'9').contains(&b) {
            saw = true;
            let d = (b - b'0') as u32;
            match v.checked_mul(10).and_then(|x| x.checked_add(d)) {
                Some(n) => v = n,
                None => return 0,
            }
        }
    }
    if saw { v } else { 0 }
}

#[inline]
pub fn xy_lengths_from_bdal(list: Option<&BinaryDataArrayList>) -> (Option<usize>, Option<usize>) {
    let Some(list) = list else {
        return (None, None);
    };

    let mut x_len = None;
    let mut y_len = None;

    for bda in &list.binary_data_arrays {
        let len = decoded_len(bda);
        if len == 0 {
            continue;
        }

        if is_y_array(bda) {
            if y_len.is_none() {
                y_len = Some(len);
            }
        } else if x_len.is_none() {
            x_len = Some(len);
        }

        if x_len.is_some() && y_len.is_some() {
            break;
        }
    }

    (x_len, y_len)
}

#[inline]
pub fn decoded_len(bda: &BinaryDataArray) -> usize {
    if !bda.decoded_binary_f64.is_empty() {
        bda.decoded_binary_f64.len()
    } else {
        bda.decoded_binary_f32.len()
    }
}

#[inline]
pub fn is_y_array(bda: &BinaryDataArray) -> bool {
    bda.cv_params.iter().any(|p| {
        matches!(
            p.accession.as_deref(),
            Some(ACC_Y_INTENSITY) | Some(ACC_Y_SNR)
        )
    })
}

#[inline]
pub fn ordered_unique_owner_ids(metadata: &[Metadatum], tag: TagId) -> Vec<u32> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();

    for m in metadata {
        if m.tag_id == tag && seen.insert(m.owner_id) {
            out.push(m.owner_id);
        }
    }

    out
}

#[inline]
pub fn collect_subtree_owner_ids(
    root_id: u32,
    children_by_parent: &HashMap<u32, Vec<u32>>,
) -> HashSet<u32> {
    let mut out = HashSet::new();
    let mut stack = vec![root_id];

    while let Some(id) = stack.pop() {
        if !out.insert(id) {
            continue;
        }
        if let Some(children) = children_by_parent.get(&id) {
            for &child_id in children {
                stack.push(child_id);
            }
        }
    }

    out
}

#[inline]
pub fn key_parent_tag(parent_id: u32, tag: TagId) -> u64 {
    ((parent_id as u64) << 8) | (tag as u8 as u64)
}

pub struct ChildIndex {
    ids_by_parent_tag: HashMap<u64, Vec<u32>>,
    children_by_parent: HashMap<u32, Vec<u32>>,
}

impl ChildIndex {
    #[inline]
    pub fn new(metadata: &[Metadatum]) -> Self {
        let mut ids_by_parent_tag: HashMap<u64, Vec<u32>> = HashMap::new();
        let mut children_by_parent: HashMap<u32, Vec<u32>> = HashMap::new();

        for m in metadata {
            ids_by_parent_tag
                .entry(key_parent_tag(m.parent_index, m.tag_id))
                .or_default()
                .push(m.owner_id);

            children_by_parent
                .entry(m.parent_index)
                .or_default()
                .push(m.owner_id);
        }

        Self {
            ids_by_parent_tag,
            children_by_parent,
        }
    }

    #[inline]
    pub fn ids(&self, parent_id: u32, tag: TagId) -> &[u32] {
        self.ids_by_parent_tag
            .get(&key_parent_tag(parent_id, tag))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    #[inline]
    pub fn first_id(&self, parent_id: u32, tag: TagId) -> Option<u32> {
        self.ids(parent_id, tag).first().copied()
    }

    #[inline]
    pub fn children(&self, parent_id: u32) -> &[u32] {
        self.children_by_parent
            .get(&parent_id)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }
}
