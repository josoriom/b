use std::collections::HashMap;

use serde::Serialize;
use serde_json::Value;

use crate::b64::decode::{Metadatum, MetadatumValue};
use crate::mzml::attr_meta::{CV_REF_ATTR, attr_key_from_tail};
use crate::mzml::schema::{SchemaNode, TagId, schema};

static XMLKEY_TO_TAIL: std::sync::OnceLock<HashMap<&'static str, u32>> = std::sync::OnceLock::new();

pub fn assign_attributes<T>(
    expected: &T,
    tag_id: TagId,
    owner_id: u32,
    parent_index: u32,
) -> Vec<Metadatum>
where
    T: Serialize,
{
    #[inline]
    fn find_node<'a>(n: &'a SchemaNode, tag_id: TagId) -> Option<&'a SchemaNode> {
        if n.self_tags.iter().any(|t| *t == tag_id) {
            return Some(n);
        }
        for c in n.children.values() {
            if let Some(x) = find_node(c, tag_id) {
                return Some(x);
            }
        }
        None
    }

    #[inline]
    fn to_metadatum_value(v: &Value) -> Option<MetadatumValue> {
        match v {
            Value::Null => None,
            Value::String(s) => Some(MetadatumValue::Text(s.clone())),
            Value::Number(n) => {
                if let Some(u) = n.as_u64() {
                    if u > (1u64 << 53) {
                        panic!("numeric attribute too large for lossless f64: {u}");
                    }
                    Some(MetadatumValue::Number(u as f64))
                } else if let Some(i) = n.as_i64() {
                    if (i.unsigned_abs() as u64) > (1u64 << 53) {
                        panic!("numeric attribute too large for lossless f64: {i}");
                    }
                    Some(MetadatumValue::Number(i as f64))
                } else if let Some(f) = n.as_f64() {
                    Some(MetadatumValue::Number(f))
                } else {
                    None
                }
            }
            Value::Bool(b) => Some(MetadatumValue::Text(b.to_string())),
            _ => panic!("non-primitive attribute value not supported: {v:?}"),
        }
    }

    #[inline]
    fn push_upper_first(dst: &mut String, s: &str) {
        let mut it = s.chars();
        let Some(c0) = it.next() else { return };
        for uc in c0.to_uppercase() {
            dst.push(uc);
        }
        dst.push_str(it.as_str());
    }

    #[inline]
    fn snake_to_camel_variants(s: &str) -> (String, Option<String>, Option<String>) {
        let mut parts = s.split('_').filter(|p| !p.is_empty());

        let first = match parts.next() {
            Some(f) => f,
            None => return (String::new(), None, None),
        };

        let mut camel = String::from(first);
        for p in parts {
            push_upper_first(&mut camel, p);
        }

        let id_variant = if camel.ends_with("Id") {
            let mut v = camel.clone();
            v.truncate(v.len() - 2);
            v.push_str("ID");
            Some(v)
        } else {
            None
        };

        let uri_variant = if camel.ends_with("Uri") {
            let mut v = camel.clone();
            v.truncate(v.len() - 3);
            v.push_str("URI");
            Some(v)
        } else {
            None
        };

        (camel, id_variant, uri_variant)
    }

    #[inline]
    fn get_expected_field<'a>(
        expected_obj: &'a serde_json::Map<String, Value>,
        field_key: &str,
    ) -> Option<&'a Value> {
        if let Some(v) = expected_obj.get(field_key) {
            return Some(v);
        }

        let (camel, id_var, uri_var) = snake_to_camel_variants(field_key);

        if !camel.is_empty() {
            if let Some(v) = expected_obj.get(&camel) {
                return Some(v);
            }
        }
        if let Some(k) = id_var {
            if let Some(v) = expected_obj.get(&k) {
                return Some(v);
            }
        }
        if let Some(k) = uri_var {
            if let Some(v) = expected_obj.get(&k) {
                return Some(v);
            }
        }

        None
    }

    let tree = schema();
    let node = tree
        .roots
        .values()
        .find_map(|root| find_node(root, tag_id))
        .unwrap_or_else(|| panic!("schema missing node for tag_id={tag_id:?}"));

    let expected_json = serde_json::to_value(expected).expect("to_value(expected)");
    let expected_obj = expected_json
        .as_object()
        .unwrap_or_else(|| panic!("expected must serialize to JSON object"));

    let xmlkey_to_tail = XMLKEY_TO_TAIL.get_or_init(|| {
        let mut m = HashMap::new();
        for tail in 9_900_000u32..=9_920_000u32 {
            if let Some(k) = attr_key_from_tail(tail) {
                m.insert(k, tail);
            }
        }
        m
    });

    let mut items: Vec<(u32, MetadatumValue)> = Vec::with_capacity(node.attributes.len());

    for (field_key, xml_keys) in node.attributes.iter() {
        let Some(v) = get_expected_field(expected_obj, field_key.as_str()) else {
            continue;
        };
        let Some(value) = to_metadatum_value(v) else {
            continue;
        };

        let xml_key = xml_keys
            .first()
            .unwrap_or_else(|| panic!("schema attribute {field_key:?} has empty xml key list"));

        let tail = *xmlkey_to_tail.get(xml_key.as_str()).unwrap_or_else(|| {
            panic!(
                "no B000 tail mapping for xml attribute key {xml_key:?} \
                 (field {field_key:?}, tag {tag_id:?})"
            )
        });

        items.push((tail, value));
    }

    items.sort_by_key(|(tail, _)| *tail);

    let mut out = Vec::with_capacity(items.len());
    let mut item_index: u32 = 0;

    for (tail, value) in items {
        out.push(Metadatum {
            item_index,
            owner_id,
            parent_index,
            tag_id,
            accession: Some(format!("{CV_REF_ATTR}:{tail:07}")),
            unit_accession: None,
            value,
        });
        item_index += 1;
    }

    out
}
