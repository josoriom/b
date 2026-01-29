use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use serde::Serialize;
use serde_json::Value;

use crate::b64::decode::{Metadatum, MetadatumValue};
use crate::mzml::attr_meta::{CV_REF_ATTR, attr_key_from_tail, attr_tail_from_key};
use crate::mzml::schema::{SchemaNode, TagId, schema};

static XMLKEY_TO_TAIL: OnceLock<HashMap<&'static str, u32>> = OnceLock::new();
static TAG_ATTR_SPECS: OnceLock<Mutex<HashMap<TagId, &'static [FieldSpec]>>> = OnceLock::new();

#[derive(Clone)]
struct FieldSpec {
    tail: u32,
    field_key: String,
    camel: String,
    id_variant: Option<String>,
    uri_variant: Option<String>,
}
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
    fn get_expected_field_from_spec<'a>(
        expected_obj: &'a serde_json::Map<String, Value>,
        spec: &FieldSpec,
    ) -> Option<&'a Value> {
        if let Some(v) = expected_obj.get(&spec.field_key) {
            return Some(v);
        }
        if !spec.camel.is_empty() {
            if let Some(v) = expected_obj.get(&spec.camel) {
                return Some(v);
            }
        }
        if let Some(k) = spec.id_variant.as_ref() {
            if let Some(v) = expected_obj.get(k) {
                return Some(v);
            }
        }
        if let Some(k) = spec.uri_variant.as_ref() {
            if let Some(v) = expected_obj.get(k) {
                return Some(v);
            }
        }
        None
    }

    #[inline]
    fn collect_xml_keys(n: &SchemaNode, out: &mut Vec<String>) {
        for xml_keys in n.attributes.values() {
            for k in xml_keys {
                out.push(k.clone());
            }
        }
        for c in n.children.values() {
            collect_xml_keys(c, out);
        }
    }

    #[inline]
    fn build_xmlkey_to_tail() -> HashMap<&'static str, u32> {
        let tree = schema();
        let mut keys: Vec<String> = Vec::new();
        for root in tree.roots.values() {
            collect_xml_keys(root, &mut keys);
        }

        keys.sort();
        keys.dedup();

        let mut m: HashMap<&'static str, u32> = HashMap::with_capacity(keys.len());

        for k in keys {
            let Some(tail) = attr_tail_from_key(k.as_str()) else {
                continue;
            };

            if let Some(&prev) = m.get(k.as_str()) {
                if prev == tail {
                    continue;
                }
                continue;
            }

            let kk: &'static str = Box::leak(k.into_boxed_str());
            m.insert(kk, tail);
        }

        m
    }

    #[inline]
    fn build_specs_for_tag(
        tag_id: TagId,
        xmlkey_to_tail: &HashMap<&'static str, u32>,
    ) -> Vec<FieldSpec> {
        let tree = schema();
        let node = tree
            .roots
            .values()
            .find_map(|root| find_node(root, tag_id))
            .unwrap_or_else(|| panic!("schema missing node for tag_id={tag_id:?}"));

        let mut specs: Vec<FieldSpec> = Vec::with_capacity(node.attributes.len());

        for (field_key, xml_keys) in node.attributes.iter() {
            let mut tail: Option<u32> = None;

            for k in xml_keys {
                if let Some(&t) = xmlkey_to_tail.get(k.as_str()) {
                    tail = Some(t);
                    break;
                }
            }

            let Some(tail) = tail else {
                continue;
            };

            let (camel, id_variant, uri_variant) = snake_to_camel_variants(field_key.as_str());

            specs.push(FieldSpec {
                tail,
                field_key: field_key.clone(),
                camel,
                id_variant,
                uri_variant,
            });
        }

        specs.sort_by_key(|s| s.tail);
        specs
    }

    let expected_json = serde_json::to_value(expected).expect("to_value(expected)");
    let expected_obj = expected_json
        .as_object()
        .unwrap_or_else(|| panic!("expected must serialize to JSON object"));

    let xmlkey_to_tail = XMLKEY_TO_TAIL.get_or_init(build_xmlkey_to_tail);

    let specs_map = TAG_ATTR_SPECS.get_or_init(|| Mutex::new(HashMap::new()));

    let specs: &'static [FieldSpec] = {
        let mut guard = specs_map.lock().expect("TAG_ATTR_SPECS lock");
        if let Some(&specs) = guard.get(&tag_id) {
            specs
        } else {
            let specs_vec = build_specs_for_tag(tag_id, xmlkey_to_tail);
            let leaked: &'static [FieldSpec] = Box::leak(specs_vec.into_boxed_slice());
            guard.insert(tag_id, leaked);
            leaked
        }
    };

    let mut out: Vec<Metadatum> = Vec::with_capacity(specs.len());
    let mut item_index: u32 = 0;

    for spec in specs {
        let Some(v) = get_expected_field_from_spec(expected_obj, spec) else {
            continue;
        };
        let Some(value) = to_metadatum_value(v) else {
            continue;
        };

        let mut accession = String::with_capacity(CV_REF_ATTR.len() + 1 + 7);
        accession.push_str(CV_REF_ATTR);
        accession.push(':');
        use core::fmt::Write;
        write!(&mut accession, "{:07}", spec.tail).expect("write accession tail");

        out.push(Metadatum {
            item_index,
            owner_id,
            parent_index,
            tag_id,
            accession: Some(accession),
            unit_accession: None,
            value,
        });
        item_index += 1;
    }

    out
}
