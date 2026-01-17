use std::collections::{HashMap, HashSet};

use crate::{
    Product, ProductList,
    b64::utilities::{
        common::{ChildIndex, child_node, find_node_by_tag, ordered_unique_owner_ids},
        parse_cv_and_user_params,
    },
    decode::Metadatum,
    mzml::{
        schema::{SchemaTree as Schema, TagId},
        structs::IsolationWindow,
    },
};

/// <productList>
#[inline]
pub fn parse_product_list(
    schema: &Schema,
    metadata: &[Metadatum],
    child_index: &ChildIndex,
) -> Option<ProductList> {
    let list_node = find_node_by_tag(schema, TagId::ProductList)?;
    let product_node = child_node(Some(list_node), TagId::Product)?;

    let allowed_isolation_window: HashSet<&str> =
        child_node(Some(product_node), TagId::IsolationWindow)
            .and_then(|n| child_node(Some(n), TagId::CvParam))
            .map(|n| n.accessions.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default();

    let mut meta_by_id: HashMap<u32, &Metadatum> = HashMap::with_capacity(metadata.len());
    for m in metadata {
        meta_by_id.insert(m.owner_id, m);
    }

    let product_list_id = metadata
        .iter()
        .find(|m| m.tag_id == TagId::ProductList)
        .map(|m| m.owner_id);

    let mut product_ids = if let Some(list_id) = product_list_id {
        child_index.ids(list_id, TagId::Product).to_vec()
    } else {
        Vec::new()
    };

    if product_ids.is_empty() {
        product_ids = ordered_unique_owner_ids(metadata, TagId::Product);
    }

    if product_ids.is_empty() {
        return None;
    }

    let mut products = Vec::with_capacity(product_ids.len());
    for product_id in product_ids {
        products.push(parse_product(
            product_id,
            &allowed_isolation_window,
            &meta_by_id,
            child_index,
        ));
    }

    Some(ProductList {
        count: Some(products.len()),
        products,
    })
}

/// <product>
#[inline]
fn parse_product(
    product_id: u32,
    allowed_isolation_window: &HashSet<&str>,
    meta_by_id: &HashMap<u32, &Metadatum>,
    child_index: &ChildIndex,
) -> Product {
    let product_parent = meta_by_id
        .get(&product_id)
        .map(|m| m.parent_index)
        .unwrap_or(0);

    let isolation_window_id = child_index
        .first_id(product_id, TagId::IsolationWindow)
        .or_else(|| child_index.first_id(product_parent, TagId::IsolationWindow));

    Product {
        spectrum_ref: None,
        source_file_ref: None,
        external_spectrum_id: None,
        isolation_window: parse_isolation_window(
            allowed_isolation_window,
            meta_by_id,
            child_index,
            isolation_window_id,
            product_id,
            product_parent,
        ),
    }
}

/// <isolationWindow>
#[inline]
fn parse_isolation_window(
    allowed_isolation_window: &HashSet<&str>,
    meta_by_id: &HashMap<u32, &Metadatum>,
    child_index: &ChildIndex,
    isolation_window_id: Option<u32>,
    product_id: u32,
    product_parent: u32,
) -> Option<IsolationWindow> {
    let mut meta = if let Some(iso_id) = isolation_window_id {
        params_for_parent(meta_by_id, child_index, iso_id)
    } else {
        params_for_parent(meta_by_id, child_index, product_id)
    };

    if meta.is_empty() && product_parent != product_id {
        meta = params_for_parent(meta_by_id, child_index, product_parent);
    }

    if meta.is_empty() && isolation_window_id.is_none() {
        return None;
    }

    let (cv_params, user_params) = parse_cv_and_user_params(allowed_isolation_window, &meta);

    Some(IsolationWindow {
        referenceable_param_group_refs: Vec::new(),
        cv_params,
        user_params,
    })
}

#[inline]
fn params_for_parent<'a>(
    meta_by_id: &HashMap<u32, &'a Metadatum>,
    child_index: &ChildIndex,
    parent_id: u32,
) -> Vec<&'a Metadatum> {
    let cv_ids = child_index.ids(parent_id, TagId::CvParam);
    let up_ids = child_index.ids(parent_id, TagId::UserParam);

    let mut out = Vec::with_capacity(cv_ids.len() + up_ids.len());

    for &id in cv_ids {
        if let Some(m) = meta_by_id.get(&id) {
            out.push(*m);
        }
    }
    for &id in up_ids {
        if let Some(m) = meta_by_id.get(&id) {
            out.push(*m);
        }
    }

    out
}
