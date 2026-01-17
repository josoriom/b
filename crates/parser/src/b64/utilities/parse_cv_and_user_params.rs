use std::collections::{HashMap, HashSet};

use crate::{
    CvParam, UserParam,
    b64::utilities::common::{is_cv_prefix, unit_cv_ref, value_to_opt_string},
    decode::Metadatum,
    mzml::{attr_meta::CV_REF_ATTR, cv_table},
};

#[inline]
pub fn parse_list_grouped_by_owner_id<'m, T, F, I>(iter: I, mut parse_item: F) -> Option<Vec<T>>
where
    I: IntoIterator<Item = &'m Metadatum>,
    F: FnMut(u32, &[&'m Metadatum]) -> T,
{
    let mut groups: HashMap<u32, Vec<&'m Metadatum>> = HashMap::new();

    for m in iter {
        groups.entry(m.owner_id).or_default().push(m);
    }

    if groups.is_empty() {
        return None;
    }

    let mut entries: Vec<(u32, Vec<&'m Metadatum>)> = groups.into_iter().collect();
    entries.sort_unstable_by_key(|(id, _)| *id);

    let mut out = Vec::with_capacity(entries.len());
    for (id, group) in entries {
        out.push(parse_item(id, group.as_slice()));
    }
    Some(out)
}

#[inline]
pub fn parse_cv_and_user_params(
    allowed: &HashSet<&str>,
    metadata: &[&Metadatum],
) -> (Vec<CvParam>, Vec<UserParam>) {
    let mut cv_params = Vec::with_capacity(metadata.len());
    let mut user_params = Vec::new();

    let allow_all = allowed.is_empty();

    for m in metadata {
        let Some(acc) = m.accession.as_deref() else {
            continue;
        };
        let Some((prefix, _)) = acc.split_once(':') else {
            continue;
        };
        if prefix == CV_REF_ATTR {
            continue;
        }

        let value = value_to_opt_string(&m.value);
        let unit_accession = m.unit_accession.clone();
        let unit_cv_ref = unit_cv_ref(&unit_accession);

        if is_cv_prefix(prefix) {
            // TODO: Check if possible allowed accessions per tag
            let _ = (allow_all, allowed);

            let name = cv_table::get(acc)
                .and_then(|v| v.as_str())
                .unwrap_or(acc)
                .to_string();

            let unit_name = unit_accession
                .as_deref()
                .and_then(|ua| cv_table::get(ua).and_then(|v| v.as_str()))
                .map(|s| s.to_string());

            cv_params.push(CvParam {
                cv_ref: Some(prefix.to_string()),
                accession: Some(acc.to_string()),
                name,
                value,
                unit_cv_ref,
                unit_name,
                unit_accession,
            });
        } else {
            let name = cv_table::get(acc)
                .and_then(|v| v.as_str())
                .unwrap_or(acc)
                .to_string();

            user_params.push(UserParam {
                name,
                r#type: None,
                unit_accession,
                unit_cv_ref,
                unit_name: None,
                value,
            });
        }
    }

    (cv_params, user_params)
}
