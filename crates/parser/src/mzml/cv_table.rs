use once_cell::sync::Lazy;
use serde_json::Value;
use std::collections::HashMap;

static RAW_JSON: &str = include_str!("cv_table.json");

pub static TABLE: Lazy<HashMap<String, Value>> = Lazy::new(|| {
    let v: Value = serde_json::from_str(RAW_JSON).unwrap();
    let mut map = HashMap::new();
    if let Value::Object(obj) = v {
        for (k, val) in obj {
            map.insert(k, val);
        }
    }
    map
});

pub fn get(key: &str) -> Option<&Value> {
    TABLE.get(key)
}

pub const CV_CODE_MS: u8 = 0;
pub const CV_CODE_UO: u8 = 1;
pub const CV_CODE_NCIT: u8 = 2;
pub const CV_CODE_PEFF: u8 = 3;
pub const CV_CODE_ATTR: u8 = 4;
pub const CV_CODE_OTHER: u8 = 255;
