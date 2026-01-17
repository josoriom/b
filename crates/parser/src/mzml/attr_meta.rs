use crate::mzml::structs::CvParam;

pub const CV_REF_ATTR: &str = "B000";

pub const CV_CODE_MS: u8 = 0;
pub const CV_CODE_UO: u8 = 1;
pub const CV_CODE_NCIT: u8 = 2;
pub const CV_CODE_PEFF: u8 = 3;
pub const CV_CODE_B000: u8 = 4;
pub const CV_CODE_UNKNOWN: u8 = 255;

/// Encode: string -> code
#[inline]
pub fn cv_ref_code_from_str(cv_ref: Option<&str>) -> u8 {
    match cv_ref {
        Some("MS") => CV_CODE_MS,
        Some("UO") => CV_CODE_UO,
        Some("NCIT") => CV_CODE_NCIT,
        Some("PEFF") => CV_CODE_PEFF,
        Some(CV_REF_ATTR) => CV_CODE_B000,
        _ => CV_CODE_UNKNOWN,
    }
}

/// Decode: code -> string
#[inline]
pub fn cv_ref_prefix_from_code(code: u8) -> Option<&'static str> {
    match code {
        CV_CODE_MS => Some("MS"),
        CV_CODE_UO => Some("UO"),
        CV_CODE_NCIT => Some("NCIT"),
        CV_CODE_PEFF => Some("PEFF"),
        CV_CODE_B000 => Some(CV_REF_ATTR), // "B000"
        _ => None,
    }
}

#[inline]
pub fn format_accession(cv_ref_code: u8, tail_raw: u32) -> Option<String> {
    let pref = cv_ref_prefix_from_code(cv_ref_code)?;

    match pref {
        "MS" => {
            let tail = normalize_ms_accession_tail(cv_ref_code, tail_raw);
            Some(format!("MS:{tail:07}"))
        }
        "UO" => Some(format!("UO:{tail_raw:07}")),
        "NCIT" => Some(format!("NCIT:C{tail_raw}")),
        x if x == CV_REF_ATTR => Some(format!("{CV_REF_ATTR}:{tail_raw}")),
        _ => Some(format!("{pref}:{tail_raw}")),
    }
}

const MS_ACCESSION_BASE: u32 = 1_000_000;
#[inline]
pub fn normalize_ms_accession_tail(cv_ref_code: u8, tail: u32) -> u32 {
    if cv_ref_code == CV_CODE_MS && tail != 0 && tail < MS_ACCESSION_BASE {
        MS_ACCESSION_BASE + tail
    } else {
        tail
    }
}

// Strings
pub const ACC_ATTR_ID: u32 = 9_910_001;
pub const ACC_ATTR_REF: u32 = 9_910_002;
pub const ACC_ATTR_NAME: u32 = 9_910_003;
pub const ACC_ATTR_LOCATION: u32 = 9_910_004;

// <cv> attributes
pub const ACC_ATTR_CV_ID: u32 = 9_900_001;
pub const ACC_ATTR_CV_FULL_NAME: u32 = 9_900_002;
pub const ACC_ATTR_CV_VERSION: u32 = 9_900_003;
pub const ACC_ATTR_CV_URI: u32 = 9_900_004;
pub const ACC_ATTR_LABEL: u32 = 9_910_020;

pub const ACC_ATTR_START_TIME_STAMP: u32 = 9_910_005;
pub const ACC_ATTR_DEFAULT_INSTRUMENT_CONFIGURATION_REF: u32 = 9_910_006;
pub const ACC_ATTR_DEFAULT_SOURCE_FILE_REF: u32 = 9_910_007;
pub const ACC_ATTR_SAMPLE_REF: u32 = 9_910_008;

pub const ACC_ATTR_DEFAULT_DATA_PROCESSING_REF: u32 = 9_910_009;
pub const ACC_ATTR_DATA_PROCESSING_REF: u32 = 9_910_010;
pub const ACC_ATTR_SOURCE_FILE_REF: u32 = 9_910_011;

pub const ACC_ATTR_NATIVE_ID: u32 = 9_910_012;
pub const ACC_ATTR_SPOT_ID: u32 = 9_910_013;
pub const ACC_ATTR_EXTERNAL_SPECTRUM_ID: u32 = 9_910_014;
pub const ACC_ATTR_SPECTRUM_REF: u32 = 9_910_015;

pub const ACC_ATTR_SCAN_SETTINGS_REF: u32 = 9_910_016;
pub const ACC_ATTR_INSTRUMENT_CONFIGURATION_REF: u32 = 9_910_017;

pub const ACC_ATTR_SOFTWARE_REF: u32 = 9_910_018;
pub const ACC_ATTR_VERSION: u32 = 9_910_019;

// Numbers
pub const ACC_ATTR_COUNT: u32 = 9_910_100;
pub const ACC_ATTR_ORDER: u32 = 9_910_101;
pub const ACC_ATTR_INDEX: u32 = 9_910_102;
pub const ACC_ATTR_SCAN_NUMBER: u32 = 9_910_103;
pub const ACC_ATTR_DEFAULT_ARRAY_LENGTH: u32 = 9_910_104;
pub const ACC_ATTR_ARRAY_LENGTH: u32 = 9_910_105;
pub const ACC_ATTR_ENCODED_LENGTH: u32 = 9_910_106;
pub const ACC_ATTR_MS_LEVEL: u32 = 9_910_107;

#[inline]
pub fn attr_cv_param(accession_tail: u32, value: &str) -> CvParam {
    CvParam {
        cv_ref: Some(CV_REF_ATTR.to_string()),
        accession: Some(format!("{}:{:07}", CV_REF_ATTR, accession_tail)),
        name: String::new(),
        value: Some(value.to_string()),
        unit_cv_ref: None,
        unit_name: None,
        unit_accession: None,
    }
}

#[inline]
pub fn push_attr_str(out: &mut Vec<CvParam>, accession_tail: u32, value: Option<&str>) {
    if let Some(v) = value {
        if !v.is_empty() {
            out.push(attr_cv_param(accession_tail, v));
        }
    }
}

#[inline]
pub fn push_attr_string(out: &mut Vec<CvParam>, accession_tail: u32, value: &str) {
    if !value.is_empty() {
        out.push(attr_cv_param(accession_tail, value));
    }
}

#[inline]
pub fn push_attr_u32(out: &mut Vec<CvParam>, accession_tail: u32, value: Option<u32>) {
    if let Some(v) = value {
        out.push(attr_cv_param(accession_tail, &v.to_string()));
    }
}

#[inline]
pub fn push_attr_usize(out: &mut Vec<CvParam>, accession_tail: u32, value: Option<usize>) {
    if let Some(v) = value {
        out.push(attr_cv_param(accession_tail, &v.to_string()));
    }
}

#[inline]
pub fn attr_key_from_tail(accession_tail: u32) -> Option<&'static str> {
    let key = match accession_tail {
        ACC_ATTR_ID => "id",
        ACC_ATTR_REF => "ref",
        ACC_ATTR_LOCATION => "location",
        ACC_ATTR_START_TIME_STAMP => "startTimeStamp",
        ACC_ATTR_DEFAULT_INSTRUMENT_CONFIGURATION_REF => "defaultInstrumentConfigurationRef",
        ACC_ATTR_DEFAULT_SOURCE_FILE_REF => "defaultSourceFileRef",
        ACC_ATTR_SAMPLE_REF => "sampleRef",
        ACC_ATTR_DEFAULT_DATA_PROCESSING_REF => "defaultDataProcessingRef",
        ACC_ATTR_DATA_PROCESSING_REF => "dataProcessingRef",
        ACC_ATTR_SOURCE_FILE_REF => "sourceFileRef",
        ACC_ATTR_NATIVE_ID => "nativeID",
        ACC_ATTR_SPOT_ID => "spotID",
        ACC_ATTR_EXTERNAL_SPECTRUM_ID => "externalSpectrumID",
        ACC_ATTR_SPECTRUM_REF => "spectrumRef",
        ACC_ATTR_SCAN_SETTINGS_REF => "scanSettingsRef",
        ACC_ATTR_INSTRUMENT_CONFIGURATION_REF => "instrumentConfigurationRef",
        ACC_ATTR_SOFTWARE_REF => "softwareRef",
        ACC_ATTR_VERSION => "version",
        ACC_ATTR_COUNT => "count",
        ACC_ATTR_ORDER => "order",
        ACC_ATTR_INDEX => "index",
        ACC_ATTR_SCAN_NUMBER => "scanNumber",
        ACC_ATTR_DEFAULT_ARRAY_LENGTH => "defaultArrayLength",
        ACC_ATTR_ARRAY_LENGTH => "arrayLength",
        ACC_ATTR_ENCODED_LENGTH => "encodedLength",
        ACC_ATTR_MS_LEVEL => "msLevel",
        _ => return None,
    };
    Some(key)
}
