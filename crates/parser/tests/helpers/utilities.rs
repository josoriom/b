use std::{fs, path::PathBuf, sync::OnceLock};

use b::utilities::{
    decode::decode,
    mzml::{
        Chromatogram, ChromatogramList, CvParam, MzML, PrecursorList, Run, ScanList, Software,
        SoftwareParam, Spectrum, SpectrumDescription,
    },
    parse_mzml::parse_mzml,
};

#[derive(Debug, Clone, Copy)]
pub enum CvRefMode {
    #[allow(dead_code)]
    Strict,
    #[allow(dead_code)]
    AllowMissingMs,
}

#[allow(dead_code)]
pub fn mzml(cache: &'static OnceLock<MzML>, path: &str) -> &'static MzML {
    cache.get_or_init(|| {
        let bytes = load_mzml_bytes(path);
        parse_mzml(&bytes, false).unwrap_or_else(|e| panic!("parse_mzml failed: {e}"))
    })
}

#[allow(dead_code)]
pub fn parse_b(cache: &'static OnceLock<MzML>, path: &str) -> &'static MzML {
    cache.get_or_init(|| {
        let bytes = load_mzml_bytes(path);
        decode(&bytes).unwrap_or_else(|e| panic!("parse_mzml failed: {e}"))
    })
}

pub fn load_mzml_bytes(path: &str) -> Vec<u8> {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(path);
    fs::read(&full).unwrap_or_else(|e| panic!("cannot read {:?}: {}", full, e))
}

#[allow(dead_code)]
pub fn spectrum_by_id<'a>(mzml: &'a MzML, id: &str) -> &'a Spectrum {
    let sl = mzml
        .run
        .spectrum_list
        .as_ref()
        .expect("spectrumList parsed");
    sl.spectra
        .iter()
        .find(|s| s.id == id)
        .unwrap_or_else(|| panic!("spectrum {id} not found"))
}

#[allow(dead_code)]
pub fn spectrum_by_index<'a>(mzml: &'a MzML, idx: usize) -> &'a Spectrum {
    let sl = mzml
        .run
        .spectrum_list
        .as_ref()
        .expect("spectrumList parsed");
    sl.spectra
        .get(idx)
        .unwrap_or_else(|| panic!("spectrum index {idx} not found"))
}

pub fn spectrum_description(s: &Spectrum) -> &SpectrumDescription {
    s.spectrum_description
        .as_ref()
        .expect("spectrumDescription parsed")
}

pub fn spectrum_scan_list(s: &Spectrum) -> &ScanList {
    if let Some(sd) = s.spectrum_description.as_ref() {
        if let Some(sl) = sd.scan_list.as_ref() {
            return sl;
        }
    }
    s.scan_list.as_ref().expect("scanList parsed")
}

pub fn spectrum_precursor_list(s: &Spectrum) -> Option<&PrecursorList> {
    if let Some(sd) = s.spectrum_description.as_ref() {
        if sd.precursor_list.is_some() {
            return sd.precursor_list.as_ref();
        }
    }
    s.precursor_list.as_ref()
}

#[allow(dead_code)]
pub fn chromatogram_list(run: &Run) -> &ChromatogramList {
    run.chromatogram_list
        .as_ref()
        .expect("chromatogramList parsed")
}

#[allow(dead_code)]
pub fn chromatogram<'a>(cl: &'a ChromatogramList, id: &str) -> &'a Chromatogram {
    cl.chromatograms
        .iter()
        .find(|c| c.id == id)
        .unwrap_or_else(|| panic!("chromatogram {id} not found"))
}

#[allow(dead_code)]
pub fn cv_by_name<'a>(cv_params: &'a [CvParam], name: &str) -> Option<&'a CvParam> {
    cv_params.iter().find(|cv| cv.name == name)
}

#[allow(dead_code)]
pub fn assert_cv_absent(cv_params: &[CvParam], name: &str) {
    let got = cv_by_name(cv_params, name);
    assert!(
        got.is_none(),
        "expected cvParam {name:?} to be absent, got: {got:?}"
    );
}

pub fn assert_cv_ref(policy: CvRefMode, got: Option<&str>, expected: &str, ctx: &str) {
    match (policy, got) {
        (_, Some(v)) if v == expected => {}
        (CvRefMode::AllowMissingMs, Some("")) if expected == "MS" => {}
        (CvRefMode::AllowMissingMs, None) if expected == "MS" => {}
        _ => panic!("wrong cv_ref for {ctx}: got {got:?}, expected {expected:?}"),
    }
}

pub fn assert_software(
    policy: CvRefMode,
    sw: &Software,
    cv_ref: &str,
    accession: &str,
    name: &str,
    version: Option<&str>,
) {
    if let Some(p) = sw.software_param.get(0) {
        assert_software_param(policy, p, cv_ref, accession, name, version);
        return;
    }

    assert_cv(policy, &sw.cv_param, name, accession, cv_ref, None, None);
    assert_eq!(sw.version.as_deref(), version, "wrong version for {name}");
}

pub fn assert_software_param(
    policy: CvRefMode,
    p: &SoftwareParam,
    cv_ref: &str,
    accession: &str,
    name: &str,
    version: Option<&str>,
) {
    assert_cv_ref(policy, p.cv_ref.as_deref(), cv_ref, name);
    assert_eq!(
        p.accession.as_str(),
        accession,
        "wrong accession for {name}"
    );
    assert_eq!(p.name.as_str(), name, "wrong name for {name}");
    assert_eq!(p.version.as_deref(), version, "wrong version for {name}");
}

pub fn assert_cv(
    policy: CvRefMode,
    cv_params: &[CvParam],
    name: &str,
    accession: &str,
    cv_ref: &str,
    value: Option<&str>,
    unit_name: Option<&str>,
) {
    let cv = cv_params
        .iter()
        .find(|cv| cv.name == name)
        .unwrap_or_else(|| panic!("cvParam with name {name} not found"));

    assert_eq!(
        cv.accession.as_deref(),
        Some(accession),
        "wrong accession for {name}"
    );

    assert_cv_ref(policy, cv.cv_ref.as_deref(), cv_ref, name);

    match value {
        Some(v) if v.is_empty() => {
            assert!(
                cv.value.as_deref().unwrap_or("").is_empty(),
                "wrong value for {name}: {:?}",
                cv.value
            );
        }
        Some(v) => assert_eq!(cv.value.as_deref(), Some(v), "wrong value for {name}"),
        None => assert!(
            cv.value.is_none(),
            "expected no value for {name}, got {:?}",
            cv.value
        ),
    }

    match unit_name {
        Some(u) => assert_eq!(
            cv.unit_name.as_deref(),
            Some(u),
            "wrong unit_name for {name}"
        ),
        None => assert!(
            cv.unit_name.is_none(),
            "expected no unit_name for {name}, got {:?}",
            cv.unit_name
        ),
    }
}

#[allow(dead_code)]
pub fn assert_cv_f64(
    policy: CvRefMode,
    cv_params: &[CvParam],
    name: &str,
    accession: &str,
    cv_ref: &str,
    expected: f64,
    unit_name: Option<&str>,
) {
    let cv = cv_params
        .iter()
        .find(|cv| cv.name == name)
        .unwrap_or_else(|| panic!("cvParam with name {name} not found"));

    assert_eq!(
        cv.accession.as_deref(),
        Some(accession),
        "wrong accession for {name}"
    );

    assert_cv_ref(policy, cv.cv_ref.as_deref(), cv_ref, name);

    let s = cv
        .value
        .as_deref()
        .unwrap_or_else(|| panic!("expected numeric value for {name}, got None"));

    let got: f64 = s
        .parse()
        .unwrap_or_else(|_| panic!("failed to parse numeric value for {name}: {s:?}"));

    let diff = (got - expected).abs();
    let tol = expected.abs().max(1.0) * 1e-6;
    assert!(
        diff <= tol,
        "wrong numeric value for {name}: got {got}, expected {expected} (tol {tol})"
    );

    match unit_name {
        Some(u) => assert_eq!(
            cv.unit_name.as_deref(),
            Some(u),
            "wrong unit_name for {name}"
        ),
        None => assert!(
            cv.unit_name.is_none(),
            "expected no unit_name for {name}, got {:?}",
            cv.unit_name
        ),
    }
}
