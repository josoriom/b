use std::{fs, path::PathBuf};

use b::utilities::{
    mzml::{CvParam, MzML},
    parse_mzml,
};

fn load_mzml_bytes() -> Vec<u8> {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/mzml/tiny.msdata.mzML0.99.9.mzML");

    fs::read(&path).unwrap_or_else(|e| panic!("cannot read {:?}: {}", path, e))
}

fn parse_1min_full() -> MzML {
    let bytes = load_mzml_bytes();
    parse_mzml(&bytes, false).unwrap_or_else(|e| panic!("parse_mzml failed: {e}"))
}

#[test]
fn tiny_msdata_mzml0_99_9_header_sections() {
    let mzml = parse_1min_full();

    let cv_list = mzml.cv_list.as_ref().expect("cvList parsed");
    assert_eq!(cv_list.cv.len(), 1);
    let cv0 = &cv_list.cv[0];
    assert_eq!(cv0.id, "MS");
    assert_eq!(
        cv0.full_name.as_deref(),
        Some("Proteomics Standards Initiative Mass Spectrometry Ontology")
    );
    assert_eq!(cv0.version.as_deref(), Some("2.0.2"));
    assert_eq!(
        cv0.uri.as_deref(),
        Some("http://psidev.sourceforge.net/ms/xml/mzdata/psi-ms.2.0.2.obo")
    );

    let file_desc = &mzml.file_description;

    assert_eq!(file_desc.file_content.cv_params.len(), 1);
    assert_cv(
        &file_desc.file_content.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );

    assert_eq!(file_desc.source_file_list.source_file.len(), 1);
    let sf0 = &file_desc.source_file_list.source_file[0];
    assert_eq!(sf0.id, "sf1");
    assert_eq!(sf0.name, "tiny1.RAW");
    assert_eq!(sf0.location, "file://F:/data/Exp01");
    assert_cv(
        &sf0.cv_param,
        "Xcalibur RAW file",
        "MS:1000563",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &sf0.cv_param,
        "SHA-1",
        "MS:1000569",
        "MS",
        Some("71be39fb2700ab2f3c8b2234b91274968b6899b1"),
        None,
    );

    let rpgl = mzml
        .referenceable_param_group_list
        .as_ref()
        .expect("referenceableParamGroupList parsed");
    assert_eq!(rpgl.referenceable_param_groups.len(), 2);

    let r0 = &rpgl.referenceable_param_groups[0];
    assert_eq!(r0.id, "CommonMS1SpectrumParams");
    assert_cv(
        &r0.cv_params,
        "positive scan",
        "MS:1000130",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &r0.cv_params,
        "full scan",
        "MS:1000498",
        "MS",
        Some(""),
        None,
    );

    let r1 = &rpgl.referenceable_param_groups[1];
    assert_eq!(r1.id, "CommonMS2SpectrumParams");
    assert_cv(
        &r1.cv_params,
        "positive scan",
        "MS:1000130",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &r1.cv_params,
        "full scan",
        "MS:1000498",
        "MS",
        Some(""),
        None,
    );

    let sample_list = mzml.sample_list.as_ref().expect("sampleList parsed");
    assert_eq!(sample_list.samples.len(), 1);
    let sample0 = &sample_list.samples[0];
    assert_eq!(sample0.id, "sp1");
    assert_eq!(sample0.name, "Sample1");

    let inst_list = mzml
        .instrument_list
        .as_ref()
        .expect("instrumentList parsed");
    assert_eq!(inst_list.instrument.len(), 1);
    let inst0 = &inst_list.instrument[0];
    assert_eq!(inst0.id, "LCQDeca");
    assert_cv(
        &inst0.cv_param,
        "LCQ Deca",
        "MS:1000554",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &inst0.cv_param,
        "instrument serial number",
        "MS:1000529",
        "MS",
        Some("23433"),
        None,
    );

    let cl0 = inst0.component_list.as_ref().expect("componentList parsed");
    assert_eq!(cl0.source.len(), 1);
    assert_eq!(cl0.analyzer.len(), 1);
    assert_eq!(cl0.detector.len(), 1);

    let src = &cl0.source[0];
    assert_eq!(src.order, Some(1));
    assert_cv(
        &src.cv_param,
        "nanoelectrospray",
        "MS:1000398",
        "MS",
        Some(""),
        None,
    );

    let an = &cl0.analyzer[0];
    assert_eq!(an.order, Some(2));
    assert_cv(
        &an.cv_param,
        "quadrupole ion trap",
        "MS:1000082",
        "MS",
        Some(""),
        None,
    );

    let det = &cl0.detector[0];
    assert_eq!(det.order, Some(3));
    assert_cv(
        &det.cv_param,
        "electron multiplier",
        "MS:1000253",
        "MS",
        Some(""),
        None,
    );

    let sw_list = mzml.software_list.as_ref().expect("softwareList parsed");
    assert_eq!(sw_list.software.len(), 3);

    let sw0 = &sw_list.software[0];
    assert_eq!(sw0.id, "Bioworks");
    assert_eq!(sw0.version.as_deref(), Some("3.3.1 sp1"));
    assert_cv(&sw0.cv_param, "Bioworks", "MS:1000533", "MS", None, None);

    let sw1 = &sw_list.software[1];
    assert_eq!(sw1.id, "ReAdW");
    assert_eq!(sw1.version.as_deref(), Some("1.0"));
    assert_cv(&sw1.cv_param, "ReAdW", "MS:1000541", "MS", None, None);

    let sw2 = &sw_list.software[2];
    assert_eq!(sw2.id, "Xcalibur");
    assert_eq!(sw2.version.as_deref(), Some("2.0.5"));
    assert_cv(&sw2.cv_param, "Xcalibur", "MS:1000532", "MS", None, None);

    let dp_list = mzml
        .data_processing_list
        .as_ref()
        .expect("dataProcessingList parsed");
    assert_eq!(dp_list.data_processing.len(), 2);

    let dp0 = &dp_list.data_processing[0];
    assert_eq!(dp0.id, "XcaliburProcessing");
    assert_eq!(dp0.processing_method.len(), 1);
    let m0 = &dp0.processing_method[0];
    assert_cv(
        &m0.cv_param,
        "deisotoping",
        "MS:1000033",
        "MS",
        Some("false"),
        None,
    );
    assert_cv(
        &m0.cv_param,
        "charge deconvolution",
        "MS:1000034",
        "MS",
        Some("false"),
        None,
    );
    assert_cv(
        &m0.cv_param,
        "peak picking",
        "MS:1000035",
        "MS",
        Some("true"),
        None,
    );

    let dp1 = &dp_list.data_processing[1];
    assert_eq!(dp1.id, "ReAdWConversion");
    assert_eq!(dp1.processing_method.len(), 1);
    let m1 = &dp1.processing_method[0];
    assert_cv(
        &m1.cv_param,
        "Conversion to mzML",
        "MS:1000544",
        "MS",
        Some(""),
        None,
    );

    let ss_list = mzml
        .scan_settings_list
        .as_ref()
        .expect("scanSettingsList parsed");
    assert_eq!(ss_list.scan_settings.len(), 1);

    let acq0 = &ss_list.scan_settings[0];
    assert_eq!(acq0.id.as_deref(), Some("aS1"));
    assert_eq!(
        acq0.instrument_configuration_ref.as_deref(),
        Some("LCQDeca")
    );

    let sfrefl = acq0
        .source_file_ref_list
        .as_ref()
        .expect("sourceFileRefList parsed");
    assert_eq!(sfrefl.source_file_refs.len(), 1);
    let sfref0 = &sfrefl.source_file_refs[0];
    assert_eq!(sfref0.r#ref, "SF2");

    let asf = file_desc
        .source_file_list
        .source_file
        .iter()
        .find(|sf| sf.id == "SF2")
        .unwrap_or_else(|| panic!("sourceFile SF2 not found"));
    assert_eq!(asf.id, "SF2");
    assert_eq!(asf.name, "parameters.par");
    assert_eq!(asf.location, "file:///C:/settings/");

    let tl = acq0.target_list.as_ref().expect("targetList parsed");
    assert_eq!(tl.targets.len(), 2);

    let t0 = &tl.targets[0];
    assert_cv(
        &t0.cv_params,
        "precursorMz",
        "MS:1000xxx",
        "MS",
        Some("123.456"),
        None,
    );
    assert_cv(
        &t0.cv_params,
        "fragmentMz",
        "MS:1000xxx",
        "MS",
        Some("456.789"),
        None,
    );
    assert_cv(
        &t0.cv_params,
        "dwell time",
        "MS:1000xxx",
        "MS",
        Some("1"),
        Some("seconds"),
    );
    assert_cv(
        &t0.cv_params,
        "active time",
        "MS:1000xxx",
        "MS",
        Some("0.5"),
        Some("seconds"),
    );
}

#[test]
fn tiny_msdata_mzml0_99_9_first_spectrum() {
    let mzml = parse_1min_full();
    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    let s0 = &sl.spectra[0];
    assert_eq!(s0.index, Some(0));
    assert_eq!(s0.id, "S19");

    assert_eq!(s0.cv_params.len(), 8);

    assert_cv(
        &s0.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "ms level",
        "MS:1000511",
        "MS",
        Some("1"),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "centroid mass spectrum",
        "MS:1000127",
        "MS",
        Some(""),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "lowest m/z value",
        "MS:1000528",
        "MS",
        Some("400.39"),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "highest m/z value",
        "MS:1000527",
        "MS",
        Some("1795.56"),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "base peak m/z",
        "MS:1000504",
        "MS",
        Some("445.347"),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        Some("120053"),
        None,
    );

    assert_cv(
        &s0.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        Some("16675500"),
        None,
    );

    let pl = s0
        .precursor_list
        .as_ref()
        .map(|p| p.precursors.len())
        .unwrap_or(0);
    assert_eq!(pl, 0);

    let scl = s0.scan_list.as_ref().expect("scanList parsed");
    assert_eq!(scl.scans.len(), 1);
    let scan0 = &scl.scans[0];

    assert_cv(
        &scan0.cv_params,
        "scan time",
        "MS:1000016",
        "MS",
        Some("5.8905"),
        Some("minute"),
    );

    assert_cv(
        &scan0.cv_params,
        "filter string",
        "MS:1000512",
        "MS",
        Some("+ c NSI Full ms [ 400.00-1800.00]"),
        None,
    );

    let swl = scan0
        .scan_window_list
        .as_ref()
        .expect("scanWindowList parsed");
    assert_eq!(swl.scan_windows.len(), 1);
    let win0 = &swl.scan_windows[0];

    assert_cv(
        &win0.cv_params,
        "scan m/z lower limit",
        "MS:1000501",
        "MS",
        Some("400"),
        None,
    );

    assert_cv(
        &win0.cv_params,
        "scan m/z upper limit",
        "MS:1000500",
        "MS",
        Some("1800"),
        None,
    );

    let bal = s0
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    let ba0 = &bal.binary_data_arrays[0];
    assert_eq!(ba0.cv_params.len(), 3);
    assert_cv(
        &ba0.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba0.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba0.cv_params,
        "m/z array",
        "MS:1000514",
        "MS",
        Some(""),
        None,
    );

    let ba1 = &bal.binary_data_arrays[1];
    assert_eq!(ba1.cv_params.len(), 3);
    assert_cv(
        &ba1.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba1.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba1.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        None,
    );
}

#[test]
fn tiny_msdata_mzml0_99_9_second_spectrum() {
    let mzml = parse_1min_full();
    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    let s1 = &sl.spectra[1];
    assert_eq!(s1.index, Some(1));
    assert_eq!(s1.id, "S20");

    assert_eq!(s1.cv_params.len(), 8);

    assert_cv(
        &s1.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "ms level",
        "MS:1000511",
        "MS",
        Some("2"),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "centroid mass spectrum",
        "MS:1000127",
        "MS",
        Some(""),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "lowest m/z value",
        "MS:1000528",
        "MS",
        Some("320.39"),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "highest m/z value",
        "MS:1000527",
        "MS",
        Some("1003.56"),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "base peak m/z",
        "MS:1000504",
        "MS",
        Some("456.347"),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        Some("23433"),
        None,
    );

    assert_cv(
        &s1.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        Some("16675500"),
        None,
    );

    let pl = s1.precursor_list.as_ref().expect("precursorList parsed");
    assert_eq!(pl.precursors.len(), 1);

    let p0 = &pl.precursors[0];
    assert_eq!(p0.spectrum_ref.as_deref(), Some("S19"));

    let iw = p0
        .isolation_window
        .as_ref()
        .expect("isolationWindow parsed");
    assert_cv(
        &iw.cv_params,
        "isolation center m/z",
        "MS:1000xxx",
        "MS",
        Some("445.34"),
        None,
    );
    assert_cv(
        &iw.cv_params,
        "isolation half width",
        "MS:1000xxx",
        "MS",
        Some("2.0"),
        None,
    );

    let sil = p0
        .selected_ion_list
        .as_ref()
        .expect("selectedIonList parsed");
    assert_eq!(sil.selected_ions.len(), 1);
    let ion0 = &sil.selected_ions[0];
    assert_eq!(ion0.cv_params.len(), 2);
    assert_cv(
        &ion0.cv_params,
        "m/z",
        "MS:1000040",
        "MS",
        Some("445.34"),
        None,
    );
    assert_cv(
        &ion0.cv_params,
        "charge state",
        "MS:1000041",
        "MS",
        Some("2"),
        None,
    );

    let act = p0.activation.as_ref().expect("activation parsed");
    assert_cv(
        &act.cv_params,
        "collision-induced dissociation",
        "MS:1000133",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &act.cv_params,
        "collision energy",
        "MS:1000045",
        "MS",
        Some("35"),
        Some("electron volt"),
    );

    let scl = s1.scan_list.as_ref().expect("scanList parsed");
    assert_eq!(scl.scans.len(), 1);
    let scan1 = &scl.scans[0];

    assert_cv(
        &scan1.cv_params,
        "scan time",
        "MS:1000016",
        "MS",
        Some("5.9905"),
        Some("minute"),
    );

    assert_cv(
        &scan1.cv_params,
        "filter string",
        "MS:1000512",
        "MS",
        Some("+ c d Full ms2  445.35@cid35.00 [ 110.00-905.00]"),
        None,
    );

    let swl = scan1
        .scan_window_list
        .as_ref()
        .expect("scanWindowList parsed");
    assert_eq!(swl.scan_windows.len(), 1);
    let win1 = &swl.scan_windows[0];

    assert_cv(
        &win1.cv_params,
        "scan m/z lower limit",
        "MS:1000501",
        "MS",
        Some("110"),
        None,
    );

    assert_cv(
        &win1.cv_params,
        "scan m/z upper limit",
        "MS:1000500",
        "MS",
        Some("905"),
        None,
    );

    let bal = s1
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    let ba0 = &bal.binary_data_arrays[0];
    assert_eq!(ba0.cv_params.len(), 3);
    assert_cv(
        &ba0.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba0.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba0.cv_params,
        "m/z array",
        "MS:1000514",
        "MS",
        Some(""),
        None,
    );
    assert_eq!(ba0.array_length, Some(20));
    assert_eq!(ba0.encoded_length, Some(216));

    let ba1 = &bal.binary_data_arrays[1];
    assert_eq!(ba1.cv_params.len(), 3);
    assert_cv(
        &ba1.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba1.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        &ba1.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        None,
    );
    assert_eq!(ba1.array_length, Some(20));
    assert_eq!(ba1.encoded_length, Some(216));
}

fn assert_cv(
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
    assert_eq!(
        cv.cv_ref.as_deref(),
        Some(cv_ref),
        "wrong cv_ref for {name}"
    );

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
