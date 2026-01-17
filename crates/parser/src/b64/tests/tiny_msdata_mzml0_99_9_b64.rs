use crate::{
    mzml::structs::MzML,
    utilities::test::{
        CvRefMode, assert_cv, assert_software_param, parse_b, spectrum_description,
        spectrum_precursor_list, spectrum_scan_list,
    },
};

use std::sync::OnceLock;

static MZML_CACHE: OnceLock<MzML> = OnceLock::new();

const PATH: &str = "data/b64/tiny.msdata.mzML0.99.9.b64";
const CV_REF_MODE: CvRefMode = CvRefMode::Strict;

fn mzml_cached() -> &'static MzML {
    parse_b(&MZML_CACHE, PATH)
}

#[test]
fn tiny_msdata_mzml0_99_9_header_sections() {
    let mzml = parse_b(&MZML_CACHE, PATH);

    // cvList
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

    // fileDescription
    let file_desc = &mzml.file_description;

    // fileContent
    let cv_params = &file_desc.file_content.cv_params;
    assert_eq!(cv_params.len(), 1);
    assert_cv(
        CV_REF_MODE,
        cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );

    // sourceFileList

    assert_eq!(file_desc.source_file_list.source_file.len(), 1);
    let sf0 = &file_desc.source_file_list.source_file[0];
    assert_eq!(sf0.id, "sf1");
    assert_eq!(sf0.name, "tiny1.RAW");
    assert_eq!(sf0.location, "file://F:/data/Exp01");
    assert_cv(
        CV_REF_MODE,
        &sf0.cv_param,
        "Thermo RAW format",
        "MS:1000563",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sf0.cv_param,
        "SHA-1",
        "MS:1000569",
        "MS",
        Some("71be39fb2700ab2f3c8b2234b91274968b6899b1"),
        None,
    );

    // referenceableParamGroupList
    let rpgl = mzml
        .referenceable_param_group_list
        .as_ref()
        .expect("referenceableParamGroupList parsed");
    assert_eq!(rpgl.referenceable_param_groups.len(), 2);
    for (idx, id) in [
        (0usize, "CommonMS1SpectrumParams"),
        (1usize, "CommonMS2SpectrumParams"),
    ] {
        let g = &rpgl.referenceable_param_groups[idx];
        assert_eq!(g.id, id);
        assert_cv(
            CV_REF_MODE,
            &g.cv_params,
            "positive scan",
            "MS:1000130",
            "MS",
            Some(""),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &g.cv_params,
            "full scan",
            "MS:1000498",
            "MS",
            Some(""),
            None,
        );
    }

    // sampleList
    let sample_list = mzml.sample_list.as_ref().expect("sampleList parsed");
    assert_eq!(sample_list.samples.len(), 1);
    let sample0 = &sample_list.samples[0];
    assert_eq!(sample0.id, "sp1");
    assert_eq!(sample0.name, "Sample1");

    // instrumentList
    let inst_list = mzml
        .instrument_list
        .as_ref()
        .expect("instrumentList parsed");
    assert_eq!(inst_list.instrument.len(), 1);
    let inst0 = &inst_list.instrument[0];
    assert_eq!(inst0.id, "LCQDeca");
    assert_cv(
        CV_REF_MODE,
        &inst0.cv_param,
        "LCQ Deca",
        "MS:1000554",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &inst0.cv_param,
        "instrument serial number",
        "MS:1000529",
        "MS",
        Some("23433"),
        None,
    );

    // componentList
    let cl0 = inst0.component_list.as_ref().expect("componentList parsed");
    assert_eq!(cl0.source.len(), 1);
    assert_eq!(cl0.analyzer.len(), 1);
    assert_eq!(cl0.detector.len(), 1);

    let src = &cl0.source[0];
    assert_eq!(src.order, Some(1));
    assert_cv(
        CV_REF_MODE,
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
        CV_REF_MODE,
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
        CV_REF_MODE,
        &det.cv_param,
        "electron multiplier",
        "MS:1000253",
        "MS",
        Some(""),
        None,
    );

    // softwareList
    let sw_list = mzml.software_list.as_ref().expect("softwareList parsed");
    assert_eq!(sw_list.software.len(), 3);

    let sw0 = &sw_list.software[0];
    assert_eq!(sw0.id, "Bioworks");
    assert_eq!(sw0.cv_param.len(), 0);
    assert_eq!(sw0.software_param.len(), 1);
    assert_software_param(
        CV_REF_MODE,
        &sw0.software_param[0],
        "MS",
        "MS:1000533",
        "Bioworks",
        Some("3.3.1 sp1"),
    );

    let sw1 = &sw_list.software[1];
    assert_eq!(sw1.id, "ReAdW");
    assert_eq!(sw1.cv_param.len(), 0);
    assert_eq!(sw1.software_param.len(), 1);
    assert_software_param(
        CV_REF_MODE,
        &sw1.software_param[0],
        "MS",
        "MS:1000541",
        "ReAdW",
        Some("1"),
    );

    let sw2 = &sw_list.software[2];
    assert_eq!(sw2.id, "Xcalibur");
    assert_eq!(sw2.cv_param.len(), 0);
    assert_eq!(sw2.software_param.len(), 1);
    assert_software_param(
        CV_REF_MODE,
        &sw2.software_param[0],
        "MS",
        "MS:1000532",
        "Xcalibur",
        Some("2.0.5"),
    );

    // dataProcessingList
    let dp_list = mzml
        .data_processing_list
        .as_ref()
        .expect("dataProcessingList parsed");
    assert_eq!(dp_list.data_processing.len(), 2);

    let dp0 = &dp_list.data_processing[0];
    assert_eq!(dp0.id, "XcaliburProcessing");
    assert_eq!(dp0.processing_method.len(), 1);
    let m0 = &dp0.processing_method[0];
    for (name, accession, value) in [
        ("deisotoping", "MS:1000033", "false"),
        ("charge deconvolution", "MS:1000034", "false"),
        ("peak picking", "MS:1000035", "true"),
    ] {
        assert_cv(
            CV_REF_MODE,
            &m0.cv_param,
            name,
            accession,
            "MS",
            Some(value),
            None,
        );
    }

    let dp1 = &dp_list.data_processing[1];
    assert_eq!(dp1.id, "ReAdWConversion");
    assert_eq!(dp1.processing_method.len(), 1);
    let m1 = &dp1.processing_method[0];
    assert_cv(
        CV_REF_MODE,
        &m1.cv_param,
        "Conversion to mzML",
        "MS:1000544",
        "MS",
        Some(""),
        None,
    );

    // scanSettingsList
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

    println!("---::>>{:#?}", ss_list);

    if let Some(sfrefl) = acq0.source_file_ref_list.as_ref() {
        assert_eq!(sfrefl.source_file_refs.len(), 1);
        assert!(!sfrefl.source_file_refs[0].r#ref.is_empty());
    }

    let tl = acq0.target_list.as_ref().expect("targetList parsed");
    assert_eq!(tl.targets.len(), 2);

    let precursor_vals = ["123.456", "231.673"];
    let fragment_vals = ["456.789", "566.328"];

    for (i, t) in tl.targets.iter().enumerate() {
        assert_cv(
            CV_REF_MODE,
            &t.cv_params,
            "isolation window target m/z",
            "MS:1000827",
            "MS",
            Some(precursor_vals[i]),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &t.cv_params,
            "product ion m/z",
            "MS:1001225",
            "MS",
            Some(fragment_vals[i]),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &t.cv_params,
            "dwell time",
            "MS:1000502",
            "MS",
            Some("1"),
            Some("second"),
        );
        assert_cv(
            CV_REF_MODE,
            &t.cv_params,
            "completion time",
            "MS:1000747",
            "MS",
            Some("0.5"),
            Some("second"),
        );
    }
}

#[test]
fn tiny_msdata_mzml0_99_9_first_spectrum() {
    let mzml = mzml_cached();
    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);
    let s0 = &sl.spectra[0];

    // spectrum
    assert_eq!(s0.index, Some(0));
    assert_eq!(s0.id, "S19");
    assert_eq!(s0.cv_params.len(), 2);
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "ms level",
        "MS:1000511",
        "MS",
        Some("1"),
        None,
    );

    // spectrumDescription
    let sd = spectrum_description(s0);
    assert_eq!(sd.cv_params.len(), 6);
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "centroid spectrum",
        "MS:1000127",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "lowest observed m/z",
        "MS:1000528",
        "MS",
        Some("400.39"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "highest observed m/z",
        "MS:1000527",
        "MS",
        Some("1795.56"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak m/z",
        "MS:1000504",
        "MS",
        Some("445.347"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        Some("120053"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        Some("16675500"),
        None,
    );

    // precursorList
    let pl = spectrum_precursor_list(s0)
        .map(|p| p.precursors.len())
        .unwrap_or(0);
    assert_eq!(pl, 0);

    // scan
    let scl = spectrum_scan_list(s0);
    assert_eq!(scl.scans.len(), 1);
    let scan0 = &scl.scans[0];
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "scan start time",
        "MS:1000016",
        "MS",
        Some("5.8905"),
        Some("minute"),
    );
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "filter string",
        "MS:1000512",
        "MS",
        Some("+ c NSI Full ms [ 400.00-1800.00]"),
        None,
    );

    // scanWindowList
    let swl = scan0
        .scan_window_list
        .as_ref()
        .expect("scanWindowList parsed");
    assert_eq!(swl.scan_windows.len(), 1);
    let win0 = &swl.scan_windows[0];
    assert_cv(
        CV_REF_MODE,
        &win0.cv_params,
        "scan window lower limit",
        "MS:1000501",
        "MS",
        Some("400"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &win0.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        Some("1800"),
        None,
    );

    // binaryDataArrayList
    let bal = s0
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    for (i, name, accession) in [
        (0usize, "m/z array", "MS:1000514"),
        (1usize, "intensity array", "MS:1000515"),
    ] {
        let ba = &bal.binary_data_arrays[i];
        assert_eq!(ba.cv_params.len(), 3);
        assert_cv(
            CV_REF_MODE,
            &ba.cv_params,
            "64-bit float",
            "MS:1000523",
            "MS",
            Some(""),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &ba.cv_params,
            "no compression",
            "MS:1000576",
            "MS",
            Some(""),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &ba.cv_params,
            name,
            accession,
            "MS",
            Some(""),
            None,
        );
    }
}

#[test]
fn tiny_msdata_mzml0_99_9_second_spectrum() {
    let mzml = mzml_cached();
    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);
    let s1 = &sl.spectra[1];

    // spectrum
    assert_eq!(s1.index, Some(1));
    assert_eq!(s1.id, "S20");
    assert_eq!(s1.cv_params.len(), 2);
    assert_cv(
        CV_REF_MODE,
        &s1.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s1.cv_params,
        "ms level",
        "MS:1000511",
        "MS",
        Some("2"),
        None,
    );

    // spectrumDescription
    let sd = spectrum_description(s1);
    assert_eq!(sd.cv_params.len(), 6);
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "centroid spectrum",
        "MS:1000127",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "lowest observed m/z",
        "MS:1000528",
        "MS",
        Some("320.39"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "highest observed m/z",
        "MS:1000527",
        "MS",
        Some("1003.56"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak m/z",
        "MS:1000504",
        "MS",
        Some("456.347"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        Some("23433"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        Some("16675500"),
        None,
    );

    // precursorList
    let pl = spectrum_precursor_list(s1).expect("precursorList parsed");
    assert_eq!(pl.precursors.len(), 1);
    let p0 = &pl.precursors[0];
    assert_eq!(p0.spectrum_ref.as_deref(), Some("S19"));

    // isolationWindow
    let iw = p0
        .isolation_window
        .as_ref()
        .expect("isolationWindow parsed");
    assert_cv(
        CV_REF_MODE,
        &iw.cv_params,
        "isolation window target m/z",
        "MS:1000827",
        "MS",
        Some("445.34"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &iw.cv_params,
        "isolation window lower offset",
        "MS:1000828",
        "MS",
        Some("2"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &iw.cv_params,
        "isolation window upper offset",
        "MS:1000829",
        "MS",
        Some("2"),
        None,
    );

    // selectedIonList
    let sil = p0
        .selected_ion_list
        .as_ref()
        .expect("selectedIonList parsed");
    assert_eq!(sil.selected_ions.len(), 1);
    let ion0 = &sil.selected_ions[0];
    assert_eq!(ion0.cv_params.len(), 2);
    assert_cv(
        CV_REF_MODE,
        &ion0.cv_params,
        "m/z",
        "MS:1000040",
        "MS",
        Some("445.34"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &ion0.cv_params,
        "charge state",
        "MS:1000041",
        "MS",
        Some("2"),
        None,
    );

    // activation
    let act = p0.activation.as_ref().expect("activation parsed");
    assert_cv(
        CV_REF_MODE,
        &act.cv_params,
        "collision-induced dissociation",
        "MS:1000133",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &act.cv_params,
        "collision energy",
        "MS:1000045",
        "MS",
        Some("35"),
        Some("electron volt"),
    );

    // scan
    let scl = spectrum_scan_list(s1);
    assert_eq!(scl.scans.len(), 1);
    let scan1 = &scl.scans[0];
    assert_cv(
        CV_REF_MODE,
        &scan1.cv_params,
        "scan start time",
        "MS:1000016",
        "MS",
        Some("5.9905"),
        Some("minute"),
    );
    assert_cv(
        CV_REF_MODE,
        &scan1.cv_params,
        "filter string",
        "MS:1000512",
        "MS",
        Some("+ c d Full ms2  445.35@cid35.00 [ 110.00-905.00]"),
        None,
    );

    // scanWindowList
    let swl = scan1
        .scan_window_list
        .as_ref()
        .expect("scanWindowList parsed");
    assert_eq!(swl.scan_windows.len(), 1);
    let win1 = &swl.scan_windows[0];
    assert_cv(
        CV_REF_MODE,
        &win1.cv_params,
        "scan window lower limit",
        "MS:1000501",
        "MS",
        Some("110"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &win1.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        Some("905"),
        None,
    );

    // binaryDataArrayList
    let bal = s1
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    for (i, name, accession, expect_len, expect_enc) in [
        (
            0usize,
            "m/z array",
            "MS:1000514",
            Some(20usize),
            Some(216usize),
        ),
        (
            1usize,
            "intensity array",
            "MS:1000515",
            Some(20usize),
            Some(216usize),
        ),
    ] {
        let ba = &bal.binary_data_arrays[i];
        assert_eq!(ba.cv_params.len(), 3);
        assert_cv(
            CV_REF_MODE,
            &ba.cv_params,
            "64-bit float",
            "MS:1000523",
            "MS",
            Some(""),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &ba.cv_params,
            "no compression",
            "MS:1000576",
            "MS",
            Some(""),
            None,
        );
        assert_cv(
            CV_REF_MODE,
            &ba.cv_params,
            name,
            accession,
            "MS",
            Some(""),
            None,
        );
        assert_eq!(ba.array_length, expect_len);
        assert_eq!(ba.encoded_length, expect_enc);
    }
}
