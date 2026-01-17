use std::sync::OnceLock;

use crate::{
    mzml::structs::MzML,
    utilities::test::{
        CvRefMode, assert_cv, assert_software_param, parse_b, spectrum_by_id, spectrum_description,
        spectrum_precursor_list, spectrum_scan_list,
    },
};

static MZML_CACHE: OnceLock<MzML> = OnceLock::new();

const PATH: &str = "data/b64/tiny1.mzML0.99.0.b64";
const CV_REF_MODE: CvRefMode = CvRefMode::AllowMissingMs;

#[test]
fn tiny1_mzml0_99_0_header_sections() {
    let mzml = parse_b(&MZML_CACHE, PATH);
    // cvList
    // let cv_list = mzml.cv_list.as_ref().expect("cvList parsed");
    // assert_eq!(cv_list.cv.len(), 1);
    // let cv0 = &cv_list.cv[0];
    // assert!(cv0.id.is_empty() || cv0.id == "MS");
    // assert_eq!(
    //     cv0.full_name.as_deref(),
    //     Some("Proteomics Standards Initiative Mass Spectrometry Ontology")
    // );
    // assert_eq!(cv0.version.as_deref(), Some("2.0.2"));
    // assert_eq!(
    //     cv0.uri.as_deref(),
    //     Some("http://psidev.sourceforge.net/ms/xml/mzdata/psi-ms.2.0.2.obo")
    // );

    // fileDescription
    let file_desc = &mzml.file_description;

    // fileContent
    assert_eq!(file_desc.file_content.cv_params.len(), 1);
    assert_cv(
        CV_REF_MODE,
        &file_desc.file_content.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );

    // sourceFileList
    assert_eq!(file_desc.source_file_list.source_file.len(), 1);
    let sf0 = &file_desc.source_file_list.source_file[0];
    assert_eq!(sf0.id, "1");
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

    // contact
    let contact = file_desc.contacts.first().expect("contacts parsed");
    assert_cv(
        CV_REF_MODE,
        &contact.cv_params,
        "contact name",
        "MS:1000586",
        "MS",
        Some("William Pennington"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &contact.cv_params,
        "contact address",
        "MS:1000587",
        "MS",
        Some("Higglesworth University, 12 Higglesworth Avenue, 12045, HI, USA"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &contact.cv_params,
        "contact URL",
        "MS:1000588",
        "MS",
        Some("http://www.higglesworth.edu/"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &contact.cv_params,
        "contact email",
        "MS:1000589",
        "MS",
        Some("wpennington@higglesworth.edu"),
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
    assert_eq!(sample0.id, "1");
    assert_eq!(sample0.name, "Sample1");

    // instrumentList
    let inst_list = mzml
        .instrument_list
        .as_ref()
        .expect("instrumentList parsed");
    assert_eq!(inst_list.instrument.len(), 1);
    let inst0 = &inst_list.instrument[0];
    assert_eq!(inst0.id, "LCQ Deca");
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
    assert_eq!(dp0.id, "Xcalibur Processing");
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
    assert_eq!(dp1.id, "ReAdW Conversion");
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
}

#[test]
fn tiny1_mzml0_99_0_spectrum_s19() {
    let mzml = parse_b(&MZML_CACHE, PATH);

    // run / sourceFileRefList
    let run = &mzml.run;
    let sfrefl = run
        .source_file_ref_list
        .as_ref()
        .expect("sourceFileRefList parsed");
    assert_eq!(sfrefl.source_file_refs.len(), 1);
    assert_eq!(sfrefl.source_file_refs[0].r#ref, "1");

    // spectrumList
    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);

    // spectrum
    let s0 = spectrum_by_id(mzml, "S19");
    assert!(s0.cv_params.iter().any(|cv| cv.name == "MSn spectrum"));
    if s0.cv_params.iter().any(|cv| cv.name == "ms level") {
        assert_cv(
            CV_REF_MODE,
            &s0.cv_params,
            "ms level",
            "MS:1000511",
            "MS",
            Some("1"),
            None,
        );
    }

    // spectrumDescription
    let sd = spectrum_description(s0);
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
        400.39,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "highest observed m/z",
        "MS:1000527",
        "MS",
        1795.56,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak m/z",
        "MS:1000504",
        "MS",
        445.347,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        120053.0,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        16675500.0,
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
        5.8905,
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
        400.0,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &win0.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        1800.0,
        None,
    );

    // binaryDataArrayList
    let bal = s0
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    // binaryDataArray (m/z)
    let mz_ba = &bal.binary_data_arrays[0];
    assert_eq!(mz_ba.array_length, Some(1313));
    assert_eq!(mz_ba.encoded_length, Some(5000));
    assert_cv(
        CV_REF_MODE,
        &mz_ba.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz_ba.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz_ba.cv_params,
        "m/z array",
        "MS:1000514",
        "MS",
        Some(""),
        None,
    );

    // binaryDataArray (intensity)
    let int_ba = &bal.binary_data_arrays[1];
    assert_eq!(int_ba.array_length, Some(1313));
    assert_eq!(int_ba.encoded_length, Some(5000));
    assert_cv(
        CV_REF_MODE,
        &int_ba.cv_params,
        "32-bit float",
        "MS:1000521",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &int_ba.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &int_ba.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        None,
    );
}

#[test]
fn tiny1_mzml0_99_0_spectrum_s20() {
    let mzml = parse_b(&MZML_CACHE, PATH);

    // spectrumList
    let sl = mzml
        .run
        .spectrum_list
        .as_ref()
        .expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);

    // spectrum
    let s1 = spectrum_by_id(mzml, "S20");
    assert!(s1.cv_params.iter().any(|cv| cv.name == "MSn spectrum"));
    if s1.cv_params.iter().any(|cv| cv.name == "ms level") {
        assert_cv(
            CV_REF_MODE,
            &s1.cv_params,
            "ms level",
            "MS:1000511",
            "MS",
            Some("2"),
            None,
        );
    }

    // spectrumDescription
    let sd = spectrum_description(s1);
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
        320.39,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "highest observed m/z",
        "MS:1000527",
        "MS",
        1003.56,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak m/z",
        "MS:1000504",
        "MS",
        456.347,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        23433.0,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sd.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        16675500.0,
        None,
    );

    // precursorList
    let pl = spectrum_precursor_list(s1).expect("precursorList parsed");
    assert_eq!(pl.precursors.len(), 1);
    let p0 = &pl.precursors[0];
    assert_eq!(p0.spectrum_ref.as_deref(), Some("19"));

    if let Some(sil) = p0.selected_ion_list.as_ref() {
        assert_eq!(sil.selected_ions.len(), 1);
        let ion0 = &sil.selected_ions[0];
        assert_cv(
            CV_REF_MODE,
            &ion0.cv_params,
            "m/z",
            "MS:1000040",
            "MS",
            445.34,
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
    }

    if let Some(act) = p0.activation.as_ref() {
        assert_cv(
            CV_REF_MODE,
            &act.cv_params,
            "collision-induced dissociation",
            "MS:1000133",
            "MS",
            Some(""),
            None,
        );
        // assert_cv(
        //     CV_REF_MODE,
        //     &act.cv_params,
        //     "collision energy",
        //     "MS:1000045",
        //     "MS",
        //     35.0,
        //     Some("electron volt"),
        // );
    }

    // scan
    let scl = spectrum_scan_list(s1);
    assert_eq!(scl.scans.len(), 1);
    let scan1 = &scl.scans[0];
    // assert_cv(
    //     CV_REF_MODE,
    //     &scan1.cv_params,
    //     "scan start time",
    //     "MS:1000016",
    //     "MS",
    //     5.9905,
    //     Some("minute"),
    // );
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
        110.0,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &win1.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        905.0,
        None,
    );

    // binaryDataArrayList
    let bal = s1
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    // binaryDataArray (m/z)
    let mz_ba = &bal.binary_data_arrays[0];
    assert_eq!(mz_ba.array_length, Some(43));
    assert_eq!(mz_ba.encoded_length, Some(5000));
    assert_cv(
        CV_REF_MODE,
        &mz_ba.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz_ba.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz_ba.cv_params,
        "m/z array",
        "MS:1000514",
        "MS",
        Some(""),
        None,
    );

    // binaryDataArray (intensity)
    let int_ba = &bal.binary_data_arrays[1];
    assert_eq!(int_ba.array_length, Some(43));
    assert_eq!(int_ba.encoded_length, Some(2500));
    assert_cv(
        CV_REF_MODE,
        &int_ba.cv_params,
        "32-bit float",
        "MS:1000521",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &int_ba.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &int_ba.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        None,
    );
}
