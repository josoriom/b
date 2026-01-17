use std::sync::OnceLock;

use crate::{
    mzml::structs::MzML,
    utilities::test::{
        CvRefMode, assert_cv, assert_cv_ref, assert_software_param, parse_b, spectrum_by_id,
        spectrum_description, spectrum_precursor_list, spectrum_scan_list,
    },
};

static MZML_CACHE: OnceLock<MzML> = OnceLock::new();

const PATH: &str = "data/b64/tiny4_LTQ-FT.mzML0.99.0.b64";
const CV_REF_MODE: CvRefMode = CvRefMode::AllowMissingMs;

#[test]
fn tiny1_mzml0_99_0_header_sections() {
    let mzml = parse_b(&MZML_CACHE, PATH);

    // cvList
    let cv_list = mzml.cv_list.as_ref().expect("cvList parsed");
    assert_eq!(cv_list.cv.len(), 1);
    let cv0 = &cv_list.cv[0];
    assert!(cv0.id.is_empty() || cv0.id == "MS");
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

    let g0 = &rpgl.referenceable_param_groups[0];
    assert_eq!(g0.id, "CommonMS1SpectrumParams");
    assert_cv(
        CV_REF_MODE,
        &g0.cv_params,
        "positive scan",
        "MS:1000130",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &g0.cv_params,
        "full scan",
        "MS:1000498",
        "MS",
        Some(""),
        None,
    );

    let g1 = &rpgl.referenceable_param_groups[1];
    assert_eq!(g1.id, "CommonMS2SpectrumParams");
    assert_cv(
        CV_REF_MODE,
        &g1.cv_params,
        "positive scan",
        "MS:1000130",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &g1.cv_params,
        "full scan",
        "MS:1000498",
        "MS",
        Some(""),
        None,
    );

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
    assert_eq!(inst_list.instrument.len(), 2);

    let ltq = &inst_list.instrument[0];
    assert_eq!(ltq.id, "LTQ");
    assert_cv(
        CV_REF_MODE,
        &ltq.cv_param,
        "LTQ",
        "MS:1000447",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &ltq.cv_param,
        "instrument serial number",
        "MS:1000529",
        "MS",
        Some("34454"),
        None,
    );
    let cl = ltq.component_list.as_ref().expect("componentList parsed");
    assert_eq!(cl.source.len(), 1);
    assert_eq!(cl.analyzer.len(), 1);
    assert_eq!(cl.detector.len(), 1);
    assert_eq!(cl.source[0].order, Some(1));
    assert_cv(
        CV_REF_MODE,
        &cl.source[0].cv_param,
        "nanoelectrospray",
        "MS:1000398",
        "MS",
        Some(""),
        None,
    );
    assert_eq!(cl.analyzer[0].order, Some(2));
    assert_cv(
        CV_REF_MODE,
        &cl.analyzer[0].cv_param,
        "linear ion trap",
        "MS:1000291",
        "MS",
        Some(""),
        None,
    );
    assert_eq!(cl.detector[0].order, Some(3));
    assert_cv(
        CV_REF_MODE,
        &cl.detector[0].cv_param,
        "electron multiplier",
        "MS:1000253",
        "MS",
        Some(""),
        None,
    );

    let ltqft = &inst_list.instrument[1];
    assert_eq!(ltqft.id, "LTQ FT");
    assert_cv(
        CV_REF_MODE,
        &ltqft.cv_param,
        "LTQ FT",
        "MS:1000448",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &ltqft.cv_param,
        "instrument serial number",
        "MS:1000529",
        "MS",
        Some("34454"),
        None,
    );
    let cl2 = ltqft.component_list.as_ref().expect("componentList parsed");
    assert_eq!(cl2.source.len(), 1);
    assert_eq!(cl2.analyzer.len(), 2);
    assert_eq!(cl2.detector.len(), 1);

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

    // run
    let run = &mzml.run;
    assert_eq!(run.id.as_str(), "Exp01");
    assert_eq!(run.sample_ref.as_deref(), Some("1"));
    // TODO: Fix default_instrument_configuration_ref parsing
    // assert_eq!(
    //     run.default_instrument_configuration_ref.as_deref(),
    //     Some("LTQ")
    // );

    // spectrumList
    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);

    let s = spectrum_by_id(mzml, "S19");
    assert!(s.cv_params.iter().any(|cv| cv.name == "MSn spectrum"));

    if let Some(mslvl) = s.cv_params.iter().find(|cv| cv.name == "ms level") {
        assert_eq!(mslvl.accession.as_deref(), Some("MS:1000511"));
        assert_cv_ref(CV_REF_MODE, mslvl.cv_ref.as_deref(), "MS", "ms level");
        assert_eq!(mslvl.value.as_deref(), Some("1"));
    }

    // spectrumDescription
    let sd = spectrum_description(s);
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
        1.66755e+007,
        None,
    );

    let scl = spectrum_scan_list(s);
    assert_eq!(scl.scans.len(), 1);
    let scan0 = &scl.scans[0];
    // TODO: Fix instrument_configuration_ref parsing
    // assert_eq!(
    //     scan0.instrument_configuration_ref.as_deref(),
    //     Some("LCQ Deca")
    // );
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "scan start time",
        "MS:1000016",
        "MS",
        5.890500,
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
    let bal = s
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

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

    let s = spectrum_by_id(mzml, "S20");
    assert!(s.cv_params.iter().any(|cv| cv.name == "MSn spectrum"));

    if let Some(mslvl) = s.cv_params.iter().find(|cv| cv.name == "ms level") {
        assert_eq!(mslvl.accession.as_deref(), Some("MS:1000511"));
        assert_cv_ref(CV_REF_MODE, mslvl.cv_ref.as_deref(), "MS", "ms level");
        assert_eq!(mslvl.value.as_deref(), Some("2"));
    }

    let sd = spectrum_description(s);
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
        1.66755e+007,
        None,
    );

    // precursorList
    let pl = spectrum_precursor_list(s).expect("precursorList parsed");
    assert_eq!(pl.precursors.len(), 1);
    let p0 = &pl.precursors[0];
    assert_eq!(p0.spectrum_ref.as_deref(), Some("19"));

    let sil = p0
        .selected_ion_list
        .as_ref()
        .expect("selectedIonList parsed");
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
        35.0,
        Some("electron volt"),
    );

    // scan
    let scl = spectrum_scan_list(s);
    assert_eq!(scl.scans.len(), 1);
    let scan0 = &scl.scans[0];
    // TODO: Fix instrument_configuration_ref parsing
    // assert_eq!(scan0.instrument_configuration_ref.as_deref(), Some("LTQ"));
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "scan start time",
        "MS:1000016",
        "MS",
        5.990500,
        Some("minute"),
    );
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "filter string",
        "MS:1000512",
        "MS",
        Some("+ c d Full ms2  445.35@cid35.00 [ 110.00-905.00]"),
        None,
    );

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
        110.0,
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &win0.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        905.0,
        None,
    );

    // binaryDataArrayList
    let bal = s
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

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
