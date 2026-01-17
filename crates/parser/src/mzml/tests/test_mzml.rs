use std::sync::OnceLock;

use crate::{
    mzml::structs::MzML,
    utilities::test::{CvRefMode, assert_cv, mzml, spectrum_precursor_list, spectrum_scan_list},
};

static MZML_CACHE: OnceLock<MzML> = OnceLock::new();

const PATH: &str = "data/mzml/test.mzML";
const CV_REF_MODE: CvRefMode = CvRefMode::Strict;

#[test]
fn anpc_mzml1_1_0_header_sections() {
    let mzml = mzml(&MZML_CACHE, PATH);

    let cv_list = mzml.cv_list.as_ref().expect("cvList parsed");
    assert_eq!(cv_list.cv.len(), 2);

    let ms = cv_list.cv.iter().find(|c| c.id == "MS").expect("MS cv");
    assert_eq!(
        ms.full_name.as_deref(),
        Some("Proteomics Standards Initiative Mass Spectrometry Ontology")
    );
    assert_eq!(ms.version.as_deref(), Some("4.1.182"));
    assert_eq!(
        ms.uri.as_deref(),
        Some("https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo")
    );

    let uo = cv_list.cv.iter().find(|c| c.id == "UO").expect("UO cv");
    assert_eq!(uo.full_name.as_deref(), Some("Unit Ontology"));
    assert_eq!(uo.version.as_deref(), Some("09:04:2014"));
    assert_eq!(
        uo.uri.as_deref(),
        Some(
            "https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"
        )
    );

    let file_desc = &mzml.file_description;

    // fileContent
    assert_eq!(file_desc.file_content.cv_params.len(), 2);
    assert_cv(
        CV_REF_MODE,
        &file_desc.file_content.cv_params,
        "MS1 spectrum",
        "MS:1000579",
        "MS",
        Some(""),
        None,
    );
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
    assert_eq!(sf0.id, "anpc_file.d_x005c_Analysis.baf");
    assert_eq!(sf0.name, "Analysis.baf");
    assert_eq!(sf0.location, r"file://Z:\inputDirectory\anpc_file.d");

    assert_eq!(sf0.cv_param.len(), 3);
    assert_cv(
        CV_REF_MODE,
        &sf0.cv_param,
        "Bruker BAF nativeID format",
        "MS:1000772",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &sf0.cv_param,
        "Bruker BAF format",
        "MS:1000815",
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
        Some("36a9346b9d32b3ef5b30e48d1a20cf1515232083"),
        None,
    );

    let rpgl = mzml
        .referenceable_param_group_list
        .as_ref()
        .expect("referenceableParamGroupList parsed");
    assert_eq!(rpgl.referenceable_param_groups.len(), 1);

    let g0 = &rpgl.referenceable_param_groups[0];
    assert_eq!(g0.id, "CommonInstrumentParams");
    assert_eq!(g0.cv_params.len(), 1);
    assert_cv(
        CV_REF_MODE,
        &g0.cv_params,
        "Bruker Daltonics maXis series",
        "MS:1001547",
        "MS",
        Some(""),
        None,
    );

    let sw_list = mzml.software_list.as_ref().expect("softwareList parsed");
    assert_eq!(sw_list.software.len(), 3);

    let baf2sql = sw_list.software.iter().find(|s| s.id == "BAF2SQL").unwrap();
    assert_eq!(baf2sql.version.as_deref(), Some("2.7.300.20-112"));
    assert_cv(
        CV_REF_MODE,
        &baf2sql.cv_param,
        "Bruker software",
        "MS:1000692",
        "MS",
        Some(""),
        None,
    );

    let micro = sw_list
        .software
        .iter()
        .find(|s| s.id == "micrOTOFcontrol")
        .unwrap();
    assert_eq!(micro.version.as_deref(), Some("5.2.109.779-16393"));
    assert_cv(
        CV_REF_MODE,
        &micro.cv_param,
        "micrOTOFcontrol",
        "MS:1000726",
        "MS",
        Some(""),
        None,
    );

    let pwiz = sw_list
        .software
        .iter()
        .find(|s| s.id == "pwiz_Reader_Bruker")
        .unwrap();
    assert_eq!(pwiz.version.as_deref(), Some("3.0.25114"));
    assert_cv(
        CV_REF_MODE,
        &pwiz.cv_param,
        "ProteoWizard software",
        "MS:1000615",
        "MS",
        Some(""),
        None,
    );

    let inst_list = mzml
        .instrument_list
        .as_ref()
        .expect("instrumentConfigurationList parsed");
    assert_eq!(inst_list.instrument.len(), 1);

    let ic1 = &inst_list.instrument[0];
    assert_eq!(ic1.id, "IC1");
    assert_cv(
        CV_REF_MODE,
        &ic1.cv_param,
        "instrument serial number",
        "MS:1000529",
        "MS",
        Some("1825265.10271"),
        None,
    );

    let cl = ic1.component_list.as_ref().expect("componentList parsed");
    assert_eq!(cl.source.len(), 1);
    assert_eq!(cl.analyzer.len(), 2);
    assert_eq!(cl.detector.len(), 2);

    let src = &cl.source[0];
    assert_eq!(src.order, Some(1));
    assert_cv(
        CV_REF_MODE,
        &src.cv_param,
        "electrospray ionization",
        "MS:1000073",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &src.cv_param,
        "electrospray inlet",
        "MS:1000057",
        "MS",
        Some(""),
        None,
    );

    let an2 = &cl.analyzer[0];
    assert_eq!(an2.order, Some(2));
    assert_cv(
        CV_REF_MODE,
        &an2.cv_param,
        "quadrupole",
        "MS:1000081",
        "MS",
        Some(""),
        None,
    );

    let an3 = &cl.analyzer[1];
    assert_eq!(an3.order, Some(3));
    assert_cv(
        CV_REF_MODE,
        &an3.cv_param,
        "time-of-flight",
        "MS:1000084",
        "MS",
        Some(""),
        None,
    );

    let det4 = &cl.detector[0];
    assert_eq!(det4.order, Some(4));
    assert_cv(
        CV_REF_MODE,
        &det4.cv_param,
        "microchannel plate detector",
        "MS:1000114",
        "MS",
        Some(""),
        None,
    );

    let det5 = &cl.detector[1];
    assert_eq!(det5.order, Some(5));
    assert_cv(
        CV_REF_MODE,
        &det5.cv_param,
        "photomultiplier",
        "MS:1000116",
        "MS",
        Some(""),
        None,
    );

    let dp_list = mzml
        .data_processing_list
        .as_ref()
        .expect("dataProcessingList parsed");
    assert_eq!(dp_list.data_processing.len(), 1);

    let dp0 = &dp_list.data_processing[0];
    assert_eq!(dp0.id, "pwiz_Reader_Bruker_conversion");
    assert_eq!(dp0.processing_method.len(), 1);

    let pm0 = &dp0.processing_method[0];
    assert_cv(
        CV_REF_MODE,
        &pm0.cv_param,
        "Conversion to mzML",
        "MS:1000544",
        "MS",
        Some(""),
        None,
    );

    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);

    let cl = run
        .chromatogram_list
        .as_ref()
        .expect("chromatogramList parsed");
    assert_eq!(cl.chromatograms.len(), 2);
}

#[test]
fn anpc_mzml1_1_0_first_spectrum() {
    let mzml = mzml(&MZML_CACHE, PATH);
    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);

    let s0 = &sl.spectra[0];
    assert_eq!(s0.index, Some(0));
    assert_eq!(s0.id, "scan=1");

    assert_eq!(s0.cv_params.len(), 6);
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "ms level",
        "MS:1000511",
        "MS",
        Some("1"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "MS1 spectrum",
        "MS:1000579",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "positive scan",
        "MS:1000130",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        Some("24998.0"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        Some("4.40132e05"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s0.cv_params,
        "profile spectrum",
        "MS:1000128",
        "MS",
        Some(""),
        None,
    );

    // scanList/scan
    let scl = spectrum_scan_list(s0);
    assert_eq!(scl.scans.len(), 1);
    // TODO: re-enable once ScanList.cv_params is parsed
    // assert_cv(
    //     CV_REF_MODE,
    //     &scl.cv_params,
    //     "no combination",
    //     "MS:1000795",
    //     "MS",
    //     Some(""),
    //     None,
    // );

    let scan0 = &scl.scans[0];
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "scan start time",
        "MS:1000016",
        "MS",
        Some("0.191"),
        Some("second"),
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
        Some("30.0"),
        Some("m/z"),
    );
    assert_cv(
        CV_REF_MODE,
        &win0.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        Some("1000.0"),
        Some("m/z"),
    );

    // precursorList must be absent for MS1
    assert!(spectrum_precursor_list(s0).is_none());

    // binaryDataArrayList
    let bal = s0
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    // m/z array
    let mz = &bal.binary_data_arrays[0];
    assert_eq!(mz.encoded_length, Some(3627008));
    assert_cv(
        CV_REF_MODE,
        &mz.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz.cv_params,
        "m/z array",
        "MS:1000514",
        "MS",
        Some(""),
        Some("m/z"),
    );

    // intensity array (32-bit float)
    let it = &bal.binary_data_arrays[1];
    assert_eq!(it.encoded_length, Some(1813504));
    assert_cv(
        CV_REF_MODE,
        &it.cv_params,
        "32-bit float",
        "MS:1000521",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &it.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &it.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        Some("number of detector counts"),
    );
}

#[test]
fn anpc_mzml1_1_0_last_spectrum() {
    let mzml = mzml(&MZML_CACHE, PATH);
    let run = &mzml.run;

    let sl = run.spectrum_list.as_ref().expect("spectrumList parsed");
    assert_eq!(sl.spectra.len(), 2);

    let s_last = sl.spectra.last().expect("last spectrum");
    assert_eq!(s_last.index, Some(3475));
    assert_eq!(s_last.id, "scan=3476");

    assert_eq!(s_last.cv_params.len(), 6);
    assert_cv(
        CV_REF_MODE,
        &s_last.cv_params,
        "ms level",
        "MS:1000511",
        "MS",
        Some("2"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s_last.cv_params,
        "MSn spectrum",
        "MS:1000580",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s_last.cv_params,
        "positive scan",
        "MS:1000130",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s_last.cv_params,
        "base peak intensity",
        "MS:1000505",
        "MS",
        Some("20032.0"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s_last.cv_params,
        "total ion current",
        "MS:1000285",
        "MS",
        Some("3.59026e05"),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &s_last.cv_params,
        "centroid spectrum",
        "MS:1000127",
        "MS",
        Some(""),
        None,
    );

    // scanList/scan
    // TODO: re-enable once ScanList.cv_params is parsed
    let scl = spectrum_scan_list(s_last);
    assert_eq!(scl.scans.len(), 1);
    // assert_cv(
    //     CV_REF_MODE,
    //     &scl.cv_params,
    //     "no combination",
    //     "MS:1000795",
    //     "MS",
    //     Some(""),
    //     None,
    // );

    let scan0 = &scl.scans[0];
    assert_cv(
        CV_REF_MODE,
        &scan0.cv_params,
        "scan start time",
        "MS:1000016",
        "MS",
        Some("452.262"),
        Some("second"),
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
        Some("30.0"),
        Some("m/z"),
    );
    assert_cv(
        CV_REF_MODE,
        &win0.cv_params,
        "scan window upper limit",
        "MS:1000500",
        "MS",
        Some("1000.0"),
        Some("m/z"),
    );

    // precursorList
    let pl = spectrum_precursor_list(s_last).expect("precursorList parsed");
    assert_eq!(pl.precursors.len(), 1);

    let p0 = &pl.precursors[0];

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
        Some("515.0"),
        Some("m/z"),
    );
    assert_cv(
        CV_REF_MODE,
        &iw.cv_params,
        "isolation window lower offset",
        "MS:1000828",
        "MS",
        Some("485.0"),
        Some("m/z"),
    );
    assert_cv(
        CV_REF_MODE,
        &iw.cv_params,
        "isolation window upper offset",
        "MS:1000829",
        "MS",
        Some("485.0"),
        Some("m/z"),
    );

    let sil = p0
        .selected_ion_list
        .as_ref()
        .expect("selectedIonList parsed");
    assert_eq!(sil.selected_ions.len(), 1);

    let ion0 = &sil.selected_ions[0];
    assert_eq!(ion0.cv_params.len(), 1);
    assert_cv(
        CV_REF_MODE,
        &ion0.cv_params,
        "selected ion m/z",
        "MS:1000744",
        "MS",
        Some("515.0"),
        Some("m/z"),
    );

    let act = p0.activation.as_ref().expect("activation parsed");
    assert_cv(
        CV_REF_MODE,
        &act.cv_params,
        "in-source collision-induced dissociation",
        "MS:1001880",
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
        Some("20.0"),
        None,
    );

    // binaryDataArrayList
    let bal = s_last
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 2);

    let mz = &bal.binary_data_arrays[0];
    assert_eq!(mz.encoded_length, Some(46296));
    assert_cv(
        CV_REF_MODE,
        &mz.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &mz.cv_params,
        "m/z array",
        "MS:1000514",
        "MS",
        Some(""),
        Some("m/z"),
    );

    let it = &bal.binary_data_arrays[1];
    assert_eq!(it.encoded_length, Some(23148));
    assert_cv(
        CV_REF_MODE,
        &it.cv_params,
        "32-bit float",
        "MS:1000521",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &it.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &it.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        Some("number of detector counts"),
    );
}

#[test]
fn anpc_mzml1_1_0_chromatograms() {
    let mzml = mzml(&MZML_CACHE, PATH);
    let run = &mzml.run;

    let cl = run
        .chromatogram_list
        .as_ref()
        .expect("chromatogramList parsed");
    assert_eq!(cl.chromatograms.len(), 2);

    // TIC
    let tic = &cl.chromatograms[0];
    assert_eq!(tic.index, Some(0));
    assert_eq!(tic.id, "TIC");

    assert_eq!(tic.cv_params.len(), 1);
    assert_cv(
        CV_REF_MODE,
        &tic.cv_params,
        "total ion current chromatogram",
        "MS:1000235",
        "MS",
        Some(""),
        None,
    );

    let bal = tic
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 3);

    // time array (64-bit float)
    let t = &bal.binary_data_arrays[0];
    assert_eq!(t.encoded_length, Some(37080));
    assert_cv(
        CV_REF_MODE,
        &t.cv_params,
        "64-bit float",
        "MS:1000523",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &t.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &t.cv_params,
        "time array",
        "MS:1000595",
        "MS",
        Some(""),
        Some("second"),
    );

    // intensity array (32-bit float)
    let i = &bal.binary_data_arrays[1];
    assert_eq!(i.encoded_length, Some(18540));
    assert_cv(
        CV_REF_MODE,
        &i.cv_params,
        "32-bit float",
        "MS:1000521",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &i.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &i.cv_params,
        "intensity array",
        "MS:1000515",
        "MS",
        Some(""),
        Some("number of detector counts"),
    );

    let ms_level = &bal.binary_data_arrays[2];
    assert_eq!(ms_level.array_length, Some(3476));
    assert_eq!(ms_level.encoded_length, Some(37080));
    assert_cv(
        CV_REF_MODE,
        &ms_level.cv_params,
        "64-bit integer",
        "MS:1000522",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &ms_level.cv_params,
        "no compression",
        "MS:1000576",
        "MS",
        Some(""),
        None,
    );
    assert_cv(
        CV_REF_MODE,
        &ms_level.cv_params,
        "non-standard data array",
        "MS:1000786",
        "MS",
        Some("ms level"),
        Some("dimensionless unit"),
    );

    // BPC
    let bpc = &cl.chromatograms[1];
    assert_eq!(bpc.index, Some(1));
    assert_eq!(bpc.id, "BPC");

    assert_eq!(bpc.cv_params.len(), 1);
    assert_cv(
        CV_REF_MODE,
        &bpc.cv_params,
        "basepeak chromatogram",
        "MS:1000628",
        "MS",
        Some(""),
        None,
    );

    let bal = bpc
        .binary_data_array_list
        .as_ref()
        .expect("binaryDataArrayList parsed");
    assert_eq!(bal.binary_data_arrays.len(), 3);

    let ms_level = &bal.binary_data_arrays[2];
    assert_eq!(ms_level.array_length, Some(3476));
    assert_eq!(ms_level.encoded_length, Some(37080));
    assert_cv(
        CV_REF_MODE,
        &ms_level.cv_params,
        "non-standard data array",
        "MS:1000786",
        "MS",
        Some("ms level"),
        Some("dimensionless unit"),
    );
}
