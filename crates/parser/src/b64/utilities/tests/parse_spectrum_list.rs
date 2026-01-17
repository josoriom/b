use std::{fs, path::PathBuf};

use crate::b64::decode::Metadatum;
use crate::b64::utilities::common::ChildIndex;
use crate::b64::utilities::{parse_header, parse_metadata, parse_spectrum_list};
use crate::mzml::schema::{TagId, schema};
use crate::{CvParam, SpectrumList};

const PATH: &str = "data/b64/test.b64";

fn read_bytes(path: &str) -> Vec<u8> {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(path);
    fs::read(&full).unwrap_or_else(|e| panic!("cannot read {:?}: {}", full, e))
}

fn parse_metadata_section_from_test_file(
    start_off: u64,
    end_off: u64,
    item_count: u32,
    expected_item_count: u32,
    meta_count: u32,
    num_count: u32,
    str_count: u32,
    compression_flag_bit: u8,
    expected_total_meta_len: usize,
    section_name: &str,
) -> Vec<Metadatum> {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let c0 = start_off as usize;
    let c1 = end_off as usize;

    assert!(
        c0 < c1,
        "invalid metadata offsets for {section_name}: start >= end"
    );
    assert!(
        c1 <= bytes.len(),
        "invalid metadata offsets for {section_name}: end out of bounds"
    );

    assert_eq!(
        item_count, expected_item_count,
        "test.b64 should contain {expected_item_count} {section_name} items"
    );

    let compressed = (header.reserved_flags & (1u8 << compression_flag_bit)) != 0;
    let slice = &bytes[c0..c1];

    let meta = parse_metadata(
        slice,
        item_count,
        meta_count,
        num_count,
        str_count,
        compressed,
        header.reserved_flags,
    )
    .expect("parse_metadata failed");

    assert_eq!(
        meta.len(),
        expected_total_meta_len,
        "unexpected {section_name} metadata count (expected {expected_total_meta_len} total items)"
    );

    meta
}

fn parse_spectrum_list_from_test_file() -> SpectrumList {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    let meta = parse_metadata_section_from_test_file(
        header.off_spec_meta,
        header.off_chrom_meta,
        header.spectrum_count,
        2,
        header.spec_meta_count,
        header.spec_num_count,
        header.spec_str_count,
        4,
        header.spec_meta_count as usize,
        "spectra",
    );

    let child_index = ChildIndex::new(&meta);

    let spectrum_list = parse_spectrum_list(schema(), &meta, &child_index)
        .expect("parse_spectrum_list returned None");
    spectrum_list
}

#[derive(Clone, Copy, Debug)]
enum ExpectedValue<'a> {
    None,
    Str(&'a str),
    Num(f64),
}

#[derive(Clone, Copy, Debug)]
struct ExpectedCv<'a> {
    cv_ref: Option<&'a str>,
    accession: Option<&'a str>,
    name: &'a str,
    value: ExpectedValue<'a>,
    unit_cv_ref: Option<&'a str>,
    unit_name: Option<&'a str>,
    unit_accession: Option<&'a str>,
}

fn assert_f64_close(got: f64, expected: f64) {
    let diff = (got - expected).abs();
    let scale = expected.abs().max(1.0);
    let tol = 1e-9 * scale;
    assert!(
        diff <= tol,
        "numeric mismatch: got={got} expected={expected} diff={diff} tol={tol}"
    );
}

fn assert_cv_param_strict(p: &CvParam, e: ExpectedCv<'_>) {
    assert_eq!(p.cv_ref.as_deref(), e.cv_ref, "cv_ref mismatch");
    assert_eq!(p.accession.as_deref(), e.accession, "accession mismatch");
    assert_eq!(p.name.as_str(), e.name, "name mismatch");
    assert_eq!(
        p.unit_cv_ref.as_deref(),
        e.unit_cv_ref,
        "unit_cv_ref mismatch"
    );
    assert_eq!(p.unit_name.as_deref(), e.unit_name, "unit_name mismatch");
    assert_eq!(
        p.unit_accession.as_deref(),
        e.unit_accession,
        "unit_accession mismatch"
    );

    match e.value {
        ExpectedValue::None => {
            let v = p.value.as_deref();
            assert!(
                v.is_none() || v == Some(""),
                "value mismatch: expected None/empty, got {:?}",
                p.value
            );
        }
        ExpectedValue::Str(s) => {
            assert_eq!(
                p.value.as_deref(),
                Some(s),
                "value mismatch: expected {:?}, got {:?}",
                s,
                p.value
            );
        }
        ExpectedValue::Num(x) => {
            let raw = p
                .value
                .as_deref()
                .unwrap_or_else(|| panic!("value mismatch: expected numeric {x}, got None"));

            let got = raw
                .parse::<f64>()
                .unwrap_or_else(|_| panic!("expected numeric value, got {:?}", raw));

            assert_f64_close(got, x);
        }
    }
}

fn fmt_expected_value(ev: ExpectedValue<'_>) -> String {
    match ev {
        ExpectedValue::None => "None".to_string(),
        ExpectedValue::Str(s) => format!("Str({s:?})"),
        ExpectedValue::Num(x) => format!("Num({x:?})"),
    }
}

fn assert_cv_params_exact(actual: &[CvParam], expected: &[ExpectedCv<'_>]) {
    let fmt_actual = |p: &CvParam| {
        format!(
            "cvRef={:?} acc={:?} name={:?} value={:?} unit=({:?},{:?},{:?})",
            p.cv_ref.as_deref(),
            p.accession.as_deref(),
            p.name,
            p.value.as_deref(),
            p.unit_cv_ref.as_deref(),
            p.unit_name.as_deref(),
            p.unit_accession.as_deref(),
        )
    };

    let fmt_expected = |e: &ExpectedCv<'_>| {
        format!(
            "cvRef={:?} acc={:?} name={:?} value={} unit=({:?},{:?},{:?})",
            e.cv_ref,
            e.accession,
            e.name,
            fmt_expected_value(e.value),
            e.unit_cv_ref,
            e.unit_name,
            e.unit_accession
        )
    };

    if actual.len() != expected.len() {
        use std::collections::BTreeMap;

        let mut a_map: BTreeMap<String, usize> = BTreeMap::new();
        let mut e_map: BTreeMap<String, usize> = BTreeMap::new();

        for p in actual {
            let k = format!("{:?}|{}", p.accession.as_deref(), p.name);
            *a_map.entry(k).or_insert(0) += 1;
        }
        for e in expected {
            let k = format!("{:?}|{}", e.accession, e.name);
            *e_map.entry(k).or_insert(0) += 1;
        }

        let mut missing: Vec<String> = Vec::new();
        let mut extra: Vec<String> = Vec::new();

        for (k, c) in &e_map {
            let ac = a_map.get(k).copied().unwrap_or(0);
            if ac < *c {
                missing.push(format!("{k} x{}", c - ac));
            }
        }
        for (k, c) in &a_map {
            let ec = e_map.get(k).copied().unwrap_or(0);
            if ec < *c {
                extra.push(format!("{k} x{}", c - ec));
            }
        }

        panic!(
            "cv_params length mismatch: actual={} expected={}\nmissing={:#?}\nextra={:#?}\n\
             --- expected ---\n{:#?}\n--- actual ---\n{:#?}\n",
            actual.len(),
            expected.len(),
            missing,
            extra,
            expected.iter().map(fmt_expected).collect::<Vec<_>>(),
            actual.iter().map(fmt_actual).collect::<Vec<_>>(),
        );
    }

    for (i, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        if a.cv_ref.as_deref() != e.cv_ref
            || a.accession.as_deref() != e.accession
            || a.name.as_str() != e.name
            || a.unit_cv_ref.as_deref() != e.unit_cv_ref
            || a.unit_name.as_deref() != e.unit_name
            || a.unit_accession.as_deref() != e.unit_accession
        {
            panic!(
                "cv_param mismatch at index {i}\nexpected: {}\nactual:   {}\nfull actual list:\n{:#?}",
                fmt_expected(e),
                fmt_actual(a),
                actual.iter().map(fmt_actual).collect::<Vec<_>>()
            );
        }

        assert_cv_param_strict(a, *e);
    }
}

#[test]
fn spectrum_list_strict_attributes_and_count() {
    let sl = parse_spectrum_list_from_test_file();

    assert_eq!(sl.count, Some(3476));
    assert_eq!(
        sl.default_data_processing_ref.as_deref(),
        Some("pwiz_Reader_Bruker_conversion")
    );

    assert_eq!(sl.spectra.len(), 2);
}

#[test]
fn spectrum0_strict_full_structure() {
    let sl = parse_spectrum_list_from_test_file();
    let s = &sl.spectra[0];

    assert_eq!(s.id.as_str(), "scan=1");
    assert_eq!(s.index, Some(0));
    assert_eq!(s.default_array_length, Some(340032));

    assert_eq!(s.native_id.as_deref(), None);
    assert_eq!(s.source_file_ref.as_deref(), None);
    assert_eq!(s.spot_id.as_deref(), None);
    assert_eq!(s.scan_number, None);

    assert_eq!(
        s.data_processing_ref.as_deref(),
        Some("pwiz_Reader_Bruker_conversion")
    );

    let expected_spec = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000511"),
            name: "ms level",
            value: ExpectedValue::Num(1.0),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000579"),
            name: "MS1 spectrum",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000130"),
            name: "positive scan",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000505"),
            name: "base peak intensity",
            value: ExpectedValue::Num(24998.0),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000285"),
            name: "total ion current",
            value: ExpectedValue::Num(4.40132e05),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000128"),
            name: "profile spectrum",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
    ];
    assert_cv_params_exact(&s.cv_params, expected_spec);

    let scan_list = s.scan_list.as_ref().expect("missing <scanList>");
    assert_eq!(scan_list.count, Some(1));
    assert_eq!(scan_list.scans.len(), 1);

    // IMPORTANT: raw XML includes <cvParam ... MS:1000795 no combination/> under <scanList>
    // let expected_scan_list_cvs = &[ExpectedCv {
    //     cv_ref: Some("MS"),
    //     accession: Some("MS:1000795"),
    //     name: "no combination",
    //     value: ExpectedValue::None,
    //     unit_cv_ref: None,
    //     unit_name: None,
    //     unit_accession: None,
    // }];
    // assert_cv_params_exact(&scan_list.cv_params, expected_scan_list_cvs);

    let scan = &scan_list.scans[0];

    assert_eq!(scan.instrument_configuration_ref.as_deref(), None);
    assert_eq!(scan.external_spectrum_id.as_deref(), None);
    assert_eq!(scan.source_file_ref.as_deref(), None);
    assert_eq!(scan.spectrum_ref.as_deref(), None);

    let expected_scan_cvs = &[ExpectedCv {
        cv_ref: Some("MS"),
        accession: Some("MS:1000016"),
        name: "scan start time",
        value: ExpectedValue::Num(0.191),
        unit_cv_ref: Some("UO"),
        unit_name: Some("second"),
        unit_accession: Some("UO:0000010"),
    }];
    assert_cv_params_exact(&scan.cv_params, expected_scan_cvs);

    // scanWindowList
    let swl = scan
        .scan_window_list
        .as_ref()
        .expect("missing <scanWindowList>");
    assert_eq!(swl.count, Some(1));
    assert_eq!(swl.scan_windows.len(), 1);

    let sw = &swl.scan_windows[0];
    let expected_sw_cvs = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000501"),
            name: "scan window lower limit",
            value: ExpectedValue::Num(30.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000500"),
            name: "scan window upper limit",
            value: ExpectedValue::Num(1000.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
    ];
    assert_cv_params_exact(&sw.cv_params, expected_sw_cvs);

    let bal = s
        .binary_data_array_list
        .as_ref()
        .expect("missing <binaryDataArrayList>");
    assert_eq!(bal.count, Some(2));
    assert_eq!(bal.binary_data_arrays.len(), 2);

    let ba0 = &bal.binary_data_arrays[0];
    let expected_ba0 = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000523"),
            name: "64-bit float",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000576"),
            name: "no compression",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000514"),
            name: "m/z array",
            value: ExpectedValue::None,
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
    ];
    assert_cv_params_exact(&ba0.cv_params, expected_ba0);

    let ba1 = &bal.binary_data_arrays[1];
    let expected_ba1 = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000521"),
            name: "32-bit float",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000576"),
            name: "no compression",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000515"),
            name: "intensity array",
            value: ExpectedValue::None,
            unit_cv_ref: Some("MS"),
            unit_name: Some("number of detector counts"),
            unit_accession: Some("MS:1000131"),
        },
    ];
    assert_cv_params_exact(&ba1.cv_params, expected_ba1);

    assert!(s.precursor_list.is_none());
    assert!(s.product_list.is_none());
}

#[test]
fn spectrum1_strict_full_structure_with_precursor_list() {
    let sl = parse_spectrum_list_from_test_file();
    let s = &sl.spectra[1];

    assert_eq!(s.id.as_str(), "scan=3476");
    assert_eq!(s.index, Some(3475));
    assert_eq!(s.default_array_length, Some(4340));

    assert_eq!(s.native_id.as_deref(), None);
    assert_eq!(s.source_file_ref.as_deref(), None);
    assert_eq!(s.spot_id.as_deref(), None);
    assert_eq!(s.scan_number, None);

    assert_eq!(
        s.data_processing_ref.as_deref(),
        Some("pwiz_Reader_Bruker_conversion")
    );

    let expected_spec = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000511"),
            name: "ms level",
            value: ExpectedValue::Num(2.0),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000580"),
            name: "MSn spectrum",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000130"),
            name: "positive scan",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000505"),
            name: "base peak intensity",
            value: ExpectedValue::Num(20032.0),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000285"),
            name: "total ion current",
            value: ExpectedValue::Num(3.59026e05),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000127"),
            name: "centroid spectrum",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
    ];
    assert_cv_params_exact(&s.cv_params, expected_spec);

    let scan_list = s.scan_list.as_ref().expect("missing <scanList>");
    assert_eq!(scan_list.count, Some(1));
    assert_eq!(scan_list.scans.len(), 1);

    // let expected_scan_list_cvs = &[ExpectedCv {
    //     cv_ref: Some("MS"),
    //     accession: Some("MS:1000795"),
    //     name: "no combination",
    //     value: ExpectedValue::None,
    //     unit_cv_ref: None,
    //     unit_name: None,
    //     unit_accession: None,
    // }];
    // assert_cv_params_exact(&scan_list.cv_params, expected_scan_list_cvs);

    let scan = &scan_list.scans[0];
    let expected_scan_cvs = &[ExpectedCv {
        cv_ref: Some("MS"),
        accession: Some("MS:1000016"),
        name: "scan start time",
        value: ExpectedValue::Num(452.262),
        unit_cv_ref: Some("UO"),
        unit_name: Some("second"),
        unit_accession: Some("UO:0000010"),
    }];
    assert_cv_params_exact(&scan.cv_params, expected_scan_cvs);

    let swl = scan
        .scan_window_list
        .as_ref()
        .expect("missing <scanWindowList>");
    assert_eq!(swl.count, Some(1));
    assert_eq!(swl.scan_windows.len(), 1);

    let sw = &swl.scan_windows[0];
    let expected_sw_cvs = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000501"),
            name: "scan window lower limit",
            value: ExpectedValue::Num(30.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000500"),
            name: "scan window upper limit",
            value: ExpectedValue::Num(1000.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
    ];
    assert_cv_params_exact(&sw.cv_params, expected_sw_cvs);

    let precursor_list = s.precursor_list.as_ref().expect("missing <precursorList>");
    assert_eq!(precursor_list.count, Some(1));
    assert_eq!(precursor_list.precursors.len(), 1);
    assert!(
        precursor_list.cv_params.is_empty(),
        "precursorList should have no cv_params in this XML"
    );

    let precursor = &precursor_list.precursors[0];

    let iw = precursor
        .isolation_window
        .as_ref()
        .expect("missing <isolationWindow>");
    let expected_iw = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000827"),
            name: "isolation window target m/z",
            value: ExpectedValue::Num(515.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000828"),
            name: "isolation window lower offset",
            value: ExpectedValue::Num(485.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000829"),
            name: "isolation window upper offset",
            value: ExpectedValue::Num(485.0),
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
    ];
    assert_cv_params_exact(&iw.cv_params, expected_iw);

    let sil = precursor
        .selected_ion_list
        .as_ref()
        .expect("missing <selectedIonList>");
    assert_eq!(sil.count, Some(1));
    assert_eq!(sil.selected_ions.len(), 1);

    let ion = &sil.selected_ions[0];
    let expected_ion = &[ExpectedCv {
        cv_ref: Some("MS"),
        accession: Some("MS:1000744"),
        name: "selected ion m/z",
        value: ExpectedValue::Num(515.0),
        unit_cv_ref: Some("MS"),
        unit_name: Some("m/z"),
        unit_accession: Some("MS:1000040"),
    }];
    assert_cv_params_exact(&ion.cv_params, expected_ion);

    let act = precursor.activation.as_ref().expect("missing <activation>");
    let expected_act = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1001880"),
            name: "in-source collision-induced dissociation",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000045"),
            name: "collision energy",
            value: ExpectedValue::Num(20.0),
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
    ];
    assert_cv_params_exact(&act.cv_params, expected_act);

    let bal = s
        .binary_data_array_list
        .as_ref()
        .expect("missing <binaryDataArrayList>");
    assert_eq!(bal.count, Some(2));
    assert_eq!(bal.binary_data_arrays.len(), 2);

    let ba0 = &bal.binary_data_arrays[0];
    let expected_ba0 = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000523"),
            name: "64-bit float",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000576"),
            name: "no compression",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000514"),
            name: "m/z array",
            value: ExpectedValue::None,
            unit_cv_ref: Some("MS"),
            unit_name: Some("m/z"),
            unit_accession: Some("MS:1000040"),
        },
    ];
    assert_cv_params_exact(&ba0.cv_params, expected_ba0);

    let ba1 = &bal.binary_data_arrays[1];
    let expected_ba1 = &[
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000521"),
            name: "32-bit float",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000576"),
            name: "no compression",
            value: ExpectedValue::None,
            unit_cv_ref: None,
            unit_name: None,
            unit_accession: None,
        },
        ExpectedCv {
            cv_ref: Some("MS"),
            accession: Some("MS:1000515"),
            name: "intensity array",
            value: ExpectedValue::None,
            unit_cv_ref: Some("MS"),
            unit_name: Some("number of detector counts"),
            unit_accession: Some("MS:1000131"),
        },
    ];
    assert_cv_params_exact(&ba1.cv_params, expected_ba1);

    assert!(s.product_list.is_none());
}
