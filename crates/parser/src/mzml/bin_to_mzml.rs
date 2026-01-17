use base64::Engine;
use base64::engine::general_purpose::STANDARD;

use quick_xml::Writer;
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};

use crate::mzml::structs::*;

#[derive(Default)]
struct IndexAcc {
    spectrum: Vec<IndexOffsetAcc>,
    chromatogram: Vec<IndexOffsetAcc>,
}

struct IndexOffsetAcc {
    id_ref: String,
    offset: u64,
}

fn nonempty<'a>(s: Option<&'a str>) -> Option<&'a str> {
    match s {
        Some(v) if !v.is_empty() => Some(v),
        _ => None,
    }
}

fn write_start_capture_offset(
    writer: &mut Writer<Vec<u8>>,
    tag: BytesStart<'_>,
) -> Result<u64, String> {
    let before = writer.get_ref().len();
    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;
    let after = writer.get_ref().len();

    let buf = writer.get_ref();
    let rel = buf[before..after]
        .iter()
        .position(|&b| b == b'<')
        .ok_or_else(|| "could not find '<' for start tag".to_string())?;

    Ok((before + rel) as u64)
}

pub fn bin_to_mzml(mzml: &MzML) -> Result<String, String> {
    let bytes = convert_bin_to_mzml_bytes(mzml)?;
    String::from_utf8(bytes).map_err(|e| e.to_string())
}

pub fn convert_bin_to_mzml_bytes(mzml: &MzML) -> Result<Vec<u8>, String> {
    let mut writer = Writer::new_with_indent(Vec::new(), b' ', 2);

    writer
        .write_event(Event::Decl(BytesDecl::new("1.0", Some("utf-8"), None)))
        .map_err(|e| e.to_string())?;

    let mut idx_tag = BytesStart::new("indexedmzML");
    idx_tag.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
    idx_tag.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
    idx_tag.push_attribute((
        "xsi:schemaLocation",
        "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd",
    ));
    writer
        .write_event(Event::Start(idx_tag))
        .map_err(|e| e.to_string())?;

    let mut mzml_tag = BytesStart::new("mzML");
    mzml_tag.push_attribute(("xmlns", "http://psi.hupo.org/ms/mzml"));
    mzml_tag.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
    mzml_tag.push_attribute((
        "xsi:schemaLocation",
        "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd",
    ));
    mzml_tag.push_attribute(("id", mzml.run.id.as_str()));
    mzml_tag.push_attribute(("version", "1.1.0"));

    writer
        .write_event(Event::Start(mzml_tag))
        .map_err(|e| e.to_string())?;

    let mut fallback_cvl = CvList::default();
    let cvl = if let Some(cvl) = &mzml.cv_list {
        cvl
    } else {
        fallback_cvl = default_cv_list();
        &fallback_cvl
    };
    write_cv_list(&mut writer, cvl)?;

    write_file_description(&mut writer, &mzml.file_description)?;

    if let Some(rpgl) = &mzml.referenceable_param_group_list {
        write_referenceable_param_group_list(&mut writer, rpgl)?;
    }

    if let Some(sl) = &mzml.sample_list {
        write_sample_list(&mut writer, sl)?;
    }

    if let Some(il) = &mzml.instrument_list {
        write_instrument_list(&mut writer, il)?;
    }

    if let Some(sw) = &mzml.software_list {
        write_software_list(&mut writer, sw)?;
    }

    if let Some(dpl) = &mzml.data_processing_list {
        write_data_processing_list(&mut writer, dpl)?;
    }

    if let Some(ssl) = &mzml.scan_settings_list {
        write_scan_settings_list(&mut writer, ssl)?;
    }

    let fallback_default_dp = mzml
        .data_processing_list
        .as_ref()
        .and_then(|dpl| dpl.data_processing.first())
        .map(|dp| dp.id.as_str());

    let mut idx = IndexAcc::default();
    write_run(&mut writer, &mzml.run, fallback_default_dp, &mut idx)?;

    writer
        .write_event(Event::End(BytesEnd::new("mzML")))
        .map_err(|e| e.to_string())?;

    let index_list_offset = write_index_list_with_offset(&mut writer, &idx)?;
    write_index_list_offset(&mut writer, index_list_offset)?;

    writer
        .write_event(Event::End(BytesEnd::new("indexedmzML")))
        .map_err(|e| e.to_string())?;

    Ok(writer.into_inner())
}

fn default_cv_list() -> CvList {
    CvList {
        count: Some(2),
        cv: vec![
            Cv {
                id: "MS".to_string(),
                full_name: Some(
                    "Proteomics Standards Initiative Mass Spectrometry Ontology".to_string(),
                ),
                version: Some("4.1.182".to_string()),
                uri: Some(
                    "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"
                        .to_string(),
                ),
            },
            Cv {
                id: "UO".to_string(),
                full_name: Some("Unit Ontology".to_string()),
                version: Some("09:04:2014".to_string()),
                uri: Some(
                    "https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"
                        .to_string(),
                ),
            },
        ],
    }
}

pub fn write_cv_list(writer: &mut Writer<Vec<u8>>, cvl: &CvList) -> Result<(), String> {
    let count = cvl.count.unwrap_or(cvl.cv.len());
    let mut tag = BytesStart::new("cvList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for cv in &cvl.cv {
        let mut cv_tag = BytesStart::new("cv");
        cv_tag.push_attribute(("id", cv.id.as_str()));
        if let Some(v) = &cv.full_name {
            cv_tag.push_attribute(("fullName", v.as_str()));
        }
        if let Some(v) = &cv.version {
            cv_tag.push_attribute(("version", v.as_str()));
        }
        if let Some(v) = &cv.uri {
            cv_tag.push_attribute(("URI", v.as_str()));
        }

        writer
            .write_event(Event::Empty(cv_tag))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("cvList")))
        .map_err(|e| e.to_string())?;

    Ok(())
}

fn write_file_description(
    writer: &mut Writer<Vec<u8>>,
    fd: &FileDescription,
) -> Result<(), String> {
    writer
        .write_event(Event::Start(BytesStart::new("fileDescription")))
        .map_err(|e| e.to_string())?;

    writer
        .write_event(Event::Start(BytesStart::new("fileContent")))
        .map_err(|e| e.to_string())?;

    write_referenceable_param_group_refs(writer, &fd.file_content.referenceable_param_group_refs)?;
    write_cv_params(writer, &fd.file_content.cv_params)?;
    write_user_params(writer, &fd.file_content.user_params)?;

    writer
        .write_event(Event::End(BytesEnd::new("fileContent")))
        .map_err(|e| e.to_string())?;

    write_source_file_list(writer, &fd.source_file_list)?;

    for c in &fd.contacts {
        writer
            .write_event(Event::Start(BytesStart::new("contact")))
            .map_err(|e| e.to_string())?;
        write_referenceable_param_group_refs(writer, &c.referenceable_param_group_refs)?;
        write_cv_params(writer, &c.cv_params)?;
        write_user_params(writer, &c.user_params)?;
        writer
            .write_event(Event::End(BytesEnd::new("contact")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("fileDescription")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_source_file_list(
    writer: &mut Writer<Vec<u8>>,
    sfl: &SourceFileList,
) -> Result<(), String> {
    let count = sfl.count.unwrap_or(sfl.source_file.len());
    let mut tag = BytesStart::new("sourceFileList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for sf in &sfl.source_file {
        let mut sf_tag = BytesStart::new("sourceFile");
        sf_tag.push_attribute(("id", sf.id.as_str()));
        if !sf.name.is_empty() {
            sf_tag.push_attribute(("name", sf.name.as_str()));
        }
        if !sf.location.is_empty() {
            sf_tag.push_attribute(("location", sf.location.as_str()));
        }

        writer
            .write_event(Event::Start(sf_tag))
            .map_err(|e| e.to_string())?;

        write_referenceable_param_group_refs(writer, &sf.referenceable_param_group_ref)?;
        write_cv_params(writer, &sf.cv_param)?;
        write_user_params(writer, &sf.user_param)?;

        writer
            .write_event(Event::End(BytesEnd::new("sourceFile")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("sourceFileList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_referenceable_param_group_list(
    writer: &mut Writer<Vec<u8>>,
    list: &ReferenceableParamGroupList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.referenceable_param_groups.len());
    let mut tag = BytesStart::new("referenceableParamGroupList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for g in &list.referenceable_param_groups {
        let mut g_tag = BytesStart::new("referenceableParamGroup");
        g_tag.push_attribute(("id", g.id.as_str()));
        writer
            .write_event(Event::Start(g_tag))
            .map_err(|e| e.to_string())?;

        write_cv_params(writer, &g.cv_params)?;
        write_user_params(writer, &g.user_params)?;

        writer
            .write_event(Event::End(BytesEnd::new("referenceableParamGroup")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("referenceableParamGroupList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_sample_list(writer: &mut Writer<Vec<u8>>, list: &SampleList) -> Result<(), String> {
    let count = list.count.unwrap_or(list.samples.len() as u32) as usize;
    let mut tag = BytesStart::new("sampleList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for s in &list.samples {
        let mut s_tag = BytesStart::new("sample");
        s_tag.push_attribute(("id", s.id.as_str()));
        if !s.name.is_empty() {
            s_tag.push_attribute(("name", s.name.as_str()));
        }
        writer
            .write_event(Event::Start(s_tag))
            .map_err(|e| e.to_string())?;

        if let Some(r) = &s.referenceable_param_group_ref {
            write_referenceable_param_group_ref(writer, r)?;
        }
        // write_cv_params(writer, &s.cv_params)?;
        // write_user_params(writer, &s.user_params)?;

        writer
            .write_event(Event::End(BytesEnd::new("sample")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("sampleList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_instrument_list(
    writer: &mut Writer<Vec<u8>>,
    list: &InstrumentList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.instrument.len());
    let mut tag = BytesStart::new("instrumentConfigurationList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for ic in &list.instrument {
        let mut ic_tag = BytesStart::new("instrumentConfiguration");
        ic_tag.push_attribute(("id", ic.id.as_str()));
        if let Some(ssr) = &ic.scan_settings_ref {
            if let Some(v) = nonempty(Some(ssr.r#ref.as_str())) {
                ic_tag.push_attribute(("scanSettingsRef", v));
            }
        }

        writer
            .write_event(Event::Start(ic_tag))
            .map_err(|e| e.to_string())?;

        write_referenceable_param_group_refs(writer, &ic.referenceable_param_group_ref)?;

        let mut top_level = Vec::new();
        let mut src = Vec::new();
        let mut an = Vec::new();
        let mut det = Vec::new();

        for p in &ic.cv_param {
            let acc = p.accession.as_deref().unwrap_or("");
            match acc {
                "MS:1000073" | "MS:1000057" => src.push(p.clone()),
                "MS:1000081" | "MS:1000084" => an.push(p.clone()),
                "MS:1000114" | "MS:1000116" => det.push(p.clone()),
                _ => top_level.push(p.clone()),
            }
        }

        write_cv_params(writer, &top_level)?;
        write_user_params(writer, &ic.user_param)?;

        if let Some(cl) = &ic.component_list {
            write_component_list(writer, cl)?;
        } else if !src.is_empty() || !an.is_empty() || !det.is_empty() {
            write_component_list_fallback(writer, &src, &an, &det)?;
        }

        if let Some(sw) = &ic.software_ref {
            let mut sw_tag = BytesStart::new("softwareRef");
            sw_tag.push_attribute(("ref", sw.r#ref.as_str()));
            writer
                .write_event(Event::Empty(sw_tag))
                .map_err(|e| e.to_string())?;
        }

        writer
            .write_event(Event::End(BytesEnd::new("instrumentConfiguration")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("instrumentConfigurationList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_component_list(writer: &mut Writer<Vec<u8>>, cl: &ComponentList) -> Result<(), String> {
    let count = cl
        .count
        .unwrap_or(cl.source.len() + cl.analyzer.len() + cl.detector.len());
    let mut tag = BytesStart::new("componentList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for s in &cl.source {
        write_component(
            writer,
            "source",
            s.order,
            &s.referenceable_param_group_ref,
            &s.cv_param,
            &s.user_param,
        )?;
    }
    for a in &cl.analyzer {
        write_component(
            writer,
            "analyzer",
            a.order,
            &a.referenceable_param_group_ref,
            &a.cv_param,
            &a.user_param,
        )?;
    }
    for d in &cl.detector {
        write_component(
            writer,
            "detector",
            d.order,
            &d.referenceable_param_group_ref,
            &d.cv_param,
            &d.user_param,
        )?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("componentList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_component_list_fallback(
    writer: &mut Writer<Vec<u8>>,
    src: &Vec<CvParam>,
    an: &Vec<CvParam>,
    det: &Vec<CvParam>,
) -> Result<(), String> {
    let mut tag = BytesStart::new("componentList");
    tag.push_attribute(("count", "3"));
    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    write_component(writer, "source", Some(1), &Vec::new(), src, &Vec::new())?;
    write_component(writer, "analyzer", Some(2), &Vec::new(), an, &Vec::new())?;
    write_component(writer, "detector", Some(3), &Vec::new(), det, &Vec::new())?;

    writer
        .write_event(Event::End(BytesEnd::new("componentList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_component(
    writer: &mut Writer<Vec<u8>>,
    name: &str,
    order: Option<u32>,
    refs: &Vec<ReferenceableParamGroupRef>,
    cvs: &Vec<CvParam>,
    ups: &Vec<UserParam>,
) -> Result<(), String> {
    let mut tag = BytesStart::new(name);
    if let Some(o) = order {
        let o_s = o.to_string();
        tag.push_attribute(("order", o_s.as_str()));
    }

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    write_referenceable_param_group_refs(writer, refs)?;
    write_cv_params(writer, cvs)?;
    write_user_params(writer, ups)?;

    writer
        .write_event(Event::End(BytesEnd::new(name)))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_software_list(writer: &mut Writer<Vec<u8>>, list: &SoftwareList) -> Result<(), String> {
    let count = list.count.unwrap_or(list.software.len());
    let mut tag = BytesStart::new("softwareList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for sw in &list.software {
        let mut sw_tag = BytesStart::new("software");
        sw_tag.push_attribute(("id", sw.id.as_str()));
        if let Some(v) = &sw.version {
            sw_tag.push_attribute(("version", v.as_str()));
        }

        writer
            .write_event(Event::Start(sw_tag))
            .map_err(|e| e.to_string())?;

        for sp in &sw.software_param {
            let mut sp_tag = BytesStart::new("softwareParam");
            if let Some(v) = &sp.cv_ref {
                sp_tag.push_attribute(("cvRef", v.as_str()));
            }
            sp_tag.push_attribute(("accession", sp.accession.as_str()));
            sp_tag.push_attribute(("name", sp.name.as_str()));
            if let Some(v) = &sp.version {
                sp_tag.push_attribute(("version", v.as_str()));
            }
            writer
                .write_event(Event::Empty(sp_tag))
                .map_err(|e| e.to_string())?;
        }

        write_cv_params(writer, &sw.cv_param)?;

        writer
            .write_event(Event::End(BytesEnd::new("software")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("softwareList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_data_processing_list(
    writer: &mut Writer<Vec<u8>>,
    list: &DataProcessingList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.data_processing.len());
    let mut tag = BytesStart::new("dataProcessingList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for dp in &list.data_processing {
        let mut dp_tag = BytesStart::new("dataProcessing");
        dp_tag.push_attribute(("id", dp.id.as_str()));
        if let Some(sw) = nonempty(dp.software_ref.as_deref()) {
            dp_tag.push_attribute(("softwareRef", sw));
        }

        writer
            .write_event(Event::Start(dp_tag))
            .map_err(|e| e.to_string())?;

        for m in &dp.processing_method {
            let mut pm = BytesStart::new("processingMethod");
            if let Some(order) = m.order {
                let s = order.to_string();
                pm.push_attribute(("order", s.as_str()));
            }
            if let Some(sw) = nonempty(m.software_ref.as_deref()) {
                pm.push_attribute(("softwareRef", sw));
            }

            writer
                .write_event(Event::Start(pm))
                .map_err(|e| e.to_string())?;

            write_referenceable_param_group_refs(writer, &m.referenceable_param_group_ref)?;
            write_cv_params(writer, &m.cv_param)?;
            write_user_params(writer, &m.user_param)?;

            writer
                .write_event(Event::End(BytesEnd::new("processingMethod")))
                .map_err(|e| e.to_string())?;
        }

        writer
            .write_event(Event::End(BytesEnd::new("dataProcessing")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("dataProcessingList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_scan_settings_list(
    writer: &mut Writer<Vec<u8>>,
    list: &ScanSettingsList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.scan_settings.len());
    let mut tag = BytesStart::new("scanSettingsList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for ss in &list.scan_settings {
        let mut ss_tag = BytesStart::new("scanSettings");
        if let Some(id) = ss.id.as_deref() {
            if !id.is_empty() {
                ss_tag.push_attribute(("id", id));
            }
        }
        if let Some(r) = ss.instrument_configuration_ref.as_deref() {
            if !r.is_empty() {
                ss_tag.push_attribute(("instrumentConfigurationRef", r));
            }
        }

        writer
            .write_event(Event::Start(ss_tag))
            .map_err(|e| e.to_string())?;

        write_referenceable_param_group_refs(writer, &ss.referenceable_param_group_refs)?;
        write_cv_params(writer, &ss.cv_params)?;
        write_user_params(writer, &ss.user_params)?;

        if let Some(sfrl) = &ss.source_file_ref_list {
            write_source_file_ref_list(writer, sfrl)?;
        }
        if let Some(tl) = &ss.target_list {
            write_target_list(writer, tl)?;
        }

        writer
            .write_event(Event::End(BytesEnd::new("scanSettings")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("scanSettingsList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_run(
    writer: &mut Writer<Vec<u8>>,
    run: &Run,
    fallback_default_dp: Option<&str>,
    idx: &mut IndexAcc,
) -> Result<(), String> {
    let mut run_tag = BytesStart::new("run");
    run_tag.push_attribute(("id", run.id.as_str()));
    if let Some(ts) = nonempty(run.start_time_stamp.as_deref()) {
        run_tag.push_attribute(("startTimeStamp", ts));
    }
    if let Some(ic) = nonempty(run.default_instrument_configuration_ref.as_deref()) {
        run_tag.push_attribute(("defaultInstrumentConfigurationRef", ic));
    }
    if let Some(sf) = nonempty(run.default_source_file_ref.as_deref()) {
        run_tag.push_attribute(("defaultSourceFileRef", sf));
    }
    if let Some(samp) = nonempty(run.sample_ref.as_deref()) {
        run_tag.push_attribute(("sampleRef", samp));
    }

    writer
        .write_event(Event::Start(run_tag))
        .map_err(|e| e.to_string())?;

    write_referenceable_param_group_refs(writer, &run.referenceable_param_group_refs)?;
    write_cv_params(writer, &run.cv_params)?;
    write_user_params(writer, &run.user_params)?;

    if let Some(sfrl) = &run.source_file_ref_list {
        write_source_file_ref_list(writer, sfrl)?;
    }

    if let Some(sl) = &run.spectrum_list {
        write_spectrum_list(writer, sl, fallback_default_dp, idx)?;
    }

    if let Some(cl) = &run.chromatogram_list {
        write_chromatogram_list(writer, cl, fallback_default_dp, idx)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("run")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn is_scanlist_level_cv(acc: &str) -> bool {
    matches!(acc, "MS:1000795")
}

fn is_scan_level_cv(acc: &str) -> bool {
    matches!(acc, "MS:1000016")
}

fn is_scanwindow_level_cv(acc: &str) -> bool {
    matches!(acc, "MS:1000501" | "MS:1000500")
}

fn is_array_meta_cv(acc: &str) -> bool {
    matches!(
        acc,
        "MS:1000523" | "MS:1000521" | "MS:1000576" | "MS:1000514" | "MS:1000515" | "MS:1000595"
    )
}

fn write_spectrum_list(
    writer: &mut Writer<Vec<u8>>,
    list: &SpectrumList,
    fallback_default_dp: Option<&str>,
    idx: &mut IndexAcc,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.spectra.len());
    let mut tag = BytesStart::new("spectrumList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    let dpr = nonempty(list.default_data_processing_ref.as_deref()).or(fallback_default_dp);
    if let Some(dp) = nonempty(dpr) {
        tag.push_attribute(("defaultDataProcessingRef", dp));
    }

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for s in &list.spectra {
        write_spectrum(writer, s, fallback_default_dp, idx)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("spectrumList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_spectrum(
    writer: &mut Writer<Vec<u8>>,
    s: &Spectrum,
    fallback_default_dp: Option<&str>,
    idx: &mut IndexAcc,
) -> Result<(), String> {
    let default_len = s
        .default_array_length
        .or_else(|| {
            s.binary_data_array_list.as_ref().and_then(|l| {
                l.binary_data_arrays.first().map(|b| {
                    b.array_length
                        .or_else(|| {
                            if !b.decoded_binary_f64.is_empty() {
                                Some(b.decoded_binary_f64.len())
                            } else if !b.decoded_binary_f32.is_empty() {
                                Some(b.decoded_binary_f32.len())
                            } else {
                                None
                            }
                        })
                        .unwrap_or(0)
                })
            })
        })
        .unwrap_or(0);

    let mut tag = BytesStart::new("spectrum");

    if let Some(idx0) = s.index {
        let idx_s = idx0.to_string();
        tag.push_attribute(("index", idx_s.as_str()));
    }

    let id_to_write = nonempty(Some(s.id.as_str()))
        .or_else(|| nonempty(s.native_id.as_deref()))
        .unwrap_or(s.id.as_str());
    tag.push_attribute(("id", id_to_write));

    let len_s = default_len.to_string();
    tag.push_attribute(("defaultArrayLength", len_s.as_str()));

    if let Some(sn) = s.scan_number {
        let sn_s = sn.to_string();
        tag.push_attribute(("scanNumber", sn_s.as_str()));
    }

    if let Some(v) = nonempty(s.native_id.as_deref()) {
        tag.push_attribute(("nativeID", v));
    }

    let dpr = nonempty(s.data_processing_ref.as_deref()).or(fallback_default_dp);
    if let Some(v) = nonempty(dpr) {
        tag.push_attribute(("dataProcessingRef", v));
    }

    if let Some(v) = nonempty(s.source_file_ref.as_deref()) {
        tag.push_attribute(("sourceFileRef", v));
    }
    if let Some(v) = nonempty(s.spot_id.as_deref()) {
        tag.push_attribute(("spotID", v));
    }

    let off = write_start_capture_offset(writer, tag)?;
    idx.spectrum.push(IndexOffsetAcc {
        id_ref: id_to_write.to_string(),
        offset: off,
    });

    write_referenceable_param_group_refs(writer, &s.referenceable_param_group_refs)?;

    let mut has_scan_related = false;
    for p in &s.cv_params {
        let acc = p.accession.as_deref().unwrap_or("");
        if is_scanlist_level_cv(acc) || is_scan_level_cv(acc) || is_scanwindow_level_cv(acc) {
            has_scan_related = true;
            break;
        }
    }

    for p in &s.cv_params {
        let acc = p.accession.as_deref().unwrap_or("");
        if is_scanlist_level_cv(acc) || is_scan_level_cv(acc) || is_scanwindow_level_cv(acc) {
            continue;
        }
        if is_array_meta_cv(acc) {
            continue;
        }
        write_cv_params(writer, &vec![p.clone()])?;
    }

    write_user_params(writer, &s.user_params)?;

    if let Some(sd) = &s.spectrum_description {
        write_spectrum_description(writer, sd)?;
    }

    if let Some(sl) = &s.scan_list {
        write_scan_list(writer, sl)?;
    } else if has_scan_related {
        write_scan_list_from_spectrum_cv(writer, &s.cv_params)?;
    }

    if let Some(pl) = &s.precursor_list {
        write_precursor_list(writer, pl)?;
    }
    if let Some(pr) = &s.product_list {
        write_product_list(writer, pr)?;
    }
    if let Some(bdal) = &s.binary_data_array_list {
        write_binary_data_array_list(writer, bdal, fallback_default_dp)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("spectrum")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_scan_list_from_spectrum_cv(
    writer: &mut Writer<Vec<u8>>,
    params: &Vec<CvParam>,
) -> Result<(), String> {
    let mut sl_tag = BytesStart::new("scanList");
    sl_tag.push_attribute(("count", "1"));
    writer
        .write_event(Event::Start(sl_tag))
        .map_err(|e| e.to_string())?;

    for p in params {
        let acc = p.accession.as_deref().unwrap_or("");
        if is_scanlist_level_cv(acc) {
            write_cv_params(writer, &vec![p.clone()])?;
        }
    }

    writer
        .write_event(Event::Start(BytesStart::new("scan")))
        .map_err(|e| e.to_string())?;

    for p in params {
        let acc = p.accession.as_deref().unwrap_or("");
        if is_scan_level_cv(acc) {
            write_cv_params(writer, &vec![p.clone()])?;
        }
    }

    let mut has_win = false;
    for p in params {
        let acc = p.accession.as_deref().unwrap_or("");
        if is_scanwindow_level_cv(acc) {
            has_win = true;
            break;
        }
    }

    if has_win {
        let mut swl = BytesStart::new("scanWindowList");
        swl.push_attribute(("count", "1"));
        writer
            .write_event(Event::Start(swl))
            .map_err(|e| e.to_string())?;

        writer
            .write_event(Event::Start(BytesStart::new("scanWindow")))
            .map_err(|e| e.to_string())?;

        for p in params {
            let acc = p.accession.as_deref().unwrap_or("");
            if is_scanwindow_level_cv(acc) {
                write_cv_params(writer, &vec![p.clone()])?;
            }
        }

        writer
            .write_event(Event::End(BytesEnd::new("scanWindow")))
            .map_err(|e| e.to_string())?;

        writer
            .write_event(Event::End(BytesEnd::new("scanWindowList")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("scan")))
        .map_err(|e| e.to_string())?;

    writer
        .write_event(Event::End(BytesEnd::new("scanList")))
        .map_err(|e| e.to_string())?;

    Ok(())
}

fn write_spectrum_description(
    writer: &mut Writer<Vec<u8>>,
    sd: &SpectrumDescription,
) -> Result<(), String> {
    writer
        .write_event(Event::Start(BytesStart::new("spectrumDescription")))
        .map_err(|e| e.to_string())?;

    write_referenceable_param_group_refs(writer, &sd.referenceable_param_group_refs)?;
    write_cv_params(writer, &sd.cv_params)?;
    write_user_params(writer, &sd.user_params)?;

    if let Some(sl) = &sd.scan_list {
        write_scan_list(writer, sl)?;
    }
    if let Some(pl) = &sd.precursor_list {
        write_precursor_list(writer, pl)?;
    }
    if let Some(pr) = &sd.product_list {
        write_product_list(writer, pr)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("spectrumDescription")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_scan_list(writer: &mut Writer<Vec<u8>>, list: &ScanList) -> Result<(), String> {
    let count = list.count.unwrap_or(list.scans.len());
    let mut tag = BytesStart::new("scanList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for s in &list.scans {
        let mut st = BytesStart::new("scan");
        if let Some(v) = nonempty(s.instrument_configuration_ref.as_deref()) {
            st.push_attribute(("instrumentConfigurationRef", v));
        }
        if let Some(v) = nonempty(s.external_spectrum_id.as_deref()) {
            st.push_attribute(("externalSpectrumID", v));
        }
        if let Some(v) = nonempty(s.source_file_ref.as_deref()) {
            st.push_attribute(("sourceFileRef", v));
        }
        if let Some(v) = nonempty(s.spectrum_ref.as_deref()) {
            st.push_attribute(("spectrumRef", v));
        }

        writer
            .write_event(Event::Start(st))
            .map_err(|e| e.to_string())?;

        write_referenceable_param_group_refs(writer, &s.referenceable_param_group_refs)?;
        write_cv_params(writer, &s.cv_params)?;
        write_user_params(writer, &s.user_params)?;

        if let Some(swl) = &s.scan_window_list {
            write_scan_window_list(writer, swl)?;
        }

        writer
            .write_event(Event::End(BytesEnd::new("scan")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("scanList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_scan_window_list(
    writer: &mut Writer<Vec<u8>>,
    list: &ScanWindowList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.scan_windows.len());
    let mut tag = BytesStart::new("scanWindowList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for w in &list.scan_windows {
        writer
            .write_event(Event::Start(BytesStart::new("scanWindow")))
            .map_err(|e| e.to_string())?;
        write_cv_params(writer, &w.cv_params)?;
        write_user_params(writer, &w.user_params)?;
        writer
            .write_event(Event::End(BytesEnd::new("scanWindow")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("scanWindowList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_precursor_list(writer: &mut Writer<Vec<u8>>, list: &PrecursorList) -> Result<(), String> {
    let count = list.count.unwrap_or(list.precursors.len());
    let mut tag = BytesStart::new("precursorList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for p in &list.precursors {
        write_precursor(writer, p)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("precursorList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_precursor(writer: &mut Writer<Vec<u8>>, p: &Precursor) -> Result<(), String> {
    let mut pt = BytesStart::new("precursor");
    if let Some(v) = nonempty(p.spectrum_ref.as_deref()) {
        pt.push_attribute(("spectrumRef", v));
    }
    if let Some(v) = nonempty(p.source_file_ref.as_deref()) {
        pt.push_attribute(("sourceFileRef", v));
    }
    if let Some(v) = nonempty(p.external_spectrum_id.as_deref()) {
        pt.push_attribute(("externalSpectrumID", v));
    }

    writer
        .write_event(Event::Start(pt))
        .map_err(|e| e.to_string())?;

    if let Some(iw) = &p.isolation_window {
        write_cv_container(
            writer,
            "isolationWindow",
            &iw.referenceable_param_group_refs,
            &iw.cv_params,
            &iw.user_params,
        )?;
    }
    if let Some(sil) = &p.selected_ion_list {
        write_selected_ion_list(writer, sil)?;
    }
    if let Some(act) = &p.activation {
        write_cv_container(
            writer,
            "activation",
            &act.referenceable_param_group_refs,
            &act.cv_params,
            &act.user_params,
        )?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("precursor")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_selected_ion_list(
    writer: &mut Writer<Vec<u8>>,
    list: &SelectedIonList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.selected_ions.len());
    let mut tag = BytesStart::new("selectedIonList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for si in &list.selected_ions {
        writer
            .write_event(Event::Start(BytesStart::new("selectedIon")))
            .map_err(|e| e.to_string())?;

        write_referenceable_param_group_refs(writer, &si.referenceable_param_group_refs)?;
        write_cv_params(writer, &si.cv_params)?;
        write_user_params(writer, &si.user_params)?;

        writer
            .write_event(Event::End(BytesEnd::new("selectedIon")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("selectedIonList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_product_list(writer: &mut Writer<Vec<u8>>, list: &ProductList) -> Result<(), String> {
    let count = list.count.unwrap_or(list.products.len());
    let mut tag = BytesStart::new("productList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for p in &list.products {
        write_product(writer, p)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("productList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_product(writer: &mut Writer<Vec<u8>>, p: &Product) -> Result<(), String> {
    let mut pt = BytesStart::new("product");
    if let Some(v) = nonempty(p.spectrum_ref.as_deref()) {
        pt.push_attribute(("spectrumRef", v));
    }
    if let Some(v) = nonempty(p.source_file_ref.as_deref()) {
        pt.push_attribute(("sourceFileRef", v));
    }
    if let Some(v) = nonempty(p.external_spectrum_id.as_deref()) {
        pt.push_attribute(("externalSpectrumID", v));
    }

    writer
        .write_event(Event::Start(pt))
        .map_err(|e| e.to_string())?;

    if let Some(iw) = &p.isolation_window {
        write_cv_container(
            writer,
            "isolationWindow",
            &iw.referenceable_param_group_refs,
            &iw.cv_params,
            &iw.user_params,
        )?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("product")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_chromatogram_list(
    writer: &mut Writer<Vec<u8>>,
    list: &ChromatogramList,
    fallback_default_dp: Option<&str>,
    idx: &mut IndexAcc,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.chromatograms.len());
    let mut tag = BytesStart::new("chromatogramList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    let dpr = nonempty(list.default_data_processing_ref.as_deref()).or(fallback_default_dp);
    if let Some(dp) = nonempty(dpr) {
        tag.push_attribute(("defaultDataProcessingRef", dp));
    }

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for c in &list.chromatograms {
        write_chromatogram(writer, c, fallback_default_dp, idx)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("chromatogramList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_chromatogram(
    writer: &mut Writer<Vec<u8>>,
    c: &Chromatogram,
    fallback_default_dp: Option<&str>,
    idx: &mut IndexAcc,
) -> Result<(), String> {
    let default_len = c
        .default_array_length
        .or_else(|| {
            c.binary_data_array_list.as_ref().and_then(|l| {
                l.binary_data_arrays.first().map(|b| {
                    b.array_length
                        .or_else(|| {
                            if !b.decoded_binary_f64.is_empty() {
                                Some(b.decoded_binary_f64.len())
                            } else if !b.decoded_binary_f32.is_empty() {
                                Some(b.decoded_binary_f32.len())
                            } else {
                                None
                            }
                        })
                        .unwrap_or(0)
                })
            })
        })
        .unwrap_or(0);

    let mut tag = BytesStart::new("chromatogram");
    tag.push_attribute(("id", c.id.as_str()));
    if let Some(v) = nonempty(c.native_id.as_deref()) {
        tag.push_attribute(("nativeID", v));
    }
    if let Some(idx0) = c.index {
        let idx_s = idx0.to_string();
        tag.push_attribute(("index", idx_s.as_str()));
    }
    let len_s = default_len.to_string();
    tag.push_attribute(("defaultArrayLength", len_s.as_str()));

    let dpr = nonempty(c.data_processing_ref.as_deref()).or(fallback_default_dp);
    if let Some(v) = nonempty(dpr) {
        tag.push_attribute(("dataProcessingRef", v));
    }

    let off = write_start_capture_offset(writer, tag)?;
    idx.chromatogram.push(IndexOffsetAcc {
        id_ref: c.id.clone(),
        offset: off,
    });

    write_referenceable_param_group_refs(writer, &c.referenceable_param_group_refs)?;
    write_cv_params(writer, &c.cv_params)?;
    write_user_params(writer, &c.user_params)?;

    if let Some(p) = &c.precursor {
        write_precursor(writer, p)?;
    }
    if let Some(p) = &c.product {
        write_product(writer, p)?;
    }

    if let Some(bdal) = &c.binary_data_array_list {
        write_binary_data_array_list(writer, bdal, fallback_default_dp)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("chromatogram")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_binary_data_array_list(
    writer: &mut Writer<Vec<u8>>,
    list: &BinaryDataArrayList,
    fallback_default_dp: Option<&str>,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.binary_data_arrays.len());
    let mut tag = BytesStart::new("binaryDataArrayList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for bda in &list.binary_data_arrays {
        write_binary_data_array(writer, bda, fallback_default_dp)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("binaryDataArrayList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn has_accession(params: &Vec<CvParam>, acc: &str) -> bool {
    params.iter().any(|p| p.accession.as_deref() == Some(acc))
}

fn write_binary_data_array(
    writer: &mut Writer<Vec<u8>>,
    bda: &BinaryDataArray,
    fallback_default_dp: Option<&str>,
) -> Result<(), String> {
    let cv_has_f64 = has_accession(&bda.cv_params, "MS:1000523");
    let cv_has_f32 = has_accession(&bda.cv_params, "MS:1000521");

    let is_f64 = bda
        .is_f64
        .unwrap_or(cv_has_f64 || (!cv_has_f32 && !bda.decoded_binary_f64.is_empty()));
    let is_f32 = bda
        .is_f32
        .unwrap_or(cv_has_f32 || (!cv_has_f64 && !bda.decoded_binary_f32.is_empty()));

    let encoded_from_field: Option<&str> = None;

    let (raw_bytes, array_len) = if encoded_from_field.is_some() {
        (Vec::new(), bda.array_length.unwrap_or(0))
    } else if cv_has_f64 || (is_f64 && !bda.decoded_binary_f64.is_empty()) {
        let mut bytes = Vec::with_capacity(bda.decoded_binary_f64.len() * 8);
        for v in &bda.decoded_binary_f64 {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        (bytes, bda.decoded_binary_f64.len())
    } else if cv_has_f32 || (is_f32 && !bda.decoded_binary_f32.is_empty()) {
        let mut bytes = Vec::with_capacity(bda.decoded_binary_f32.len() * 4);
        for v in &bda.decoded_binary_f32 {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        (bytes, bda.decoded_binary_f32.len())
    } else {
        (Vec::new(), bda.array_length.unwrap_or(0))
    };

    let encoded = if let Some(s) = encoded_from_field {
        s.to_string()
    } else if raw_bytes.is_empty() {
        String::new()
    } else {
        STANDARD.encode(&raw_bytes)
    };

    let mut tag = BytesStart::new("binaryDataArray");

    if let Some(al) = bda.array_length.or(Some(array_len)).filter(|&v| v > 0) {
        let al_s = al.to_string();
        tag.push_attribute(("arrayLength", al_s.as_str()));
    }

    let el_s = encoded.len().to_string();
    tag.push_attribute(("encodedLength", el_s.as_str()));

    let dpr = nonempty(bda.data_processing_ref.as_deref()).or(fallback_default_dp);
    if let Some(dp) = nonempty(dpr) {
        tag.push_attribute(("dataProcessingRef", dp));
    }

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    write_referenceable_param_group_refs(writer, &bda.referenceable_param_group_refs)?;
    write_cv_params(writer, &bda.cv_params)?;

    if is_f64 && !cv_has_f64 && !cv_has_f32 {
        write_simple_cv(writer, "64-bit float")?;
    } else if is_f32 && !cv_has_f32 && !cv_has_f64 {
        write_simple_cv(writer, "32-bit float")?;
    }

    if !has_accession(&bda.cv_params, "MS:1000576") {
        write_simple_cv(writer, "no compression")?;
    }

    write_user_params(writer, &bda.user_params)?;

    writer
        .write_event(Event::Start(BytesStart::new("binary")))
        .map_err(|e| e.to_string())?;
    if !encoded.is_empty() {
        writer
            .write_event(Event::Text(BytesText::new(encoded.as_str())))
            .map_err(|e| e.to_string())?;
    }
    writer
        .write_event(Event::End(BytesEnd::new("binary")))
        .map_err(|e| e.to_string())?;

    writer
        .write_event(Event::End(BytesEnd::new("binaryDataArray")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_target_list(writer: &mut Writer<Vec<u8>>, list: &TargetList) -> Result<(), String> {
    let count = list.count.unwrap_or(list.targets.len());
    let mut tag = BytesStart::new("targetList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for t in &list.targets {
        writer
            .write_event(Event::Start(BytesStart::new("target")))
            .map_err(|e| e.to_string())?;

        write_referenceable_param_group_refs(writer, &t.referenceable_param_group_refs)?;
        write_cv_params(writer, &t.cv_params)?;
        write_user_params(writer, &t.user_params)?;

        writer
            .write_event(Event::End(BytesEnd::new("target")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("targetList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_source_file_ref_list(
    writer: &mut Writer<Vec<u8>>,
    list: &SourceFileRefList,
) -> Result<(), String> {
    let count = list.count.unwrap_or(list.source_file_refs.len());
    let mut tag = BytesStart::new("sourceFileRefList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for r in &list.source_file_refs {
        let mut rf = BytesStart::new("sourceFileRef");
        rf.push_attribute(("ref", r.r#ref.as_str()));
        writer
            .write_event(Event::Empty(rf))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("sourceFileRefList")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_referenceable_param_group_refs(
    writer: &mut Writer<Vec<u8>>,
    refs: &Vec<ReferenceableParamGroupRef>,
) -> Result<(), String> {
    for r in refs {
        write_referenceable_param_group_ref(writer, r)?;
    }
    Ok(())
}

fn write_referenceable_param_group_ref(
    writer: &mut Writer<Vec<u8>>,
    r: &ReferenceableParamGroupRef,
) -> Result<(), String> {
    let mut tag = BytesStart::new("referenceableParamGroupRef");
    tag.push_attribute(("ref", r.r#ref.as_str()));
    writer
        .write_event(Event::Empty(tag))
        .map_err(|e| e.to_string())
}

fn write_cv_params(writer: &mut Writer<Vec<u8>>, params: &Vec<CvParam>) -> Result<(), String> {
    for cv in params {
        let mut tag = BytesStart::new("cvParam");

        if let Some(v) = cv.cv_ref.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("cvRef", v));
        }
        if let Some(v) = cv.accession.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("accession", v));
        }
        tag.push_attribute(("name", cv.name.as_str()));

        let value_s = cv.value.as_deref().unwrap_or("");
        tag.push_attribute(("value", value_s));

        if let Some(v) = cv.unit_cv_ref.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("unitCvRef", v));
        }
        if let Some(v) = cv.unit_accession.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("unitAccession", v));
        }
        if let Some(v) = cv.unit_name.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("unitName", v));
        }

        writer
            .write_event(Event::Empty(tag))
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn write_user_params(writer: &mut Writer<Vec<u8>>, params: &Vec<UserParam>) -> Result<(), String> {
    for up in params {
        let mut tag = BytesStart::new("userParam");
        tag.push_attribute(("name", up.name.as_str()));

        if let Some(v) = up.r#type.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("type", v));
        }

        let value_s = up.value.as_deref().unwrap_or("");
        tag.push_attribute(("value", value_s));

        if let Some(v) = up.unit_cv_ref.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("unitCvRef", v));
        }
        if let Some(v) = up.unit_accession.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("unitAccession", v));
        }
        if let Some(v) = up.unit_name.as_deref().and_then(|s| nonempty(Some(s))) {
            tag.push_attribute(("unitName", v));
        }

        writer
            .write_event(Event::Empty(tag))
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn write_cv_container(
    writer: &mut Writer<Vec<u8>>,
    tag_name: &str,
    refs: &Vec<ReferenceableParamGroupRef>,
    cvs: &Vec<CvParam>,
    ups: &Vec<UserParam>,
) -> Result<(), String> {
    writer
        .write_event(Event::Start(BytesStart::new(tag_name)))
        .map_err(|e| e.to_string())?;

    write_referenceable_param_group_refs(writer, refs)?;
    write_cv_params(writer, cvs)?;
    write_user_params(writer, ups)?;

    writer
        .write_event(Event::End(BytesEnd::new(tag_name)))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_simple_cv(writer: &mut Writer<Vec<u8>>, name: &str) -> Result<(), String> {
    let mut tag = BytesStart::new("cvParam");

    match name {
        "m/z array" => {
            tag.push_attribute(("cvRef", "MS"));
            tag.push_attribute(("accession", "MS:1000514"));
            tag.push_attribute(("name", "m/z array"));
            tag.push_attribute(("value", ""));
        }
        "intensity array" => {
            tag.push_attribute(("cvRef", "MS"));
            tag.push_attribute(("accession", "MS:1000515"));
            tag.push_attribute(("name", "intensity array"));
            tag.push_attribute(("value", ""));
        }
        "time array" => {
            tag.push_attribute(("cvRef", "MS"));
            tag.push_attribute(("accession", "MS:1000595"));
            tag.push_attribute(("name", "time array"));
            tag.push_attribute(("value", ""));
        }
        "32-bit float" => {
            tag.push_attribute(("cvRef", "MS"));
            tag.push_attribute(("accession", "MS:1000521"));
            tag.push_attribute(("name", "32-bit float"));
            tag.push_attribute(("value", ""));
        }
        "64-bit float" => {
            tag.push_attribute(("cvRef", "MS"));
            tag.push_attribute(("accession", "MS:1000523"));
            tag.push_attribute(("name", "64-bit float"));
            tag.push_attribute(("value", ""));
        }
        "no compression" => {
            tag.push_attribute(("cvRef", "MS"));
            tag.push_attribute(("accession", "MS:1000576"));
            tag.push_attribute(("name", "no compression"));
            tag.push_attribute(("value", ""));
        }
        _ => {
            tag.push_attribute(("name", name));
            tag.push_attribute(("value", ""));
        }
    }

    writer
        .write_event(Event::Empty(tag))
        .map_err(|e| e.to_string())
}

fn write_index_list_with_offset(
    writer: &mut Writer<Vec<u8>>,
    idx: &IndexAcc,
) -> Result<u64, String> {
    let mut count = 0usize;
    if !idx.spectrum.is_empty() {
        count += 1;
    }
    if !idx.chromatogram.is_empty() {
        count += 1;
    }

    let mut tag = BytesStart::new("indexList");
    let count_s = count.to_string();
    tag.push_attribute(("count", count_s.as_str()));

    let off = write_start_capture_offset(writer, tag)?;

    if !idx.spectrum.is_empty() {
        write_index(writer, "spectrum", &idx.spectrum)?;
    }
    if !idx.chromatogram.is_empty() {
        write_index(writer, "chromatogram", &idx.chromatogram)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("indexList")))
        .map_err(|e| e.to_string())?;

    Ok(off)
}

fn write_index(
    writer: &mut Writer<Vec<u8>>,
    name: &str,
    offsets: &Vec<IndexOffsetAcc>,
) -> Result<(), String> {
    let mut tag = BytesStart::new("index");
    tag.push_attribute(("name", name));

    writer
        .write_event(Event::Start(tag))
        .map_err(|e| e.to_string())?;

    for o in offsets {
        let mut ot = BytesStart::new("offset");
        ot.push_attribute(("idRef", o.id_ref.as_str()));

        writer
            .write_event(Event::Start(ot))
            .map_err(|e| e.to_string())?;

        let s = o.offset.to_string();
        writer
            .write_event(Event::Text(BytesText::new(s.as_str())))
            .map_err(|e| e.to_string())?;

        writer
            .write_event(Event::End(BytesEnd::new("offset")))
            .map_err(|e| e.to_string())?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("index")))
        .map_err(|e| e.to_string())?;
    Ok(())
}

fn write_index_list_offset(writer: &mut Writer<Vec<u8>>, off: u64) -> Result<(), String> {
    writer
        .write_event(Event::Start(BytesStart::new("indexListOffset")))
        .map_err(|e| e.to_string())?;

    let s = off.to_string();
    writer
        .write_event(Event::Text(BytesText::new(s.as_str())))
        .map_err(|e| e.to_string())?;

    writer
        .write_event(Event::End(BytesEnd::new("indexListOffset")))
        .map_err(|e| e.to_string())?;
    Ok(())
}
