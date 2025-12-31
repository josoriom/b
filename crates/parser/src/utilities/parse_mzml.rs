use base64::Engine;
use base64::engine::general_purpose::STANDARD;
use miniz_oxide::inflate::decompress_to_vec_zlib;
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use serde::{Deserialize, Serialize};
use std::io::{BufRead, Cursor};

use crate::utilities::mzml::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CvEntry {
    pub id: String,
    pub full_name: Option<String>,
    pub version: Option<String>,
    pub uri: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexOffset {
    pub id_ref: Option<String>,
    pub offset: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexList {
    pub spectrum: Vec<IndexOffset>,
    pub chromatogram: Vec<IndexOffset>,
    pub index_list_offset: Option<u64>,
    pub file_checksum: Option<String>,
}

fn get_attr_any(start: &BytesStart, names: &[&[u8]]) -> Option<String> {
    for a in start.attributes().with_checks(false).flatten() {
        let key = a.key.as_ref();
        if names.iter().any(|n| *n == key) {
            return a.unescape_value().ok().map(|v| v.to_string());
        }
    }
    None
}

fn get_attr(start: &BytesStart, name: &[u8]) -> Option<String> {
    get_attr_any(start, &[name])
}

fn get_attr_u32(start: &BytesStart, name: &[u8]) -> Option<u32> {
    get_attr(start, name).and_then(|s| s.parse().ok())
}

fn get_attr_usize(start: &BytesStart, name: &[u8]) -> Option<usize> {
    get_attr(start, name).and_then(|s| s.parse().ok())
}

#[inline]
fn local_name(mut raw: &[u8]) -> &[u8] {
    if raw.first() == Some(&b'{') {
        if let Some(end) = raw.iter().position(|&b| b == b'}') {
            raw = &raw[end + 1..];
        }
    }
    if let Some(colon) = raw.iter().rposition(|&b| b == b':') {
        &raw[colon + 1..]
    } else {
        raw
    }
}

fn parse_referenceable_param_group_ref(start: &BytesStart) -> ReferenceableParamGroupRef {
    ReferenceableParamGroupRef {
        r#ref: get_attr(start, b"ref").unwrap_or_default(),
    }
}

fn parse_cv_param(start: &BytesStart) -> CvParam {
    CvParam {
        cv_ref: get_attr_any(start, &[b"cvRef", b"cvLabel"]),
        accession: get_attr(start, b"accession"),
        name: get_attr(start, b"name").unwrap_or_default(),
        value: get_attr(start, b"value"),
        unit_cv_ref: get_attr_any(start, &[b"unitCvRef", b"unitCvLabel"]),
        unit_name: get_attr(start, b"unitName"),
        unit_accession: get_attr(start, b"unitAccession"),
    }
}

fn parse_user_param(start: &BytesStart) -> UserParam {
    UserParam {
        name: get_attr(start, b"name").unwrap_or_default(),
        r#type: get_attr(start, b"type"),
        unit_accession: get_attr(start, b"unitAccession"),
        unit_cv_ref: get_attr_any(start, &[b"unitCvRef", b"unitCvLabel"]),
        unit_name: get_attr(start, b"unitName"),
        value: get_attr(start, b"value"),
    }
}

fn parse_software_param(start: &BytesStart) -> SoftwareParam {
    SoftwareParam {
        cv_ref: get_attr_any(start, &[b"cvRef", b"cvLabel"]),
        accession: get_attr(start, b"accession").unwrap_or_default(),
        name: get_attr(start, b"name").unwrap_or_default(),
        version: get_attr(start, b"version"),
    }
}

fn parse_source_file_ref(start: &BytesStart) -> SourceFileRef {
    SourceFileRef {
        r#ref: get_attr(start, b"ref").unwrap_or_default(),
    }
}

fn skip_element<R: BufRead>(reader: &mut Reader<R>, end: &[u8]) -> Result<(), String> {
    let mut depth = 1usize;
    let mut buf = Vec::with_capacity(512);

    while depth != 0 {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(_) => depth += 1,
            Event::End(e) => {
                if depth == 1 && e.name().as_ref() == end {
                    break;
                }
                depth = depth.saturating_sub(1);
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(())
}

fn read_text_content<R: BufRead>(reader: &mut Reader<R>, end: &[u8]) -> Result<String, String> {
    let mut buf = Vec::with_capacity(512);
    let mut out = String::new();

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Text(t) => out.push_str(&t.decode().map_err(|e| e.to_string())?),
            Event::CData(t) => out.push_str(&String::from_utf8_lossy(&t.into_inner())),
            Event::End(e) if e.name().as_ref() == end => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(out)
}

fn push_params_empty(
    e: &BytesStart,
    referenceable_param_group_refs: &mut Vec<ReferenceableParamGroupRef>,
    cv_params: &mut Vec<CvParam>,
    user_params: &mut Vec<UserParam>,
) -> bool {
    match e.name().as_ref() {
        b"referenceableParamGroupRef" => {
            referenceable_param_group_refs.push(parse_referenceable_param_group_ref(e));
            true
        }
        b"cvParam" => {
            cv_params.push(parse_cv_param(e));
            true
        }
        b"userParam" => {
            user_params.push(parse_user_param(e));
            true
        }
        _ => false,
    }
}

fn push_params_start<R: BufRead>(
    reader: &mut Reader<R>,
    e: &BytesStart,
    referenceable_param_group_refs: &mut Vec<ReferenceableParamGroupRef>,
    cv_params: &mut Vec<CvParam>,
    user_params: &mut Vec<UserParam>,
) -> Result<bool, String> {
    match e.name().as_ref() {
        b"referenceableParamGroupRef" => {
            referenceable_param_group_refs.push(parse_referenceable_param_group_ref(e));
            skip_element(reader, b"referenceableParamGroupRef")?;
            Ok(true)
        }
        b"cvParam" => {
            cv_params.push(parse_cv_param(e));
            skip_element(reader, b"cvParam")?;
            Ok(true)
        }
        b"userParam" => {
            user_params.push(parse_user_param(e));
            skip_element(reader, b"userParam")?;
            Ok(true)
        }
        _ => Ok(false),
    }
}

fn maybe_set_ms_level(spectrum: &mut Spectrum, p: &CvParam) {
    if spectrum.ms_level.is_some() || p.name != "ms level" {
        return;
    }
    let Some(v) = p.value.as_deref() else { return };
    let Ok(n) = v.parse::<u32>() else { return };
    spectrum.ms_level = Some(n);
}

/// <mzML>
pub fn parse_mzml(bytes: &[u8], slim: bool) -> Result<MzML, String> {
    let mut reader = Reader::from_reader(Cursor::new(bytes));
    reader.config_mut().trim_text(true);

    let mut buf = Vec::with_capacity(1024);
    let mut mzml = MzML::default();
    let mut in_mzml = false;

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"mzML" {
                    in_mzml = true;
                    buf.clear();
                    continue;
                }
                if !in_mzml {
                    buf.clear();
                    continue;
                }
                match e.name().as_ref() {
                    b"cvList" => mzml.cv_list = Some(parse_cv_list(&mut reader, &e)?),
                    b"fileDescription" => {
                        mzml.file_description = parse_file_description(&mut reader, &e)?
                    }
                    b"referenceableParamGroupList" => {
                        mzml.referenceable_param_group_list =
                            Some(parse_referenceable_param_group_list(&mut reader, &e)?);
                    }
                    b"sampleList" => mzml.sample_list = Some(parse_sample_list(&mut reader, &e)?),
                    b"instrumentList" | b"instrumentConfigurationList" => {
                        mzml.instrument_list = parse_instrument_list(&mut reader, &e)?;
                    }
                    b"softwareList" => {
                        mzml.software_list = Some(parse_software_list(&mut reader, &e)?)
                    }
                    b"dataProcessingList" => {
                        mzml.data_processing_list =
                            Some(parse_data_processing_list(&mut reader, &e)?);
                    }
                    b"scanSettingsList" | b"acquisitionSettingsList" => {
                        mzml.scan_settings_list = parse_scan_settings_list(&mut reader, &e)?;
                    }
                    b"run" => mzml.run = parse_run(&mut reader, &e, slim)?,
                    _ => skip_element(&mut reader, e.name().as_ref())?,
                }
            }
            Event::Empty(e) => {
                if !in_mzml {
                    buf.clear();
                    continue;
                }
                match e.name().as_ref() {
                    b"referenceableParamGroupList" => {
                        let mut list = ReferenceableParamGroupList::default();
                        list.count = get_attr_usize(&e, b"count");
                        mzml.referenceable_param_group_list = Some(list);
                    }
                    b"sampleList" => {
                        let mut list = SampleList::default();
                        list.count = get_attr_u32(&e, b"count");
                        mzml.sample_list = Some(list);
                    }
                    b"instrumentList" | b"instrumentConfigurationList" => {
                        let mut list = InstrumentList::default();
                        list.count = get_attr_usize(&e, b"count");
                        mzml.instrument_list = Some(list);
                    }
                    b"softwareList" => {
                        let mut list = SoftwareList::default();
                        list.count = get_attr_usize(&e, b"count");
                        mzml.software_list = Some(list);
                    }
                    b"dataProcessingList" => {
                        let mut list = DataProcessingList::default();
                        list.count = get_attr_usize(&e, b"count");
                        mzml.data_processing_list = Some(list);
                    }
                    b"scanSettingsList" | b"acquisitionSettingsList" => {
                        let mut list = ScanSettingsList::default();
                        list.count = get_attr_usize(&e, b"count");
                        mzml.scan_settings_list = Some(list);
                    }
                    _ => {}
                }
            }
            Event::End(e) if e.name().as_ref() == b"mzML" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(mzml)
}

/// <cvList>
pub fn parse_cv_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<CvList, String> {
    let mut list = CvList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) if e.name().as_ref() == b"cv" => list.cv.push(parse_cv_tag(&e)?),
            Event::Start(e) => {
                if e.name().as_ref() == b"cv" {
                    list.cv.push(parse_cv_tag(&e)?);
                    skip_element(reader, b"cv")?;
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"cvList" => break,
            Event::Eof => return Err("Unexpected EOF while parsing <cvList>".into()),
            _ => {}
        }
        buf.clear();
    }
    Ok(list)
}

fn parse_cv_tag(start: &BytesStart) -> Result<Cv, String> {
    Ok(Cv {
        id: get_attr_any(start, &[b"id", b"cvLabel"]).unwrap_or_default(),
        full_name: get_attr(start, b"fullName"),
        version: get_attr(start, b"version"),
        uri: get_attr(start, b"URI"),
    })
}

/// <fileDescription>
fn parse_file_description<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<FileDescription, String> {
    let mut file_content = FileContent::default();
    let mut source_file_list = SourceFileList::default();
    let mut contacts = Vec::new();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => match e.name().as_ref() {
                b"fileContent" => file_content = parse_file_content(reader, &e)?,
                b"sourceFileList" => source_file_list = parse_source_file_list(reader, &e)?,
                b"contact" => contacts.push(parse_contact(reader, &e)?),
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::Empty(e) if e.name().as_ref() == b"contact" => contacts.push(Contact::default()),
            Event::End(e) if e.name().as_ref() == b"fileDescription" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(FileDescription {
        file_content,
        source_file_list,
        contacts,
    })
}

/// <fileContent>
fn parse_file_content<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<FileContent, String> {
    let mut fc = FileContent::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut fc.referenceable_param_group_refs,
                    &mut fc.cv_params,
                    &mut fc.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut fc.referenceable_param_group_refs,
                    &mut fc.cv_params,
                    &mut fc.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"fileContent" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(fc)
}

/// <sourceFileList>
fn parse_source_file_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SourceFileList, String> {
    let mut list = SourceFileList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"sourceFile" {
                    list.source_file.push(parse_source_file(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"sourceFile" => {
                list.source_file.push(SourceFile {
                    id: get_attr(&e, b"id").unwrap_or_default(),
                    name: get_attr_any(&e, &[b"name", b"sourceFileName"]).unwrap_or_default(),
                    location: get_attr_any(&e, &[b"location", b"sourceFileLocation"])
                        .unwrap_or_default(),
                    ..Default::default()
                });
            }
            Event::End(e) if e.name().as_ref() == b"sourceFileList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <sourceFile>
fn parse_source_file<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SourceFile, String> {
    let mut sf = SourceFile {
        id: get_attr(start, b"id").unwrap_or_default(),
        name: get_attr_any(start, &[b"name", b"sourceFileName"]).unwrap_or_default(),
        location: get_attr_any(start, &[b"location", b"sourceFileLocation"]).unwrap_or_default(),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut sf.referenceable_param_group_ref,
                    &mut sf.cv_param,
                    &mut sf.user_param,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut sf.referenceable_param_group_ref,
                    &mut sf.cv_param,
                    &mut sf.user_param,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"sourceFile" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(sf)
}

/// <contact>
fn parse_contact<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<Contact, String> {
    let mut c = Contact::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut c.referenceable_param_group_refs,
                    &mut c.cv_params,
                    &mut c.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut c.referenceable_param_group_refs,
                    &mut c.cv_params,
                    &mut c.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"contact" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(c)
}

/// <referenceableParamGroupList>
fn parse_referenceable_param_group_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ReferenceableParamGroupList, String> {
    let mut list = ReferenceableParamGroupList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"referenceableParamGroup" {
                    list.referenceable_param_groups
                        .push(parse_referenceable_param_group(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"referenceableParamGroup" => {
                list.referenceable_param_groups
                    .push(ReferenceableParamGroup {
                        id: get_attr(&e, b"id").unwrap_or_default(),
                        ..Default::default()
                    });
            }
            Event::End(e) if e.name().as_ref() == b"referenceableParamGroupList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <referenceableParamGroup>
fn parse_referenceable_param_group<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ReferenceableParamGroup, String> {
    let mut group = ReferenceableParamGroup {
        id: get_attr(start, b"id").unwrap_or_default(),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => match e.name().as_ref() {
                b"cvParam" => group.cv_params.push(parse_cv_param(&e)),
                b"userParam" => group.user_params.push(parse_user_param(&e)),
                _ => {}
            },
            Event::Start(e) => match e.name().as_ref() {
                b"cvParam" => {
                    group.cv_params.push(parse_cv_param(&e));
                    skip_element(reader, b"cvParam")?;
                }
                b"userParam" => {
                    group.user_params.push(parse_user_param(&e));
                    skip_element(reader, b"userParam")?;
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::End(e) if e.name().as_ref() == b"referenceableParamGroup" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(group)
}

/// <sampleList>
fn parse_sample_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SampleList, String> {
    let mut list = SampleList::default();
    list.count = get_attr_u32(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"sample" {
                    list.samples.push(parse_sample(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"sample" => {
                list.samples.push(Sample {
                    id: get_attr(&e, b"id").unwrap_or_default(),
                    name: get_attr(&e, b"name").unwrap_or_default(),
                    ..Default::default()
                });
            }
            Event::End(e) if e.name().as_ref() == b"sampleList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(list)
}

/// <sample>
fn parse_sample<R: BufRead>(reader: &mut Reader<R>, start: &BytesStart) -> Result<Sample, String> {
    let mut sample = Sample {
        id: get_attr(start, b"id").unwrap_or_default(),
        name: get_attr(start, b"name").unwrap_or_default(),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => match e.name().as_ref() {
                b"referenceableParamGroupRef" => {
                    if sample.referenceable_param_group_ref.is_none() {
                        sample.referenceable_param_group_ref =
                            Some(parse_referenceable_param_group_ref(&e));
                    }
                }
                b"cvParam" => sample.cv_params.push(parse_cv_param(&e)),
                b"userParam" => sample.user_params.push(parse_user_param(&e)),
                _ => {}
            },
            Event::Start(e) => match e.name().as_ref() {
                b"referenceableParamGroupRef" => {
                    if sample.referenceable_param_group_ref.is_none() {
                        sample.referenceable_param_group_ref =
                            Some(parse_referenceable_param_group_ref(&e));
                    }
                    skip_element(reader, b"referenceableParamGroupRef")?;
                }
                b"cvParam" => {
                    sample.cv_params.push(parse_cv_param(&e));
                    skip_element(reader, b"cvParam")?;
                }
                b"userParam" => {
                    sample.user_params.push(parse_user_param(&e));
                    skip_element(reader, b"userParam")?;
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::End(e) if e.name().as_ref() == b"sample" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(sample)
}

/// <instrumentList> / <instrumentConfigurationList>
fn parse_instrument_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Option<InstrumentList>, String> {
    let mut list = InstrumentList::default();
    list.count = get_attr_usize(start, b"count");
    let end_tag = start.name();

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => match e.name().as_ref() {
                b"instrument" | b"instrumentConfiguration" => {
                    list.instrument.push(parse_instrument(reader, &e, true)?)
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::Empty(e) => match e.name().as_ref() {
                b"instrument" | b"instrumentConfiguration" => {
                    list.instrument.push(parse_instrument(reader, &e, false)?)
                }
                _ => {}
            },
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(if list.instrument.is_empty() && list.count.is_none() {
        None
    } else {
        Some(list)
    })
}

/// <instrument> / <instrumentConfiguration>
fn parse_instrument<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
    has_body: bool,
) -> Result<Instrument, String> {
    let scan_settings_ref =
        get_attr(start, b"scanSettingsRef").map(|r| ScanSettingsRef { r#ref: r });
    let software_ref = get_attr(start, b"softwareRef").map(|r| InstrumentSoftwareRef { r#ref: r });

    let mut instrument = Instrument {
        id: get_attr(start, b"id").unwrap_or_default(),
        scan_settings_ref,
        software_ref,
        ..Default::default()
    };

    if !has_body {
        return Ok(instrument);
    }

    let end_tag = start.name();
    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                if e.name().as_ref() == b"componentList" {
                    instrument.component_list = Some(ComponentList {
                        count: get_attr_usize(&e, b"count"),
                        source: Vec::new(),
                        analyzer: Vec::new(),
                        detector: Vec::new(),
                    });
                } else {
                    push_params_empty(
                        &e,
                        &mut instrument.referenceable_param_group_ref,
                        &mut instrument.cv_param,
                        &mut instrument.user_param,
                    );
                }
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut instrument.referenceable_param_group_ref,
                    &mut instrument.cv_param,
                    &mut instrument.user_param,
                )? {
                    if e.name().as_ref() == b"componentList" {
                        instrument.component_list = Some(parse_component_list(reader, &e)?);
                    } else {
                        skip_element(reader, e.name().as_ref())?;
                    }
                }
            }
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(instrument)
}

/// <componentList>
fn parse_component_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ComponentList, String> {
    let mut list = ComponentList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => match e.name().as_ref() {
                b"source" => {
                    let s = parse_component(reader, &e)?;
                    list.source.push(Source {
                        order: s.order,
                        referenceable_param_group_ref: s.referenceable_param_group_ref,
                        cv_param: s.cv_param,
                        user_param: s.user_param,
                    });
                }
                b"analyzer" => {
                    let a = parse_component(reader, &e)?;
                    list.analyzer.push(Analyzer {
                        order: a.order,
                        referenceable_param_group_ref: a.referenceable_param_group_ref,
                        cv_param: a.cv_param,
                        user_param: a.user_param,
                    });
                }
                b"detector" => {
                    let d = parse_component(reader, &e)?;
                    list.detector.push(Detector {
                        order: d.order,
                        referenceable_param_group_ref: d.referenceable_param_group_ref,
                        cv_param: d.cv_param,
                        user_param: d.user_param,
                    });
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::Empty(e) => match e.name().as_ref() {
                b"source" => list.source.push(Source {
                    order: get_attr_u32(&e, b"order"),
                    ..Default::default()
                }),
                b"analyzer" => list.analyzer.push(Analyzer {
                    order: get_attr_u32(&e, b"order"),
                    ..Default::default()
                }),
                b"detector" => list.detector.push(Detector {
                    order: get_attr_u32(&e, b"order"),
                    ..Default::default()
                }),
                _ => {}
            },
            Event::End(e) if e.name().as_ref() == b"componentList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

#[derive(Debug, Clone)]
struct ComponentTmp {
    order: Option<u32>,
    referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    cv_param: Vec<CvParam>,
    user_param: Vec<UserParam>,
}

fn parse_component<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ComponentTmp, String> {
    let mut tmp = ComponentTmp {
        order: get_attr_u32(start, b"order"),
        referenceable_param_group_ref: Vec::new(),
        cv_param: Vec::new(),
        user_param: Vec::new(),
    };

    let end_tag = start.name();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut tmp.referenceable_param_group_ref,
                    &mut tmp.cv_param,
                    &mut tmp.user_param,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut tmp.referenceable_param_group_ref,
                    &mut tmp.cv_param,
                    &mut tmp.user_param,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(tmp)
}

/// <scanSettingsList> / <acquisitionSettingsList>
fn parse_scan_settings_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Option<ScanSettingsList>, String> {
    let mut list = ScanSettingsList::default();
    list.count = get_attr_usize(start, b"count");
    let end_tag = start.name();

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => match e.name().as_ref() {
                b"scanSettings" | b"acquisitionSettings" => {
                    list.scan_settings.push(parse_scan_settings(reader, &e)?)
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::Empty(e) => match e.name().as_ref() {
                b"scanSettings" | b"acquisitionSettings" => list.scan_settings.push(ScanSettings {
                    id: get_attr(&e, b"id"),
                    instrument_configuration_ref: get_attr(&e, b"instrumentConfigurationRef"),
                    ..Default::default()
                }),
                _ => {}
            },
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(if list.scan_settings.is_empty() && list.count.is_none() {
        None
    } else {
        Some(list)
    })
}

/// <scanSettings>
fn parse_scan_settings<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ScanSettings, String> {
    let mut s = ScanSettings {
        id: get_attr(start, b"id"),
        instrument_configuration_ref: get_attr(start, b"instrumentConfigurationRef"),
        ..Default::default()
    };

    let end_tag = start.name();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut s.referenceable_param_group_refs,
                    &mut s.cv_params,
                    &mut s.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut s.referenceable_param_group_refs,
                    &mut s.cv_params,
                    &mut s.user_params,
                )? {
                    match e.name().as_ref() {
                        b"sourceFileRefList" => {
                            s.source_file_ref_list = Some(parse_source_file_ref_list(reader, &e)?)
                        }
                        b"targetList" => s.target_list = Some(parse_target_list(reader, &e)?),
                        _ => skip_element(reader, e.name().as_ref())?,
                    }
                }
            }
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(s)
}

/// <sourceFileRefList>
fn parse_source_file_ref_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SourceFileRefList, String> {
    let mut list = SourceFileRefList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) if e.name().as_ref() == b"sourceFileRef" => {
                list.source_file_refs.push(parse_source_file_ref(&e))
            }
            Event::Start(e) => {
                if e.name().as_ref() == b"sourceFileRef" {
                    list.source_file_refs.push(parse_source_file_ref(&e));
                    skip_element(reader, b"sourceFileRef")?;
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"sourceFileRefList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <targetList>
fn parse_target_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<TargetList, String> {
    let mut list = TargetList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"target" {
                    list.targets.push(parse_target(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"target" => {
                list.targets.push(Target::default())
            }
            Event::End(e) if e.name().as_ref() == b"targetList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <target>
fn parse_target<R: BufRead>(reader: &mut Reader<R>, _start: &BytesStart) -> Result<Target, String> {
    let mut target = Target::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut target.referenceable_param_group_refs,
                    &mut target.cv_params,
                    &mut target.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut target.referenceable_param_group_refs,
                    &mut target.cv_params,
                    &mut target.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"target" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(target)
}

/// <softwareList>
fn parse_software_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SoftwareList, String> {
    let mut list = SoftwareList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) if e.name().as_ref() == b"software" => list.software.push(Software {
                id: get_attr(&e, b"id").unwrap_or_default(),
                version: get_attr(&e, b"version"),
                ..Default::default()
            }),
            Event::Start(e) => {
                if e.name().as_ref() == b"software" {
                    list.software.push(parse_software(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"softwareList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <software>
fn parse_software<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Software, String> {
    let mut s = Software {
        id: get_attr(start, b"id").unwrap_or_default(),
        version: get_attr(start, b"version"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => match local_name(e.name().as_ref()) {
                b"softwareParam" => s.software_param.push(parse_software_param(&e)),
                b"cvParam" => s.cv_param.push(parse_cv_param(&e)),
                _ => {}
            },
            Event::Start(e) => match local_name(e.name().as_ref()) {
                b"softwareParam" => {
                    s.software_param.push(parse_software_param(&e));
                    skip_element(reader, e.name().as_ref())?;
                }
                b"cvParam" => {
                    s.cv_param.push(parse_cv_param(&e));
                    skip_element(reader, e.name().as_ref())?;
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::End(e) if local_name(e.name().as_ref()) == b"software" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(s)
}

/// <dataProcessingList>
fn parse_data_processing_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<DataProcessingList, String> {
    let mut list = DataProcessingList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"dataProcessing" {
                    list.data_processing
                        .push(parse_data_processing(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"dataProcessing" => {
                list.data_processing.push(DataProcessing {
                    id: get_attr(&e, b"id").unwrap_or_default(),
                    ..Default::default()
                })
            }
            Event::End(e) if e.name().as_ref() == b"dataProcessingList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(list)
}

/// <dataProcessing>
fn parse_data_processing<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<DataProcessing, String> {
    let mut dp = DataProcessing {
        id: get_attr(start, b"id").unwrap_or_default(),
        software_ref: get_attr(start, b"softwareRef"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"processingMethod" {
                    dp.processing_method
                        .push(parse_processing_method(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"processingMethod" => {
                dp.processing_method.push(ProcessingMethod {
                    order: get_attr_u32(&e, b"order"),
                    software_ref: get_attr(&e, b"softwareRef"),
                    ..Default::default()
                })
            }
            Event::End(e) if e.name().as_ref() == b"dataProcessing" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(dp)
}

/// <processingMethod>
fn parse_processing_method<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ProcessingMethod, String> {
    let mut pm = ProcessingMethod {
        order: get_attr_u32(start, b"order"),
        software_ref: get_attr(start, b"softwareRef"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut pm.referenceable_param_group_ref,
                    &mut pm.cv_param,
                    &mut pm.user_param,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut pm.referenceable_param_group_ref,
                    &mut pm.cv_param,
                    &mut pm.user_param,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"processingMethod" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(pm)
}

/// <run>
fn parse_run<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
    slim: bool,
) -> Result<Run, String> {
    let mut run = Run {
        id: get_attr(start, b"id").unwrap_or_default(),
        start_time_stamp: get_attr(start, b"startTimeStamp"),
        default_instrument_configuration_ref: get_attr(start, b"defaultInstrumentConfigurationRef")
            .or_else(|| get_attr(start, b"instrumentRef")),
        default_source_file_ref: get_attr(start, b"defaultSourceFileRef"),
        sample_ref: get_attr(start, b"sampleRef"),
        ..Default::default()
    };

    let end_tag = start.name();
    let mut buf = Vec::with_capacity(1024);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut run.referenceable_param_group_refs,
                    &mut run.cv_params,
                    &mut run.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut run.referenceable_param_group_refs,
                    &mut run.cv_params,
                    &mut run.user_params,
                )? {
                    match e.name().as_ref() {
                        b"sourceFileRefList" => {
                            run.source_file_ref_list = Some(parse_source_file_ref_list(reader, &e)?)
                        }
                        b"spectrumList" => {
                            if slim {
                                skip_element(reader, b"spectrumList")?;
                            } else {
                                run.spectrum_list = Some(parse_spectrum_list(reader, &e)?);
                            }
                        }
                        b"chromatogramList" => {
                            if slim {
                                skip_element(reader, b"chromatogramList")?;
                            } else {
                                run.chromatogram_list = Some(parse_chromatogram_list(reader, &e)?);
                            }
                        }
                        _ => skip_element(reader, e.name().as_ref())?,
                    }
                }
            }
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(run)
}

/// <spectrumList>
fn parse_spectrum_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SpectrumList, String> {
    let mut list = SpectrumList::default();
    list.count = get_attr_usize(start, b"count");
    list.default_data_processing_ref = get_attr(start, b"defaultDataProcessingRef");

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"spectrum" {
                    list.spectra.push(parse_spectrum(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"spectrum" => list.spectra.push(Spectrum {
                id: get_attr(&e, b"id").unwrap_or_default(),
                index: get_attr_u32(&e, b"index"),
                scan_number: get_attr_u32(&e, b"scanNumber"),
                ms_level: get_attr_u32(&e, b"msLevel"),
                default_array_length: get_attr_usize(&e, b"defaultArrayLength"),
                native_id: get_attr(&e, b"nativeID"),
                data_processing_ref: get_attr(&e, b"dataProcessingRef"),
                source_file_ref: get_attr(&e, b"sourceFileRef"),
                spot_id: get_attr(&e, b"spotID"),
                ..Default::default()
            }),
            Event::End(e) if e.name().as_ref() == b"spectrumList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <spectrum>
fn parse_spectrum<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Spectrum, String> {
    let mut spectrum = Spectrum {
        id: get_attr(start, b"id").unwrap_or_default(),
        index: get_attr_u32(start, b"index"),
        scan_number: get_attr_u32(start, b"scanNumber"),
        ms_level: get_attr_u32(start, b"msLevel"),
        default_array_length: get_attr_usize(start, b"defaultArrayLength"),
        native_id: get_attr(start, b"nativeID"),
        data_processing_ref: get_attr(start, b"dataProcessingRef"),
        source_file_ref: get_attr(start, b"sourceFileRef"),
        spot_id: get_attr(start, b"spotID"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(2048);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => match e.name().as_ref() {
                b"referenceableParamGroupRef" => spectrum
                    .referenceable_param_group_refs
                    .push(parse_referenceable_param_group_ref(&e)),
                b"cvParam" => {
                    let p = parse_cv_param(&e);
                    maybe_set_ms_level(&mut spectrum, &p);
                    spectrum.cv_params.push(p);
                }
                b"userParam" => spectrum.user_params.push(parse_user_param(&e)),
                _ => {}
            },
            Event::Start(e) => match e.name().as_ref() {
                b"referenceableParamGroupRef" => {
                    spectrum
                        .referenceable_param_group_refs
                        .push(parse_referenceable_param_group_ref(&e));
                    skip_element(reader, b"referenceableParamGroupRef")?;
                }
                b"cvParam" => {
                    let p = parse_cv_param(&e);
                    maybe_set_ms_level(&mut spectrum, &p);
                    spectrum.cv_params.push(p);
                    skip_element(reader, b"cvParam")?;
                }
                b"userParam" => {
                    spectrum.user_params.push(parse_user_param(&e));
                    skip_element(reader, b"userParam")?;
                }
                b"spectrumDescription" => {
                    spectrum.spectrum_description = Some(parse_spectrum_description(reader, &e)?)
                }
                b"scanList" => spectrum.scan_list = Some(parse_scan_list(reader, &e)?),
                b"precursorList" => {
                    spectrum.precursor_list = Some(parse_precursor_list(reader, &e)?)
                }
                b"productList" => spectrum.product_list = Some(parse_product_list(reader, &e)?),
                b"binaryDataArrayList" => {
                    spectrum.binary_data_array_list =
                        Some(parse_binary_data_array_list(reader, &e)?)
                }
                b"binaryDataArray" => {
                    spectrum
                        .binary_data_array_list
                        .get_or_insert_with(Default::default)
                        .binary_data_arrays
                        .push(parse_binary_data_array(reader, &e)?);
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::End(e) if e.name().as_ref() == b"spectrum" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(spectrum)
}

/// <spectrumDescription>
fn parse_spectrum_description<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<SpectrumDescription, String> {
    let mut sd = SpectrumDescription::default();
    let mut buf = Vec::with_capacity(1024);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut sd.referenceable_param_group_refs,
                    &mut sd.cv_params,
                    &mut sd.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut sd.referenceable_param_group_refs,
                    &mut sd.cv_params,
                    &mut sd.user_params,
                )? {
                    match e.name().as_ref() {
                        b"scanList" => sd.scan_list = Some(parse_scan_list(reader, &e)?),
                        b"scan" => {
                            let scan = parse_scan(reader, &e)?;
                            let list = sd.scan_list.get_or_insert_with(Default::default);
                            list.scans.push(scan);
                            list.count = Some(list.scans.len());
                        }
                        b"precursorList" => {
                            sd.precursor_list = Some(parse_precursor_list(reader, &e)?)
                        }
                        b"productList" => sd.product_list = Some(parse_product_list(reader, &e)?),
                        _ => skip_element(reader, e.name().as_ref())?,
                    }
                }
            }
            Event::End(e) if e.name().as_ref() == b"spectrumDescription" => break,
            Event::Eof => {
                return Err("unexpected EOF while parsing <spectrumDescription>".to_string());
            }
            _ => {}
        }
        buf.clear();
    }

    Ok(sd)
}

/// <scanList>
fn parse_scan_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ScanList, String> {
    let mut list = ScanList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"scan" {
                    list.scans.push(parse_scan(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"scan" => list.scans.push(Scan {
                instrument_configuration_ref: get_attr(&e, b"instrumentConfigurationRef")
                    .or_else(|| get_attr(&e, b"instrumentRef")),
                external_spectrum_id: get_attr(&e, b"externalSpectrumID"),
                source_file_ref: get_attr(&e, b"sourceFileRef"),
                spectrum_ref: get_attr(&e, b"spectrumRef"),
                ..Default::default()
            }),
            Event::End(e) if e.name().as_ref() == b"scanList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <scan>
fn parse_scan<R: BufRead>(reader: &mut Reader<R>, start: &BytesStart) -> Result<Scan, String> {
    let mut scan = Scan {
        instrument_configuration_ref: get_attr(start, b"instrumentConfigurationRef")
            .or_else(|| get_attr(start, b"instrumentRef")),
        external_spectrum_id: get_attr(start, b"externalSpectrumID"),
        source_file_ref: get_attr(start, b"sourceFileRef"),
        spectrum_ref: get_attr(start, b"spectrumRef"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut scan.referenceable_param_group_refs,
                    &mut scan.cv_params,
                    &mut scan.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut scan.referenceable_param_group_refs,
                    &mut scan.cv_params,
                    &mut scan.user_params,
                )? {
                    match e.name().as_ref() {
                        b"scanWindowList" | b"selectionWindowList" => {
                            scan.scan_window_list = Some(parse_scan_window_list(reader, &e)?)
                        }
                        _ => skip_element(reader, e.name().as_ref())?,
                    }
                }
            }
            Event::End(e) if e.name().as_ref() == b"scan" => break,
            Event::Eof => return Err("unexpected EOF while parsing <scan>".to_string()),
            _ => {}
        }
        buf.clear();
    }

    Ok(scan)
}

/// <scanWindowList> / <selectionWindowList>
fn parse_scan_window_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ScanWindowList, String> {
    let mut list = ScanWindowList::default();
    list.count = get_attr_usize(start, b"count");

    let end_tag = start.name();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                let qname = e.name();
                let n = local_name(qname.as_ref());
                if n == b"scanWindow" || n == b"selectionWindow" {
                    list.scan_windows.push(parse_scan_window(reader, &e)?);
                } else {
                    skip_element(reader, n)?;
                }
            }
            Event::Empty(e) => {
                let qname = e.name();
                let n = local_name(qname.as_ref());
                if n == b"scanWindow" || n == b"selectionWindow" {
                    list.scan_windows.push(ScanWindow::default());
                }
            }
            Event::End(e) if e.name() == end_tag => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <scanWindow> / <selectionWindow>
fn parse_scan_window<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<ScanWindow, String> {
    let mut w = ScanWindow::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => match e.name().as_ref() {
                b"cvParam" => w.cv_params.push(parse_cv_param(&e)),
                b"userParam" => w.user_params.push(parse_user_param(&e)),
                _ => {}
            },
            Event::Start(e) => match e.name().as_ref() {
                b"cvParam" => {
                    w.cv_params.push(parse_cv_param(&e));
                    skip_element(reader, b"cvParam")?;
                }
                b"userParam" => {
                    w.user_params.push(parse_user_param(&e));
                    skip_element(reader, b"userParam")?;
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::End(e)
                if {
                    let qname = e.name();
                    let n = local_name(qname.as_ref());
                    n == b"scanWindow" || n == b"selectionWindow"
                } =>
            {
                break;
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(w)
}

/// <precursorList>
fn parse_precursor_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<PrecursorList, String> {
    let mut list = PrecursorList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"precursor" {
                    list.precursors.push(parse_precursor(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"precursor" => {
                list.precursors.push(Precursor {
                    spectrum_ref: get_attr(&e, b"spectrumRef"),
                    source_file_ref: get_attr(&e, b"sourceFileRef"),
                    external_spectrum_id: get_attr(&e, b"externalSpectrumID"),
                    ..Default::default()
                })
            }
            Event::End(e) if e.name().as_ref() == b"precursorList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(list)
}

/// <ionSelection>
fn parse_ion_selection<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<SelectedIon, String> {
    let mut s = SelectedIon::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut s.referenceable_param_group_refs,
                    &mut s.cv_params,
                    &mut s.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut s.referenceable_param_group_refs,
                    &mut s.cv_params,
                    &mut s.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if local_name(e.name().as_ref()) == b"ionSelection" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(s)
}

/// <precursor>
fn parse_precursor<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Precursor, String> {
    let mut p = Precursor {
        spectrum_ref: get_attr(start, b"spectrumRef"),
        source_file_ref: get_attr(start, b"sourceFileRef"),
        external_spectrum_id: get_attr(start, b"externalSpectrumID"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                let qname = e.name();
                let raw = qname.as_ref();
                match local_name(raw) {
                    b"isolationWindow" => {
                        p.isolation_window = Some(parse_isolation_window(reader, &e)?)
                    }
                    b"selectedIonList" => {
                        p.selected_ion_list = Some(parse_selected_ion_list(reader, &e)?)
                    }
                    b"ionSelection" => {
                        let ion = parse_ion_selection(reader, &e)?;
                        let list = p.selected_ion_list.get_or_insert_with(Default::default);
                        list.selected_ions.push(ion);
                        list.count = Some(list.selected_ions.len());
                    }
                    b"activation" => p.activation = Some(parse_activation(reader, &e)?),
                    b"referenceableParamGroupRef" | b"cvParam" | b"userParam" => {
                        let w = p.isolation_window.get_or_insert_with(Default::default);
                        if !push_params_start(
                            reader,
                            &e,
                            &mut w.referenceable_param_group_refs,
                            &mut w.cv_params,
                            &mut w.user_params,
                        )? {
                            skip_element(reader, raw)?;
                        }
                    }
                    _ => skip_element(reader, raw)?,
                }
            }
            Event::Empty(e) => {
                let qname = e.name();
                let raw = qname.as_ref();
                match local_name(raw) {
                    b"ionSelection" => {
                        let list = p.selected_ion_list.get_or_insert_with(Default::default);
                        list.selected_ions.push(SelectedIon::default());
                        list.count = Some(list.selected_ions.len());
                    }
                    b"isolationWindow" => p.isolation_window = Some(IsolationWindow::default()),
                    b"selectedIonList" => p.selected_ion_list = Some(SelectedIonList::default()),
                    b"referenceableParamGroupRef" | b"cvParam" | b"userParam" => {
                        let w = p.isolation_window.get_or_insert_with(Default::default);
                        push_params_empty(
                            &e,
                            &mut w.referenceable_param_group_refs,
                            &mut w.cv_params,
                            &mut w.user_params,
                        );
                    }
                    _ => {}
                }
            }
            Event::End(e) if local_name(e.name().as_ref()) == b"precursor" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(p)
}

/// <isolationWindow>
fn parse_isolation_window<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<IsolationWindow, String> {
    let mut w = IsolationWindow::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut w.referenceable_param_group_refs,
                    &mut w.cv_params,
                    &mut w.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut w.referenceable_param_group_refs,
                    &mut w.cv_params,
                    &mut w.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if local_name(e.name().as_ref()) == b"isolationWindow" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(w)
}

/// <selectedIonList>
fn parse_selected_ion_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<SelectedIonList, String> {
    let mut list = SelectedIonList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"selectedIon" {
                    list.selected_ions.push(parse_selected_ion(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"selectedIon" => {
                list.selected_ions.push(SelectedIon::default())
            }
            Event::End(e) if e.name().as_ref() == b"selectedIonList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <selectedIon>
fn parse_selected_ion<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<SelectedIon, String> {
    let mut s = SelectedIon::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut s.referenceable_param_group_refs,
                    &mut s.cv_params,
                    &mut s.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut s.referenceable_param_group_refs,
                    &mut s.cv_params,
                    &mut s.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if e.name().as_ref() == b"selectedIon" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(s)
}

/// <activation>
fn parse_activation<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<Activation, String> {
    let mut a = Activation::default();
    let mut buf = Vec::with_capacity(512);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut a.referenceable_param_group_refs,
                    &mut a.cv_params,
                    &mut a.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut a.referenceable_param_group_refs,
                    &mut a.cv_params,
                    &mut a.user_params,
                )? {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if local_name(e.name().as_ref()) == b"activation" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(a)
}

/// <productList>
fn parse_product_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ProductList, String> {
    let mut list = ProductList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"product" {
                    list.products.push(parse_product(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"product" => list.products.push(Product {
                spectrum_ref: get_attr(&e, b"spectrumRef"),
                source_file_ref: get_attr(&e, b"sourceFileRef"),
                external_spectrum_id: get_attr(&e, b"externalSpectrumID"),
                ..Default::default()
            }),
            Event::End(e) if e.name().as_ref() == b"productList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <product>
fn parse_product<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Product, String> {
    let mut p = Product {
        spectrum_ref: get_attr(start, b"spectrumRef"),
        source_file_ref: get_attr(start, b"sourceFileRef"),
        external_spectrum_id: get_attr(start, b"externalSpectrumID"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if local_name(e.name().as_ref()) == b"isolationWindow" {
                    p.isolation_window = Some(parse_isolation_window(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::End(e) if local_name(e.name().as_ref()) == b"product" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(p)
}

/// <binaryDataArrayList>
fn parse_binary_data_array_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<BinaryDataArrayList, String> {
    let mut list = BinaryDataArrayList::default();
    list.count = get_attr_usize(start, b"count");

    let mut buf = Vec::with_capacity(512);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"binaryDataArray" {
                    list.binary_data_arrays
                        .push(parse_binary_data_array(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"binaryDataArray" => {
                list.binary_data_arrays.push(BinaryDataArray {
                    array_length: get_attr_usize(&e, b"arrayLength"),
                    encoded_length: get_attr_usize(&e, b"encodedLength"),
                    data_processing_ref: get_attr(&e, b"dataProcessingRef"),
                    ..Default::default()
                })
            }
            Event::End(e) if e.name().as_ref() == b"binaryDataArrayList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <binaryDataArray>
fn parse_binary_data_array<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<BinaryDataArray, String> {
    let mut a = BinaryDataArray {
        array_length: get_attr_usize(start, b"arrayLength"),
        encoded_length: get_attr_usize(start, b"encodedLength"),
        data_processing_ref: get_attr(start, b"dataProcessingRef"),
        ..Default::default()
    };

    let mut binary_b64: Option<String> = None;
    let mut buf = Vec::with_capacity(1024);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => {
                push_params_empty(
                    &e,
                    &mut a.referenceable_param_group_refs,
                    &mut a.cv_params,
                    &mut a.user_params,
                );
            }
            Event::Start(e) => {
                if !push_params_start(
                    reader,
                    &e,
                    &mut a.referenceable_param_group_refs,
                    &mut a.cv_params,
                    &mut a.user_params,
                )? {
                    if e.name().as_ref() == b"binary" {
                        binary_b64 = Some(read_text_content(reader, b"binary")?);
                    } else {
                        skip_element(reader, e.name().as_ref())?;
                    }
                }
            }
            Event::End(e) if e.name().as_ref() == b"binaryDataArray" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    let flags = binary_array_flags(&a);
    a.is_f32 = Some(flags.is_f32);
    a.is_f64 = Some(flags.is_f64);

    let Some(b64) = binary_b64.as_deref() else {
        return Ok(a);
    };
    let b64 = b64.trim();
    if b64.is_empty() {
        return Ok(a);
    }

    let mut bytes = Vec::new();
    if let Some(n) = a.encoded_length {
        bytes.reserve(n.saturating_mul(3) / 4 + 8);
    } else {
        bytes.reserve(b64.len().saturating_mul(3) / 4 + 8);
    }

    STANDARD
        .decode_vec(b64.as_bytes(), &mut bytes)
        .map_err(|e| format!("base64 decode failed: {e}"))?;

    if flags.is_zlib {
        bytes =
            decompress_to_vec_zlib(&bytes).map_err(|e| format!("zlib decompress failed: {e:?}"))?;
    }

    let want = match a.array_length {
        Some(n) => n,
        None => {
            if flags.is_f64 {
                bytes.len() / 8
            } else if flags.is_f32 {
                bytes.len() / 4
            } else {
                0
            }
        }
    };

    if flags.is_f64 {
        a.decoded_binary_f64.clear();
        decode_f64_into(&bytes, true, want, &mut a.decoded_binary_f64);
    } else if flags.is_f32 {
        a.decoded_binary_f32.clear();
        decode_f32_into(&bytes, true, want, &mut a.decoded_binary_f32);
    }

    Ok(a)
}

#[derive(Debug, Clone, Copy)]
struct BinaryArrayFlags {
    is_zlib: bool,
    is_f64: bool,
    is_f32: bool,
}

fn has_acc(cv_params: &[CvParam], acc: &str) -> bool {
    cv_params
        .iter()
        .any(|p| p.accession.as_deref() == Some(acc))
}

fn binary_array_flags(binary_data_array: &BinaryDataArray) -> BinaryArrayFlags {
    let cv = &binary_data_array.cv_params;

    let is_zlib = has_acc(cv, "MS:1000574");
    let is_f64 = has_acc(cv, "MS:1000523");
    let is_f32 = has_acc(cv, "MS:1000521");

    let (is_f64, is_f32) = match (is_f64, is_f32) {
        (true, false) => (true, false),
        (false, true) => (false, true),
        _ => (true, false),
    };

    BinaryArrayFlags {
        is_zlib,
        is_f64,
        is_f32,
    }
}

fn decode_f64_into(bytes: &[u8], little: bool, want: usize, out: &mut Vec<f64>) {
    let n = want.min(bytes.len() / 8);
    out.reserve(n);
    let slice = &bytes[..n * 8];

    if little {
        for c in slice.chunks_exact(8) {
            out.push(f64::from_bits(u64::from_le_bytes([
                c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7],
            ])));
        }
    } else {
        for c in slice.chunks_exact(8) {
            out.push(f64::from_bits(u64::from_be_bytes([
                c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7],
            ])));
        }
    }
}

fn decode_f32_into(bytes: &[u8], little: bool, want: usize, out: &mut Vec<f32>) {
    let n = want.min(bytes.len() / 4);
    out.reserve(n);
    let slice = &bytes[..n * 4];

    if little {
        for c in slice.chunks_exact(4) {
            out.push(f32::from_bits(u32::from_le_bytes([c[0], c[1], c[2], c[3]])));
        }
    } else {
        for c in slice.chunks_exact(4) {
            out.push(f32::from_bits(u32::from_be_bytes([c[0], c[1], c[2], c[3]])));
        }
    }
}

/// <chromatogramList>
fn parse_chromatogram_list<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<ChromatogramList, String> {
    let mut list = ChromatogramList::default();
    list.count = get_attr_usize(start, b"count");
    list.default_data_processing_ref = get_attr(start, b"defaultDataProcessingRef");

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => {
                if e.name().as_ref() == b"chromatogram" {
                    list.chromatograms.push(parse_chromatogram(reader, &e)?);
                } else {
                    skip_element(reader, e.name().as_ref())?;
                }
            }
            Event::Empty(e) if e.name().as_ref() == b"chromatogram" => {
                list.chromatograms.push(Chromatogram {
                    id: get_attr(&e, b"id").unwrap_or_default(),
                    native_id: get_attr(&e, b"nativeID"),
                    index: get_attr_u32(&e, b"index"),
                    default_array_length: get_attr_usize(&e, b"defaultArrayLength"),
                    data_processing_ref: get_attr(&e, b"dataProcessingRef"),
                    ..Default::default()
                })
            }
            Event::End(e) if e.name().as_ref() == b"chromatogramList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(list)
}

/// <chromatogram>
fn parse_chromatogram<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<Chromatogram, String> {
    let mut c = Chromatogram {
        id: get_attr(start, b"id").unwrap_or_default(),
        native_id: get_attr(start, b"nativeID"),
        index: get_attr_u32(start, b"index"),
        default_array_length: get_attr_usize(start, b"defaultArrayLength"),
        data_processing_ref: get_attr(start, b"dataProcessingRef"),
        ..Default::default()
    };

    let mut buf = Vec::with_capacity(1024);
    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Empty(e) => match e.name().as_ref() {
                b"referenceableParamGroupRef" => c
                    .referenceable_param_group_refs
                    .push(parse_referenceable_param_group_ref(&e)),
                b"cvParam" => c.cv_params.push(parse_cv_param(&e)),
                b"userParam" => c.user_params.push(parse_user_param(&e)),
                b"binaryDataArray" => {
                    let a = parse_binary_data_array(reader, &e)?;
                    c.binary_data_array_list
                        .get_or_insert_with(Default::default)
                        .binary_data_arrays
                        .push(a);
                }
                _ => {}
            },
            Event::Start(e) => match e.name().as_ref() {
                b"referenceableParamGroupRef" => {
                    c.referenceable_param_group_refs
                        .push(parse_referenceable_param_group_ref(&e));
                    skip_element(reader, b"referenceableParamGroupRef")?;
                }
                b"cvParam" => {
                    c.cv_params.push(parse_cv_param(&e));
                    skip_element(reader, b"cvParam")?;
                }
                b"userParam" => {
                    c.user_params.push(parse_user_param(&e));
                    skip_element(reader, b"userParam")?;
                }
                b"precursor" => c.precursor = Some(parse_precursor(reader, &e)?),
                b"product" => c.product = Some(parse_product(reader, &e)?),
                b"binaryDataArrayList" => {
                    c.binary_data_array_list = Some(parse_binary_data_array_list(reader, &e)?)
                }
                b"binaryDataArray" => {
                    let a = parse_binary_data_array(reader, &e)?;
                    c.binary_data_array_list
                        .get_or_insert_with(Default::default)
                        .binary_data_arrays
                        .push(a);
                }
                _ => skip_element(reader, e.name().as_ref())?,
            },
            Event::End(e) if e.name().as_ref() == b"chromatogram" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(c)
}

/// <indexList> / <indexListOffset> / <fileChecksum>
pub fn parse_index_list(bytes: &[u8]) -> Result<Option<IndexList>, String> {
    let mut reader = Reader::from_reader(Cursor::new(bytes));
    reader.config_mut().trim_text(true);

    let mut out = IndexList {
        spectrum: Vec::new(),
        chromatogram: Vec::new(),
        index_list_offset: None,
        file_checksum: None,
    };

    let mut saw_any = false;
    let mut buf = Vec::with_capacity(1024);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) => match e.name().as_ref() {
                b"indexList" => {
                    let (s, c) = parse_index_list_tag(&mut reader, &e)?;
                    if !s.is_empty() || !c.is_empty() {
                        out.spectrum = s;
                        out.chromatogram = c;
                        saw_any = true;
                    }
                }
                b"indexListOffset" => {
                    let t = read_text_content(&mut reader, b"indexListOffset")?;
                    out.index_list_offset = t.trim().parse::<u64>().ok();
                    if out.index_list_offset.is_some() {
                        saw_any = true;
                    }
                }
                b"fileChecksum" => {
                    out.file_checksum = Some(read_text_content(&mut reader, b"fileChecksum")?);
                    if out
                        .file_checksum
                        .as_deref()
                        .map(|s| !s.is_empty())
                        .unwrap_or(false)
                    {
                        saw_any = true;
                    }
                }
                _ => {}
            },
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(if saw_any { Some(out) } else { None })
}

/// <indexList>
fn parse_index_list_tag<R: BufRead>(
    reader: &mut Reader<R>,
    _start: &BytesStart,
) -> Result<(Vec<IndexOffset>, Vec<IndexOffset>), String> {
    let mut spectrum = Vec::new();
    let mut chromatogram = Vec::new();
    let mut buf = Vec::with_capacity(1024);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) if e.name().as_ref() == b"index" => {
                let (name, offsets) = parse_index_tag(reader, &e)?;
                if name == "spectrum" {
                    spectrum = offsets;
                } else if name == "chromatogram" {
                    chromatogram = offsets;
                }
            }
            Event::End(e) if e.name().as_ref() == b"indexList" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok((spectrum, chromatogram))
}

/// <index>
fn parse_index_tag<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<(String, Vec<IndexOffset>), String> {
    let name = get_attr(start, b"name").unwrap_or_default();
    let mut offsets = Vec::new();
    let mut buf = Vec::with_capacity(1024);

    loop {
        match reader
            .read_event_into(&mut buf)
            .map_err(|e| e.to_string())?
        {
            Event::Start(e) if e.name().as_ref() == b"offset" => {
                offsets.push(parse_offset_tag(reader, &e)?)
            }
            Event::Empty(e) if e.name().as_ref() == b"offset" => offsets.push(IndexOffset {
                id_ref: get_attr_any(&e, &[b"idRef", b"idref"]),
                offset: 0,
            }),
            Event::End(e) if e.name().as_ref() == b"index" => break,
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }

    Ok((name, offsets))
}

/// <offset>
fn parse_offset_tag<R: BufRead>(
    reader: &mut Reader<R>,
    start: &BytesStart,
) -> Result<IndexOffset, String> {
    let id_ref = get_attr_any(start, &[b"idRef", b"idref"]);
    let t = read_text_content(reader, b"offset")?;
    let offset = t.trim().parse::<u64>().unwrap_or(0);
    Ok(IndexOffset { id_ref, offset })
}
