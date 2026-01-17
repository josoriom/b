use serde::{Deserialize, Serialize};

/// <mzML>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MzML {
    pub cv_list: Option<CvList>,
    pub file_description: FileDescription,
    pub referenceable_param_group_list: Option<ReferenceableParamGroupList>,
    pub sample_list: Option<SampleList>,
    pub instrument_list: Option<InstrumentList>,
    pub software_list: Option<SoftwareList>,
    pub data_processing_list: Option<DataProcessingList>,
    pub scan_settings_list: Option<ScanSettingsList>,
    pub run: Run,
}

/// <cvList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct CvList {
    pub count: Option<usize>,
    pub cv: Vec<Cv>,
}

/// <cv>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Cv {
    pub id: String,
    pub full_name: Option<String>,
    pub version: Option<String>,
    pub uri: Option<String>,
}

/// <cvParam>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct CvParam {
    pub cv_ref: Option<String>,
    pub accession: Option<String>,
    pub name: String,
    pub value: Option<String>,
    pub unit_cv_ref: Option<String>,
    pub unit_name: Option<String>,
    pub unit_accession: Option<String>,
}

/// <userParam>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct UserParam {
    pub name: String,
    pub r#type: Option<String>,
    pub unit_accession: Option<String>,
    pub unit_cv_ref: Option<String>,
    pub unit_name: Option<String>,
    pub value: Option<String>,
}

/// <referenceableParamGroupRef>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ReferenceableParamGroupRef {
    pub r#ref: String,
}

/// <dataProcessingList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DataProcessingList {
    pub count: Option<usize>,
    pub data_processing: Vec<DataProcessing>,
}

/// <dataProcessing>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DataProcessing {
    pub id: String,
    pub software_ref: Option<String>,
    pub processing_method: Vec<ProcessingMethod>,
}

/// <processingMethod>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ProcessingMethod {
    pub order: Option<u32>,
    pub software_ref: Option<String>,
    pub referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    pub cv_param: Vec<CvParam>,
    pub user_param: Vec<UserParam>,
}

/// <fileDescription>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FileDescription {
    pub file_content: FileContent,
    pub source_file_list: SourceFileList,
    pub contacts: Vec<Contact>,
}

/// <sourceFileList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SourceFileList {
    pub count: Option<usize>,
    pub source_file: Vec<SourceFile>,
}

/// <sourceFile>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SourceFile {
    pub id: String,
    pub name: String,
    pub location: String,
    pub referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    pub cv_param: Vec<CvParam>,
    pub user_param: Vec<UserParam>,
}

/// <fileContent>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FileContent {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <contact>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Contact {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <instrumentList> / <instrumentConfigurationList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct InstrumentList {
    pub count: Option<usize>,
    pub instrument: Vec<Instrument>,
}

/// <instrument> / <instrumentConfiguration>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Instrument {
    pub id: String,
    pub scan_settings_ref: Option<ScanSettingsRef>,
    pub cv_param: Vec<CvParam>,
    pub user_param: Vec<UserParam>,
    pub referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    pub component_list: Option<ComponentList>,
    pub software_ref: Option<InstrumentSoftwareRef>,
}

/// attribute scanSettingsRef
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScanSettingsRef {
    pub r#ref: String,
}

/// <componentList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ComponentList {
    pub count: Option<usize>,
    pub source: Vec<Source>,
    pub analyzer: Vec<Analyzer>,
    pub detector: Vec<Detector>,
}

/// <source>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Source {
    pub order: Option<u32>,
    pub referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    pub cv_param: Vec<CvParam>,
    pub user_param: Vec<UserParam>,
}

/// <analyzer>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Analyzer {
    pub order: Option<u32>,
    pub referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    pub cv_param: Vec<CvParam>,
    pub user_param: Vec<UserParam>,
}

/// <detector>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Detector {
    pub order: Option<u32>,
    pub referenceable_param_group_ref: Vec<ReferenceableParamGroupRef>,
    pub cv_param: Vec<CvParam>,
    pub user_param: Vec<UserParam>,
}

/// <instrumentSoftwareRef> / <softwareRef>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct InstrumentSoftwareRef {
    pub r#ref: String,
}

/// <referenceableParamGroupList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ReferenceableParamGroupList {
    pub count: Option<usize>,
    pub referenceable_param_groups: Vec<ReferenceableParamGroup>,
}

/// <referenceableParamGroup>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ReferenceableParamGroup {
    pub id: String,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <sampleList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SampleList {
    pub count: Option<u32>,
    pub samples: Vec<Sample>,
}

/// <sample>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Sample {
    pub id: String,
    pub name: String,
    pub referenceable_param_group_ref: Option<ReferenceableParamGroupRef>,
}

/// <scanSettingsList> / <acquisitionSettingsList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScanSettingsList {
    pub count: Option<usize>,
    pub scan_settings: Vec<ScanSettings>,
}

/// <scanSettings> / <acquisitionSettings>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScanSettings {
    pub id: Option<String>,
    pub instrument_configuration_ref: Option<String>,
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
    pub source_file_ref_list: Option<SourceFileRefList>,
    pub target_list: Option<TargetList>,
}

/// <sourceFileRefList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SourceFileRefList {
    pub count: Option<usize>,
    pub source_file_refs: Vec<SourceFileRef>,
}

/// <sourceFileRef>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SourceFileRef {
    pub r#ref: String,
}

/// <targetList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct TargetList {
    pub count: Option<usize>,
    pub targets: Vec<Target>,
}

/// <target>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Target {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <softwareList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SoftwareList {
    pub count: Option<usize>,
    pub software: Vec<Software>,
}

/// <software>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Software {
    pub id: String,
    pub version: Option<String>,
    pub software_param: Vec<SoftwareParam>,
    pub cv_param: Vec<CvParam>,
}

/// <softwareParam>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SoftwareParam {
    pub cv_ref: Option<String>,
    pub accession: String,
    pub name: String,
    pub version: Option<String>,
}

/// <run>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Run {
    pub id: String,
    pub start_time_stamp: Option<String>,
    pub default_instrument_configuration_ref: Option<String>,
    pub default_source_file_ref: Option<String>,
    pub sample_ref: Option<String>,

    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,

    pub source_file_ref_list: Option<SourceFileRefList>,
    pub spectrum_list: Option<SpectrumList>,
    pub chromatogram_list: Option<ChromatogramList>,
}

/// <spectrumList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SpectrumList {
    pub count: Option<usize>,
    pub default_data_processing_ref: Option<String>,
    pub spectra: Vec<Spectrum>,
}

/// <spectrumDescription> (1.0.0)
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SpectrumDescription {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,

    pub scan_list: Option<ScanList>,
    pub precursor_list: Option<PrecursorList>,
    pub product_list: Option<ProductList>,
}

/// <scanList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScanList {
    pub count: Option<usize>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
    pub scans: Vec<Scan>,
}

/// <scan>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Scan {
    pub instrument_configuration_ref: Option<String>,
    pub external_spectrum_id: Option<String>,
    pub source_file_ref: Option<String>,
    pub spectrum_ref: Option<String>,

    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,

    pub scan_window_list: Option<ScanWindowList>,
}

/// <scanWindowList> / <selectionWindowList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScanWindowList {
    pub count: Option<usize>,
    pub scan_windows: Vec<ScanWindow>,
}

/// <scanWindow> / <selectionWindow>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScanWindow {
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <precursorList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PrecursorList {
    pub count: Option<usize>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
    pub precursors: Vec<Precursor>,
}

/// <precursor>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Precursor {
    pub spectrum_ref: Option<String>,
    pub source_file_ref: Option<String>,
    pub external_spectrum_id: Option<String>,

    pub isolation_window: Option<IsolationWindow>,
    pub selected_ion_list: Option<SelectedIonList>,
    pub activation: Option<Activation>,
}

/// <isolationWindow>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct IsolationWindow {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <selectedIonList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SelectedIonList {
    pub count: Option<usize>,
    pub selected_ions: Vec<SelectedIon>,
}

/// <selectedIon>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SelectedIon {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <activation>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Activation {
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,
}

/// <productList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ProductList {
    pub count: Option<usize>,
    pub products: Vec<Product>,
}

/// <product>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Product {
    pub spectrum_ref: Option<String>,
    pub source_file_ref: Option<String>,
    pub external_spectrum_id: Option<String>,

    pub isolation_window: Option<IsolationWindow>,
}

/// <binaryDataArrayList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct BinaryDataArrayList {
    pub count: Option<usize>,
    pub binary_data_arrays: Vec<BinaryDataArray>,
}

/// <binaryDataArray>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct BinaryDataArray {
    pub array_length: Option<usize>,
    pub encoded_length: Option<usize>,
    pub data_processing_ref: Option<String>,

    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,

    pub is_f32: Option<bool>,
    pub is_f64: Option<bool>,

    pub decoded_binary_f32: Vec<f32>,
    pub decoded_binary_f64: Vec<f64>,
}

/// <chromatogramList>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ChromatogramList {
    pub count: Option<usize>,
    pub default_data_processing_ref: Option<String>,
    pub chromatograms: Vec<Chromatogram>,
}

/// <chromatogram>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Chromatogram {
    pub id: String,
    pub native_id: Option<String>,
    pub index: Option<u32>,
    pub default_array_length: Option<usize>,
    pub data_processing_ref: Option<String>,

    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,

    pub precursor: Option<Precursor>,
    pub product: Option<Product>,

    pub binary_data_array_list: Option<BinaryDataArrayList>,
}

/// <spectrum>
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Spectrum {
    // Attributes
    pub id: String,
    pub index: Option<u32>,
    pub scan_number: Option<u32>,
    pub default_array_length: Option<usize>,
    pub native_id: Option<String>,
    pub data_processing_ref: Option<String>,
    pub source_file_ref: Option<String>,
    pub spot_id: Option<String>,
    pub ms_level: Option<u32>,

    // Children
    pub referenceable_param_group_refs: Vec<ReferenceableParamGroupRef>,
    pub cv_params: Vec<CvParam>,
    pub user_params: Vec<UserParam>,

    /// 1.0.0: <spectrumDescription>
    pub spectrum_description: Option<SpectrumDescription>,

    /// 1.1.x: <spectrum>
    pub scan_list: Option<ScanList>,
    pub precursor_list: Option<PrecursorList>,
    pub product_list: Option<ProductList>,
    pub binary_data_array_list: Option<BinaryDataArrayList>,
}
