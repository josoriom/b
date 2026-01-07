use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::default::Default;
use std::sync::OnceLock;

static SCHEMA: OnceLock<SchemaTree> = OnceLock::new();

pub fn schema() -> &'static SchemaTree {
    SCHEMA.get_or_init(|| {
        let mut tree: SchemaTree =
            serde_json::from_str(include_str!("schema.json")).expect("schema.json (tree)");
        tree.build_index();
        tree
    })
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum TagId {
    FileContent = 0,
    SourceFile = 1,
    Contact = 2,
    ReferenceableParamGroup = 3,
    Sample = 4,

    #[serde(rename = "instrumentConfiguration")]
    Instrument = 5,

    #[serde(rename = "source")]
    ComponentSource = 6,
    #[serde(rename = "analyzer")]
    ComponentAnalyzer = 7,
    #[serde(rename = "detector")]
    ComponentDetector = 8,

    Software = 9,
    ProcessingMethod = 10,
    ScanSettings = 11,
    Target = 12,
    Run = 13,

    Spectrum = 14,
    SpectrumDescription = 15,
    Scan = 16,
    ScanWindow = 17,
    Precursor = 18,
    IsolationWindow = 19,
    SelectedIon = 20,
    Activation = 21,
    Product = 22,
    BinaryDataArray = 23,

    Chromatogram = 24,

    FileDescription = 25,

    SourceFileList = 26,
    SourceFileRef = 27,
    SourceFileRefList = 28,

    ReferenceableParamGroupList = 29,
    ReferenceableParamGroupRef = 30,

    SampleList = 31,

    InstrumentConfigurationList = 32,
    ComponentList = 33,

    SoftwareList = 34,
    SoftwareParam = 35,
    SoftwareRef = 36,

    DataProcessing = 37,
    DataProcessingList = 38,

    ScanSettingsList = 39,
    AcquisitionSettings = 40,
    AcquisitionSettingsList = 41,

    TargetList = 42,

    SpectrumList = 43,
    ScanList = 44,
    ScanWindowList = 45,

    PrecursorList = 46,
    SelectedIonList = 47,
    ProductList = 48,

    BinaryDataArrayList = 49,
    Binary = 50,

    ChromatogramList = 51,

    CvParam = 52,
    UserParam = 53,

    Unknown = 255,
}

impl TagId {
    #[inline]
    pub fn from_xml_tag(s: &str) -> TagId {
        match s {
            "fileContent" => TagId::FileContent,
            "sourceFile" => TagId::SourceFile,
            "contact" => TagId::Contact,
            "referenceableParamGroup" => TagId::ReferenceableParamGroup,
            "sample" => TagId::Sample,

            "instrumentConfiguration" => TagId::Instrument,
            "source" => TagId::ComponentSource,
            "analyzer" => TagId::ComponentAnalyzer,
            "detector" => TagId::ComponentDetector,

            "software" => TagId::Software,
            "processingMethod" => TagId::ProcessingMethod,
            "scanSettings" => TagId::ScanSettings,
            "target" => TagId::Target,
            "run" => TagId::Run,

            "spectrum" => TagId::Spectrum,
            "spectrumDescription" => TagId::SpectrumDescription,
            "scan" => TagId::Scan,
            "scanWindow" => TagId::ScanWindow,
            "precursor" => TagId::Precursor,
            "isolationWindow" => TagId::IsolationWindow,
            "selectedIon" => TagId::SelectedIon,
            "activation" => TagId::Activation,
            "product" => TagId::Product,
            "binaryDataArray" => TagId::BinaryDataArray,

            "chromatogram" => TagId::Chromatogram,
            "fileDescription" => TagId::FileDescription,

            "sourceFileList" => TagId::SourceFileList,
            "sourceFileRef" => TagId::SourceFileRef,
            "sourceFileRefList" => TagId::SourceFileRefList,

            "referenceableParamGroupList" => TagId::ReferenceableParamGroupList,
            "referenceableParamGroupRef" => TagId::ReferenceableParamGroupRef,

            "sampleList" => TagId::SampleList,

            "instrumentConfigurationList" => TagId::InstrumentConfigurationList,
            "componentList" => TagId::ComponentList,

            "softwareList" => TagId::SoftwareList,
            "softwareParam" => TagId::SoftwareParam,
            "softwareRef" => TagId::SoftwareRef,

            "dataProcessing" => TagId::DataProcessing,
            "dataProcessingList" => TagId::DataProcessingList,

            "scanSettingsList" => TagId::ScanSettingsList,
            "acquisitionSettings" => TagId::AcquisitionSettings,
            "acquisitionSettingsList" => TagId::AcquisitionSettingsList,

            "targetList" => TagId::TargetList,

            "spectrumList" => TagId::SpectrumList,
            "scanList" => TagId::ScanList,
            "scanWindowList" => TagId::ScanWindowList,

            "precursorList" => TagId::PrecursorList,
            "selectedIonList" => TagId::SelectedIonList,
            "productList" => TagId::ProductList,

            "binaryDataArrayList" => TagId::BinaryDataArrayList,
            "binary" => TagId::Binary,

            "chromatogramList" => TagId::ChromatogramList,

            "cvParam" => TagId::CvParam,
            "userParam" => TagId::UserParam,

            _ => TagId::Unknown,
        }
    }

    #[inline]
    pub fn from_u8(b: u8) -> Option<TagId> {
        if b <= TagId::UserParam as u8 || b == TagId::Unknown as u8 {
            Some(TagId::from(b))
        } else {
            None
        }
    }

    #[inline]
    pub fn as_u8(self) -> u8 {
        self as u8
    }
}

impl From<u8> for TagId {
    #[inline]
    fn from(b: u8) -> Self {
        match b {
            0 => TagId::FileContent,
            1 => TagId::SourceFile,
            2 => TagId::Contact,
            3 => TagId::ReferenceableParamGroup,
            4 => TagId::Sample,
            5 => TagId::Instrument,
            6 => TagId::ComponentSource,
            7 => TagId::ComponentAnalyzer,
            8 => TagId::ComponentDetector,
            9 => TagId::Software,
            10 => TagId::ProcessingMethod,
            11 => TagId::ScanSettings,
            12 => TagId::Target,
            13 => TagId::Run,
            14 => TagId::Spectrum,
            15 => TagId::SpectrumDescription,
            16 => TagId::Scan,
            17 => TagId::ScanWindow,
            18 => TagId::Precursor,
            19 => TagId::IsolationWindow,
            20 => TagId::SelectedIon,
            21 => TagId::Activation,
            22 => TagId::Product,
            23 => TagId::BinaryDataArray,
            24 => TagId::Chromatogram,
            25 => TagId::FileDescription,
            26 => TagId::SourceFileList,
            27 => TagId::SourceFileRef,
            28 => TagId::SourceFileRefList,
            29 => TagId::ReferenceableParamGroupList,
            30 => TagId::ReferenceableParamGroupRef,
            31 => TagId::SampleList,
            32 => TagId::InstrumentConfigurationList,
            33 => TagId::ComponentList,
            34 => TagId::SoftwareList,
            35 => TagId::SoftwareParam,
            36 => TagId::SoftwareRef,
            37 => TagId::DataProcessing,
            38 => TagId::DataProcessingList,
            39 => TagId::ScanSettingsList,
            40 => TagId::AcquisitionSettings,
            41 => TagId::AcquisitionSettingsList,
            42 => TagId::TargetList,
            43 => TagId::SpectrumList,
            44 => TagId::ScanList,
            45 => TagId::ScanWindowList,
            46 => TagId::PrecursorList,
            47 => TagId::SelectedIonList,
            48 => TagId::ProductList,
            49 => TagId::BinaryDataArrayList,
            50 => TagId::Binary,
            51 => TagId::ChromatogramList,
            52 => TagId::CvParam,
            53 => TagId::UserParam,
            255 => TagId::Unknown,
            _ => TagId::Unknown,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Use {
    Required,
    Optional,
}

impl Default for Use {
    #[inline]
    fn default() -> Self {
        Use::Optional
    }
}

fn default_child_key_by_tag() -> Vec<Option<String>> {
    vec![None; 256]
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SchemaNode {
    #[serde(rename = "self", default)]
    pub self_tags: Vec<TagId>,

    #[serde(default)]
    pub attributes: HashMap<String, Vec<String>>,

    #[serde(default)]
    pub children: HashMap<String, SchemaNode>,

    #[serde(rename = "use", default)]
    pub use_: Use,

    #[serde(default)]
    pub accessions: Vec<String>,

    #[serde(default)]
    pub attributes_use: HashMap<String, Use>,

    #[serde(skip, default = "default_child_key_by_tag")]
    pub child_key_by_tag: Vec<Option<String>>,
}

impl SchemaNode {
    #[inline]
    pub fn build_child_index(&mut self) {
        if self.child_key_by_tag.len() != 256 {
            self.child_key_by_tag = default_child_key_by_tag();
        } else {
            for x in &mut self.child_key_by_tag {
                *x = None;
            }
        }

        for (child_key, child) in &mut self.children {
            child.build_child_index();

            for &tag in &child.self_tags {
                let idx = tag as usize;
                if idx < 256 {
                    self.child_key_by_tag[idx] = Some(child_key.clone());
                }
            }
        }
    }

    #[inline]
    pub fn child_key_for_tag(&self, tag: TagId) -> Option<&str> {
        self.child_key_by_tag
            .get(tag as usize)
            .and_then(|x| x.as_deref())
    }
}

fn default_key_by_tag() -> Vec<Option<String>> {
    vec![None; 256]
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SchemaTree {
    #[serde(flatten)]
    pub roots: HashMap<String, SchemaNode>,

    #[serde(skip, default = "default_key_by_tag")]
    root_key_by_tag: Vec<Option<String>>,
}

impl SchemaTree {
    pub fn build_index(&mut self) {
        if self.root_key_by_tag.len() != 256 {
            self.root_key_by_tag = default_key_by_tag();
        } else {
            for x in &mut self.root_key_by_tag {
                *x = None;
            }
        }

        for (root_key, root_node) in self.roots.iter_mut() {
            root_node.build_child_index();

            for &tag in &root_node.self_tags {
                let idx = tag.as_u8() as usize;
                if idx < 256 {
                    self.root_key_by_tag[idx] = Some(root_key.clone());
                }
            }
        }
    }

    #[inline]
    pub fn root_key_for_tag(&self, tag: TagId) -> Option<&str> {
        self.root_key_by_tag
            .get(tag.as_u8() as usize)
            .and_then(|x| x.as_deref())
    }

    #[inline]
    pub fn root_by_key(&self, root_key: &str) -> Option<&SchemaNode> {
        self.roots.get(root_key)
    }

    #[inline]
    pub fn root_by_tag(&self, tag: TagId) -> Option<&SchemaNode> {
        self.root_key_for_tag(tag).and_then(|k| self.roots.get(k))
    }

    #[inline]
    pub fn root_by_xml_tag(&self, xml_tag: &str) -> Option<&SchemaNode> {
        self.root_by_tag(TagId::from_xml_tag(xml_tag))
    }
}
