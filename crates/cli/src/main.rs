use clap::{ArgAction, ArgGroup, Args, Parser, Subcommand};
use serde::Serialize;
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
};

use octo::{
    b64::{decode, encode},
    mzml::{bin_to_mzml::bin_to_mzml, parse_mzml::parse_mzml, structs::*},
};

const VERSION: &str = "0.0.0";

const HELP_TXT: &str = r#"octo v0.0.0

USAGE:
  octo -h | --help
  octo -v | --version

  octo convert (--mzml-to-b64 | --mzml-to-b32 | --b64-to-mzml) [--input-path DIR] [--output-path DIR] [--level 0..22]
  octo show --file-path PATH (--general | --run | --spectrum-list | --chromatogram-list | --spectrum | --chromatogram) [--items SPEC] [--binary]

CONVERT FLAGS:
  --mzml-to-b64        .mzML -> .b64
  --mzml-to-b32        .mzML -> .b32
  --b64-to-mzml        .b64/.b32 -> .mzML
  --input-path DIR     default: crates/parser/data/mzml
  --output-path DIR    default: crates/parser/data/b64
  --level 0..22        default: 12
  --overwrite          default: false (skip if output already exists)

SHOW FLAGS:
  --file-path PATH     input file (.mzML/.b64/.b32)
  --general
  --run
  --spectrum-list
  --chromatogram-list
  --spectrum
  --chromatogram
  --items SPEC         only with --spectrum/--chromatogram, default: 0-100 (END is exclusive)
  --binary             only with --spectrum/--chromatogram, include decoded arrays

EXAMPLES:
  octo convert --mzml-to-b64 --input-path crates/parser/data/mzml --output-path crates/parser/data/b64
  octo convert --b64-to-mzml --input-path crates/parser/data/b64 --output-path crates/parser/data/mzml_out

  octo show --file-path crates/parser/data/mzml/tiny.msdata.mzML0.99.9.mzML --general
  octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --run
  octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum --items 5
  octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum --items 0-10 --binary
"#;

#[derive(Parser)]
#[command(
    name = "b",
    about = "octo CLI",
    version = VERSION,
    disable_help_flag = true,
    disable_version_flag = true,
    disable_help_subcommand = true
)]
struct Cli {
    #[arg(short = 'h', long = "help", action = ArgAction::SetTrue)]
    help: bool,

    #[arg(short = 'v', long = "version", action = ArgAction::SetTrue)]
    version: bool,

    #[command(subcommand)]
    cmd: Option<Cmd>,
}

#[derive(Subcommand)]
enum Cmd {
    Convert(ConvertArgs),
    Show(ShowArgs),
}

#[derive(Args)]
#[command(
    group(
        ArgGroup::new("convert_mode")
            .args(["mzml_to_b64", "mzml_to_b32", "b64_to_mzml"])
            .required(true)
            .multiple(false)
    )
)]
struct ConvertArgs {
    #[arg(long)]
    input_path: Option<PathBuf>,

    #[arg(long)]
    output_path: Option<PathBuf>,

    #[arg(long = "level", default_value_t = 12, value_parser = clap::value_parser!(u8).range(0..=22))]
    compression_level: u8,

    #[arg(long, default_value_t = false, action = ArgAction::SetTrue)]
    overwrite: bool,

    #[command(flatten)]
    which: ConvertWhich,
}

#[derive(Args)]
struct ConvertWhich {
    #[arg(long = "mzml-to-b64")]
    mzml_to_b64: bool,

    #[arg(long = "mzml-to-b32")]
    mzml_to_b32: bool,

    #[arg(long = "b64-to-mzml")]
    b64_to_mzml: bool,
}

#[derive(Args)]
#[command(
    group(
        ArgGroup::new("items_scope")
            .args(["spectrum", "chromatogram"])
            .multiple(false)
    )
)]
struct ShowArgs {
    #[arg(long = "file-path")]
    file_path: PathBuf,

    #[command(flatten)]
    which: ShowWhich,

    #[arg(long, default_value_t = false, action = ArgAction::Set, requires = "items_scope")]
    binary: bool,

    #[arg(long, default_value = "0-100", requires = "items_scope")]
    items: String,
}

#[derive(Args)]
#[group(required = true, multiple = false)]
struct ShowWhich {
    #[arg(long)]
    general: bool,

    #[arg(long)]
    run: bool,

    #[arg(long = "spectrum-list")]
    spectrum_list: bool,

    #[arg(long)]
    spectrum: bool,

    #[arg(long = "chromatogram-list")]
    chromatogram_list: bool,

    #[arg(long)]
    chromatogram: bool,
}

enum ItemsSpec {
    One(usize),
    Range { start: usize, end_exclusive: usize },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    if cli.version {
        println!("{VERSION}");
        return Ok(());
    }

    if cli.help || cli.cmd.is_none() {
        print!("{HELP_TXT}");
        return Ok(());
    }

    match cli.cmd.unwrap() {
        Cmd::Convert(cmd) => convert(cmd).map_err(|e| e.into()),
        Cmd::Show(cmd) => show(cmd).map_err(|e| e.into()),
    }
}

fn parse_items_spec(s: &str) -> Result<ItemsSpec, String> {
    let s = s.trim();
    if let Some((a, b)) = s.split_once('-') {
        let start: usize = a.trim().parse().map_err(|_| "bad items start")?;
        let end_exclusive: usize = b.trim().parse().map_err(|_| "bad items end")?;
        if end_exclusive < start {
            return Err("items end must be >= start".to_string());
        }
        return Ok(ItemsSpec::Range {
            start,
            end_exclusive,
        });
    }
    let idx: usize = s.parse().map_err(|_| "bad items index")?;
    Ok(ItemsSpec::One(idx))
}

fn slice_indices(len: usize, spec: &ItemsSpec) -> (usize, usize, bool) {
    match *spec {
        ItemsSpec::One(i) => {
            if i >= len {
                (len, len, true)
            } else {
                (i, i + 1, true)
            }
        }
        ItemsSpec::Range {
            start,
            end_exclusive,
        } => {
            let s = start.min(len);
            let e = end_exclusive.min(len);
            (s, e, false)
        }
    }
}

fn workspace_root() -> PathBuf {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    here.ancestors()
        .nth(2)
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

fn file_ext_lower(path: &Path) -> String {
    path.extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase()
}

fn out_name_for_mzml_file(path: &Path, out_ext: &str) -> Option<String> {
    if file_ext_lower(path) != "mzml" {
        return None;
    }
    let stem = path.file_stem()?.to_string_lossy();
    Some(format!("{stem}.{out_ext}"))
}

fn out_name_for_bin_file_as_mzml(path: &Path) -> Option<String> {
    let ext = file_ext_lower(path);
    if ext != "b64" && ext != "b32" {
        return None;
    }
    let stem = path.file_stem()?.to_string_lossy();
    Some(format!("{stem}.mzML"))
}

fn read_mzml_or_b64(file_path: &Path) -> Result<MzML, String> {
    let bytes = fs::read(file_path).map_err(|e| format!("read failed: {e}"))?;
    let ext = file_ext_lower(file_path);

    if ext == "b64" || ext == "b32" {
        return decode(&bytes).map_err(|e| format!("decode failed: {e}"));
    }
    if ext == "mzml" {
        return parse_mzml(&bytes, false).map_err(|e| format!("parse_mzml failed: {e}"));
    }

    Err(format!(
        "unsupported file extension: {ext:?} (expected .mzML or .b64/.b32)"
    ))
}

fn collect_files_with_exts(input_root: &Path, exts: &[&str]) -> Result<Vec<PathBuf>, String> {
    let mut out = Vec::new();
    let mut stack = vec![input_root.to_path_buf()];

    while let Some(dir) = stack.pop() {
        let entries = fs::read_dir(&dir).map_err(|e| format!("read dir failed: {e}"))?;
        for entry in entries {
            let entry = entry.map_err(|e| format!("read dir entry failed: {e}"))?;
            let p = entry.path();
            if p.is_dir() {
                stack.push(p);
                continue;
            }
            if !p.is_file() {
                continue;
            }
            let ext = file_ext_lower(&p);
            let mut ok = false;
            for want in exts {
                if ext == *want {
                    ok = true;
                    break;
                }
            }
            if ok {
                out.push(p);
            }
        }
    }

    out.sort();
    Ok(out)
}

fn print_json<T: Serialize>(v: &T) -> Result<(), String> {
    let mut val = serde_json::to_value(v).map_err(|e| format!("json failed: {e}"))?;
    prune_json(&mut val);
    let s = serde_json::to_string_pretty(&val).map_err(|e| format!("json failed: {e}"))?;
    println!("{s}");
    Ok(())
}

#[derive(Serialize)]
struct GeneralOut<'a> {
    cv_list: &'a Option<CvList>,
    file_description: &'a FileDescription,
    referenceable_param_group_list: &'a Option<ReferenceableParamGroupList>,
    sample_list: &'a Option<SampleList>,
    instrument_list: &'a Option<InstrumentList>,
    software_list: &'a Option<SoftwareList>,
    data_processing_list: &'a Option<DataProcessingList>,
    scan_settings_list: &'a Option<ScanSettingsList>,
}

#[derive(Serialize)]
struct RunOut<'a> {
    id: &'a str,
    start_time_stamp: &'a Option<String>,
    default_instrument_configuration_ref: &'a Option<String>,
    default_source_file_ref: &'a Option<String>,
    sample_ref: &'a Option<String>,
    referenceable_param_group_refs: &'a Vec<ReferenceableParamGroupRef>,
    cv_params: &'a Vec<CvParam>,
    user_params: &'a Vec<UserParam>,
    source_file_ref_list: &'a Option<SourceFileRefList>,
}

#[derive(Serialize)]
struct SpectrumListOut<'a> {
    count: &'a Option<usize>,
    default_data_processing_ref: &'a Option<String>,
}

#[derive(Serialize)]
struct ChromatogramListOut<'a> {
    count: &'a Option<usize>,
    default_data_processing_ref: &'a Option<String>,
}

#[derive(Serialize)]
struct BinaryDataArrayOut<'a> {
    array_length: Option<usize>,
    encoded_length: Option<usize>,
    data_processing_ref: &'a Option<String>,
    referenceable_param_group_refs: &'a Vec<ReferenceableParamGroupRef>,
    cv_params: &'a Vec<CvParam>,
    user_params: &'a Vec<UserParam>,
    is_f32: Option<bool>,
    is_f64: Option<bool>,

    #[serde(skip_serializing_if = "Option::is_none")]
    decoded_binary_f32: Option<&'a [f32]>,

    #[serde(skip_serializing_if = "Option::is_none")]
    decoded_binary_f64: Option<&'a [f64]>,
}

#[derive(Serialize)]
struct BinaryDataArrayListOut<'a> {
    count: &'a Option<usize>,
    binary_data_arrays: Vec<BinaryDataArrayOut<'a>>,
}

fn bda_out<'a>(ba: &'a BinaryDataArray, binary: bool) -> BinaryDataArrayOut<'a> {
    BinaryDataArrayOut {
        array_length: ba.array_length,
        encoded_length: ba.encoded_length,
        data_processing_ref: &ba.data_processing_ref,
        referenceable_param_group_refs: &ba.referenceable_param_group_refs,
        cv_params: &ba.cv_params,
        user_params: &ba.user_params,
        is_f32: ba.is_f32,
        is_f64: ba.is_f64,
        decoded_binary_f32: if binary {
            Some(ba.decoded_binary_f32.as_slice())
        } else {
            None
        },
        decoded_binary_f64: if binary {
            Some(ba.decoded_binary_f64.as_slice())
        } else {
            None
        },
    }
}

fn bda_list_out<'a>(
    list: Option<&'a BinaryDataArrayList>,
    binary: bool,
) -> Option<BinaryDataArrayListOut<'a>> {
    let list = list?;
    let mut out = Vec::with_capacity(list.binary_data_arrays.len());
    for ba in &list.binary_data_arrays {
        out.push(bda_out(ba, binary));
    }
    Some(BinaryDataArrayListOut {
        count: &list.count,
        binary_data_arrays: out,
    })
}

#[derive(Serialize)]
struct SpectrumOut<'a> {
    id: &'a str,
    index: &'a Option<u32>,
    scan_number: &'a Option<u32>,
    default_array_length: &'a Option<usize>,
    native_id: &'a Option<String>,
    data_processing_ref: &'a Option<String>,
    source_file_ref: &'a Option<String>,
    spot_id: &'a Option<String>,
    ms_level: &'a Option<u32>,

    referenceable_param_group_refs: &'a Vec<ReferenceableParamGroupRef>,
    cv_params: &'a Vec<CvParam>,
    user_params: &'a Vec<UserParam>,

    spectrum_description: &'a Option<SpectrumDescription>,
    scan_list: &'a Option<ScanList>,
    precursor_list: &'a Option<PrecursorList>,
    product_list: &'a Option<ProductList>,

    binary_data_array_list: Option<BinaryDataArrayListOut<'a>>,
}

fn spectrum_out<'a>(s: &'a Spectrum, binary: bool) -> SpectrumOut<'a> {
    SpectrumOut {
        id: s.id.as_str(),
        index: &s.index,
        scan_number: &s.scan_number,
        default_array_length: &s.default_array_length,
        native_id: &s.native_id,
        data_processing_ref: &s.data_processing_ref,
        source_file_ref: &s.source_file_ref,
        spot_id: &s.spot_id,
        ms_level: &s.ms_level,
        referenceable_param_group_refs: &s.referenceable_param_group_refs,
        cv_params: &s.cv_params,
        user_params: &s.user_params,
        spectrum_description: &s.spectrum_description,
        scan_list: &s.scan_list,
        precursor_list: &s.precursor_list,
        product_list: &s.product_list,
        binary_data_array_list: bda_list_out(s.binary_data_array_list.as_ref(), binary),
    }
}

#[derive(Serialize)]
struct ChromatogramOut<'a> {
    id: &'a str,
    native_id: &'a Option<String>,
    index: &'a Option<u32>,
    default_array_length: &'a Option<usize>,
    data_processing_ref: &'a Option<String>,

    referenceable_param_group_refs: &'a Vec<ReferenceableParamGroupRef>,
    cv_params: &'a Vec<CvParam>,
    user_params: &'a Vec<UserParam>,

    precursor: &'a Option<Precursor>,
    product: &'a Option<Product>,

    binary_data_array_list: Option<BinaryDataArrayListOut<'a>>,
}

fn chromatogram_out<'a>(c: &'a Chromatogram, binary: bool) -> ChromatogramOut<'a> {
    ChromatogramOut {
        id: c.id.as_str(),
        native_id: &c.native_id,
        index: &c.index,
        default_array_length: &c.default_array_length,
        data_processing_ref: &c.data_processing_ref,
        referenceable_param_group_refs: &c.referenceable_param_group_refs,
        cv_params: &c.cv_params,
        user_params: &c.user_params,
        precursor: &c.precursor,
        product: &c.product,
        binary_data_array_list: bda_list_out(c.binary_data_array_list.as_ref(), binary),
    }
}

fn show(cmd: ShowArgs) -> Result<(), String> {
    let cwd = std::env::current_dir().map_err(|e| format!("get current dir failed: {e}"))?;
    let file_path = resolve_user_path(&cwd, &cmd.file_path);

    let mzml = read_mzml_or_b64(&file_path)?;

    if cmd.which.general {
        let out = GeneralOut {
            cv_list: &mzml.cv_list,
            file_description: &mzml.file_description,
            referenceable_param_group_list: &mzml.referenceable_param_group_list,
            sample_list: &mzml.sample_list,
            instrument_list: &mzml.instrument_list,
            software_list: &mzml.software_list,
            data_processing_list: &mzml.data_processing_list,
            scan_settings_list: &mzml.scan_settings_list,
        };
        return print_json(&out);
    }

    if cmd.which.run {
        let r: &Run = &mzml.run;
        let out = RunOut {
            id: r.id.as_str(),
            start_time_stamp: &r.start_time_stamp,
            default_instrument_configuration_ref: &r.default_instrument_configuration_ref,
            default_source_file_ref: &r.default_source_file_ref,
            sample_ref: &r.sample_ref,
            referenceable_param_group_refs: &r.referenceable_param_group_refs,
            cv_params: &r.cv_params,
            user_params: &r.user_params,
            source_file_ref_list: &r.source_file_ref_list,
        };
        return print_json(&out);
    }

    if cmd.which.spectrum_list {
        let out = mzml
            .run
            .spectrum_list
            .as_ref()
            .map(|sl: &SpectrumList| SpectrumListOut {
                count: &sl.count,
                default_data_processing_ref: &sl.default_data_processing_ref,
            });
        return print_json(&out);
    }

    if cmd.which.chromatogram_list {
        let out = mzml
            .run
            .chromatogram_list
            .as_ref()
            .map(|cl: &ChromatogramList| ChromatogramListOut {
                count: &cl.count,
                default_data_processing_ref: &cl.default_data_processing_ref,
            });
        return print_json(&out);
    }

    let items_spec = parse_items_spec(&cmd.items)?;

    if cmd.which.spectrum {
        let spectra: &[Spectrum] = mzml
            .run
            .spectrum_list
            .as_ref()
            .map(|sl| sl.spectra.as_slice())
            .unwrap_or(&[]);
        let (s, e, single) = slice_indices(spectra.len(), &items_spec);
        if s == e {
            return Err("items out of bounds".to_string());
        }
        if single {
            return print_json(&spectrum_out(&spectra[s], cmd.binary));
        }
        let mut out = Vec::with_capacity(e - s);
        for i in s..e {
            out.push(spectrum_out(&spectra[i], cmd.binary));
        }
        return print_json(&out);
    }

    if cmd.which.chromatogram {
        let chromatograms: &[Chromatogram] = mzml
            .run
            .chromatogram_list
            .as_ref()
            .map(|cl| cl.chromatograms.as_slice())
            .unwrap_or(&[]);
        let (s, e, single) = slice_indices(chromatograms.len(), &items_spec);
        if s == e {
            return Err("items out of bounds".to_string());
        }
        if single {
            return print_json(&chromatogram_out(&chromatograms[s], cmd.binary));
        }
        let mut out = Vec::with_capacity(e - s);
        for i in s..e {
            out.push(chromatogram_out(&chromatograms[i], cmd.binary));
        }
        return print_json(&out);
    }

    Err("no show mode selected".to_string())
}

fn convert(cmd: ConvertArgs) -> Result<(), String> {
    let workspace = workspace_root();
    let cwd = std::env::current_dir().map_err(|e| format!("get current dir failed: {e}"))?;

    let input_root = match cmd.input_path.as_deref() {
        Some(p) => resolve_user_path(&cwd, p),
        None => workspace.join("crates/parser/data/mzml"),
    };

    let output_root = match cmd.output_path.as_deref() {
        Some(p) => resolve_user_path(&cwd, p),
        None => workspace.join("crates/parser/data/b64"),
    };

    fs::create_dir_all(&output_root).map_err(|e| format!("create output dir failed: {e}"))?;

    const MB: f64 = 1024.0 * 1024.0;

    if cmd.which.mzml_to_b64 || cmd.which.mzml_to_b32 {
        let out_ext = if cmd.which.mzml_to_b32 { "b32" } else { "b64" };
        let f32_compress = cmd.which.mzml_to_b32;

        let files = collect_files_with_exts(&input_root, &["mzml"])?;
        if files.is_empty() {
            return Err(format!(
                "no .mzML files found under {}",
                input_root.display()
            ));
        }

        let mut ok = 0u32;
        let mut failed = 0u32;
        let mut skipped = 0u32;

        let total = files.len();
        for (i, in_path) in files.into_iter().enumerate() {
            let idx = i + 1;

            let rel = match in_path.strip_prefix(&input_root) {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("{}: cannot make relative path", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let out_name = match out_name_for_mzml_file(&in_path, out_ext) {
                Some(v) => v,
                None => continue,
            };

            let parent_rel = rel.parent().unwrap_or_else(|| Path::new(""));
            let out_dir = output_root.join(parent_rel);
            let out_path = out_dir.join(out_name);

            if !cmd.overwrite {
                if let Ok(m) = fs::metadata(&out_path) {
                    if m.is_file() && m.len() > 0 {
                        let in_mb = fs::metadata(&in_path)
                            .map(|m| m.len() as f64 / MB)
                            .unwrap_or(0.0);
                        let out_mb = m.len() as f64 / MB;

                        println!(
                            "[{}/{}] skip: {}  input={:.2} MB, output={:.2} MB",
                            idx,
                            total,
                            out_path.display(),
                            in_mb,
                            out_mb
                        );

                        skipped += 1;
                        continue;
                    }
                }
            }

            if let Err(e) = fs::create_dir_all(&out_dir) {
                eprintln!("{}: create output dir failed: {e}", out_dir.display());
                failed += 1;
                continue;
            }

            let bytes = match fs::read(&in_path) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{}: read failed: {e}", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let mzml = match parse_mzml(&bytes, false) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{}: parse_mzml failed: {e}", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let encoded = encode(&mzml, cmd.compression_level, f32_compress);

            let in_mb = bytes.len() as f64 / MB;
            let out_mb = encoded.len() as f64 / MB;

            println!(
                "[{}/{}] output: {}  input={:.2} MB, output={:.2} MB",
                idx,
                total,
                out_path.display(),
                in_mb,
                out_mb
            );

            if let Err(e) = fs::write(&out_path, encoded) {
                eprintln!("{}: write failed: {e}", out_path.display());
                failed += 1;
                continue;
            }

            ok += 1;
        }

        println!("converted_ok={ok} converted_failed={failed} converted_skipped={skipped}");
        if failed != 0 {
            return Err("some files failed".to_string());
        }
        return Ok(());
    }

    if cmd.which.b64_to_mzml {
        let files = collect_files_with_exts(&input_root, &["b64", "b32"])?;
        if files.is_empty() {
            return Err(format!(
                "no .b64/.b32 files found under {}",
                input_root.display()
            ));
        }

        let mut ok = 0u32;
        let mut failed = 0u32;
        let mut skipped = 0u32;

        let total = files.len();
        for (i, in_path) in files.into_iter().enumerate() {
            let idx = i + 1;

            let rel = match in_path.strip_prefix(&input_root) {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("{}: cannot make relative path", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let out_name = match out_name_for_bin_file_as_mzml(&in_path) {
                Some(v) => v,
                None => continue,
            };

            let parent_rel = rel.parent().unwrap_or_else(|| Path::new(""));
            let out_dir = output_root.join(parent_rel);
            let out_path = out_dir.join(out_name);

            if !cmd.overwrite {
                if let Ok(m) = fs::metadata(&out_path) {
                    if m.is_file() && m.len() > 0 {
                        let in_mb = fs::metadata(&in_path)
                            .map(|m| m.len() as f64 / MB)
                            .unwrap_or(0.0);
                        let out_mb = m.len() as f64 / MB;

                        println!(
                            "[{}/{}] skip: {}  input={:.2} MB, output={:.2} MB",
                            idx,
                            total,
                            out_path.display(),
                            in_mb,
                            out_mb
                        );

                        skipped += 1;
                        continue;
                    }
                }
            }

            if let Err(e) = fs::create_dir_all(&out_dir) {
                eprintln!("{}: create output dir failed: {e}", out_dir.display());
                failed += 1;
                continue;
            }

            let in_bytes = match fs::read(&in_path) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{}: read failed: {e}", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let mzml = match read_mzml_or_b64_from_bytes(&in_path, &in_bytes) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{}: {e}", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let xml = match bin_to_mzml(&mzml) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{}: bin_to_mzml failed: {e}", in_path.display());
                    failed += 1;
                    continue;
                }
            };

            let in_mb = in_bytes.len() as f64 / MB;
            let out_mb = xml.len() as f64 / MB;

            println!(
                "[{}/{}] output: {}  input={:.2} MB, output={:.2} MB",
                idx,
                total,
                out_path.display(),
                in_mb,
                out_mb
            );

            if let Err(e) = fs::write(&out_path, xml.as_bytes()) {
                eprintln!("{}: write failed: {e}", out_path.display());
                failed += 1;
                continue;
            }

            ok += 1;
        }

        println!("converted_ok={ok} converted_failed={failed} converted_skipped={skipped}");
        if failed != 0 {
            return Err("some files failed".to_string());
        }
        return Ok(());
    }

    Err("no convert mode selected".to_string())
}

fn prune_json(v: &mut Value) -> bool {
    match v {
        Value::Null => false,

        Value::String(s) => !s.is_empty(),

        Value::Array(xs) => {
            xs.retain_mut(|x| prune_json(x));
            !xs.is_empty()
        }

        Value::Object(m) => {
            let keys: Vec<String> = m.keys().cloned().collect();
            for k in keys {
                let keep = m.get_mut(&k).map(prune_json).unwrap_or(false);
                if !keep {
                    m.remove(&k);
                }
            }
            !m.is_empty()
        }

        Value::Bool(_) | Value::Number(_) => true,
    }
}

fn read_mzml_or_b64_from_bytes(file_path: &Path, bytes: &[u8]) -> Result<MzML, String> {
    let ext = file_ext_lower(file_path);

    if ext == "b64" || ext == "b32" {
        return decode(bytes).map_err(|e| format!("decode failed: {e}"));
    }
    if ext == "mzml" {
        return parse_mzml(bytes, false).map_err(|e| format!("parse_mzml failed: {e}"));
    }

    Err(format!(
        "unsupported file extension: {ext:?} (expected .mzML or .b64/.b32)"
    ))
}

fn resolve_user_path(cwd: &Path, p: &Path) -> PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        cwd.join(p)
    }
}
