# `b` CLI — Commands & Flags (with examples)

This CLI has **two subcommands**:

- `convert` — batch convert folders (recursive)
- `show` — print JSON metadata for a single file (`.mzML`, `.b64`, or `.b32`)

Examples below use:

- `tiny.msdata.mzML0.99.9.mzML`
- `tiny.msdata.mzML0.99.9.b64`

---

## 1) `octo convert`

Batch converts files under a folder **recursively**, preserving the relative folder structure under `--output-path`.

### Common flags

- `--input-path <DIR>`

  - Directory to scan recursively (required)

- `--output-path <DIR>`

  - Output directory root (subfolders are created) (required)

- `--level <0..22>`
  - Default: `12`
  - Compression level (only relevant for mzML → b64/b32)

### Choose exactly one conversion mode (required)

- `--mzml-to-b64`

  - Converts `.mzML` → `.b64`
  - Uses `encode(mzml, level, f32_compress=false)`

- `--mzml-to-b32`

  - Converts `.mzML` → `.b32`
  - Uses `encode(mzml, level, f32_compress=true)`

- `--b64-to-mzml`
  - Converts `.b64` and `.b32` → `.mzML`
  - Uses `decode(...)` into `MzML`, then `bin_to_mzml(&mzml)` to produce mzML XML

---

### 1.1) Convert `.mzML` → `.b64/.b32`

Defaults:

- `--level 12`

Run (defaults):

```bash
octo convert --mzml-to-b64 --input-path crates/parser/data/mzml --output-path crates/parser/data/b64
```

Custom compression level:

```bash
octo convert --mzml-to-b64 --input-path crates/parser/data/mzml --output-path crates/parser/data/b64 --level 12
```

---

### 1.2) Convert `.b64/.b32` → `.mzML`

### Example:

```bash
octo convert --b64-to-mzml --input-path crates/parser/data/b64 --output-path crates/parser/data/mzml_out
```

(If the input tree contains `.b32` too, it is also converted by this mode.)

---

## 2) `octo show`

Prints JSON to stdout for a **single file**.

### Common flags

- `--file-path <PATH>`
  - Path to the file you want to inspect.
  - Supported file types: `.mzML`, `.b64`, `.b32`
  - If relative, it is resolved under the workspace root (per your `workspace_root + resolve_under_root` logic).

### Choose exactly one output target (required)

You must pass **one** of:

- `--general`
- `--run`
- `--spectrum-list`
- `--chromatogram-list`
- `--spectrum`
- `--chromatogram`

### Extra flags (only for `--spectrum` and `--chromatogram`)

These flags are valid **only when** selecting `--spectrum` or `--chromatogram`:

- `--show-numbers`

  - Default: `false`
  - `false`: prints spectrum/chromatogram metadata **without** decoded arrays
  - `true`: includes decoded arrays (`decoded_binary_f32/decoded_binary_f64`) in JSON

- `--items <SPEC>`
  - Default: `0-100`
  - Selects which items to print from:
    - `SpectrumList.spectra` (when using `--spectrum`)
    - `ChromatogramList.chromatograms` (when using `--chromatogram`)
  - Supported formats:
    - `<N>` (single item): `--items 5`
    - `<START>-<END>` (range, end is **exclusive**): `--items 0-100` prints items `0..99`

Notes:

- `--items` and `--show-numbers` are rejected/meaningless unless `--spectrum` or `--chromatogram` is selected.
- Output JSON is pruned:
  - `null` fields removed
  - `""` (empty string) fields removed
  - empty arrays/objects removed after pruning

---

### 2.1) `octo show --general`

Shows the general metadata (top-level fields in `MzML`, excluding `run`):

- `cv_list`
- `file_description`
- `referenceable_param_group_list`
- `sample_list`
- `instrument_list`
- `software_list`
- `data_processing_list`
- `scan_settings_list`

### Example:

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --general
```

---

### 2.2) `octo show --run`

Shows `Run` metadata **excluding** `spectrum_list` and `chromatogram_list`:

Includes:

- `id`
- `start_time_stamp`
- `default_instrument_configuration_ref`
- `default_source_file_ref`
- `sample_ref`
- `referenceable_param_group_refs`
- `cv_params`
- `user_params`
- `source_file_ref_list`

### Example:

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --run
```

---

### 2.3) `octo show --spectrum-list`

Shows only `SpectrumList` list metadata (NOT the `spectra` array):

- `count`
- `default_data_processing_ref`

### Example:

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum-list
```

---

### 2.4) `octo show --chromatogram-list`

Shows only `ChromatogramList` list metadata (NOT the `chromatograms` array):

- `count`
- `default_data_processing_ref`

### Example:

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --chromatogram-list
```

---

### 2.5) `octo show --spectrum`

Prints one or more `Spectrum` items from:

- `mzml.run.spectrum_list.spectra`

### Supports:

- `--items` (which spectra)
- `--show-numbers` (include decoded arrays)
  Same examples (b64 and b32):

### Formats:

- `.mzML`
- `.b32`
- `.b64`

### Examples:

First 100 (0..99):

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum
```

One item (index 5):

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum --items 5
```

Range (10..99):

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum --items 10-100
```

Include decoded arrays:

```bash
octo show --file-path crates/parser/data/b64/tiny.msdata.mzML0.99.9.b64 --spectrum --items 0-10 --show-numbers
```

---

### 2.6) `octo show --chromatogram`

Same as `2.5` but with the flag `--chromatogram`

---
