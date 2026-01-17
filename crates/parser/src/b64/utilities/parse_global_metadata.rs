use crate::{
    b64::{
        encode::{HDR_CODEC_MASK, HDR_CODEC_ZSTD},
        utilities::{common::decompress_zstd_allow_aligned_padding, parse_metadata},
    },
    decode::Metadatum,
};

const GLOBAL_HEADER_SIZE_32: usize = 32;
const GLOBAL_HEADER_SIZE_36: usize = 36;

#[inline]
fn u32_at(bytes: &[u8], off: usize) -> Result<u32, String> {
    let end = off
        .checked_add(4)
        .ok_or_else(|| "global metadata: section too small".to_string())?;
    if end > bytes.len() {
        return Err("global metadata: section too small".to_string());
    }
    Ok(u32::from_le_bytes(bytes[off..end].try_into().unwrap()))
}

pub fn parse_global_metadata(
    bytes: &[u8],
    item_count: u32,
    meta_count: u32,
    num_count: u32,
    str_count: u32,
    compressed: bool,
    reserved_flags: u8,
) -> Result<Vec<Metadatum>, String> {
    let codec_id = reserved_flags & HDR_CODEC_MASK;

    let owned;
    let bytes = if compressed {
        if codec_id != HDR_CODEC_ZSTD {
            return Err(format!("unsupported metadata codec_id={codec_id}"));
        }
        owned = decompress_zstd_allow_aligned_padding(bytes)?;
        owned.as_slice()
    } else {
        bytes
    };

    if bytes.len() < GLOBAL_HEADER_SIZE_32 + 4 {
        return Err("global metadata: section too small".to_string());
    }

    let header_size = if bytes.len() >= GLOBAL_HEADER_SIZE_36 + 4
        && u32_at(bytes, GLOBAL_HEADER_SIZE_36)? == 0
    {
        GLOBAL_HEADER_SIZE_36
    } else if bytes.len() >= GLOBAL_HEADER_SIZE_32 + 4 && u32_at(bytes, GLOBAL_HEADER_SIZE_32)? == 0
    {
        GLOBAL_HEADER_SIZE_32
    } else {
        return Err("global metadata: missing header or corrupted CI".to_string());
    };

    let mut pos = 0usize;

    let n_file_description = u32_at(bytes, pos)?;
    pos += 4;

    let n_run = if header_size == GLOBAL_HEADER_SIZE_36 {
        let v = u32_at(bytes, pos)?;
        pos += 4;
        v
    } else {
        0
    };

    let n_ref_param_groups = u32_at(bytes, pos)?;
    pos += 4;
    let n_samples = u32_at(bytes, pos)?;
    pos += 4;
    let n_instrument_configs = u32_at(bytes, pos)?;
    pos += 4;
    let n_software = u32_at(bytes, pos)?;
    pos += 4;
    let n_data_processing = u32_at(bytes, pos)?;
    pos += 4;
    let n_acquisition_settings = u32_at(bytes, pos)?;
    pos += 4;
    let n_cvs = u32_at(bytes, pos)?;

    let derived_item_count = n_file_description
        .wrapping_add(n_run)
        .wrapping_add(n_ref_param_groups)
        .wrapping_add(n_samples)
        .wrapping_add(n_instrument_configs)
        .wrapping_add(n_software)
        .wrapping_add(n_data_processing)
        .wrapping_add(n_acquisition_settings)
        .wrapping_add(n_cvs);

    let item_count = if derived_item_count != 0 {
        derived_item_count
    } else {
        item_count
    };

    parse_metadata(
        &bytes[header_size..],
        item_count,
        meta_count,
        num_count,
        str_count,
        false,
        reserved_flags,
    )
}
