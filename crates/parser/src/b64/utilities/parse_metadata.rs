use crate::{
    b64::utilities::common::*,
    decode::{Metadatum, MetadatumValue},
    mzml::{attr_meta::format_accession, schema::TagId},
};

const HDR_CODEC_MASK: u8 = 0x0F;
const HDR_CODEC_ZSTD: u8 = 1;

pub fn parse_metadata(
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

    let item_count_usize = item_count as usize;
    let meta_count_usize = meta_count as usize;
    let num_count_usize = num_count as usize;
    let str_count_usize = str_count as usize;

    let mut pos = 0usize;

    let ci = read_u32_vec(bytes, &mut pos, item_count_usize + 1)?;
    let moi = read_u32_vec(bytes, &mut pos, meta_count_usize)?; // MOI
    let mpi = read_u32_vec(bytes, &mut pos, meta_count_usize)?; // MPI

    let mti = take(bytes, &mut pos, meta_count_usize, "metadatum tag id")?;
    let mri = take(bytes, &mut pos, meta_count_usize, "metadatum ref id")?;
    let man = read_u32_vec(bytes, &mut pos, meta_count_usize)?;
    let muri = take(bytes, &mut pos, meta_count_usize, "metadatum unit ref id")?;
    let muan = read_u32_vec(bytes, &mut pos, meta_count_usize)?;
    let vk = take(bytes, &mut pos, meta_count_usize, "metadatum value kind")?;
    let vi = read_u32_vec(bytes, &mut pos, meta_count_usize)?;

    let vn = read_f64_vec(bytes, &mut pos, num_count_usize)?;
    let voff = read_u32_vec(bytes, &mut pos, str_count_usize)?;
    let vlen = read_u32_vec(bytes, &mut pos, str_count_usize)?;

    let vs_needed = vs_len_bytes(vk, &vi, &voff, &vlen)?;
    let vs = take(bytes, &mut pos, vs_needed, "string values")?;

    if !compressed {
        let trailing = &bytes[pos..];
        if trailing.len() > 7 || trailing.iter().any(|&b| b != 0) {
            return Err("trailing bytes in metadata section".to_string());
        }
    } else if pos != bytes.len() {
        return Err("trailing bytes in decompressed metadata section".to_string());
    }

    if ci.is_empty() || ci[0] != 0 {
        return Err("CI[0] must be 0".to_string());
    }
    if ci[item_count_usize] as usize != meta_count_usize {
        return Err("CI[last] must equal meta_count".to_string());
    }

    let mut prev = 0u32;
    for &x in &ci {
        if x < prev || (x as usize) > meta_count_usize {
            return Err("CI is not monotonic or out of range".to_string());
        }
        prev = x;
    }

    let mut out = Vec::with_capacity(meta_count_usize);

    for item_index in 0..item_count_usize {
        let start = ci[item_index] as usize;
        let end = ci[item_index + 1] as usize;

        for j in start..end {
            let tag_id = TagId::from_u8(mti[j]).unwrap_or(TagId::Unknown);

            let value = match vk[j] {
                0 => {
                    let idx = vi[j] as usize;
                    if idx >= vn.len() {
                        return Err("numeric VI out of range".to_string());
                    }
                    MetadatumValue::Number(vn[idx])
                }
                1 => {
                    let idx = vi[j] as usize;
                    if idx >= voff.len() || idx >= vlen.len() {
                        return Err("string VI out of range".to_string());
                    }
                    let off = voff[idx] as usize;
                    let len = vlen[idx] as usize;
                    if off.checked_add(len).map_or(true, |e| e > vs.len()) {
                        return Err("string slice out of bounds".to_string());
                    }
                    MetadatumValue::Text(String::from_utf8_lossy(&vs[off..off + len]).into_owned())
                }
                2 => MetadatumValue::Empty,
                _ => MetadatumValue::Empty,
            };

            let accession = format_accession(mri[j], man[j]);
            let unit_accession = format_accession(muri[j], muan[j]);

            out.push(Metadatum {
                item_index: item_index as u32,
                owner_id: moi[j],
                parent_index: mpi[j],
                tag_id,
                accession,
                unit_accession,
                value,
            });
        }
    }

    Ok(out)
}
