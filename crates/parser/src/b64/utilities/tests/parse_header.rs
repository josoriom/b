use std::{fs, path::PathBuf};

use crate::b64::utilities::parse_header;

const PATH: &str = "data/b64/test.b64";

fn read_bytes(path: &str) -> Vec<u8> {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(path);
    fs::read(&full).unwrap_or_else(|e| panic!("cannot read {:?}: {}", full, e))
}

#[test]
fn check_header() {
    let bytes = read_bytes(PATH);
    let header = parse_header(&bytes).expect("parse_header failed");

    assert_eq!(header.file_signature, [66, 48, 48, 48]);
    assert_eq!(header.endianness_flag, 0);
    assert_eq!(header.reserved_alignment, [1, 0, 0]);

    assert_eq!(header.off_spec_index, 192);
    assert_eq!(header.off_chrom_index, 256);
    assert_eq!(header.off_spec_meta, 320);
    assert_eq!(header.off_chrom_meta, 824);
    assert_eq!(header.off_global_meta, 1128);

    assert_eq!(header.size_container_spect_x, 67);
    assert_eq!(header.off_container_spect_x, 1952);
    assert_eq!(header.size_container_spect_y, 76);
    assert_eq!(header.off_container_spect_y, 2024);
    assert_eq!(header.size_container_chrom_x, 67);
    assert_eq!(header.off_container_chrom_x, 2104);
    assert_eq!(header.size_container_chrom_y, 76);
    assert_eq!(header.off_container_chrom_y, 2176);

    assert_eq!(header.spectrum_count, 2);
    assert_eq!(header.chrom_count, 2);

    assert_eq!(header.spec_meta_count, 54);
    assert_eq!(header.spec_num_count, 32);
    assert_eq!(header.spec_str_count, 3);

    assert_eq!(header.chrom_meta_count, 38);
    assert_eq!(header.chrom_num_count, 15);
    assert_eq!(header.chrom_str_count, 5);

    assert_eq!(header.global_meta_count, 49);
    assert_eq!(header.global_num_count, 8);
    assert_eq!(header.global_str_count, 25);

    assert_eq!(header.block_count_spect_x, 1);
    assert_eq!(header.block_count_spect_y, 1);
    assert_eq!(header.block_count_chrom_x, 1);
    assert_eq!(header.block_count_chrom_y, 1);

    assert_eq!(header.reserved_flags, 113);
    assert_eq!(header.chrom_x_format, 2);
    assert_eq!(header.chrom_y_format, 1);
    assert_eq!(header.spect_x_format, 2);
    assert_eq!(header.spect_y_format, 1);
    assert_eq!(header.compression_level, 12);
    assert_eq!(header.array_filter, 1);

    assert_eq!(header.reserved, [0; 13]);

    let len = bytes.len() as u64;
    for &off in &[
        header.off_spec_index,
        header.off_chrom_index,
        header.off_spec_meta,
        header.off_chrom_meta,
        header.off_global_meta,
        header.off_container_spect_x,
        header.off_container_spect_y,
        header.off_container_chrom_x,
        header.off_container_chrom_y,
    ] {
        assert!(off < len, "offset {off} out of bounds (len={len})");
    }
}
