const HDR_FLAG_GLOBAL_META_COMP: u8 = 1 << 6;
const HEADER_SIZE: usize = 192;

pub fn parse_header(bytes: &[u8]) -> Result<Header, String> {
    if bytes.len() < HEADER_SIZE {
        return Err("header: file too small".to_string());
    }

    let mut r = Reader::new(bytes);

    let file_signature = r.read_arr::<4>("file_signature")?;
    let endianness_flag = r.read_u8("endianness_flag")?;
    let reserved_alignment = r.read_arr::<3>("reserved_alignment")?;

    if &file_signature != b"B000" {
        return Err("header: invalid file_signature (expected \"B000\")".into());
    }
    if endianness_flag != 0 {
        return Err("header: expected little-endian endianness_flag=0".into());
    }

    let off_spec_index = r.read_u64_le("off_spec_index")?;
    let off_chrom_index = r.read_u64_le("off_chrom_index")?;
    let off_spec_meta = r.read_u64_le("off_spec_meta")?;
    let off_chrom_meta = r.read_u64_le("off_chrom_meta")?;
    let off_global_meta = r.read_u64_le("off_global_meta")?;

    let size_container_spect_x = r.read_u64_le("size_container_spect_x")?;
    let off_container_spect_x = r.read_u64_le("off_container_spect_x")?;
    let size_container_spect_y = r.read_u64_le("size_container_spect_y")?;
    let off_container_spect_y = r.read_u64_le("off_container_spect_y")?;
    let size_container_chrom_x = r.read_u64_le("size_container_chrom_x")?;
    let off_container_chrom_x = r.read_u64_le("off_container_chrom_x")?;
    let size_container_chrom_y = r.read_u64_le("size_container_chrom_y")?;
    let off_container_chrom_y = r.read_u64_le("off_container_chrom_y")?;

    let spectrum_count = r.read_u32_le("spectrum_count")?;
    let chrom_count = r.read_u32_le("chrom_count")?;

    let spec_meta_count = r.read_u32_le("spec_meta_count")?;
    let spec_num_count = r.read_u32_le("spec_num_count")?;
    let spec_str_count = r.read_u32_le("spec_str_count")?;

    let chrom_meta_count = r.read_u32_le("chrom_meta_count")?;
    let chrom_num_count = r.read_u32_le("chrom_num_count")?;
    let chrom_str_count = r.read_u32_le("chrom_str_count")?;

    let global_meta_count = r.read_u32_le("global_meta_count")?;
    let global_num_count = r.read_u32_le("global_num_count")?;
    let global_str_count = r.read_u32_le("global_str_count")?;

    let block_count_spect_x = r.read_u32_le("block_count_spect_x")?;
    let block_count_spect_y = r.read_u32_le("block_count_spect_y")?;
    let block_count_chrom_x = r.read_u32_le("block_count_chrom_x")?;
    let block_count_chrom_y = r.read_u32_le("block_count_chrom_y")?;

    let reserved_flags = r.read_u8("reserved_flags")?;
    let chrom_x_format = r.read_u8("chrom_x_format")?;
    let chrom_y_format = r.read_u8("chrom_y_format")?;
    let spect_x_format = r.read_u8("spect_x_format")?;
    let spect_y_format = r.read_u8("spect_y_format")?;
    let compression_level = r.read_u8("compression_level")?;
    let array_filter = r.read_u8("array_filter")?;
    let reserved = r.read_arr::<13>("reserved")?;

    let _ = reserved_flags & HDR_FLAG_GLOBAL_META_COMP;

    Ok(Header {
        file_signature,
        endianness_flag,
        reserved_alignment,

        off_spec_index,
        off_chrom_index,
        off_spec_meta,
        off_chrom_meta,
        off_global_meta,

        size_container_spect_x,
        off_container_spect_x,
        size_container_spect_y,
        off_container_spect_y,
        size_container_chrom_x,
        off_container_chrom_x,
        size_container_chrom_y,
        off_container_chrom_y,

        spectrum_count,
        chrom_count,

        spec_meta_count,
        spec_num_count,
        spec_str_count,

        chrom_meta_count,
        chrom_num_count,
        chrom_str_count,

        global_meta_count,
        global_num_count,
        global_str_count,

        block_count_spect_x,
        block_count_spect_y,
        block_count_chrom_x,
        block_count_chrom_y,

        reserved_flags,
        chrom_x_format,
        chrom_y_format,
        spect_x_format,
        spect_y_format,
        compression_level,
        array_filter,

        reserved,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    pub file_signature: [u8; 4],
    pub endianness_flag: u8,
    pub reserved_alignment: [u8; 3],

    pub off_spec_index: u64,
    pub off_chrom_index: u64,
    pub off_spec_meta: u64,
    pub off_chrom_meta: u64,
    pub off_global_meta: u64,

    pub size_container_spect_x: u64,
    pub off_container_spect_x: u64,
    pub size_container_spect_y: u64,
    pub off_container_spect_y: u64,
    pub size_container_chrom_x: u64,
    pub off_container_chrom_x: u64,
    pub size_container_chrom_y: u64,
    pub off_container_chrom_y: u64,

    pub spectrum_count: u32,
    pub chrom_count: u32,

    pub spec_meta_count: u32,
    pub spec_num_count: u32,
    pub spec_str_count: u32,

    pub chrom_meta_count: u32,
    pub chrom_num_count: u32,
    pub chrom_str_count: u32,

    pub global_meta_count: u32,
    pub global_num_count: u32,
    pub global_str_count: u32,

    pub block_count_spect_x: u32,
    pub block_count_spect_y: u32,
    pub block_count_chrom_x: u32,
    pub block_count_chrom_y: u32,

    pub reserved_flags: u8,
    pub chrom_x_format: u8,
    pub chrom_y_format: u8,
    pub spect_x_format: u8,
    pub spect_y_format: u8,
    pub compression_level: u8,
    pub array_filter: u8,

    pub reserved: [u8; 13],
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    #[inline]
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    #[inline]
    fn need(&self, n: usize, field: &str) -> Result<(), String> {
        if self.pos + n <= self.bytes.len() {
            Ok(())
        } else {
            Err(format!(
                "header: not enough bytes for {field} at offset {} (need {n}, have {})",
                self.pos,
                self.bytes.len().saturating_sub(self.pos)
            ))
        }
    }

    #[inline]
    fn read_u8(&mut self, field: &str) -> Result<u8, String> {
        self.need(1, field)?;
        let v = self.bytes[self.pos];
        self.pos += 1;
        Ok(v)
    }

    #[inline]
    fn read_u32_le(&mut self, field: &str) -> Result<u32, String> {
        self.need(4, field)?;
        let v = u32::from_le_bytes(self.bytes[self.pos..self.pos + 4].try_into().unwrap());
        self.pos += 4;
        Ok(v)
    }

    #[inline]
    fn read_u64_le(&mut self, field: &str) -> Result<u64, String> {
        self.need(8, field)?;
        let v = u64::from_le_bytes(self.bytes[self.pos..self.pos + 8].try_into().unwrap());
        self.pos += 8;
        Ok(v)
    }

    #[inline]
    fn read_arr<const N: usize>(&mut self, field: &str) -> Result<[u8; N], String> {
        self.need(N, field)?;
        let v: [u8; N] = self.bytes[self.pos..self.pos + N].try_into().unwrap();
        self.pos += N;
        Ok(v)
    }
}
