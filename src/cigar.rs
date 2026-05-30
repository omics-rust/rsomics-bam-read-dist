// BAM flag bits (SAM spec).
pub(crate) const FLAG_QCFAIL: u16 = 0x0200;
pub(crate) const FLAG_DUPLICATE: u16 = 0x0400;
pub(crate) const FLAG_SECONDARY: u16 = 0x0100;
pub(crate) const FLAG_UNMAPPED: u16 = 0x0004;

// CIGAR op codes (BAM spec).
const CIGAR_MATCH: u32 = 0; // M
const CIGAR_DEL: u32 = 2; // D
const CIGAR_REF_SKIP: u32 = 3; // N
const CIGAR_SOFT_CLIP: u32 = 4; // S

/// Exonic blocks for one read from the CIGAR string.
///
/// M advances position and emits a block; D/N/S advance without emitting (matching RSeQC fetch_exon).
pub(crate) fn fetch_exon_blocks(
    start: i32,
    cigar_ops: impl Iterator<Item = (u8, u32)>,
) -> Vec<(i32, i32)> {
    let mut blocks = Vec::new();
    let mut pos = start;
    for (op, len) in cigar_ops {
        #[allow(clippy::cast_possible_wrap)]
        let len = len as i32;
        match u32::from(op) {
            CIGAR_MATCH => {
                blocks.push((pos, pos + len));
                pos += len;
            }
            CIGAR_DEL | CIGAR_REF_SKIP | CIGAR_SOFT_CLIP => {
                pos += len;
            }
            _ => {}
        }
    }
    blocks
}
