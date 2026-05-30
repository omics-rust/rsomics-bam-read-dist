use std::num::NonZero;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

use crate::bed::FeatureIndex;
use crate::cigar::{FLAG_DUPLICATE, FLAG_QCFAIL, FLAG_SECONDARY, FLAG_UNMAPPED, fetch_exon_blocks};
use crate::tally::{ReadDistResult, classify_tag};

/// Run the full read distribution analysis over `bam_path` against the `bed_path` gene model.
pub fn run_read_dist(
    bam_path: &Path,
    bed_path: &Path,
    workers: NonZero<usize>,
) -> Result<ReadDistResult> {
    eprintln!("processing {} ...", bed_path.display());
    let index = FeatureIndex::from_bed12(bed_path)?;
    eprintln!("Done");

    eprintln!("processing {} ...", bam_path.display());

    let mut reader = rsomics_bamio::open_with_workers(bam_path, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    let mut result = ReadDistResult {
        cds_exons_bases: index.cds.total_bases,
        utr5_bases: index.utr5.total_bases,
        utr3_bases: index.utr3.total_bases,
        intron_bases: index.intron.total_bases,
        tss_up_1kb_bases: index.tss_up_1kb.total_bases,
        tss_up_5kb_bases: index.tss_up_5kb.total_bases,
        tss_up_10kb_bases: index.tss_up_10kb.total_bases,
        tes_down_1kb_bases: index.tes_down_1kb.total_bases,
        tes_down_5kb_bases: index.tes_down_5kb.total_bases,
        tes_down_10kb_bases: index.tes_down_10kb.total_bases,
        ..Default::default()
    };

    let mut rec = RawRecord::default();
    loop {
        let bytes_read = raw::read_record(reader.get_mut(), &mut rec)?;
        if bytes_read == 0 {
            break;
        }

        let flags = rec.flags();
        if flags & (FLAG_QCFAIL | FLAG_DUPLICATE | FLAG_SECONDARY | FLAG_UNMAPPED) != 0 {
            continue;
        }

        result.total_reads += 1;

        let tid = rec.reference_sequence_id();
        if tid < 0 {
            continue;
        }
        #[allow(clippy::cast_sign_loss)]
        let Some(chrom) = ref_names.get(tid as usize) else {
            continue;
        };
        let chrom_upper = chrom.to_uppercase();

        let read_start = rec.alignment_start();
        let blocks = fetch_exon_blocks(read_start, rec.cigar_ops());
        result.total_tags += blocks.len() as u64;

        for (bstart, bend) in blocks {
            // Midpoint: int(st) + int((int(end)-int(st))/2), matching RSeQC.
            let mid = bstart + (bend - bstart) / 2;
            classify_tag(mid, &chrom_upper, &index, &mut result);
        }
    }

    eprintln!("Finished");

    Ok(result)
}
