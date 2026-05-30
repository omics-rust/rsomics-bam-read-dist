use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::index::RegionIndex;
use crate::intervals::{Iv, subtract_sorted, union_merge};

/// All feature-region indexes built from a BED12 gene model.
pub struct FeatureIndex {
    pub(crate) cds: RegionIndex,
    pub(crate) utr5: RegionIndex,
    pub(crate) utr3: RegionIndex,
    pub(crate) intron: RegionIndex,
    pub(crate) tss_up_1kb: RegionIndex,
    pub(crate) tss_up_5kb: RegionIndex,
    pub(crate) tss_up_10kb: RegionIndex,
    pub(crate) tes_down_1kb: RegionIndex,
    pub(crate) tes_down_5kb: RegionIndex,
    pub(crate) tes_down_10kb: RegionIndex,
}

impl FeatureIndex {
    /// Parse a BED12 file and build all feature indexes.
    #[allow(clippy::similar_names, clippy::too_many_lines)]
    pub fn from_bed12(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| RsomicsError::Io(std::io::Error::other(format!("reading BED12: {e}"))))?;

        let mut raw_cds: Vec<Iv> = Vec::new();
        let mut raw_utr5: Vec<Iv> = Vec::new();
        let mut raw_utr3: Vec<Iv> = Vec::new();
        let mut raw_intron: Vec<Iv> = Vec::new();
        let mut raw_tss_up_1kb: Vec<Iv> = Vec::new();
        let mut raw_tss_up_5kb: Vec<Iv> = Vec::new();
        let mut raw_tss_up_10kb: Vec<Iv> = Vec::new();
        let mut raw_tes_down_1kb: Vec<Iv> = Vec::new();
        let mut raw_tes_down_5kb: Vec<Iv> = Vec::new();
        let mut raw_tes_down_10kb: Vec<Iv> = Vec::new();

        for line in content.lines() {
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                eprintln!("[NOTE:input bed must be 12-column] skipped this line: {line}");
                continue;
            }
            let chrom = fields[0];
            let Ok(tx_start) = fields[1].parse::<i32>() else {
                eprintln!("[NOTE:input bed must be 12-column] skipped this line: {line}");
                continue;
            };
            let Ok(tx_end) = fields[2].parse::<i32>() else {
                eprintln!("[NOTE:input bed must be 12-column] skipped this line: {line}");
                continue;
            };
            let strand = fields[5];
            let Ok(cds_start) = fields[6].parse::<i32>() else {
                continue;
            };
            let Ok(cds_end) = fields[7].parse::<i32>() else {
                continue;
            };
            let Ok(block_count) = fields[9].parse::<usize>() else {
                continue;
            };

            let block_sizes: Vec<i32> = fields[10]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse().ok())
                .collect();
            let block_starts: Vec<i32> = fields[11]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse::<i32>().ok())
                .map(|o| o + tx_start)
                .collect();

            if block_sizes.len() < block_count || block_starts.len() < block_count {
                continue;
            }

            let exon_starts = &block_starts[..block_count];
            let exon_ends: Vec<i32> = exon_starts
                .iter()
                .zip(block_sizes.iter())
                .map(|(&s, &sz)| s + sz)
                .collect();

            // CDS exons: intersection of each exon with [cds_start, cds_end).
            for (&es, &ee) in exon_starts.iter().zip(exon_ends.iter()) {
                if ee <= cds_start || es >= cds_end {
                    continue;
                }
                raw_cds.push(Iv {
                    chrom: chrom.to_uppercase(),
                    start: es.max(cds_start),
                    end: ee.min(cds_end),
                });
            }

            // UTRs: strand-aware, matching RSeQC getUTR.
            // '+': 5'UTR = exon parts < cdsStart; 3'UTR = exon parts > cdsEnd.
            // '-': reversed.
            for (&es, &ee) in exon_starts.iter().zip(exon_ends.iter()) {
                if es < cds_start {
                    let iv = Iv {
                        chrom: chrom.to_uppercase(),
                        start: es,
                        end: ee.min(cds_start),
                    };
                    if strand == "+" {
                        raw_utr5.push(iv);
                    } else {
                        raw_utr3.push(iv);
                    }
                }
                if ee > cds_end {
                    let iv = Iv {
                        chrom: chrom.to_uppercase(),
                        start: es.max(cds_end),
                        end: ee,
                    };
                    if strand == "+" {
                        raw_utr3.push(iv);
                    } else {
                        raw_utr5.push(iv);
                    }
                }
            }

            // Introns: gaps between consecutive exons.
            if block_count > 1 {
                for i in 0..(block_count - 1) {
                    let intron_st = exon_ends[i];
                    let intron_end = exon_starts[i + 1];
                    if intron_st < intron_end {
                        raw_intron.push(Iv {
                            chrom: chrom.to_uppercase(),
                            start: intron_st,
                            end: intron_end,
                        });
                    }
                }
            }

            // TSS upstream / TES downstream windows (1/5/10 kb), matching RSeQC getIntergenic.
            // '+': upstream of tx_start, downstream of tx_end.
            // '-': upstream of tx_end, downstream of tx_start.
            for &size in &[1000i32, 5000, 10000] {
                let (up_st, up_end, down_st, down_end) = if strand == "-" {
                    (tx_end, tx_end + size, 0i32.max(tx_start - size), tx_start)
                } else {
                    (0i32.max(tx_start - size), tx_start, tx_end, tx_end + size)
                };
                let up_iv = Iv {
                    chrom: chrom.to_uppercase(),
                    start: up_st,
                    end: up_end,
                };
                let down_iv = Iv {
                    chrom: chrom.to_uppercase(),
                    start: down_st,
                    end: down_end,
                };
                match size {
                    1000 => {
                        raw_tss_up_1kb.push(up_iv);
                        raw_tes_down_1kb.push(down_iv);
                    }
                    5000 => {
                        raw_tss_up_5kb.push(up_iv);
                        raw_tes_down_5kb.push(down_iv);
                    }
                    _ => {
                        raw_tss_up_10kb.push(up_iv);
                        raw_tes_down_10kb.push(down_iv);
                    }
                }
            }
        }

        let cds_merged = union_merge(&mut raw_cds);
        let utr5_merged = union_merge(&mut raw_utr5);
        let utr3_merged = union_merge(&mut raw_utr3);
        let intron_merged = union_merge(&mut raw_intron);
        let tss_up_1kb_merged = union_merge(&mut raw_tss_up_1kb);
        let tss_up_5kb_merged = union_merge(&mut raw_tss_up_5kb);
        let tss_up_10kb_merged = union_merge(&mut raw_tss_up_10kb);
        let tes_down_1kb_merged = union_merge(&mut raw_tes_down_1kb);
        let tes_down_5kb_merged = union_merge(&mut raw_tes_down_5kb);
        let tes_down_10kb_merged = union_merge(&mut raw_tes_down_10kb);

        // Priority: CDS > 5'UTR > 3'UTR > Intron > TSS/TES windows.
        let utr5_clean = subtract_sorted(utr5_merged.clone(), &cds_merged);
        let utr3_clean = subtract_sorted(utr3_merged.clone(), &cds_merged);
        let intron_clean = {
            let tmp = subtract_sorted(intron_merged.clone(), &cds_merged);
            let tmp = subtract_sorted(tmp, &utr5_merged);
            subtract_sorted(tmp, &utr3_merged)
        };

        let mut all_genic: Vec<Iv> = Vec::new();
        all_genic.extend_from_slice(&cds_merged);
        all_genic.extend_from_slice(&utr5_merged);
        all_genic.extend_from_slice(&utr3_merged);
        all_genic.extend_from_slice(&intron_merged);
        let all_genic_merged = union_merge(&mut all_genic);

        let clean_window = |mut ivs: Vec<Iv>| -> Vec<Iv> {
            let merged = union_merge(&mut ivs);
            subtract_sorted(merged, &all_genic_merged)
        };

        Ok(Self {
            cds: RegionIndex::from_intervals(cds_merged),
            utr5: RegionIndex::from_intervals(utr5_clean),
            utr3: RegionIndex::from_intervals(utr3_clean),
            intron: RegionIndex::from_intervals(intron_clean),
            tss_up_1kb: RegionIndex::from_intervals(clean_window(tss_up_1kb_merged)),
            tss_up_5kb: RegionIndex::from_intervals(clean_window(tss_up_5kb_merged)),
            tss_up_10kb: RegionIndex::from_intervals(clean_window(tss_up_10kb_merged)),
            tes_down_1kb: RegionIndex::from_intervals(clean_window(tes_down_1kb_merged)),
            tes_down_5kb: RegionIndex::from_intervals(clean_window(tes_down_5kb_merged)),
            tes_down_10kb: RegionIndex::from_intervals(clean_window(tes_down_10kb_merged)),
        })
    }
}
