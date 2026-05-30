use std::io::Write;

use serde::Serialize;

use crate::bed::FeatureIndex;

/// Per-feature read counters and base totals.
#[derive(Debug, Default, Serialize)]
pub struct ReadDistResult {
    pub total_reads: u64,
    pub total_tags: u64,
    pub unassigned_tags: u64,

    pub cds_exons_bases: i64,
    pub cds_exons_tags: u64,
    pub utr5_bases: i64,
    pub utr5_tags: u64,
    pub utr3_bases: i64,
    pub utr3_tags: u64,
    pub intron_bases: i64,
    pub intron_tags: u64,

    pub tss_up_1kb_bases: i64,
    pub tss_up_1kb_tags: u64,
    pub tss_up_5kb_bases: i64,
    pub tss_up_5kb_tags: u64,
    pub tss_up_10kb_bases: i64,
    pub tss_up_10kb_tags: u64,

    pub tes_down_1kb_bases: i64,
    pub tes_down_1kb_tags: u64,
    pub tes_down_5kb_bases: i64,
    pub tes_down_5kb_tags: u64,
    pub tes_down_10kb_bases: i64,
    pub tes_down_10kb_tags: u64,
}

impl ReadDistResult {
    #[must_use]
    pub fn assigned_tags(&self) -> u64 {
        self.total_tags - self.unassigned_tags
    }

    /// Emit the exact text format `RSeQC` `read_distribution.py` prints to stdout.
    ///
    /// `"%-30s%d"` header lines; `"%-20s%-20d%-20d%-18.2f"` table rows; `+1` denominator
    /// matches RSeQC source exactly.
    pub fn write_rseqc<W: Write>(&self, mut out: W) -> std::io::Result<()> {
        writeln!(out, "{:<30}{}", "Total Reads", self.total_reads)?;
        writeln!(out, "{:<30}{}", "Total Tags", self.total_tags)?;
        writeln!(out, "{:<30}{}", "Total Assigned Tags", self.assigned_tags())?;
        writeln!(
            out,
            "====================================================================="
        )?;
        writeln!(
            out,
            "{:<20}{:<20}{:<20}{:<20}",
            "Group", "Total_bases", "Tag_count", "Tags/Kb"
        )?;
        let row = |out: &mut W, name: &str, bases: i64, tags: u64| -> std::io::Result<()> {
            let tags_per_kb = tags as f64 * 1000.0 / (bases + 1) as f64;
            writeln!(out, "{name:<20}{bases:<20}{tags:<20}{tags_per_kb:<18.2}")
        };
        row(
            &mut out,
            "CDS_Exons",
            self.cds_exons_bases,
            self.cds_exons_tags,
        )?;
        row(&mut out, "5'UTR_Exons", self.utr5_bases, self.utr5_tags)?;
        row(&mut out, "3'UTR_Exons", self.utr3_bases, self.utr3_tags)?;
        row(&mut out, "Introns", self.intron_bases, self.intron_tags)?;
        row(
            &mut out,
            "TSS_up_1kb",
            self.tss_up_1kb_bases,
            self.tss_up_1kb_tags,
        )?;
        row(
            &mut out,
            "TSS_up_5kb",
            self.tss_up_5kb_bases,
            self.tss_up_5kb_tags,
        )?;
        row(
            &mut out,
            "TSS_up_10kb",
            self.tss_up_10kb_bases,
            self.tss_up_10kb_tags,
        )?;
        row(
            &mut out,
            "TES_down_1kb",
            self.tes_down_1kb_bases,
            self.tes_down_1kb_tags,
        )?;
        row(
            &mut out,
            "TES_down_5kb",
            self.tes_down_5kb_bases,
            self.tes_down_5kb_tags,
        )?;
        row(
            &mut out,
            "TES_down_10kb",
            self.tes_down_10kb_bases,
            self.tes_down_10kb_tags,
        )?;
        writeln!(
            out,
            "====================================================================="
        )?;
        Ok(())
    }
}

/// Classify one tag (exon-block midpoint) and increment the appropriate counter.
///
/// Priority: CDS → 5'UTR (exclusive) → 3'UTR (exclusive) → both UTR → intron →
/// both TSS+TES 10kb → TSS_up/TES_down 1/5/10kb → unassigned. Matches RSeQC source.
pub(crate) fn classify_tag(mid: i32, chrom: &str, idx: &FeatureIndex, res: &mut ReadDistResult) {
    if idx.cds.contains(chrom, mid) {
        res.cds_exons_tags += 1;
        return;
    }
    let in_utr5 = idx.utr5.contains(chrom, mid);
    let in_utr3 = idx.utr3.contains(chrom, mid);
    if in_utr5 && !in_utr3 {
        res.utr5_tags += 1;
        return;
    }
    if in_utr3 && !in_utr5 {
        res.utr3_tags += 1;
        return;
    }
    if in_utr5 && in_utr3 {
        // Overlapping UTR5+UTR3 — unassigned per RSeQC.
        res.unassigned_tags += 1;
        return;
    }
    if idx.intron.contains(chrom, mid) {
        res.intron_tags += 1;
        return;
    }
    if idx.tss_up_10kb.contains(chrom, mid) && idx.tes_down_10kb.contains(chrom, mid) {
        // Overlaps both 10kb windows — unassigned per RSeQC.
        res.unassigned_tags += 1;
        return;
    }
    if idx.tss_up_1kb.contains(chrom, mid) {
        res.tss_up_1kb_tags += 1;
        res.tss_up_5kb_tags += 1;
        res.tss_up_10kb_tags += 1;
        return;
    }
    if idx.tss_up_5kb.contains(chrom, mid) {
        res.tss_up_5kb_tags += 1;
        res.tss_up_10kb_tags += 1;
        return;
    }
    if idx.tss_up_10kb.contains(chrom, mid) {
        res.tss_up_10kb_tags += 1;
        return;
    }
    if idx.tes_down_1kb.contains(chrom, mid) {
        res.tes_down_1kb_tags += 1;
        res.tes_down_5kb_tags += 1;
        res.tes_down_10kb_tags += 1;
        return;
    }
    if idx.tes_down_5kb.contains(chrom, mid) {
        res.tes_down_5kb_tags += 1;
        res.tes_down_10kb_tags += 1;
        return;
    }
    if idx.tes_down_10kb.contains(chrom, mid) {
        res.tes_down_10kb_tags += 1;
        return;
    }
    res.unassigned_tags += 1;
}
