//! Classify mapped reads into genomic regions from a BAM file + BED12 gene model.
//!
//! Mirrors the algorithm of `RSeQC` `read_distribution.py` (LGPL):
//!   - for each mapped, non-QC-fail, non-duplicate, non-secondary, non-unmapped read,
//!     split the read into exonic blocks using CIGAR (M/D/N/S ops);
//!   - classify each block by the midpoint of that block;
//!   - assign to the highest-priority overlapping feature (CDS > 5'UTR > 3'UTR > Intron
//!     > TSS up > TES down > unassigned);
//!   - emit "Total Reads / Total Tags / Total Assigned Tags" + feature table.
//!
//! ## Region construction (from BED12)
//!
//! CDS exons, 5'UTR exons, 3'UTR exons, introns are extracted per-transcript and
//! then union-merged + mutually subtracted in the same priority order `RSeQC` uses:
//!   1. CDS exons = merged CDS
//!   2. 5'UTR = merged 5'UTR − CDS
//!   3. 3'UTR = merged 3'UTR − CDS
//!   4. Introns = merged introns − CDS − 5'UTR − 3'UTR
//!   5. `TSS_up`/{`TES_down`} in 1/5/10 kb windows: each window minus all exonic/intronic features
//!
//! ## Origin
//!
//! This crate is an independent Rust reimplementation based on:
//! - `RSeQC`: `read_distribution.py` (LGPL-2.1+), Wang et al. 2012
//!   <https://doi.org/10.1093/bioinformatics/bts356>
//! - The SAM/BAM format specification (MIT)
//! - BED12 format specification
//! - Black-box behaviour testing against `RSeQC` 5.0.4
//!
//! License: MIT OR Apache-2.0.
//! Upstream credit: `RSeQC` <https://rseqc.sourceforge.net/> (LGPL-2.1+).

#![allow(clippy::cast_precision_loss)]

mod bed;
mod cigar;
mod driver;
mod index;
mod intervals;
mod tally;

pub use bed::FeatureIndex;
pub use driver::run_read_dist;
pub use tally::ReadDistResult;
