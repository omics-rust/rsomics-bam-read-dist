//! Classify mapped reads into genomic regions from a BAM + BED12 gene model.
//!
//! ## Origin
//!
//! Independent Rust reimplementation based on:
//! - `RSeQC` `read_distribution.py` (LGPL-2.1+), Wang et al. 2012
//!   <https://doi.org/10.1093/bioinformatics/bts356>
//! - SAM/BAM format specification; BED12 format specification
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
