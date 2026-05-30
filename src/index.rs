use std::collections::HashMap;

use coitrees::{COITree, Interval as CoiInterval, IntervalTree};

use crate::intervals::{Iv, union_merge};

/// Per-chromosome COITree index for one feature class (CDS, UTR, intron, …).
pub(crate) struct RegionIndex {
    trees: HashMap<String, COITree<(), u32>>,
    pub(crate) total_bases: i64,
}

impl RegionIndex {
    pub(crate) fn from_intervals(mut ivs: Vec<Iv>) -> Self {
        let merged = union_merge(&mut ivs);
        let total_bases: i64 = merged.iter().map(|iv| i64::from(iv.end - iv.start)).sum();

        let mut raw: HashMap<String, Vec<CoiInterval<()>>> = HashMap::new();
        for iv in &merged {
            raw.entry(iv.chrom.clone())
                .or_default()
                // coitrees uses end-inclusive coordinates
                .push(CoiInterval::new(iv.start, iv.end - 1, ()));
        }
        let trees = raw
            .into_iter()
            .map(|(chrom, intervals)| (chrom, COITree::new(&intervals)))
            .collect();

        Self { trees, total_bases }
    }

    /// Returns `true` if `point` (0-based) overlaps any interval in this set.
    pub(crate) fn contains(&self, chrom: &str, point: i32) -> bool {
        let Some(tree) = self.trees.get(chrom) else {
            return false;
        };
        let mut found = false;
        tree.query(point, point, |_node| {
            found = true;
        });
        found
    }
}
