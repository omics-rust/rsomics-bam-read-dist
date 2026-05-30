/// One genomic interval [start, end) (half-open, 0-based).
#[derive(Clone, Debug)]
pub(crate) struct Iv {
    pub(crate) chrom: String,
    pub(crate) start: i32,
    pub(crate) end: i32,
}

/// Sort and merge overlapping/adjacent intervals; returns the merged list.
pub(crate) fn union_merge(ivs: &mut [Iv]) -> Vec<Iv> {
    if ivs.is_empty() {
        return Vec::new();
    }
    ivs.sort_unstable_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
    let mut merged: Vec<Iv> = Vec::with_capacity(ivs.len());
    for iv in ivs.iter() {
        if let Some(last) = merged.last_mut()
            && last.chrom == iv.chrom
            && iv.start <= last.end
        {
            last.end = last.end.max(iv.end);
            continue;
        }
        merged.push(iv.clone());
    }
    merged
}

/// Set-difference of two already-merged interval lists.
pub(crate) fn subtract_sorted(base: Vec<Iv>, minus: &[Iv]) -> Vec<Iv> {
    if minus.is_empty() {
        return base;
    }
    let mut result = Vec::with_capacity(base.len());
    for iv in base {
        let relevant: Vec<&Iv> = minus
            .iter()
            .filter(|m| m.chrom == iv.chrom && m.end > iv.start && m.start < iv.end)
            .collect();
        if relevant.is_empty() {
            result.push(iv);
            continue;
        }
        let mut cursor = iv.start;
        for m in &relevant {
            if m.start > cursor {
                result.push(Iv {
                    chrom: iv.chrom.clone(),
                    start: cursor,
                    end: m.start.min(iv.end),
                });
            }
            cursor = cursor.max(m.end);
            if cursor >= iv.end {
                break;
            }
        }
        if cursor < iv.end {
            result.push(Iv {
                chrom: iv.chrom.clone(),
                start: cursor,
                end: iv.end,
            });
        }
    }
    result
}
