use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bam_read_dist(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bam-read-dist");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/large_reads.bam");
    let bed = manifest.join("tests/golden/genes.bed12");
    c.bench_function("rsomics-bam-read-dist golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["-i", bam.to_str().unwrap(), "-r", bed.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bam_read_dist);
criterion_main!(benches);
