#![feature(test)]

extern crate test;

use core::hint::black_box;
use dwt::{wavelet, Operation};
use test::Bencher;

#[bench]
fn forward_haar_0004(bencher: &mut Bencher) {
    forward_haar(4, bencher);
}
#[bench]
fn forward_haar_0016(bencher: &mut Bencher) {
    forward_haar(16, bencher);
}
#[bench]
fn forward_haar_0064(bencher: &mut Bencher) {
    forward_haar(64, bencher);
}
#[bench]
fn forward_haar_0256(bencher: &mut Bencher) {
    forward_haar(256, bencher);
}
#[bench]
fn forward_haar_1024(bencher: &mut Bencher) {
    forward_haar(1024, bencher);
}
#[bench]
fn forward_haar_4096(bencher: &mut Bencher) {
    forward_haar(4096, bencher);
}

fn forward_haar(size: usize, bencher: &mut Bencher) {
    let mut data = vec![42.0; size];
    let operation = Operation::Forward;
    let wavelet = wavelet::Db4::new();
    let level = (size as f64).log2() as usize;
    bencher.iter(|| black_box(wavelet.transform(&mut data, operation, level)));
}
