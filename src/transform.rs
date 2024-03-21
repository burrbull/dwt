use num_traits::{NumAssign, Zero};

use crate::wavelet::Wavelet;

#[cfg(feature = "std")]
macro_rules! copy(
    ($source:ident, $destination:ident, $n:expr) => ({
        use core::ptr::copy_nonoverlapping as copy;
        unsafe { copy($source.as_ptr(), $destination.as_mut_ptr(), $n) };
    });
);

macro_rules! zero(
    ($buffer:expr) => ({
        use core::ptr::write_bytes as write;
        unsafe { write($buffer.as_mut_ptr(), 0, $buffer.len()) };
    });
);

impl<T, const L: usize> Wavelet<T, L>
where
    T: Zero + NumAssign + Copy,
{
    /// Perform the transform.
    ///
    /// The number of points should be divisible by `2^level`. If the operation
    /// is forward, the data are replaced by the approximation and detail
    /// coefficients stored in the first and second halves of `data`,
    /// respectively. If the operation is inverse, the data are assumed to be
    /// stored according to the above convention.
    #[cfg(feature = "std")]
    pub fn transform(&self, data: &mut [T], operation: crate::Operation, level: usize) {
        use crate::Operation;

        if level == 0 {
            return;
        }
        let n = data.len();
        assert!(n % (1 << level) == 0);
        let mut work = vec![T::zero(); n];

        match operation {
            Operation::Forward => {
                for i in 0..level {
                    let n = n >> i;
                    let (approx, detail) = work.split_at_mut(n >> 1);
                    self.forward_step(data, n, Some(approx), Some(detail));
                    copy!(work, data, n);
                }
            }
            Operation::Inverse => {
                for i in (0..level).rev() {
                    let n = n >> i;
                    let (approx, detail) = data.split_at(n >> 1);
                    self.inverse_step(approx, detail, n, &mut work);
                    copy!(work, data, n);
                }
            }
        }
    }

    #[inline(always)]
    pub fn forward_step(
        &self,
        data: &[T],
        n: usize,
        approx: Option<&mut [T]>,
        detail: Option<&mut [T]>,
    ) {
        let nm = L * n - self.offset;
        let nh = n >> 1;
        match (approx, detail) {
            (Some(approx), Some(detail)) => {
                for (i, (h, g)) in approx
                    .iter_mut()
                    .zip(detail.iter_mut())
                    .enumerate()
                    .take(nh)
                {
                    h.set_zero();
                    g.set_zero();
                    let mut k = 2 * i + nm;
                    for j in 0..L {
                        k += j;
                        let kn = k % n;
                        *h += self.dec_lo[j] * data[kn];
                        *g += self.dec_hi[j] * data[kn];
                    }
                }
            }
            (None, Some(detail)) => {
                for (i, g) in detail.iter_mut().enumerate().take(nh) {
                    g.set_zero();
                    let mut k = 2 * i + nm;
                    for j in 0..L {
                        k += j;
                        *g += self.dec_hi[j] * data[k % n];
                    }
                }
            }
            (Some(approx), None) => {
                for (i, h) in approx.iter_mut().enumerate().take(nh) {
                    h.set_zero();
                    let mut k = 2 * i + nm;
                    for j in 0..L {
                        k += j;
                        *h += self.dec_lo[j] * data[k % n];
                    }
                }
            }
            (None, None) => panic!("Pass destination buffer"),
        }
    }

    #[inline(always)]
    pub fn inverse_step(&self, approx: &[T], detail: &[T], n: usize, work: &mut [T]) {
        zero!(work);
        let nm = L * n - self.offset;
        let nh = n >> 1;
        for i in 0..nh {
            let (h, g) = (approx[i], detail[i]);
            let k = 2 * i + nm;
            for j in 0..L {
                let k = (k + j) % n;
                work[k] += self.rec_lo[j] * h + self.rec_hi[j] * g;
            }
        }
    }
}
