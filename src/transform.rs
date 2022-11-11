use num_traits::{float::FloatCore, NumAssign};

use crate::wavelet::Wavelet;

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
    T: FloatCore + NumAssign,
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
                    self.forward_step(data, n, &mut work);
                    copy!(work, data, n);
                }
            }
            Operation::Inverse => {
                for i in (0..level).rev() {
                    let n = n >> i;
                    self.inverse_step(data, n, &mut work);
                    copy!(work, data, n);
                }
            }
        }
    }

    #[inline(always)]
    pub fn forward_step(&self, data: &[T], n: usize, work: &mut [T])
    where
        T: FloatCore,
    {
        let nm = L * n - self.offset;
        let nh = n >> 1;
        for i in 0..nh {
            let (mut h, mut g) = (T::zero(), T::zero());
            let k = 2 * i + nm;
            for j in 0..L {
                let k = (k + j) % n;
                h += self.dec_lo[j] * data[k];
                g += self.dec_hi[j] * data[k];
            }
            work[i] = h;
            work[i + nh] = g;
        }
    }

    #[inline(always)]
    pub fn inverse_step(&self, data: &[T], n: usize, work: &mut [T])
    where
        T: FloatCore,
    {
        let nm = L * n - self.offset;
        let nh = n >> 1;
        for i in 0..nh {
            let (h, g) = (data[i], data[i + nh]);
            let k = 2 * i + nm;
            for j in 0..L {
                let k = (k + j) % n;
                work[k] += self.rec_lo[j] * h + self.rec_hi[j] * g;
            }
        }
    }
}
