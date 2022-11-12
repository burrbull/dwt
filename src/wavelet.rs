//! Wavelets.

use num_traits::float::{FloatConst, FloatCore};

/// A wavelet.
///
/// `T` - float type.
/// `L` - the number of coefficients.
pub struct Wavelet<T, const L: usize> {
    /// The offset of the coefficients.
    pub offset: usize,

    /// The coefficients of the decomposition low-pass filter.
    pub dec_lo: [T; L],
    /// The coefficients of the decomposition high-pass filter.
    pub dec_hi: [T; L],

    /// The coefficients of the reconstruction low-pass filter.
    pub rec_lo: [T; L],
    /// The coefficients of the reconstruction high-pass filter.
    pub rec_hi: [T; L],
}

impl<T, const L: usize> Wavelet<T, L>
where
    T: FloatCore + FloatConst,
{
    /// Create an orthogonal wavelet.
    pub fn new_orthogonal(declo: &[T], offset: usize) -> Self {
        let dec_lo = core::array::from_fn(|i| declo[i]);
        let rec_lo = core::array::from_fn(|i| dec_lo[L - i - 1]);
        let rec_hi = core::array::from_fn(|i| {
            let c = declo[i];
            if i % 2 == 0 {
                c
            } else {
                -c
            }
        });
        let dec_hi = core::array::from_fn(|i| rec_hi[L - i - 1]);

        Self {
            offset,
            dec_lo,
            rec_lo,
            rec_hi,
            dec_hi,
        }
    }
}

/// A Haar wavelet.
pub struct Haar;

impl Haar {
    /// Create a wavelet.
    pub fn new<T>() -> Wavelet<T, 2>
    where
        T: FloatCore + FloatConst,
    {
        let value = T::FRAC_1_SQRT_2();
        Wavelet {
            offset: 0,
            dec_lo: [value, value],
            dec_hi: [value, -value],
            rec_lo: [value, value],
            rec_hi: [value, -value],
        }
    }
}
