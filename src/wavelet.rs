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
