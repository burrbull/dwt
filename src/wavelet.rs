//! Wavelets.

use num_traits::float::{FloatConst, FloatCore};

mod coiflet;
mod daubechies;
mod symlet;

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

macro_rules! ortho {
    ($($struct:ident, $L:literal, $cnst:path, $doc:literal;)*) => {
        $(
            #[doc = $doc]
            pub struct $struct;
            impl $struct {
                pub fn with_offset(offset: usize) -> Wavelet<f64, $L> {
                    Wavelet::new_orthogonal(
                        $cnst.map(|b| f64::from_bits(b)).as_slice(),
                        offset,
                    )
                }
                pub fn new() -> Wavelet<f64, $L> {
                    Self::with_offset(0)
                }

                pub fn with_offset_f32(offset: usize) -> Wavelet<f32, $L> {
                    Wavelet::new_orthogonal(
                        $cnst.map(|b| f64::from_bits(b) as f32).as_slice(),
                        offset,
                    )
                }
                pub fn new_f32() -> Wavelet<f32, $L> {
                    Self::with_offset_f32(0)
                }
            }
        )*
    };
}

ortho! {
    Coif1, 6,  coiflet::COIF1, "Coiflet 1";
    Coif2, 12, coiflet::COIF2, "Coiflet 2";
    Coif3, 18, coiflet::COIF3, "Coiflet 3";
    Coif4, 24, coiflet::COIF4, "Coiflet 4";
    Coif5, 30, coiflet::COIF5, "Coiflet 5";

    Db2,  4,  daubechies::DB2, "Daubechies 2 wavelet";
    Db3,  6,  daubechies::DB3, "Daubechies 3 wavelet";
    Db4,  8,  daubechies::DB4, "Daubechies 4 wavelet";
    Db5,  10, daubechies::DB5, "Daubechies 5 wavelet";
    Db6,  12, daubechies::DB6, "Daubechies 6 wavelet";
    Db7,  14, daubechies::DB7, "Daubechies 7 wavelet";
    Db8,  16, daubechies::DB8, "Daubechies 8 wavelet";
    Db9,  18, daubechies::DB9, "Daubechies 9 wavelet";
    Db10, 20, daubechies::DB10, "Daubechies 10 wavelet";
    Db11, 22, daubechies::DB11, "Daubechies 11 wavelet";
    Db12, 24, daubechies::DB12, "Daubechies 12 wavelet";
    Db13, 26, daubechies::DB13, "Daubechies 13 wavelet";
    Db14, 28, daubechies::DB14, "Daubechies 14 wavelet";
    Db15, 30, daubechies::DB15, "Daubechies 15 wavelet";
    Db16, 32, daubechies::DB16, "Daubechies 16 wavelet";
    Db17, 34, daubechies::DB17, "Daubechies 17 wavelet";
    Db18, 36, daubechies::DB18, "Daubechies 18 wavelet";
    Db19, 38, daubechies::DB19, "Daubechies 19 wavelet";
    Db20, 40, daubechies::DB20, "Daubechies 20 wavelet";
    Db21, 42, daubechies::DB21, "Daubechies 21 wavelet";
    Db22, 44, daubechies::DB22, "Daubechies 22 wavelet";
    Db23, 46, daubechies::DB23, "Daubechies 23 wavelet";
    Db24, 48, daubechies::DB24, "Daubechies 24 wavelet";
    Db25, 50, daubechies::DB25, "Daubechies 25 wavelet";
    Db26, 52, daubechies::DB26, "Daubechies 26 wavelet";
    Db27, 54, daubechies::DB27, "Daubechies 27 wavelet";
    Db28, 56, daubechies::DB28, "Daubechies 28 wavelet";
    Db29, 58, daubechies::DB29, "Daubechies 29 wavelet";
    Db30, 60, daubechies::DB30, "Daubechies 30 wavelet";
    Db31, 62, daubechies::DB31, "Daubechies 31 wavelet";
    Db32, 64, daubechies::DB32, "Daubechies 32 wavelet";
    Db33, 66, daubechies::DB33, "Daubechies 33 wavelet";
    Db34, 68, daubechies::DB34, "Daubechies 34 wavelet";
    Db35, 70, daubechies::DB35, "Daubechies 35 wavelet";
    Db36, 72, daubechies::DB36, "Daubechies 36 wavelet";
    Db37, 74, daubechies::DB37, "Daubechies 37 wavelet";
    Db38, 76, daubechies::DB38, "Daubechies 38 wavelet";
    Db39, 78, daubechies::DB39, "Daubechies 39 wavelet";
    Db40, 80, daubechies::DB40, "Daubechies 40 wavelet";
    Db41, 82, daubechies::DB41, "Daubechies 41 wavelet";
    Db42, 84, daubechies::DB42, "Daubechies 42 wavelet";
    Db43, 86, daubechies::DB43, "Daubechies 43 wavelet";
    Db44, 88, daubechies::DB44, "Daubechies 44 wavelet";
    Db45, 90, daubechies::DB45, "Daubechies 45 wavelet";

    Sym2,  4,  symlet::SYM2, "Symlet 2";
    Sym3,  6,  symlet::SYM3, "Symlet 3";
    Sym4,  8,  symlet::SYM4, "Symlet 4";
    Sym5,  10, symlet::SYM5, "Symlet 5";
    Sym6,  12, symlet::SYM6, "Symlet 6";
    Sym7,  14, symlet::SYM7, "Symlet 7";
    Sym8,  16, symlet::SYM8, "Symlet 8";
    Sym9,  18, symlet::SYM9, "Symlet 9";
    Sym10, 20, symlet::SYM10, "Symlet 10";
    Sym11, 22, symlet::SYM11, "Symlet 11";
    Sym12, 24, symlet::SYM12, "Symlet 12";
    Sym13, 26, symlet::SYM13, "Symlet 13";
    Sym14, 28, symlet::SYM14, "Symlet 14";
    Sym15, 30, symlet::SYM15, "Symlet 15";
    Sym16, 32, symlet::SYM16, "Symlet 16";
    Sym17, 34, symlet::SYM17, "Symlet 17";
    Sym18, 36, symlet::SYM18, "Symlet 18";
    Sym19, 38, symlet::SYM19, "Symlet 19";
    Sym20, 40, symlet::SYM20, "Symlet 20";
    Sym21, 42, symlet::SYM21, "Symlet 21";
    Sym22, 44, symlet::SYM22, "Symlet 22";
    Sym23, 46, symlet::SYM23, "Symlet 23";
    Sym24, 48, symlet::SYM24, "Symlet 24";
    Sym25, 50, symlet::SYM25, "Symlet 25";
    Sym26, 52, symlet::SYM26, "Symlet 26";
    Sym27, 54, symlet::SYM27, "Symlet 27";
    Sym28, 56, symlet::SYM28, "Symlet 28";
    Sym29, 58, symlet::SYM29, "Symlet 29";
    Sym30, 60, symlet::SYM30, "Symlet 30";
    Sym31, 62, symlet::SYM31, "Symlet 31";
    Sym32, 64, symlet::SYM32, "Symlet 32";
    Sym33, 66, symlet::SYM33, "Symlet 33";
    Sym34, 68, symlet::SYM34, "Symlet 34";
}
