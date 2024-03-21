//! [Discrete wavelet transform][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Discrete_wavelet_transform

// The implementation is based on:
// http://www.gnu.org/software/gsl/manual/html_node/Wavelet-Transforms.html
#![cfg_attr(not(feature = "std"), no_std)]
#![allow(clippy::new_ret_no_self)]

mod transform;

pub mod wavelet;

/// A transform operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Operation {
    /// The forward transform.
    Forward,
    /// The inverse transform.
    Inverse,
}

/// Perform the transform.
///
/// The function is a shortcut for `Transform::transform`.
#[inline(always)]
#[cfg(feature = "std")]
pub fn transform<T, const L: usize>(
    data: &mut [T],
    operation: Operation,
    wavelet: &wavelet::Wavelet<T, L>,
    level: usize,
) where
    T: num_traits::Zero + num_traits::NumAssign + Copy,
{
    wavelet.transform(data, operation, level);
}
