//! Crate-internal input validation predicates.
//!
//! A plain `x <= 0.0` guard silently accepts NaN, because every comparison
//! involving NaN is `false`. These predicates return `true` only for finite
//! values that satisfy the bound, so `if !finite_gt(x, 0.0) { return Err(..) }`
//! rejects negative values, zero, NaN and ±∞ in one test.
//!
//! Not every feature uses every predicate, hence the module-level allow.
#![allow(dead_code)]

use crate::error::UshmaError;

/// Build an [`UshmaError::InvalidParameter`] out of line.
///
/// Keeping `format!` machinery off the hot path lets small validated
/// functions stay inlinable: `Err(invalid_param(format_args!("…")))`.
#[cold]
#[inline(never)]
pub(crate) fn invalid_param(reason: std::fmt::Arguments<'_>) -> UshmaError {
    UshmaError::InvalidParameter {
        reason: reason.to_string(),
    }
}

/// `true` when `x` is finite and strictly greater than `bound`.
#[inline]
#[must_use]
pub(crate) fn finite_gt(x: f64, bound: f64) -> bool {
    x.is_finite() && x > bound
}

/// `true` when `x` is finite and greater than or equal to `bound`.
#[inline]
#[must_use]
pub(crate) fn finite_ge(x: f64, bound: f64) -> bool {
    x.is_finite() && x >= bound
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_gt_rejects_nan_and_inf() {
        assert!(finite_gt(1.0, 0.0));
        assert!(!finite_gt(0.0, 0.0));
        assert!(!finite_gt(-1.0, 0.0));
        assert!(!finite_gt(f64::NAN, 0.0));
        assert!(!finite_gt(f64::INFINITY, 0.0));
        assert!(!finite_gt(f64::NEG_INFINITY, 0.0));
        assert!(!finite_gt(1.0, f64::NAN));
    }

    #[test]
    fn finite_ge_accepts_bound() {
        assert!(finite_ge(0.0, 0.0));
        assert!(finite_ge(2.0, 1.0));
        assert!(!finite_ge(-1e-300, 0.0));
        assert!(!finite_ge(f64::NAN, 0.0));
        assert!(!finite_ge(f64::INFINITY, 0.0));
    }
}
