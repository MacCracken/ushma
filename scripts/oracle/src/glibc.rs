//! The libm the oracle's transcendentals come from.
//!
//! ushma 1.4.0 calls the std f64 methods (`x.exp()`, `x.ln()`, `x.powf(y)`, …), which on Linux
//! call the host **glibc** libm. So every transcendental bit in the corpus is glibc's, and the
//! corpus is only reproducible on the same glibc. (`ceil`, `powi` and the inline methods are the
//! exception, see `m_math.rs`: they are exact or a fixed operation order, so their bits are
//! platform-free.) glibc also picks an implementation per CPU at load time (ifunc): on x86_64
//! its exp, log, pow, atan, tanh, expm1, … have FMA variants, chosen when the CPU has FMA and
//! AVX2, and their results can differ in rare cases, so the CPU class is part of the provenance.
//!
//! `platform()` is written into every TSV header, and `--check` refuses to compare a corpus
//! recorded on a different platform (`main.rs`). Copied from kimiya's `scripts/oracle/`.

use std::ffi::{CStr, c_char};

unsafe extern "C" {
    /// glibc's own version string, e.g. `"2.44"` (`<gnu/libc-version.h>`).
    fn gnu_get_libc_version() -> *const c_char;
}

/// The running glibc's version.
pub fn version() -> String {
    // SAFETY: gnu_get_libc_version takes no arguments and returns a pointer to a static,
    // NUL-terminated string that lives for the whole process.
    let p = unsafe { gnu_get_libc_version() };
    if p.is_null() {
        return "unknown".to_string();
    }
    // SAFETY: non-null, NUL-terminated and static, as above.
    unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned()
}

/// The ifunc class glibc's x86_64 libm selects: `fma` when the CPU has FMA and AVX2.
#[cfg(target_arch = "x86_64")]
fn cpu_class() -> &'static str {
    if std::arch::is_x86_feature_detected!("fma") && std::arch::is_x86_feature_detected!("avx2") {
        "fma"
    } else {
        "no-fma"
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn cpu_class() -> &'static str {
    "generic"
}

/// `glibc 2.44 x86_64 fma`: the platform every transcendental row was computed on.
pub fn platform() -> String {
    format!(
        "glibc {} {} {}",
        version(),
        std::env::consts::ARCH,
        cpu_class()
    )
}
