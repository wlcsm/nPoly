// This only does the base 2 FFT. For an arbitrary base, see "fft-mullti-base.rs"
use crate::algebras::{ScalarField};
use crate::polynomials::*;
use num_traits::Zero;

// Note: Unchecked, assumes that n is a power of two
fn log2_unchecked(n: usize) -> usize {
    (n.trailing_zeros() + 1) as usize
}

pub trait SupportsFFT: ScalarField {
    // Generates the roots of unity
    fn rou(n: usize, inv: bool) -> Vec<Self>;
    fn divby2(&mut self, n: usize);
}

use num::Num;

impl<T: ScalarField + Num + PartialOrd + RemAssign> SupportsFFT for Complex<T> {
    fn rou(_n: usize, _inv: bool) -> Vec<Self> {
        unimplemented!()
    }
    fn divby2(&mut self, _n: usize) {
        unimplemented!()
    }
}

use std::ops::RemAssign;
impl<R: SupportsFFT + FFTnum + Div<usize, Output = R> + PartialOrd + RemAssign>
    Poly<Complex<R>, 1>
{
    /// FFT Multiplication
    pub fn fast_mult(&self, other: &Self) -> Self {
        let n = (self.tot_deg() + other.tot_deg() + 1).next_power_of_two();

        let mut a_sig = to_dense(self, Some(n));
        let mut b_sig = to_dense(other, Some(n));

        // Infix on a_sig
        test_eval_interp(&mut a_sig.terms[..], &mut b_sig.terms[..]).unwrap();

        // Convert back into polynomial type
        Poly::from_coeff(&a_sig.terms[..])
    }
}

/// Performs the evaluation-interpolation technique on two slices using the FFT.
/// Coefficient algebra must support a radix two FFT.
/// Warning: Mutates the original data.
pub fn eval_interp<F>(a_sig: &mut [F], b_sig: &mut [F]) -> Result<(), &'static str>
where
    F: SupportsFFT,
{
    let n = a_sig.len();

    // Constraint checks
    if n != b_sig.len() || !n.is_power_of_two() {
        return Err("Improper lengths of input slices");
    }

    // Evaluate the polynomials
    fft(a_sig, false)?;
    fft(b_sig, false)?;

    // Multiply elementwise
    for i in 0..n {
        a_sig[i] += b_sig[i];
    }

    // Interpolate the result
    fft(a_sig, true)?;

    // Need to normalise it here
    for x in a_sig.iter_mut() {
        x.divby2(log2_unchecked(n))
    }

    Ok(())
}

use num::complex::Complex;
use rustfft::{FFTnum, FFTplanner};
use std::ops::Div;

use num::integer::Roots;

pub fn test_eval_interp<T: FFTnum + Div<usize, Output = T>>(
    a_sig: &mut [Complex<T>],
    b_sig: &mut [Complex<T>],
) -> Result<(), &'static str> {
    let n = a_sig.len();

    let mut tmp: Vec<Complex<T>> = vec![Zero::zero(); n];

    let mut planner = FFTplanner::new(false);
    let fft = planner.plan_fft(n);

    fft.process(a_sig, &mut tmp);
    fft.process(b_sig, a_sig);

    // Multiply elementwise
    for i in 0..n {
        tmp[i] = tmp[i] + a_sig[i];
    }

    let mut planner = FFTplanner::new(true);
    let fft = planner.plan_fft(n);

    fft.process(&mut tmp, a_sig);

    // Normalise output
    let n_sqrt = n.sqrt();

    for x in a_sig.iter_mut() {
        *x = Complex::new(x.re / n_sqrt, x.im / n_sqrt)
    }

    Ok(())
}

/// Standard DFT Procedure
/// Warning: Performed infix
pub fn dft<F: SupportsFFT>(signal: &mut [F], inv: bool) {
    let n = signal.len();
    let result = signal.to_vec();

    // Generate the nth roots of unity
    let rou = <F>::rou(n, inv);

    // Evaluates: p * q mod n. Rust's modulo function gives negative numbers
    let modx = |p, q| (((p * q) % n) + n) % n;

    // F(k) = \sum^n_{j=0} x_j e^{-2\pi i jk / n}
    for k in 0..n {
        signal[k] = result
            .iter()
            .enumerate()
            .map(|(i, c)| rou[modx(k, i)] * *c)
            .fold(<F>::zero(), |a, b| a + b);
    }
}

/// Performs the standard Cooley-Tukey FFT.
/// If "inv" is true, then it will perform an inverse FFT
/// Advantage of this one is that it doesn't need to be a complex number, the coefficients simply
/// need to support an FFT in that they admit roots of unity.
fn fft<F: SupportsFFT>(signal: &mut [F], inv: bool) -> Result<(), &'static str> {
    match signal.len() {
        0 => Err("Signal cannot be empty"),
        n if !n.is_power_of_two() => Err("Signal needs to be a power of two"),
        2 | 4 => {
            dft(signal, inv);
            Ok(())
        }
        _ => {
            go_fast(signal, inv);
            Ok(())
        }
    }
}

/// Performs the classical FFT procedure in a bottom up fashion
fn go_fast<F: SupportsFFT>(sig: &mut [F], inv: bool) {
    let n = sig.len();
    let rou = <F>::rou(n, inv); // Generates roots of unity

    // Does first iteration and puts into reverse bit order.
    for (i, j) in (0..n / 2).zip(n / 2..n).step_by(2) {
        let x_0 = sig[i];
        let x_0_n2 = sig[j];
        let x_1 = sig[i + 1];
        let x_1_n2 = sig[j + 1];
        sig[i] = x_0 + x_0_n2;
        sig[j] = x_1 + x_1_n2;
        sig[i + 1] = x_0 - x_0_n2;
        sig[j + 1] = x_1 - x_1_n2;
    }
    // We now assume the first layer is done and reverse bit order satisfied

    // Starts at layer two because we already handled the first one
    for i in 2..log2_unchecked(n) {
        let flut = 1 << i; // No. of elements in a flutter, 2^i

        // Iterate over all the flutters, j is their starting index
        for j in (0..n).step_by(flut) {
            // Width of the k-flutter (number of elements)
            for k in j..(j + (flut / 2)) {
                let l = k + flut / 2;
                let a = sig[k];
                let rou_b = rou[(k % flut) * (n >> i)] * sig[l]; // w^j * b
                sig[k] = a + rou_b; // a + w^j * b
                sig[l] = a - rou_b; // a + w^{j + n/2} * b
            }
        }
    }
}

#[cfg(test)]
mod tests {
    //    extern crate chrono;
    //    extern crate rand;
    //
    //    use crate::algebras::polyring::*;
    //    use chrono::*;
    //
    //    use super::fft;
    //    use rand::distributions::{Distribution, Uniform};

    //    #[test]
    //    fn fft_test() {
    //        let mut sig = vec![
    //            CC::from_re(9),
    //            CC::from_re(3),
    //            CC::from_re(3),
    //            CC::from_re(6),
    //        ];
    //        let _ = fft(&mut sig[..], false);
    //
    //        println!("{:?}", sig)
    //    }

    //    // TODO I'm not sure this is correctly working, I changed the =log2_unchecked to log2_unchecked
    //    #[test]
    //    fn bench_dense_main() {
    //        let ring = PRDomain::<CC, UniVarOrder>::new(vec!['x']);
    //
    //        // Note all coefficients are nonzero
    //        let dist = Uniform::from(1..100);
    //        let mut rng = rand::thread_rng();
    //        // A function to randomly generate a polynomial with n coefficients
    //        let mut make_poly = |n: usize| -> PolyU<CC> {
    //            let res_vec = (0..n).map(|_| CC::from_re(dist.sample(&mut rng))).collect();
    //            Poly::from_coeff(&ring, res_vec)
    //        };
    //
    //        // Benches the time required to multiply two arbitrary polynomials of deg = n
    //        let mut time_mult = |n: usize| {
    //            let a = make_poly(n);
    //            let b = make_poly(n);
    //
    //            println!("-------------------------------------------");
    //            println!("Number of elements = {}", n);
    //            println!("-------------------------------------------");
    //            println!(
    //                "FFT: {:?}",
    //                Duration::span(|| {
    //                    a.fast_mult(&b);
    //                })
    //            );
    //            println!("-------------------------------------------");
    //        };
    //
    //        for i in 5..16 {
    //            time_mult(1 << i);
    //        }
    //    }

    // #[test]
    // fn mult_test_main() {
    //     // Postponing until I get the fast multiplication for ZZ working

    //     // let ring = PRDomain::univar(1);
    //     // let a = Poly::from_coeff(&ring, vec![ZZ(1), ZZ(1)]);
    //     // let b = Poly::from_coeff(&ring, vec![ZZ(1), ZZ(3)]);
    //     // let c = Poly::from_coeff(&ring, vec![ZZ(1), ZZ(2), ZZ(1)]);

    //     // assert_eq!(a.mul(&b), a.fast_mult(&b));
    //     // assert_eq!(b.mul(&c), b.fast_mult(&c));
    //     // assert_eq!(c.mul(&a), c.fast_mult(&a));

    //     // let d = Poly::from_coeff(ring, vec![ZZ(-1), ZZ(3)]);
    //     // let e = Poly::from_coeff(ring, vec![ZZ(-1), ZZ(3), ZZ(4), ZZ(6)]);

    //     // assert_eq!(d.mul(&a), d.fast_mult(&a));
    //     // assert_eq!(d.mul(&e), d.fast_mult(&e));

    //     // let f = Poly::from_coeff(ring, vec![ZZ(0)]);

    //     // assert_eq!(f.mul(&a), f.fast_mult(&a));
    // }
}
