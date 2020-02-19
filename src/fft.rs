extern crate chrono;
extern crate num_complex;

use num_complex::Complex64;

use mathutils::*;
use crate::algebras::*;
use crate::algebras::complex::*;
use crate::polyu::*;
use std::f64::consts::PI;

pub trait FastMult: Group {
    fn fast_mult(self, b: Self) -> Self;
}

pub trait SupportsFFT: Ring + Copy {
    // Generates the roots of unity
    fn rou(n: usize, inv: bool) -> Vec<Self>;
}

impl FastMult for PolyU<CC> {

    fn fast_mult(self, other: Self) -> Self {

        let n = next_2pow(self.deg() + other.deg() + 1);
        let mut a_sig = to_coeffs(&self.terms[..], n);
        let mut b_sig = to_coeffs(&other.terms[..], n);

        // Infix on a_sig
        eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

        // Convert back into polynomial type
        PolyU::from_coeff(None, a_sig).unwrap()
    }
}

impl FastMult for PolyU<i32> {

    fn fast_mult(self, other: Self) -> Self {

        let n = next_2pow(self.deg() + other.deg() + 1);
        let mut a_sig = to_coeffs_complex(&self.terms[..], n);
        let mut b_sig = to_coeffs_complex(&other.terms[..], n);

        eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

        // Because we converted it to complex for the ROU
        let c_parsed = a_sig.into_iter().map(|x| x.0.re.round() as i32).collect();

        // Convert back into polynomial type
        PolyU::from_coeff(None, c_parsed).unwrap()
    }
}

impl SupportsFFT for CC {

    fn rou(n: usize, inv: bool) -> Vec<Self> {
        // Generates all the nth roots of unity
        // Changes it depending on whether computing the dft or the inverse
        let sign = if inv { 1.0 } else { -1.0 };
        let base = Complex64::new(0.0, sign * 2.0 * PI / n as f64);
        (0..n).map(|k| CC(base.scale(k as f64).exp())).collect()
    }
}

// Warning: Does it infix by default for speed, infix on a_sig
pub fn eval_interp<T>(a_sig: &mut [T], b_sig: &mut [T]) -> Result<(), &'static str> 
    where T: SupportsFFT {

    let n = a_sig.len();

    // Constraint checks
    if n != b_sig.len() || !is_2_pow(n){
        return Err("Improper lengths of input slices")
    }

    // Evaluate the polynomials
    perform_fft(a_sig, false)?;
    perform_fft(b_sig, false)?;

    // Multiply elementwise
    for i in 0..n {
        a_sig[i] *= b_sig[i];
    }

    // Interpolate the result
    perform_fft(a_sig, true)?;
    Ok(())
}

fn to_coeffs_complex(input: &[Monomial<i32>], n: usize) -> Vec<CC> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<CC> = Vec::with_capacity(n);
    for mono in input.into_iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(mono.deg, CC::zero());
        result.push(CC(Complex64::new(mono.coeff as f64, 0.0)));
    }
    // Pad the rest
    result.resize(n, CC::zero());
    result
}

fn to_coeffs<T: SupportsFFT>(input: &[Monomial<T>], n: usize) -> Vec<T> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<T> = Vec::with_capacity(n);
    for mono in input.into_iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(mono.deg, <T>::zero());
        result.push(mono.coeff);
    }
    // Pad the rest
    result.resize(n, <T>::zero());
    result
}

fn perform_fft<T: SupportsFFT>(signal: &mut [T], inv: bool) -> Result<(), &'static str> {
    // Sample is the signal. Note that is fully expanded,
    // Performs the FFT inline

    // Basic contraint checking
    match signal.len() {
        0 => Err("Signal cannot be empty"),
        1 => Ok(()),
        2 => {
            let x0 = signal[0];
            let x1 = signal[1];
            signal[0] = x0 + x1;
            signal[1] = x0 - x1;
            Ok(())
        },
        n if !is_2_pow(n) => Err("Signal needs to be a power of two"),
        _ => {
            go_fast(signal, inv);
            Ok(())
        },
    }
}

fn go_fast<T: SupportsFFT>(signal: &mut [T], inv: bool) {
    // Assumes that the length of 'signal' is >= 4

    let n   = signal.len();
    let rou = <T>::rou(n, inv); // Generates roots of unity

    // Does first iteration and puts into reverse bit order.
    for (i, j) in (0..n / 2).zip(n / 2..n).step_by(2) {
        let x_0       = signal[i];
        let x_0_n2    = signal[j];
        let x_1       = signal[i + 1];
        let x_1_n2    = signal[j + 1];
        signal[j]     = x_1 + x_1_n2;
        signal[i]     = x_0 + x_0_n2;
        signal[i + 1] = x_0 - x_0_n2;
        signal[j + 1] = x_1 - x_1_n2;
    }
    // We now assume the first layer is done and reverse bit order satisfied

    // Starts at index two because we already handled the first one
    for i in 2..=log2(n) {
        let flut = 1 << i; // No. of elements in a flutter
        // Iterate over all the flutters, j is their starting index
        for j in (0..n).step_by(flut) {
            // Width of the k-flutter (number of elements)
            for (k, l) in (j..j + (flut >> 1)).zip(j + (flut >> 1)..j + flut) {
                let a = signal[k];
                let rou_b = rou[(k % flut) * (n >> i)] * signal[l]; // w^j * b
                signal[k] = a + rou_b; // a + w^j * b
                signal[l] = a - rou_b; // a + w^{j + n/2} * b
            }
        }
    }
    // Normalize if calculating the inverse DFT
    if inv {
        for el in signal.iter_mut() {
            *el /= n as i32;
        }
    }
}

// Don't think this abstraction is really necessary but I like the organisation
mod mathutils {
    // Does log2 for integers: Assumes its a power of two, or rather it rounds down
    pub fn log2(n: usize) -> usize { 
        if n == 1 {0} else {log2(n >> 1) + 1}
    }

    // Quick test to see if a number is a power of two
    pub fn is_2_pow(n: usize) -> bool {
        (n - 1) & n == 0
    }

    // Rounds up to the nearest power of two.
    pub fn next_2pow(n: usize) -> usize {
        if is_2_pow(n) { n } else { 1 << (log2(n) + 1) }
    }
}

#[cfg(test)]
mod tests {
    extern crate chrono;
    extern crate rand;
    // use crate::polyu::*;
    use super::*;
    use chrono::*;

    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn bench_dense() {
        // Note all coefficients are nonzero
        let between = Uniform::from(1..100);
        let mut rng = rand::thread_rng();
        // A function to randomly generate a polynomial with n coefficients
        let mut make_poly = |n: usize| -> PolyU<i32> {
            let res_vec = (0..n).map(|_| between.sample(&mut rng)).collect();
            PolyU::<i32>::from_coeff(None, res_vec).unwrap()
        };

        // Benches the time required to multiply two arbitrary polynomials of deg = n
        let mut time_mult = |n: usize| {
            let a = make_poly(n);
            let b = make_poly(n);

            println!("-------------------------------------------");
            println!("Number of elements = {}", n);
            println!("-------------------------------------------");
            println!("FFT: {:?}",
                Duration::span(|| {
                    a.fast_mult(b);
                })
            );
            println!("-------------------------------------------");
        };

        time_mult(1 << 3);
        time_mult(1 << 4);
        time_mult(1 << 6);
        time_mult(1 << 8);
        time_mult(1 << 10);
        time_mult(1 << 11);
        time_mult(1 << 12);
        time_mult(1 << 13);
    }

    #[test]
    fn mult_test() {
        let a = PolyU::<i32>::from_coeff(None, vec![1, 1]).unwrap();
        let b = PolyU::<i32>::from_coeff(None, vec![1, 3]).unwrap();
        let c = PolyU::<i32>::from_coeff(None, vec![1, 2, 1]).unwrap();

        assert_eq!(a.clone() * b.clone(), a.clone().fast_mult(b.clone()));
        assert_eq!(b.clone() * c.clone(), b.clone().fast_mult(c.clone()));
        assert_eq!(c.clone() * a.clone(), c.clone().fast_mult(a.clone()));

        let d = PolyU::<i32>::from_coeff(None, vec![-1, 3]).unwrap();
        let e = PolyU::<i32>::from_coeff(None, vec![-1, 3, 4, 6]).unwrap();

        assert_eq!(d.clone() * a.clone(), d.clone().fast_mult(a.clone()));
        assert_eq!(d.clone() * e.clone(), d.clone().fast_mult(e.clone()));

        let f = PolyU::<i32>::from_coeff(None, vec![0]).unwrap();

        assert_eq!(f.clone() * a.clone(), f.clone().fast_mult(a.clone()));
    }
}
