extern crate chrono;
extern crate num_complex;
extern crate num_traits;

use num_complex::Complex64;
use num_traits::identities::Zero;

use crate::*;
use mathutils::*;

pub trait Char0: From<f64> {
}

pub fn fast_mult(polya: &PolyU<i32>, polyb: &PolyU<i32>) -> Result<PolyU<i32>, &'static str> {
    // At the moment, need to round up to the nearest power of two
    let n = next_2pow(polya.deg() + polyb.deg() + 1);

    // Expand the monomial vectors into a complex coefficient vector
    let a_sig = &mut to_coeffs_complex(&polya.terms[..], n)[..];
    let b_sig = &mut to_coeffs_complex(&polyb.terms[..], n)[..];

    // Evaluate the polynomials
    perform_fft(a_sig, false)?;
    perform_fft(b_sig, false)?;

    // Multiply elementwise
    for i in 0..n {
        a_sig[i] *= b_sig[i];
    }

    // Interpolate the result
    perform_fft(a_sig, true)?;

    // Convert back into i32
    let c_parsed = a_sig.into_iter().map(|x| x.re.round() as i32).collect();

    // Convert back into polynomial type
    Ok(PolyU::<i32>::from_coeff("x".to_string(), c_parsed).unwrap())
}

fn to_coeffs_complex(input: &[Monomial<i32>], n: usize) -> Vec<Complex64> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<Complex64> = Vec::with_capacity(n);
    for mono in input.into_iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(mono.deg, Complex64::zero());
        result.push(Complex64::from(mono.coeff as f64));
    }
    // Pad the rest
    result.resize(n, Complex64::zero());
    result
}

fn perform_fft(signal: &mut [Complex64], inv: bool) -> Result<(), &'static str> {
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

fn go_fast(signal: &mut [Complex64], inv: bool) {
    // Assumes that the length of 'signal' is >= 4

    let n   = signal.len();
    let rou = gen_rou(n, inv);

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
            *el /= n as f64;
        }
    }
}

// Don't think this abstraction is really necessary but I like the organisation
mod mathutils {
    use num_complex::Complex64;
    use std::f64::consts::PI;

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

    pub fn gen_rou(n: usize, inverse: bool) -> Vec<Complex64> {
        // Generates all the nth roots of unity
        // Changes it depending on whether computing the dft or the inverse
        let sign = if inverse { 1.0 } else { -1.0 };
        let base = Complex64::new(0.0, sign * 2.0 * PI / n as f64);
        (0..n).map(|k| {base.scale(k as f64).exp() }).collect()
    }
}

#[cfg(test)]
mod tests {
    extern crate chrono;
    extern crate rand;
    use super::fft::*;
    use chrono::*;
    use polyu::*;

    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn bench_dense() {
        // Note all coefficients are nonzero
        let between = Uniform::from(1..100);
        let mut rng = rand::thread_rng();
        // A function to randomly generate a polynomial with n coefficients
        let mut make_poly = |n: usize| -> PolyU<i32> {
            let res_vec = (0..n).map(|_| between.sample(&mut rng)).collect();
            PolyU::<i32>::from_coeff("x".to_string(), res_vec).unwrap()
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
                    fast_mult(&a, &b).unwrap();
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
        let a = PolyU::<i32>::from_coeff("x".to_string(), vec![1, 1]).unwrap();
        let b = PolyU::<i32>::from_coeff("x".to_string(), vec![1, 3]).unwrap();
        let c = PolyU::<i32>::from_coeff("x".to_string(), vec![1, 2, 1]).unwrap();

        assert_eq!(a.mul(&b), fast_mult(&a, &b).unwrap());
        assert_eq!(b.mul(&c), fast_mult(&b, &c).unwrap());
        assert_eq!(a.mul(&c), fast_mult(&a, &c).unwrap());

        let d = PolyU::<i32>::from_coeff("x".to_string(), vec![-1, 3]).unwrap();
        let e = PolyU::<i32>::from_coeff("x".to_string(), vec![-1, 3, 4, 6]).unwrap();

        assert_eq!(d.mul(&a), fast_mult(&d, &a).unwrap());
        assert_eq!(d.mul(&e), fast_mult(&d, &e).unwrap());

        let f = PolyU::<i32>::from_coeff("x".to_string(), vec![0]).unwrap();

        assert_eq!(a.mul(&f), fast_mult(&a, &f).unwrap());
    }
}
