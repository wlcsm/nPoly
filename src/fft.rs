extern crate chrono;
extern crate num_complex;

use num_complex::Complex64;

use mathutils::*;
use crate::algebras::*;
use crate::algebras::complex::*;
use crate::algebras::integers::*;
use crate::polyu::*;
use std::f64::consts::PI;

pub trait FastMult: Ring {
    fn fast_mult(&self, b: &Self, m: usize) -> Self;
}

pub trait SupportsFFT: ScalarRing {
    // Generates the roots of unity
    fn rou(n: usize, inv: bool) -> Vec<Self>;
    fn divby2(self, n: usize) -> Self;
}

// TODO these two implementations should be one macro
impl FastMult for PolyU<CC> {

    fn fast_mult(&self, other: &Self, m: usize) -> Self {

        let n = next_npow(self.deg() + other.deg() + 1, m);
        let mut a_sig = to_coeffs(&self.terms[..], n);
        let mut b_sig = to_coeffs(&other.terms[..], n);

        // Infix on a_sig
        eval_interp(&mut a_sig[..], &mut b_sig[..], m).unwrap();

        // Need to normalise it here
        for x in a_sig.iter_mut() {
            x.0 /= n as f64
        }

        // Convert back into polynomial type
        PolyU::from_coeff(None, a_sig).unwrap()
    }
}

impl FastMult for PolyU<ZZ> {

    fn fast_mult(&self, other: &Self, m: usize) -> Self {

        let n = next_npow(self.deg() + other.deg() + 1, m);
        let mut a_sig = to_coeffs_complex(&self.terms[..], n);
        let mut b_sig = to_coeffs_complex(&other.terms[..], n);

        eval_interp(&mut a_sig[..], &mut b_sig[..], m).unwrap();

        // Because we converted it to complex for the ROU
        // We also need to normalise it here
        let c_parsed = a_sig.into_iter()
                            .map(|x| 
                                ZZ((x.0 / n as f64).re.round() as i32)
                            ).collect();

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

    fn divby2(self, n: usize) -> Self {
        CC(self.0 / (1 << n) as f64)
    }
}

// Warning: Does it infix by default for speed, infix on a_sig
pub fn eval_interp<T>(a_sig: &mut [T], b_sig: &mut [T], m: usize) -> Result<(), &'static str> 
    where T: SupportsFFT {

    let n = a_sig.len();

    // Constraint checks
    if n != b_sig.len() || !is_n_pow(n, m) {
        return Err("Improper lengths of input slices")
    }

    // Evaluate the polynomials
    perform_fft(a_sig, false, m)?;
    perform_fft(b_sig, false, m)?;

    // Multiply elementwise
    for i in 0..n {
        a_sig[i].mul_ass(&b_sig[i]);
    }

    // Interpolate the result
    perform_fft(a_sig, true, m)?;
    Ok(())
}

// TODO these two should also be a macro
fn to_coeffs_complex(input: &[Monomial<ZZ>], n: usize) -> Vec<CC> {
    // Expands the input into the expanded coefficient vector (coerced into complex)
    // Then padded with zeros to length n

    let mut result: Vec<CC> = Vec::with_capacity(n);
    for mono in input.into_iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(mono.deg, CC::zero());
        result.push(CC(Complex64::new(mono.coeff.0 as f64, 0.0)));
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

pub fn dft<T: SupportsFFT>(signal: &mut [T], inv: bool) {

    let n = signal.len();
    let result = signal.to_vec();

    // Generate the nth roots of unity
    let rou = <T>::rou(n, inv);

    // Evaluates: p * q mod n.
    let modx = |p, q| (((p * q) % n) + n) % n;

    // F(k) = \sum^n_{j=0} x_j e^{-2\pi i jk / n}
    for k in 0..n {
        let term: T = result.iter().enumerate()
                                    .map(|(i,c)| rou[modx(k, i)].mul(c))
                                    .fold(<T>::zero(), |a, b| a.add(&b));
        signal[k] = term;
    }
}

pub fn perform_fft<T: SupportsFFT>(signal: &mut [T], inv: bool, base: usize) -> Result<(), &'static str> {
    // Sample is the signal. Performs the FFT inline
    // Does not normalise the inverse!

    // Basic contraint checking
    match signal.len() {
        0 => Err("Signal cannot be empty"),
        n if !is_n_pow(n, 2) => Err("Signal needs to be a power of the base"),
        p if p == base => {
            dft(signal, inv);
            Ok(())
        },
        _ => {
            base_2_fft::go_fast(signal, inv);
            Ok(())
        },
    }
}

mod base_2_fft {

    use super::*;

    pub fn go_fast<T: SupportsFFT>(signal: &mut [T], inv: bool) {
        // Assumes that the length of 'signal' is >= 4

        let n   = signal.len();
        let rou = <T>::rou(n, inv); // Generates roots of unity

        // Does first iteration and puts into reverse bit order.
        for (i, j) in (0..n / 2).zip(n / 2..n).step_by(2) {
            let x_0       = signal[i];
            let x_0_n2    = signal[j];
            let x_1       = signal[i + 1];
            let x_1_n2    = signal[j + 1];
            signal[i]     = x_0.add(&x_0_n2);
            signal[j]     = x_1.add(&x_1_n2);
            signal[i + 1] = x_0.sub(&x_0_n2);
            signal[j + 1] = x_1.sub(&x_1_n2);
        }
        // We now assume the first layer is done and reverse bit order satisfied

        // Starts at index two because we already handled the first one
        for i in 2..=logn(n, 2) {
            let flut = 1 << i; // No. of elements in a flutter
            // Iterate over all the flutters, j is their starting index
            for j in (0..n).step_by(flut) {
                // Width of the k-flutter (number of elements)
                for (k, l) in (j..j + (flut >> 1)).zip(j + (flut >> 1)..j + flut) {
                    let a = signal[k];
                    let rou_b = rou[(k % flut) * (n >> i)].mul(&signal[l]); // w^j * b
                    signal[k] = a.add(&rou_b); // a + w^j * b
                    signal[l] = a.sub(&rou_b); // a + w^{j + n/2} * b
                }
            }
        }
    }
}


pub mod base_n_fft {
    use super::*;

    pub fn go_fast_base_n<T: SupportsFFT>(sig: &mut [T], inv: bool, m: usize) {
        // Assumes that the length of 'signal' is >= 4

        let n   = sig.len();
        let rou = <T>::rou(n, inv); // Generates roots of unity

        // Does first iteration and puts into reverse bit order.
        let mut result = vec![<T>::zero(); n];

        for i in (0..n/m).step_by(m) {
            for l in 0..m {
                for p in 0..m {
                    result[l * (n/m) + i + p] = sig[(p * n/m) + i + l];
                }
                dft(&mut result[l * (n/m) + i .. l*(n/m) + i + m], inv);
            }
        }

        // We now assume the first layer is done and reverse bit order satisfied

        // Starts at index two because we already handled the first one
        for i in 2..=logn(n, m) {
            let flut = pow(m, i); // No. of elements in a flutter
            // Iterate over all the flutters, j is their starting index
            for j in (0..n).step_by(flut) {
                // Width of the k-flutter (number of elements)
                let mut tmp = vec![<T>::zero(); m];
                for ind in 0..flut/m {

                    for p in 0..m {
                        tmp[p] = rou[ind * p* (n / flut)].mul(&result[j + p * (flut/m) + ind]);
                    }

                    go_fast_base_n(&mut tmp[..], inv, m);

                    for (p, c) in tmp.iter().enumerate() {
                        result[j + p * (flut/m) + ind] = *c;
                    }
                }
            }
        }
        for (i, a) in sig.iter_mut().zip(result.into_iter()) {
            *i = a;
        }
    }
}

// Don't think this abstraction is really necessary but I like the organisation
// I didn't idiot proof these so it is easy to break them with edge cases
pub mod mathutils {
    // Rounds down
    pub fn logn(num: usize, n: usize) -> usize { 
        if num < n { 0 } else {logn(num / n, n) + 1}
    }

    pub fn is_n_pow(num: usize, m: usize) -> bool {
        if m == 2{
            (num - 1) & num == 0
        } else {
            pow(m, logn(num, m)) == num 
        }
    }

    pub fn next_npow(num: usize, n: usize) -> usize {
        if is_n_pow(num, n) { num } else { pow(n, logn(num, n) + 1) }
    }

    // Exponentiation function: n^exp
    pub fn pow(n: usize, exp: usize) -> usize {
        if n == 2 {
            1 << exp
        } else {
            match exp {
                0 => 1,
                1 => n,
                e => pow(n, e / 2) * pow(n, e / 2) * (if e % 2 == 0 {1} else {n})
            }
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate chrono;
    extern crate rand;
    use super::*;
    use chrono::*;

    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn bench_dense_main() {
        // Note all coefficients are nonzero
        let between = Uniform::from(1..100);
        let mut rng = rand::thread_rng();
        // A function to randomly generate a polynomial with n coefficients
        let mut make_poly = |n: usize| -> PolyU<ZZ> {
            let res_vec = (0..n).map(|_| ZZ(between.sample(&mut rng))).collect();
            PolyU::from_coeff(None, res_vec).unwrap()
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
                    // a.fast_mult(&b, 2);
                    a.mul(&b);
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
        time_mult(1 << 14);
        time_mult(1 << 15);
    }

    #[test]
    fn mult_test_main() {
        let a = PolyU::from_coeff(None, vec![ZZ(1), ZZ(1)]).unwrap();
        let b = PolyU::from_coeff(None, vec![ZZ(1), ZZ(3)]).unwrap();
        let c = PolyU::from_coeff(None, vec![ZZ(1), ZZ(2), ZZ(1)]).unwrap();

        assert_eq!(a.mul(&b), a.fast_mult(&b, 2));
        assert_eq!(b.mul(&c), b.fast_mult(&c, 2));
        assert_eq!(c.mul(&a), c.fast_mult(&a, 2));

        let d = PolyU::from_coeff(None, vec![ZZ(-1), ZZ(3)]).unwrap();
        let e = PolyU::from_coeff(None, vec![ZZ(-1), ZZ(3), ZZ(4), ZZ(6)]).unwrap();

        assert_eq!(d.mul(&a), d.fast_mult(&a, 2));
        assert_eq!(d.mul(&e), d.fast_mult(&e, 2));

        let f = PolyU::from_coeff(None, vec![ZZ(0)]).unwrap();

        assert_eq!(f.mul(&a), f.fast_mult(&a, 2));
    }
}
