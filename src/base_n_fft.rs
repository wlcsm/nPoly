extern crate chrono;
extern crate num_complex;

use num_complex::Complex64;

use mathutils::*;
use crate::algebras::*;
use crate::algebras::complex::*;
use crate::algebras::integers::*;
use crate::polyu::*;
use std::f64::consts::PI;


pub fn fast_mult(polya: &PolyU<ZZ>, other: &PolyU<ZZ>) -> PolyU<ZZ> {

    let n = next_3pow(polya.deg() + other.deg() + 1);
    let mut a_sig = to_coeffs_complex(&polya.terms[..], n);
    let mut b_sig = to_coeffs_complex(&other.terms[..], n);

    eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

    // Because we converted it to complex for the ROU
    // We also need to normalise it here
    let c_parsed = a_sig.into_iter()
                        .map(|x| 
                            ZZ((x.0 / n as f64).re.round() as i32)
                        ).collect();

    // Convert back into polynomial type
    PolyU::from_coeff(None, c_parsed).unwrap()
}

fn rou(n: usize, inv: bool) -> Vec<CC> {
    // Generates all the nth roots of unity
    // Changes it depending on whether computing the dft or the inverse
    let sign = if inv { 1.0 } else { -1.0 };
    let base = Complex64::new(0.0, sign * 2.0 * PI / n as f64);
    (0..n).map(|k| CC(base.scale(k as f64).exp())).collect()
}

// Warning: Does it infix by default for speed, infix on a_sig
pub fn eval_interp(a_sig: &mut [CC], b_sig: &mut [CC]) -> Result<(), &'static str> {

    let n = a_sig.len();

    // Constraint checks
    if n != b_sig.len() || !is_3_pow(n){
        return Err("Improper lengths of input slices")
    }

    // Evaluate the polynomials
    perform_base_n_fft(a_sig, false)?;
    perform_base_n_fft(b_sig, false)?;

    // Multiply elementwise
    for i in 0..n {
        a_sig[i].mul_ass(&b_sig[i]);
    }

    // Interpolate the result
    perform_base_n_fft(a_sig, true)?;
    Ok(())
}

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

// fn to_coeffs(input: &[Monomial<ZZ>], n: usize) -> Vec<CC> {
//     // Expands the input into the expanded coefficient vector (coerced into complex)
//     // Then padded with zeros to length n

//     let mut result: Vec<ZZ> = Vec::with_capacity(n);
//     for mono in input.into_iter() {
//         // Fill the gap between monomials with zeros, then add the monomial
//         result.resize(mono.deg, <ZZ>::zero());
//         result.push(mono.coeff);
//     }
//     // Pad the rest
//     result.resize(n, <ZZ>::zero());
//     result
// }

fn perform_base_n_fft(signal: &mut [CC], inv: bool) -> Result<(), &'static str> {
    // Sample is the signal. Performs the FFT inline
    // Does not normalise the inverse!

    // Basic contraint checking
    match signal.len() {
        0 => Err("Signal cannot be empty"),
        1 => Ok(()),
        2 => {
            let x0 = signal[0];
            let x1 = signal[1];
            signal[0] = x0.add(&x1);
            signal[1] = x0.sub(&x1);
            Ok(())
        },
        n if !is_3_pow(n) => Err("Signal needs to be a power of three"),
        _ => {
            go_fast(signal, inv);
            Ok(())
        },
    }
}

pub fn dft(signal: &mut [CC], inv: bool) {

    let n = signal.len();
    let result = signal.to_vec();

    // Generate the nth roots of unity
    let rou = rou(n, inv);

    // Evaluates: p * q mod n.
    let modx = |p, q| (((p * q) % n) + n) % n;

    // Must normalise the output
    let scaling = if inv {1.0 / n as f64} else {1.0};

    // F(k) = \sum^n_{j=0} x_j e^{-2\pi i jk / n}
    for k in 0..n {
        let term: CC = result.iter().enumerate()
                                    .map(|(i,c)| rou[modx(k, i)].mul(c))
                                    .fold(<CC>::zero(), |a, b| a.add(&b));
        signal[k] = CC(term.0 * scaling);
    }
}

fn go_fast(sig: &mut [CC], inv: bool) {
    // Assumes that the length of 'signal' is >= 4

    let n   = sig.len();
    let rou = rou(n, inv); // Generates roots of unity

    // Does first iteration and puts into reverse bit order.
    let mut result = vec![<CC>::zero(); n];

    // TODO The problem is here, these are meant to be dispersed around the array, not grouped together.
    let i = 0;
    let j = n/3;
    let k = 2 * n / 3;
    for ind in (0..n/3).step_by(3) {
        for l in 0..3 {
            result[l * (n/3) + ind + 0] = sig[i + ind + l];
            result[l * (n/3) + ind + 1] = sig[j + ind + l];
            result[l * (n/3) + ind + 2] = sig[k + ind + l];
            dft(&mut result[l * (n/3) + ind .. l*(n/3) + ind + 3], inv);
        }
    }

    // We now assume the first layer is done and reverse bit order satisfied

    // Starts at index two because we already handled the first one
    for i in 2..=log3(n) {
        let flut = pow(3, i); // No. of elements in a flutter
        // Iterate over all the flutters, j is their starting index
        for j in (0..n).step_by(flut) {
            // Width of the k-flutter (number of elements)
            let a = 0;
            let b = flut/3;
            let c = 2 * flut / 3;
            let mut tmp = vec![<CC>::zero(); 3];
            for ind in 0..flut/3 {
                tmp[0] = result[j + a + ind];
                tmp[1] = rou[ind * (n / flut)].mul(&result[j + b + ind]);
                tmp[2] = rou[ind * 2 * (n / flut)].mul(&result[j + c + ind]);
                dft(&mut tmp[..], inv);
                result[j + a + ind] = tmp[0];
                result[j + b + ind] = tmp[1];
                result[j + c + ind] = tmp[2];
            }
        }
    }
    for (i, a) in sig.iter_mut().zip(result.into_iter()) {
        *i = a;
    }
}

// Don't think this abstraction is really necessary but I like the organisation
mod mathutils {
    // Does log2 for integers: Assumes its a power of two, or rather it rounds up
    pub fn log3(n: usize) -> usize { 
        match n {
            0 => panic!("Can't take log of zero"),
            1 | 2 => 0,
            n => log3(n / 3) + 1,
        }
    }

    // Quick test to see if a number is a power of two
    pub fn is_3_pow(n: usize) -> bool {
        fn aux(p: f64) -> bool {
            if p == 1.0 {
                true
            } else if (p > 1.0) & (p < 3.0) {
                false
            } else {
                aux(p / 3.0)
            }
        }
        if n == 0 { false } else { aux(n as f64) }
    }

    // Rounds up to the nearest power of three.
    pub fn next_3pow(n: usize) -> usize {
        if is_3_pow(n) { n } else { pow(3, log3(n) + 1) }
    }

    pub fn pow(n: usize, exp: usize) -> usize {
        match exp {
            0 => 1,
            1 => n,
            e => pow(n, e / 2) * pow(n, e / 2) * (if e % 2 == 0 {1} else {n})
        }
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
                    fast_mult(&a, &b);
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
    fn maths() {
        assert_eq!(true, is_3_pow(9));
        assert_eq!(true, is_3_pow(81));
        assert_eq!(false, is_3_pow(28));

        assert_eq!(0, log3(1));
        assert_eq!(1, log3(3));
        assert_eq!(2, log3(24));
        assert_eq!(4, log3(81));

        assert_eq!(16, pow(2, 4));
        assert_eq!(16, pow(4, 2));
        assert_eq!(81, pow(9, 2));

        assert_eq!(27, next_3pow(24));
        assert_eq!(81, next_3pow(29));
        assert_eq!(9, next_3pow(9));
    }
    #[test]
    fn mult_test() {
        let yeet = vec![Monomial::new(ZZ(1),0), Monomial::new(ZZ(2), 1), Monomial::new(ZZ(3), 2)];
        let mut b = to_coeffs_complex(&yeet[..], 9);
        let mut c = to_coeffs_complex(&yeet[..], 9);
        c.push(CC(Complex64::new(9.0, 0.0)));
        c.resize(27, <CC>::zero());

        // It can't do 3 yet
        // perform_base_n_fft(&mut a[..], false).unwrap();
        // println!("a: base3 = {:?}", a.iter().map(|x| x.0.re.round() as i32).collect::<Vec<i32>>());
        // perform_base_n_fft(&mut a[..], true).unwrap();
        // println!("a: base3inv = {:?}", a.iter().map(|x| x.0.re.round() as i32).collect::<Vec<i32>>());

        perform_base_n_fft(&mut b[..], false).unwrap();
        println!("b: base3 = {:?}", b.iter().map(|x| x.0.re.round() as i32).collect::<Vec<i32>>());
        perform_base_n_fft(&mut b[..], true).unwrap();
        println!("b: base3inv = {:?}", b.iter().map(|x| x.0.re.round() as i32).collect::<Vec<i32>>());

        perform_base_n_fft(&mut c[..], false).unwrap();
        println!("c: base3 = {:?}", c.iter().map(|x| x.0.re.round() as i32).collect::<Vec<i32>>());
        perform_base_n_fft(&mut c[..], true).unwrap();
        println!("c: base3inv = {:?}", c.iter().map(|x| x.0.re.round() as i32).collect::<Vec<i32>>());
    }
}
