extern crate num_complex;
extern crate num_traits;
extern crate chrono;

use num_complex::Complex64;
use num_traits::identities::Zero;

use crate::*;
use mathutils::*;


// Evaluate Interpolate configuration
pub struct EIConfig<'a> {
    polya: &'a PolyU,
    polyb: &'a PolyU,
    eval   : Box<dyn Fn(&mut [Complex64])>,
    interp : Box<dyn Fn(&mut [Complex64])>,
    n : usize
}

#[derive(Copy, Clone)]
pub enum Algs {
    FFT,
    DFT,
}

impl<'a> EIConfig<'a> {
    pub fn new(polya: &'a PolyU, polyb: &'a PolyU, alg: Algs) -> EIConfig<'a> {
        // Takes the two polynomials and the type of algorithm to be used and then
        // calculates the other details
        EIConfig {
            polya,
            polyb,
            eval : match alg {
                Algs::FFT => Box::new( |a| {perform_fft(a, false).unwrap();}),
                Algs::DFT => Box::new( |a| dft(a, false)),
            },
            interp : match alg {
                Algs::FFT => Box::new( |a| {perform_fft(a, true).unwrap();}),
                Algs::DFT => Box::new( |a| dft(a, true)),
            },
            n: match alg {
                Algs::FFT => next_2pow(polya.deg() + polyb.deg() + 1),
                Algs::DFT => polya.deg() + polyb.deg() + 1,
            }
        }
    }
}

pub fn fast_mult(config: EIConfig) -> PolyU {

    // Expand the monomial vectors into a complex coefficient vector
    let a_sig = &mut to_coeffs_complex(&config.polya.terms[..], config.n)[..];
    let b_sig = &mut to_coeffs_complex(&config.polyb.terms[..], config.n)[..];

    // Evaluate the polynomials 
    (config.eval)(a_sig);
    (config.eval)(b_sig);

    // Multiply elementwise: TODO Make this infix
    let mut c_sig: Vec<Complex64> = a_sig.into_iter().zip(b_sig.into_iter())
                                       .map(|(a,b)| *a * *b).collect();
                                        
    // Interpolate the result
    (config.interp)(&mut c_sig[..]);

    // Convert back into i32
    let c_parsed = c_sig.into_iter().map(|x| x.re.round() as i32).collect();

    // Convert back into polynomial type
    PolyU::from_coeff("x".to_string(), c_parsed).unwrap()
}

fn to_coeffs_complex(input: &[Monomial], n: usize) -> Vec<Complex64> {
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

pub fn dft(signal: &mut [Complex64], inv: bool) {

    let n = signal.len();
    let result = signal.to_vec();

    // Generate the nth roots of unity
    let rou = gen_rou(n, inv);

    // Evaluates: p * q mod n.
    let modx = |p, q| (((p * q) % n) + n) % n;

    // Must normalise the output
    let scaling = if inv {1.0 / n as f64} else {1.0};

    // F(k) = \sum^n_{j=0} x_j e^{-2\pi i jk / n}
    for k in 0..n {
        let term: Complex64 = result.iter().enumerate()
                                    .map(|(i,c)| rou[modx(k, i)] * c)
                                    .sum();
        signal[k] = term * scaling;
    }
}

fn perform_fft(signal: &mut [Complex64], inv: bool) -> Result<(), &'static str> { 
    // Sample is the signal. Note that is fully expanded,
    // Performs the FFT inline

    // Basic contraint checking
    match signal.len() {
        0 => Err("Signal cannot be empty"),
        1 => Ok(()), // Already in DFT form
        2 => {
            let x0 = signal[0];
            let x1 = signal[1];
            signal[0] = x0 + x1;
            signal[1] = x0 - x1;
            Ok(())
        },
        n if !is_2_pow(n) => 
             Err("Signal needs to be a power of two"),
        _ => {
            go_fast(signal, inv); 
            Ok(())
        }

    }
}

fn go_fast(signal: &mut [Complex64], inv: bool) {
    // Assumes that the length of 'signal' is >= 4

    let n = signal.len();
    let rou = gen_rou(n, inv);

    // Does first iteration and puts into reverse bit order.
    let mut i = 0;
    while i < n / 2 {
        let x_0 = signal[i];
        let x_1 = signal[i + 1];
        let x_0_n2 = signal[i + (n / 2)];
        let x_1_n2 = signal[i + (n / 2) + 1];
        signal[i] = x_0 + x_0_n2;
        signal[i + 1] = x_0 - x_0_n2;
        signal[i + (n / 2)] = x_1 + x_1_n2;
        signal[i + (n / 2) + 1] = x_1 - x_1_n2;
        i += 2;
    }
    // We now assume the first layer is done and reverse bit order satisfied

    // Starts at index two because we already handled the first one
    for i in 2..=log2(n) {
        let flut = 1 << i; // No. of elements in a flutter
        // Iterate over all the flutters, j is their starting index
        for j in (0..n).step_by(flut) {
            // Width of the k-flutter (number of elements)
            for (k,l) in (j..j+(flut>>1)).zip(j+(flut>>1)..j+flut) {
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
    pub fn log2(n: usize) -> usize { if n != 1 {log2(n >> 1) + 1} else {0} }

    // Quick test to see if a number is a power of two
    pub fn is_2_pow(n: usize) -> bool {(n - 1) & n == 0}

    // Rounds up to the nearest power of two.
    pub fn next_2pow(n: usize) -> usize { 
        if is_2_pow(n) {n} else {1 << (log2(n) + 1)}
    }

    pub fn gen_rou(n: usize, inverse: bool) -> Vec<Complex64> {
        // Generates all the nth roots of unity 
        // Changes it depending on whether computing the dft or the inverse
        let sign = if inverse {1.0} else {-1.0};
        (0..n).map(|k| {
            Complex64::new(0.0, sign * 2.0 * PI / n as f64).scale(k as f64).exp()
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    extern crate chrono;
    extern crate rand; 
    use chrono::*;
    use super::fft::*;
    use crate::*;
    use polyu::*;


    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn bench_dense() {
        // Note all coefficients are nonzero
        let between = Uniform::from(1..100);
        let mut rng = rand::thread_rng();
        // A function to randomly generate a polynomial with n coefficients
        let mut make_poly = |n: usize| -> PolyU {
            let res_vec = (0..n).map(|_| between.sample(&mut rng)).collect();
            PolyU::from_coeff("x".to_string(), res_vec).unwrap()
        };

        // Benches the time required to multiply two arbitrary polynomials of deg = n
        let mut time_mult = |n: usize| {
            let a = make_poly(n);
            let b = make_poly(n);
            let ab_fft = EIConfig::new(&a, &b, Algs::FFT);
            // let ab_dft = EIConfig::new(&a, &b, Algs::DFT);
            println!("-------------------------------------------");
            println!("Number of elements = {}", n);
            println!("-------------------------------------------");
            println!("FFT: {:?}", Duration::span(|| {fast_mult(ab_fft);}));
            // println!("DFT: {:?}", Duration::span(|| {fast_mult(ab_dft);}));
            println!("STD: {:?}", Duration::span(|| {a.mul(&b);}));
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
    fn general_testing() {
        mult_test(Algs::DFT);
        mult_test(Algs::FFT);
    }

    fn mult_test(alg: Algs) {

        let a = PolyU::from_coeff("x".to_string(), vec![1,1]).unwrap();
        let b = PolyU::from_coeff("x".to_string(), vec![1,3]).unwrap();
        let c = PolyU::from_coeff("x".to_string(), vec![1,2,1]).unwrap();

        let opab = EIConfig::new(&a, &b, alg);
        let bench = EIConfig::new(&a, &b, alg);
        let opbc = EIConfig::new(&b, &c, alg);
        let opac = EIConfig::new(&a, &c, alg);

        println!("{:?}", Duration::span(|| {fast_mult(bench);}));
        assert_eq!(a.mul(&b), fast_mult(opab));
        assert_eq!(b.mul(&c), fast_mult(opbc));
        assert_eq!(c.mul(&a), fast_mult(opac));

        let d = PolyU::from_coeff("x".to_string(), vec![-1,3]).unwrap();
        let e = PolyU::from_coeff("x".to_string(), vec![-1,3,4,6]).unwrap();
        let opda = EIConfig::new(&d, &a, alg);
        let opde = EIConfig::new(&d, &e, alg);
        assert_eq!(d.mul(&a), fast_mult(opda));
        assert_eq!(d.mul(&e), fast_mult(opde));

        let f = PolyU::from_coeff("x".to_string(), vec![0]).unwrap();
        let opaf = EIConfig::new(&a, &f, alg);
        assert_eq!(a.mul(&f), fast_mult(opaf));

        
    }
 }