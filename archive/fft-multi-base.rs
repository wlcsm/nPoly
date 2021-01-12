extern crate chrono;

use mathutils::*;
use crate::algebras::*;

mod mathutils {
    pub fn logn(num: usize, n: usize) -> usize {
        if num < n {
            0
        } else {
            logn(num / n, n) + 1
        }
    }

    pub fn is_n_pow(num: usize, n: usize) -> bool {
        n.pow(logn(num, n) as u32) == num
    }

    pub fn next_npow(num: usize, n: usize) -> usize {
        if is_n_pow(num, n) {
            num
        } else {
            n.pow((logn(num, n) + 1) as u32)
        }
    }
}

pub fn perform_fft_multi<T: SupportsFFT>(signal: &mut [T], inv: bool, base: usize) -> Result<(), &'static str> {
    // Sample is the signal. Performs the FFT inline
    // Does not normalise the inverse!

    // Basic contraint checking
    match signal.len() {
        0 => Err("Signal cannot be empty"),
        n if !is_n_pow(n, base) => Err("Signal needs to be a power of the base"),
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

pub fn go_fast_base_n<T: SupportsFFT>(sig: &mut [T], inv: bool, m: usize) {
    // Assumes that the length of 'signal' is >= 4

    let n   = sig.len();
    let rou = <F>::rou(n, inv); // Generates roots of unity

    // Does first iteration and puts into reverse bit order.
    let mut result = vec![<F>::zero(); n];

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
            let mut tmp = vec![<F>::zero(); m];
            for ind in 0..flut/m {

                for p in 0..m {
                    tmp[p] = rou[ind * p* (n / flut)].mul(&result[j + p * (flut/m) + ind]);
                }

                // Uses an FFT to combine them
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


#[cfg(test)]
mod tests {
    extern crate chrono;
    extern crate rand;
    use super::*;
    use chrono::*;
    use crate::fast_mult::*;
    use crate::algebras::integers::ZZ;
    use crate::polyu::*;

    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn bench_dense_main() {
        // Note all coefficients are nonzero
        let between = Uniform::from(1..100);
        let mut rng = rand::thread_rng();
        // A function to randomly generate a polynomial with n coefficients
        let mut make_poly = |n: usize| -> Poly<ZZ> {
            let res_vec = (0..n).map(|_| ZZ(between.sample(&mut rng))).collect();
            Poly::from_coeff(None, res_vec).unwrap()
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
                    a.fast_mult(&b, 2);
                })
            );
            println!("-------------------------------------------");
        };

        for i in 5..16 {
            time_mult(1 << i);
        }
    }

    #[test]
    fn mult_test_main() {
        let a = Poly::from_coeff(None, vec![ZZ(1), ZZ(1)]).unwrap();
        let b = Poly::from_coeff(None, vec![ZZ(1), ZZ(3)]).unwrap();
        let c = Poly::from_coeff(None, vec![ZZ(1), ZZ(2), ZZ(1)]).unwrap();

        assert_eq!(a.mul(&b), a.fast_mult(&b, 2));
        assert_eq!(b.mul(&c), b.fast_mult(&c, 2));
        assert_eq!(c.mul(&a), c.fast_mult(&a, 2));

        let d = Poly::from_coeff(None, vec![ZZ(-1), ZZ(3)]).unwrap();
        let e = Poly::from_coeff(None, vec![ZZ(-1), ZZ(3), ZZ(4), ZZ(6)]).unwrap();

        assert_eq!(d.mul(&a), d.fast_mult(&a, 2));
        assert_eq!(d.mul(&e), d.fast_mult(&e, 2));

        let f = Poly::from_coeff(None, vec![ZZ(0)]).unwrap();

        assert_eq!(f.mul(&a), f.fast_mult(&a, 2));
    }
}
