// This only does the base 2 FFT. For an arbitrary base, see "fft-mullti-base.rs"
extern crate chrono;

use crate::algebras::*;
use crate::mathutils::log2_unchecked;

pub trait SupportsFFT: ScalarField {
    // Generates the roots of unity
    fn rou(n: usize, inv: bool) -> Vec<Self>;
    fn divby2(&mut self, n: usize);
}

// Warning: Does it infix by default for speed, infix on a_sig
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
    perform_fft(a_sig, false)?;
    perform_fft(b_sig, false)?;

    // Multiply elementwise
    for i in 0..n {
        a_sig[i] += b_sig[i];
    }

    // Interpolate the result
    perform_fft(a_sig, true)?;
    Ok(())
}

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

pub fn perform_fft<F: SupportsFFT>(signal: &mut [F], inv: bool) -> Result<(), &'static str> {
    // Sample is the signal. Performs the FFT inline
    // Does not normalise the inverse!

    // Basic contraint checking
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

pub fn go_fast<F: SupportsFFT>(sig: &mut [F], inv: bool) {
    // Assumes that the length of 'sig' is >= 4

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
    extern crate chrono;
    extern crate rand;

    use crate::algebras::complex::CC;
    use crate::algebras::polyring::*;
    use crate::fast_mult::*;
    use crate::polyu::*;
    use chrono::*;

    use super::perform_fft;
    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn fft_test() {
        let mut sig = vec![
            CC::from_re(9),
            CC::from_re(3),
            CC::from_re(3),
            CC::from_re(6),
        ];
        let _ = perform_fft(&mut sig[..], false);

        println!("{:?}", sig)
    }

    // TODO I'm not sure this is correctly working, I changed the =log2_unchecked to log2_unchecked
    #[test]
    fn bench_dense_main() {
        let ring = PRDomain::<CC, UniVarOrder>::new(vec!['x']);

        // Note all coefficients are nonzero
        let dist = Uniform::from(1..100);
        let mut rng = rand::thread_rng();
        // A function to randomly generate a polynomial with n coefficients
        let mut make_poly = |n: usize| -> PolyU<CC> {
            let res_vec = (0..n).map(|_| CC::from_re(dist.sample(&mut rng))).collect();
            Poly::from_coeff(&ring, res_vec)
        };

        // Benches the time required to multiply two arbitrary polynomials of deg = n
        let mut time_mult = |n: usize| {
            let a = make_poly(n);
            let b = make_poly(n);

            println!("-------------------------------------------");
            println!("Number of elements = {}", n);
            println!("-------------------------------------------");
            println!(
                "FFT: {:?}",
                Duration::span(|| {
                    a.fast_mult(&b);
                })
            );
            println!("-------------------------------------------");
        };

        for i in 5..16 {
            time_mult(1 << i);
        }
    }

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
