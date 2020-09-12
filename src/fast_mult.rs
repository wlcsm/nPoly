use crate::algebras::polyring::*;
use crate::fft::*;
use crate::mathutils::log2_unchecked;
use crate::polyu::*;

/// Defines a fast multiplication function that can be used to multiply polynomials asymptotically
/// quicker than the standard schoolbook algorithm that is implemented for all polynomials
pub trait FastMult {
    fn fast_mult(&self, b: &Self) -> Self;
}

impl<'a, T: SupportsFFT> FastMult for PolyU<'a, T> {
    // FFT Multiplication
    fn fast_mult(&self, other: &Self) -> Self {
        let n = (self.deg() + other.deg() + 1).next_power_of_two();

        let mut a_sig = to_coeff_vec(&self.terms, n);
        let mut b_sig = to_coeff_vec(&other.terms, n);

        // Infix on a_sig
        eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

        // Need to normalise it here
        for x in a_sig.iter_mut() {
            x.divby2(log2_unchecked(n))
        }

        // Convert back into polynomial type
        Poly::from_coeff(self.ring.unwrap(), a_sig)
    }
}

pub mod karatsuba {

    use crate::algebras::polyring::*;
    use crate::algebras::*;
    use crate::polyu::*;

    pub fn karatsuba<'a, P: PolyRingUni>(
        poly_a: &Poly<'a, P>,
        poly_b: &Poly<'a, P>,
    ) -> Poly<'a, P> {
        //  Pad to the nearest power of two
        let n = (std::cmp::max(poly_a.deg(), poly_b.deg()) + 1).next_power_of_two();

        let sig_a = super::to_coeff_vec(&poly_a.terms, n);
        let sig_b = super::to_coeff_vec(&poly_b.terms, n);

        let res = dense(&sig_a, &sig_b);
        Poly::from_coeff(poly_a.ring.unwrap(), res)
    }

    /// Assumes that the arrays are both the same power of two.]
    /// TODO I think I can relax this requirement later.
    ///
    /// I have a weird line below where I use the q1, q2, q3, q4 iterators to
    /// combine everything together at the end. This is because in Rust we cannot initialize an
    /// empty array (safely), we need to use a Vector, or pre-fill it with values. The way to get
    /// around this is to use iterators.
    ///
    /// Suppose our inputs are of length N (a power of two at the moment) and so we half them into
    /// two pieces of length N/2, then when we create z0, z1, z2, these are all of length N-1.
    ///
    /// To recombine them, we need to evaluate z2 * X^N + (z1 - z2 - z0) * X^(N/2) + z0, keep in
    /// mind that we are using coefficient vectors.
    ///
    /// Since z0 has N elements, the elements 0...N/2 will not interact with any of the other terms
    /// so we just prepend them to the vector to get q1.
    /// The elements N/2...N-1 will need to be added to the bottom half of the elements in (z1 - z2 -
    /// z0) * X^(N/2), this gives us q2.
    /// There is an extra element that is in (z1 - z2 - z0) * X^(N/2), but not z2 * X^N nor z0
    /// which is the X^(N-1) term, this gives us "extra"
    /// Similarly, the second half of the elements in (z1 - z2 - z0) * X^(N/2) will be added to the
    /// first 0...N/2 elements of z2 * X^N, this gives us q3.
    /// Then we jut append the remaining half of z2 to the end to give us q4.
    ///
    /// Then we can collect the iterator into a vector to get the desired output.
    fn dense<'a, R: ScalarRing>(arr_a: &[R], arr_b: &[R]) -> Vec<R> {
        let n = arr_a.len();

        if n == 2 {
            return vec![
                arr_a[0] * arr_b[0],
                arr_a[1] * arr_b[0] + arr_a[0] * arr_b[1],
                arr_a[1] * arr_b[1],
            ];
        }

        // Actually splits the polynomials
        let (low_a, high_a) = arr_a.split_at(n / 2);
        let (low_b, high_b) = arr_b.split_at(n / 2);

        let z0 = dense(low_a, low_b);

        let z1 = dense(
            &dense_add_slices(low_a, high_a)[..],
            &dense_add_slices(low_b, high_b)[..],
        );
        let z2 = dense(high_a, high_b);

        // Pieces together the parts manually
        let s = n / 2;

        let q1 = z0[..s].iter().map(|&x| x);

        let q2 = izip!(
            z0[s..].iter(),
            z1[..s].iter(),
            z2[..s].iter(),
            z0[..s].iter()
        )
        .map(|(&z0_h, &z1_l, &z2_l, &z0_l)| z0_h + z1_l - z2_l - z0_l);

        let extra = std::iter::once(z1[s - 1] - z2[s - 1] - z0[s - 1]);

        let q3 = izip!(
            z2[..s].iter(),
            z1[s..].iter(),
            z2[s..].iter(),
            z0[s..].iter()
        )
        .map(|(&z2_l, &z1_h, &z2_h, &z0_h)| z2_l + z1_h - z2_h - z0_h);

        let q4 = z2[s - 1..].iter().map(|&x| x);

        let res = q1.chain(q2).chain(extra).chain(q3).chain(q4).collect();

        res
    }

    fn dense_add_slices<'a, R: ScalarRing>(arr_a: &[R], arr_b: &[R]) -> Vec<R> {
        izip!(arr_a.into_iter(), arr_b.into_iter())
            .map(|(a, b)| *a + *b)
            .collect()
    }

    #[cfg(test)]
    mod test {
        extern crate test;
        use super::*;
        use crate::algebras::real::RR;
        use crate::bench::*;
        use crate::polyu::UniVarOrder;
        use chrono::*;
        use rand::distributions::uniform::UniformSampler;

        #[test]
        fn karatsuba_small_test() {
            let res = karatsuba_tester(vec![2, 10, 40]);
            for (size, time) in res {
                println!(
                    "Karatsuba's algorithm took {} to multiply polynomials of degree {}",
                    time, size
                )
            }
        }

        #[test]
        fn karatsuba_medium_test() {
            let res = karatsuba_tester(vec![1 << 7, 1 << 8, 1 << 9]);
            for (size, time) in res {
                println!(
                    "Karatsuba's algorithm took {} to multiply polynomials of degree {}",
                    time, size
                )
            }
        }

        /// Takes a vector of usize and then generates two polynomials of that size and multiplies
        /// them with Karatsuba's algorithm. Returns vector containing the original argument and
        /// the duration of the multiplication in a tuple
        fn karatsuba_tester(test_sizes: Vec<usize>) -> Vec<(usize, Duration)> {
            let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);

            let mut rng = rand::thread_rng();
            let dist = MyUniformDistRR::new(RR(-10.0), RR(10.0));

            // A function to randomly generate a polynomial with n coefficients
            let mut make_poly = |n: usize| -> PolyU<RR> {
                let res_vec = (0..n).map(|_| dist.sample(&mut rng)).collect();
                Poly::from_coeff(&ring, res_vec)
            };

            test_sizes
                .into_iter()
                .map(|size| {
                    let execution_time = {
                        let a = make_poly(size);
                        let b = make_poly(size);

                        Duration::span(|| {
                            karatsuba(&a, &b);
                        })
                    };
                    (size, execution_time)
                })
                .collect()
        }
    }
}
