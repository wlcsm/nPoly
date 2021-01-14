use crate::algebras::*;
use itertools::izip;
use crate::polynomials::*;

pub fn karatsuba<R: ScalarRing>(
    poly_a: &Poly<R, 1>,
    poly_b: &Poly<R, 1>,
) -> Poly<R, 1> {
    //  Pad to the nearest power of two
    let n = (std::cmp::max(poly_a.tot_deg(), poly_b.tot_deg()) + 1).next_power_of_two();

    let sig_a = to_dense(&poly_a, Some(n));
    let sig_b = to_dense(&poly_b, Some(n));

    let res = dense(&sig_a.terms[..], &sig_b.terms[..]);
    Poly::from_coeff(&res[..])
}

/// Assumes that the arrays are both the same power of two.]
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
/// Then we can collect the iterator into a vector to get the desired output.
fn dense<R: ScalarRing>(arr_a: &[R], arr_b: &[R]) -> Vec<R> {
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
        &*add_slice(low_a, high_a),
        &*add_slice(low_b, high_b),
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

    q1.chain(q2).chain(extra).chain(q3).chain(q4).collect()
}

fn add_slice<R: ScalarRing>(arr_a: &[R], arr_b: &[R]) -> Box<[R]> {
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

    // #[test]
    // fn karatsuba_small_test() {
    //     let res = karatsuba_tester(vec![2, 10, 40]);
    //     for (size, time) in res {
    //         println!(
    //             "Karatsuba's algorithm took {} to multiply polynomials of degree {}",
    //             time, size
    //         )
    //     }
    // }

    // #[test]
    // fn karatsuba_medium_test() {
    //     let res = karatsuba_tester(vec![1 << 7, 1 << 8, 1 << 9]);
    //     for (size, time) in res {
    //         println!(
    //             "Karatsuba's algorithm took {} to multiply polynomials of degree {}",
    //             time, size
    //         )
    //     }
    // }

    // /// Takes a vector of usize and then generates two polynomials of that size and multiplies
    // /// them with Karatsuba's algorithm. Returns vector containing the original argument and
    // /// the duration of the multiplication in a tuple
    // fn karatsuba_tester(test_sizes: Vec<usize>) -> Vec<(usize, Duration)> {
    //     let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);

    //     let mut rng = rand::thread_rng();
    //     let dist = RR::gen_sampler(RR(10.0));

    //     // A function to randomly generate a polynomial with n coefficients
    //     let mut make_poly = |n: usize| -> PolyU<RR> {
    //         let res_vec = (0..n).map(|_| dist.sample(&mut rng)).collect();
    //         Poly::from_coeff(&ring, res_vec)
    //     };

    //     test_sizes
    //         .into_iter()
    //         .map(|size| {
    //             let execution_time = {
    //                 let a = make_poly(size);
    //                 let b = make_poly(size);

    //                 Duration::span(|| {
    //                     karatsuba(&a, &b);
    //                 })
    //             };
    //             (size, execution_time)
    //         })
    //         .collect()
    // }
}
