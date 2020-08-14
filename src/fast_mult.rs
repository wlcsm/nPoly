
use crate::algebras::polyring::*;
// use crate::algebras::*;
use crate::fft::*;
use crate::mathutils::log2_unchecked;
// use crate::polyu::*;

pub trait FastMult {
    fn fast_mult(&self, b: &Self) -> Self;
}

// // TODO these two should also be a macro
// fn to_coeff_vec_complex(input: &Vec<Term<ZZ>>, n: usize) -> Vec<CC> {
//     // Expands the input into the expanded coefficient vector (coerced into complex)
//     // Then padded with zeros to length n

//     let mut result: Vec<CC> = Vec::with_capacity(n);
//     for Term { coeff, deg } in input {
//         // Fill the gap between monomials with zeros, then add the monomial
//         result.resize(deg.0 as usize, CC::zero());
//         result.push(CC::from_re(coeff.0));
//     }
//     // Pad the rest
//     result.resize(n, CC::zero());
//     result
// }

use num_traits::Zero;

/// Expands the input into a coefficient vector padded with zeros to length n
pub fn to_coeff_vec<P: PolyRing>(input: &[Term<P>], n: usize) -> Vec<P::Coeff> {

    let mut result: Vec<P::Coeff> = Vec::with_capacity(n);
    for Term { coeff, mon } in input.iter() {
        // Fill the gap between monomials with zeros, then add the monomial
        result.resize(mon.tot_deg(), <P::Coeff>::zero());
        result.push(*coeff);
    }
    // Pad the rest
    result.resize(n, <P::Coeff>::zero());
    result
}

use crate::polyu::*;

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

    use crate::polyu::*;
    use crate::algebras::polyring::*;
    use crate::algebras::*;

    pub fn karatsuba<'a, P: PolyRingUni>(poly_a: &Poly<'a, P>, poly_b: &Poly<'a, P>) -> Poly<'a, P> {

        let n_a = poly_a.deg().next_power_of_two();
        let n_b = poly_b.deg().next_power_of_two();

        let sig_a = super::to_coeff_vec(&poly_a.terms, n_a);
        let sig_b = super::to_coeff_vec(&poly_b.terms, n_b);

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
            return vec![arr_a[0]*arr_b[0], arr_a[1]*arr_b[0] + arr_a[0]*arr_b[1], arr_a[1]*arr_b[1]]
        }
    
        // Actually splits the polynomials
        let (low_a, high_a) = arr_a.split_at(n / 2);
        let (low_b, high_b) = arr_b.split_at(n / 2);
        
        let z0 = dense(low_a, low_b);
        let z1 = dense(&dense_add_slices(low_a, high_a)[..], &dense_add_slices(low_b, high_b)[..]);
        let z2 = dense(high_a, high_b);

        // Pieces together the parts manually
        let s = n / 2;

        let q1 = z0[..s].iter().map(|&x| x);

        let q2 = izip!(z0[s..].iter(), z1[..s].iter(),
                       z2[..s].iter(), z0[..s].iter())
                        .map(|(&z0_h, &z1_l, &z2_l, &z0_l)| z0_h + z1_l - z2_l - z0_l);

        let extra = std::iter::once(z1[s-1] - z2[s-1] - z0[s-1]);

        let q3 = izip!(z2[..s].iter(), z1[s..].iter(),
                       z2[s..].iter(), z0[s..].iter())
                        .map(|(&z2_l, &z1_h, &z2_h, &z0_h)| z2_l + z1_h - z2_h - z0_h);

        let q4 = z2[s-1..].iter().map(|&x| x);

        let res = q1.chain(q2).chain(extra).chain(q3).chain(q4).collect();

        res
    }

    fn dense_add_slices<'a, R: ScalarRing>(arr_a: &[R], arr_b: &[R]) -> Vec<R> {
        izip!(arr_a.into_iter(), arr_b.into_iter()).map(|(a, b)| *a + *b).collect()
    }
    // TODO Had some troubles
    // fn sparse<'a, P: PolyRingUni>(arr_a: &[Term<P>], arr_b: &[Term<P>]) -> Vec<Term<P>> {

    //     let n_a = arr_a.len();
    //     let n_b = arr_b.len();

    //     if (n_a <= 1) || (n_b <= 1) {
    //         return match (n_a, n_b) {
    //             (0, 0) => Vec::new(),
    //             (1, 0) => arr_a.to_vec(),
    //             (0, 1) => arr_b.to_vec(),
    //             (1, 1) => vec![arr_a[0] * arr_b[0]],
    //         }
    //     }
    
    //     // Calculates where to split the polynomial at
    //     let split_value = std::cmp::max(arr_a[n_a-1].mon.tot_deg(), arr_b[n_b-1].mon.tot_deg()) / 2;
        
    //     // Gets the digit to split the polynomials at
    //     let split_index_a = split_at(arr_a, split_value);
    //     let split_index_b = split_at(arr_b, split_value);

    //     // Actually splits the polynomials
    //     let (low_a, high_a) = arr_a.split_at(split_index_a);
    //     let (low_b, high_b) = arr_b.split_at(split_index_b);
        
    //     // FIXME This line is a hack, I need a ring to make a polynomial so I can add them
    //     // together. Maybe the add function should only be implemented on slices?
    //     let ring = <P>::new(vec!['x']);

    //     let z0 = Poly::from_terms( sparse(low_a, low_b), Some(&ring));
    //     let z1 = Poly::from_terms( sparse(low_a + scale_down(high_a, split_value), low_b + scale_down(high_b, split_value)), Some(&ring));
    //     let z2 = Poly::from_terms( sparse(high_a, high_b), Some(&ring));

    //     z2 + z1 - scale_down(&z2[..], split_value) + scale_up(&z0[..], split_value) + z0
            
    // }

    // fn scale_down<'a, P: PolyRingUni>(terms: &[Term<P>], deg: usize) -> Vec<Term<P>> {
    //     terms.into_iter().map(|t| Term::new(t.coeff, UniIndex(t.mon.0 - deg))).collect()
    // }
    // fn scale_up<'a, P: PolyRingUni>(terms: &[Term<P>], deg: usize) -> Vec<Term<P>> {
    //     terms.into_iter().map(|t| Term::new(t.coeff, UniIndex(t.mon.0 + deg))).collect()
    // }

    #[cfg(test)]
    mod test {
        extern crate test;
        use crate::algebras::real::RR;
        use super::*;
        use crate::polyu::UniVarOrder;
        use crate::parse::*;
        use test::Bencher;

        #[test]
        fn karatsuba_mult_test() {

            let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);
            let a = Poly::from_str(&ring, "1.0x^1 + 3.0x^2 + 5.0x^5").unwrap();
            let b = Poly::from_str(&ring, "1.0x^1 + 2.0x^3 + 2.0x^5").unwrap();

            let res = karatsuba(&a, &b);
            println!("{}", res)
        }

        #[bench]
        fn karatsuba_bench(b: &mut Bencher) {
            let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);
            let poly_a = Poly::from_str(&ring, "1.0x^1 + 3.0x^2 + 5.0x^5").unwrap();
            let poly_b = Poly::from_str(&ring, "1.0x^1 + 2.0x^3 + 2.0x^5").unwrap();

            b.iter(|| karatsuba(&poly_a, &poly_b));
        }
    }


}

// impl<'a, P: PolyRing<Var=U1>> FastMult for Poly<'a, P> {
//     fn fast_mult(&self, other: &Self) -> Self {
//         let n = (self.deg() + other.deg() + 1).next_power_of_two();
//         let mut a_sig = to_coeff_vec(&self.terms, n);
//         let mut b_sig = to_coeff_vec(&other.terms, n);

//         // Infix on a_sig
//         eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

//         // Need to normalise it here
//         for x in a_sig.iter_mut() {
//             x.divby2(log2_unchecked(n))
//         }

//         // Convert back into polynomial type
//         Poly::from_coeff(self.ring, a_sig)
//     }
// }
// Frozen this at the moment becaue it requires a change of type from CC to ZZ
// which I don't feel like doing at the moment

// impl<'a, P: PolyRing<Coeff=CC, Var=Univariate>> FastMult for Poly<'a, P> {

//     fn fast_mult(&self, other: &Self) -> Self {

//         let n = (self.deg() + other.deg() + 1).next_power_of_two();
//         let mut a_sig = to_coeff_vec_complex(&self.terms, n);
//         let mut b_sig = to_coeff_vec_complex(&other.terms, n);

//         eval_interp(&mut a_sig[..], &mut b_sig[..]).unwrap();

//         // Because we converted it to complex for the ROU
//         // We also need to normalise it here
//         let c_parsed = a_sig.into_iter()
//                             .map(|x|
//                                 ZZ((x.0 / n as f64).re.round() as i32)
//                             ).collect();

//         // Convert back into polynomial type
//         Poly::from_coeff(self.ring, c_parsed)
//     }
// }
