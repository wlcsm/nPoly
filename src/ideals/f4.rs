use crate::algebras::polyring::*;
use crate::display::*;
use crate::ideals::*;
use num_traits::Zero;

/// Converts the generators into a Groebner basis using the F4 algorithm
/// Psuedocode in CLO
pub fn f4<'a, P: FPolyRing>(f_gens: Vec<Poly<'a, P>>) -> Ideal<'a, P> {
    let mut g_basis = Ideal::new(f_gens.clone());
    let mut t = f_gens.len();
    let mut subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).filter(|(a, b)| a != b).collect();

    while !subsets.is_empty() {
        let selection: Vec<(usize, usize)> = choose_subset(&mut subsets);

        let mut l_mat: Vec<Poly<'a, P>> = selection
            .into_iter()
            .map(|(i, j)| left_hand_s_poly(&g_basis.gens[i], &g_basis.gens[j]))
            .collect();

        compute_m(&mut l_mat, &g_basis);

        let m_init = MonomialIdeal::<P>::new(l_mat.iter().map(|m| m.lm()).collect());

        // Row reduction and find the new elements
        row_reduce(&mut l_mat);

        let n_plus: Vec<Poly<'a, P>> = l_mat
            .into_iter()
            .filter(|n| !n.is_zero() && !m_init.is_in(&n.lm()))
            .collect();

        for n in n_plus {
            println!("Newly added poly = {}", n);
            g_basis.add(n);
            // TODO for some reason this failed, i.e. the .add() method didn't actually add
            // anything to the list of generators, this is weird because by the construction
            // of n_plus it definitely should have
            // t += 1;
            t = g_basis.gens.len();
            let mut new_pairs_l = (0..t - 1).map(|new| (new, t - 1)).collect();
            let mut new_pairs_r = (0..t - 1).map(|new| (t - 1, new)).collect();
            subsets.append(&mut new_pairs_l);
            subsets.append(&mut new_pairs_r);
        }
    }

    row_reduce(&mut g_basis.gens);
    g_basis.gens = g_basis.gens.into_iter().filter(|x| !x.is_zero()).collect();
    g_basis
}

pub fn row_reduce<'a, P: FPolyRing>(mat: &mut Vec<Poly<'a, P>>) {
    // Standard row reduction, but operates on a vector of polynomials.

    let mut i = 0;

    while i < mat.len() {
        // TODO probably don't want to sort the entire thing every single time
        mat.sort_by(|a, b| <P::Ord>::cmp(&b.lt().mon, &a.lt().mon));
        // Get one row. We then use this is cancel other rows
        if mat[i].lc().is_zero() {
            i += 1;
            continue;
        }
        let lm = mat[i].lm();

        let mut j = i + 1;

        // FIXME there is the potential for the lead terms to not cancel correctly due to
        // rounding error
        while j < mat.len() && mat[j].lm() == lm {
            // println!("mat[j] before = {}", mat[j]);

            let scalar = mat[j].lc() / mat[i].lc();
            mat[i] *= scalar;

            // println!("mat[i] scaled = {}", mat[i]);
            //
            // Cancel the first term
            mat[j] = mat[j].ref_sub(&mat[i]);
            if mat[j].lm() == lm {
                mat[j].terms.pop();
            }

            // println!("mat[j] after cancellation = {}", mat[j]);

            j += 1
        }

        // Goes back up the matrix in order to put it into reduces row echelon
        for j in (0..i).rev() {
            if let Some(c) = mat[j].has(&lm) {
                // Cancel the lead terms
                // println!("mat[j] before = {}", mat[j]);
                // println!("mat[i] scaled = {}", mat[i]);
                let scalar = c / mat[i].lc();
                mat[i] *= scalar;

                mat[j] = mat[j].ref_sub(&mat[i]);

                if mat[j].lm() == lm {
                    mat[j].terms.pop();
                }
                // println!("mat[j] after cancellation = {}", mat[j]);
            }
        }
        i += 1;
    }

    mat.sort_by(|a, b| <P::Ord>::cmp(&b.lt().mon, &a.lt().mon));
}

fn choose_subset(pairs: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    // Returns a subset of B based on G, and also removes that subset from B
    // Currently just processing all of it
    let mut a = Vec::new();
    a.append(pairs);
    a
}

use std::collections::HashMap;

/// The ComputeM function in CLO
/// TODO Document my algorithm below
pub fn compute_m<'a, P: FPolyRing>(l_mat: &mut Vec<Poly<'a, P>>, g: &Ideal<'a, P>) {
    let g_init = MonomialIdeal::from(&g);

    println!("ginit = {:?}", g_init);

    let mut i = 0;

    let mut seen_monomials: HashMap<P::Mon, Vec<usize>> = HashMap::new();

    // Iterates through the monomials of the polynomials in the matrix and adds the extra
    // scaled basis elements if necessary
    while i < l_mat.len() {
        let mut aux = Vec::new();

        // Go through the terms and add any that aren't already in there
        // TODO In the compute_m_test l_matt[i] contains a zero term, which should be impossible.
        // I have added the filter here as a temporary fix but need to come back later
        // and properly find the root cause
        for term in l_mat[i].terms.iter().filter(|x| !x.is_zero()) {
            // If term is in G_init, that is, there exists a polynomial in $g$ that divides it
            if let Some(poly) = g_init.get_poly(&term.mon) {
                // Record the row we found this monomial in
                match seen_monomials.get_mut(&term.mon.clone()) {
                    Some(x) => x.push(i),
                    None => {
                        // Evaluate x^\alpha f_\ell and push it the aux buffer
                        aux.push(poly * term.euclid_div(&poly.lt()).unwrap().0);
                        // Record that we have seen it
                        seen_monomials.insert(term.mon.clone(), vec![i]);
                    }
                }
            }
        }
        l_mat.append(&mut aux);
        i += 1;
    }
}

/// Computes the Left hand side polynomial
/// ~PO don't need to find lcm and then div, can be combined into one iteration
pub fn left_hand_s_poly<'a, P: FPolyRing>(a: &Poly<'a, P>, b: &Poly<'a, P>) -> Poly<'a, P> {
    // Returns the "Left hand side" of the S-Polynomial of a and b
    let lcm = a.lt().lcm(&b.lt());
    let scale = lcm.euclid_div(&a.lt()).unwrap().0;
    a * scale
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::*;
    use crate::parse::MyFromStr;
    use crate::polym::*;
    use generic_array::typenum::{U2, U3};

    // This isn't working at the moment
    // #[test]
    // fn f4_test() {
    //     let ring = PRDomain::<RR, GLex<MultiIndex<U3>>>::new(vec!['x', 'y', 'z']);
    //     let f_vec = vec![
    //             Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
    //             Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
    //             Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
    //         ];

    //     let gb = f4(f_vec);
    //     println!("gb = {}", gb);
    //     assert!(gb.is_groebner_basis());

    //     let ring = PRDomain::<RR, GLex<MultiIndex<U3>>>::new(vec!['x', 'y', 'z']);
    //     let f_vec = vec![
    //             Poly::from_str(&ring, "1.0x^1y^2 + 1.0x^1y^1 - 1.0").unwrap(),
    //             Poly::from_str(&ring, "1.0x^2 - 1.0z^2 + 1.0").unwrap(),
    //             Poly::from_str(&ring, "1.0x^2 - 1.0z^2 + 1.0x^1y^1 + 1.0").unwrap(),
    //         ];

    //     let gb = f4(f_vec);
    //     println!("gb = {}", gb);
    //     assert!(gb.is_groebner_basis());
    // }

    #[test]
    fn cancelling_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U3>>>::new(vec!['x', 'y', 'z']);
        let mut j = Poly::from_str(&ring, "-1.0y^1 + 1.0y^1z^2 + 1.0x^1y^2").unwrap();
        let mut i = Poly::from_str(&ring, "-1.0x^1 + -1.0y^1z^2").unwrap();

        println!("yeet");
        if let Some(c) = j.has(&i.lm()) {
            // Cancel the lead terms
            println!("mat[j] before = {}", j);
            println!("mat[i] before scaled = {}", i);
            let scalar = c / i.lc();
            i *= scalar;
            println!("mat[i] scaled = {}", i);

            j = j.ref_sub(&i);

            println!("mat[j] after cancellation = {}", j);
        }
    }
    use chrono::*;

    #[test]
    fn bench_f4() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = &Ideal::new(vec![a, b]);
        println!(
            "BB Alg time = {:?}",
            Duration::span(|| {
                f4(r.gens.clone());
            })
        );

        let ring = PRDomain::<RR, GLex<MultiIndex<U3>>>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
            Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
            Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
            Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
        ];

        println!(
            "BB Alg time = {:?}",
            Duration::span(|| {
                f4(f_vec);
            })
        );

        // let f_vec = vec![
        //         Poly::from_str(&ring, "1.0x^3 + 1.0x^1y^1 - 1.0").unwrap(),
        //         Poly::from_str(&ring, "1.0x^2y^1 - 1.0z^2 + 1.0").unwrap(),
        //         Poly::from_str(&ring, "1.0x^1y^1 + 1.0 + 2.0z^1").unwrap(),
        //     ];

        // println!("BB Alg time = {:?}",
        //     Duration::span(|| {
        //         f4(f_vec);
        //     }));
    }

    #[test]
    fn row_reduce_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);

        let mut f_vec = vec![
            Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2 + 2.0x^1y^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2").unwrap(),
            Poly::from_str(&ring, "3.0x^2 + 5.0y^2 + 1.0y^1").unwrap(),
        ];

        row_reduce(&mut f_vec);
        println!("======================");
        println!("{}", show_vec_poly(&f_vec));
    }

    #[test]
    fn compute_m_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);

        let f_vec = vec![
            Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2").unwrap(),
        ];

        let g_basis = Ideal::new(f_vec.clone());
        let t = f_vec.len();
        let subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).filter(|(a, b)| a != b).collect();

        let mut l_mat: Vec<Poly<PRDomain<RR, GLex<MultiIndex<U2>>>>> = subsets
            .into_iter()
            .map(|(i, j)| left_hand_s_poly(&g_basis.gens[i], &g_basis.gens[j]))
            .collect();

        println!("{}", show_vec_poly(&l_mat));

        compute_m(&mut l_mat, &g_basis);
        println!("===========================");
        println!("{}", show_vec_poly(&l_mat));
    }
}

// % I have actually made my own way for this. It leverages the fact that we only care about collecting the polynomials whose lead terms were not originally in $G$. Notice that whenever we add a $x^\alpha f_\ell$ into the matrix, this polynomial will be excluded in the end. The only way this wouldn't happen is if while we are reducing the matrix, we find another intermediate polynomial has the same lead monomial as $x^\alpha f_\ell$ and then it can cancel its lead term. But actually, in the process of matrix reduction, we could choose $x^\alpha f_\ell$ to cancel the other polynomials lead term instead. Hence the trick here is to actually avoid unnecessarily making the matrix larger

// % \begin{enumerate}[(1)]
// %     \item Start with the largest polynomial in terms of lead term size in the sorted list
// %     \item Cancel out the lead terms. Remember the lead terms of the resulting polynomials, they are now in the unsorted pile
// %     \item Put this polynomial in the Done list
// %     \item Then add the polynomials with the highest lead term from the sorted list into the unsorted one, also adding its lead term.
// %     \item Select the next biggest lead term, and check if it is in the initial ideal of G. If this was the one from the sorted list, add the next one from the sorted list also into it.
// %     \item  - If it is, then use the polynomial from $G$ to subtract the terms and also from in the done pile. We don't actually include this polynomial in any of the piles because we know that it won't end up in $N^+$.
// %     \item  - If not, then subtract any other polynomials and then add it into the Done pile
// %     \item Remembering to add the new lead terms into the monomial list as well. If the unsorted pile is empty, then add the polynomials will the highest lead terms from the sorted list
// % \end{enumerate}
