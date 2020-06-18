use crate::algebras::polyring::*;
use crate::algebras::*;
use crate::ideals::*;

fn vec_poly_str<'a, P: PolyRing>(poly: &Vec<Poly<'a, P>>) -> String {
    poly.iter()
        .map(|p| format!("{}", p))
        .collect::<Vec<String>>()
        .join("\n")
}

// The F4 algorithm
pub fn f4<'a, P: FPolyRing>(f_gens: Vec<Poly<'a, P>>) -> Ideal<'a, P> {
    let mut g_basis = Ideal::new(f_gens.clone());
    let mut t = f_gens.len();
    println!("f+gens = {}", vec_poly_str(&f_gens));
    println!("t = {}", t);
    let mut subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).filter(|(a, b)| a != b).collect();

    while !subsets.is_empty() {
        let selection: Vec<(usize, usize)> = choose_subset(&mut subsets);

        println!("selection = {:?}", selection);
        println!("g_basis = \n{}", g_basis);
        let mut l_mat: Vec<Poly<'a, P>> = selection
            .into_iter()
            .map(|(i, j)| left_hand_s_poly(&g_basis.gens[i], &g_basis.gens[j]))
            .collect();

        compute_m(&mut l_mat, &g_basis);
        let m_init = MonomialIdeal::new(l_mat.iter().map(|m| m.lt()).collect());

        println!("compute_m = \n{}", vec_poly_str(&l_mat));

        // Row reduction and find the new elements
        row_reduce(&mut l_mat);
        let n_plus: Vec<Poly<'a, P>> = l_mat
            .into_iter()
            .filter(|n| !n.is_zero() && !m_init.is_in(&n.lt()))
            .collect();

        println!("n_plus = \n{}", vec_poly_str(&n_plus));

        println!("subsets = \n{:?}", subsets);
        for n in n_plus {
            println!("next to be added to n = {}", n);
            g_basis.add(n);
            t += 1;
            let mut new_pairs_l = (0..t-1).map(|new| (new, t-1)).collect();
            let mut new_pairs_r = (0..t-1).map(|new| (t-1, new)).collect();
            subsets.append(&mut new_pairs_l);
            subsets.append(&mut new_pairs_r);
        }
        println!("len = {}", g_basis.num_gens());
        println!("t = {}", t);
        println!("subsets = \n{:?}", subsets);
    }
    g_basis
}

pub fn row_reduce<'a, P: FPolyRing>(mat: &mut Vec<Poly<'a, P>>) {
    // Standard row reduction, but operates on a vector of polynomials.

    mat.sort_by(|a, b| <P::Ord>::cmp(&b.lt().deg, &a.lt().deg));

    for i in 0..mat.len() {
        // Get one row. We then use this is cancel other rows
        let lm = mat[i].lm();
        if lm == <P::Var as Zero>::zero() {
            continue
        }
        let mut j = i + 1;

        while j < mat.len() && mat[j].lm() == lm {

            let scalar = mat[j].lc().div(&mat[i].lc()).unwrap();
            mat[i].scale_ass(scalar);

            // Cancel the first term
            mat[j] = mat[j].sub(&mat[i]);

            j += 1
        }
        // Preserve the order
        mat[i..j].sort_by(|a, b| <P::Ord>::cmp(&b.lt().deg, &a.lt().deg));

        println!("Matrix so far = [\n{}\n]", vec_poly_str(mat));

        // Goes back up the matrix in order to put it into reduces row echelon
        for j in (0..i).rev() {
            if let Some(c) = mat[j].has(&lm) {
                // Cancel the lead terms
                let scalar = c.div(&mat[i].lc()).unwrap();
                mat[i].scale_ass(scalar);

                mat[j] = mat[j].sub(&mat[i]);
            }
        }
        println!("Matrix after RREF = [\n{}\n]", vec_poly_str(mat));
    }
}

pub fn choose_subset(pairs: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    // Returns a subset of B based on G, and also removes that subset from B
    // Currently just processing all of it
    let mut a = Vec::new();
    a.append(pairs);
    a
}

use std::collections::BTreeMap;

pub fn compute_m<'a, P: FPolyRing>(l_mat: &mut Vec<Poly<'a, P>>, g: &Ideal<'a, P>) {
    // This line is a bit of a hack
    let g_init = MonomialIdeal::from(&g);

    let mut i = 0;

    let mut seen_monomials: BTreeMap<Term<P>, Vec<usize>> = BTreeMap::new();

    // Iterates through the monomials of the polynomials in the matrix and adds the extra
    // scaled basis elements if necessary
    while i < l_mat.len() {
        let mut aux = Vec::new();

        // Go through the terms and add any that aren't already in there
        for term in l_mat[i].terms.iter() {
            // If term is in G_init, then give back the scaled f
            if let Some(poly) = g_init.get_poly(&term) {
                // Record the row we found this monomial in
                match seen_monomials.get_mut(&term) {
                    Some(x) => x.push(i),
                    None => {
                        // Evaluate x^\alpha f_ell and pushes it the aux buffer
                        aux.push(poly.term_scale(&term.div(poly.lt()).unwrap()));
                        // Record that we have seen it
                        seen_monomials.insert(term.clone(), vec![i]);
                    }
                }
            }
        }
        l_mat.append(&mut aux);
        i += 1;
    }
}

pub fn left_hand_s_poly<'a, P: FPolyRing>(a: &Poly<'a, P>, b: &Poly<'a, P>) -> Poly<'a, P> {
    // Returns the "Left hand side" of the S-Polynomial of a and b
    let lcm = a.lt().lcm(&b.lt());
    let scale = lcm.div(a.lt()).unwrap();
    a.term_scale(&scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::*;
    use crate::parse::MyFromStr;
    use crate::polym::*;
    use generic_array::typenum::{U2, U3};
    use crate::ideals::groebner_basis::*;

    #[test]
    fn f4_test() {
        let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
                Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
                Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
            ];

        let gb = f4(f_vec);
        println!("The Grobner basis = \n{}", gb);
        assert!(is_groebner_basis(&gb));

        let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
                Poly::from_str(&ring, "1.0x^1y^2 + 1.0x^1y^1 - 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2 + 1.0x^1y^1 + 1.0").unwrap(),
            ];

        let gb = f4(f_vec);
        println!("The Grobner basis = \n{}", gb);
        assert!(is_groebner_basis(&gb));
    }

    #[test]
    fn row_reduce_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);

        let mut f_vec = vec![
            Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2 + 2.0x^1y^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2").unwrap(),
            Poly::from_str(&ring, "3.0x^2 + 5.0y^2 + 1.0y^1").unwrap(),
        ];

        row_reduce(&mut f_vec);
        println!("======================");
        println!("{}", vec_poly_str(&f_vec));

    }

    #[test]
    fn compute_m_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, Lex>::new(vec!['x', 'y']);

        let f_vec = vec![
            Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2").unwrap(),
        ];

        let g_basis = Ideal::new(f_vec.clone());
        let t = f_vec.len();
        let subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).filter(|(a, b)| a != b).collect();

        let mut l_mat: Vec<Poly<PRDomain<RR, MultiIndex<U2>, Lex>>> = subsets
            .into_iter()
            .map(|(i, j)| left_hand_s_poly(&g_basis.gens[i], &g_basis.gens[j]))
            .collect();

        println!("{}", vec_poly_str(&l_mat));

        compute_m(&mut l_mat, &g_basis);
        println!("===========================");
        println!("{}", vec_poly_str(&l_mat));
    }
}
