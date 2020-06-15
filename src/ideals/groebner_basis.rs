use crate::algebras::polyring::*;
use crate::ideals::*;

// Check if an ideal is a Groebner basis.
// Uses Buchberger's original criterion with no added optimisations
// That is, it checks that all S-polynomials reduce to zero from division
// by the ideal
impl<P: FPolyRing> Ideal<P> {
    pub fn is_groebner_basis(&self) -> bool {
        let n = self.gens.len();
        (0..n)
            .zip(0..n)
            .filter(|(i, j)| i >= j)
            .map(|(i, j)| s_poly(&self.gens[i], &self.gens[j]).reduce(&self))
            .all(|x| !x.unwrap().is_zero())
    }

    /// Standard implementation of Buchberger's Algorithm
    pub fn bb_algorithm(&self) -> Self {
        let mut gb = self.clone();
        let mut n = gb.gens.len();
        let mut pairs = unord_pairs_int(n);

        // Iterates through all the pairs until there are no more pairs to check
        while !pairs.is_empty() {
            let mut new_pairs = Vec::new();

            for (i, j) in pairs {
                let r = s_poly(&gb.gens[i], &gb.gens[j]).reduce(&gb).unwrap();
                if !r.is_zero() {
                    // Add the remainder to the ideal and update the pairs that we
                    // need to check
                    gb.add(r);
                    n += 1;
                    new_pairs.append(&mut (0..n - 1).map(|k| (k, n)).collect());
                }
            }
            pairs = new_pairs;
        }
        gb
    }

    /// Implementation of the Improved Buchberger's algorithm found at the end of
    /// "Ideals, Varieties, and Algorithms" by Cox, Little, and O'Shea 4th edition.
    pub fn bb_algorithm_plus(&self) -> Self {
        let mut gb = self.clone();
        let mut n = gb.gens.len();
        let mut unchecked_pairs = unord_pairs_int(n);

        // Iterates through all the pairs until there are no more pairs to check
        while !unchecked_pairs.is_empty() {
            let mut new_pairs = Vec::new();

            for (i, j) in unchecked_pairs {
                if gb.gens[i].lt().gcd(&gb.gens[j].lt()) != Term::one() {
                    let r = s_poly(&gb.gens[i], &gb.gens[j]).reduce(&gb).unwrap();
                    if !r.is_zero() {
                        // Add the remainder to the ideal and update the pairs that we
                        // need to check
                        gb.add(r);
                        n += 1;
                        new_pairs.append(&mut (0..n - 1).map(|k| (k, n)).collect());
                    }
                }
            }
            unchecked_pairs = new_pairs;
        }
        gb
    }
}

/// Gets all the unordered pairs between zero and the input n.
/// Not the most efficient implementation but the simplest
fn unord_pairs_int(n: usize) -> Vec<(usize, usize)> {
    (0..n).zip(0..n).filter(|(i, j)| i >= j).collect()
}

/// Calculates the S-polynomial of the pair
fn s_poly<P: FPolyRing>(lhs: &Poly<P>, rhs: &Poly<P>) -> Poly<P> {

    let lcm = lhs.lt().lcm(&rhs.lt());
    let a_newlead = lcm.div(lhs.lt()).unwrap();
    let b_newlead = lcm.div(rhs.lt()).unwrap();

    lhs.term_scale(&a_newlead).sub(&rhs.term_scale(&b_newlead))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::RR;
    use crate::parse::*;
    use crate::polym::*;
    use crate::polyu::*;
    use generic_array::typenum::U2;

    #[test]
    fn division_test() {
        // These tests work
        // Univariate
        let ring = PRDomain::<RR, UniIndex, UnivarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3.0x^2 + 5.0x^98").unwrap();
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        let (q, r) = a.euclid_div(&b).unwrap();
        println!("q = {:?}", q);
        println!("r = {}", r);

        // Multivariate
        let ring = PRDomain::<RR, MultiIndex<U2>, Lex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "3.0x^2y^6 + 5.0x^98y^2").unwrap();
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        println!("a = {}", a);
        println!("b = {}", b);
        let (q, r) = a.euclid_div(&b).unwrap();
        println!("q = {:?}", q);
        println!("r = {}", r);
    }

    #[test]
    fn bb_alg_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = bb_algorithm(&Ideal::new(vec![a, b]));
        println!("r = {:?}", r);
        println!("is GB? = {:?}", is_groebner_basis(&r));
    }

    #[test]
    fn is_gb_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let x = Poly::from_str(&ring, "1.0x^1").unwrap();
        let y = Poly::from_str(&ring, "1.0y^1").unwrap();
        let x_y_ideal = Ideal::new(vec![x, y]);
        assert!(x_y_ideal.is_groebner_basis());

        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let not_gb = Ideal::new(vec![a, b]);
        assert!(!not_gb.is_groebner_basis());
    }
}
