use crate::algebras::polyring::*;
use crate::ideals::*;

use crate::algebras::EuclideanDomain;

pub fn is_groebner_basis<'a, P: FPolyRing>(g: &Ideal<'a, P>) -> bool {
    let n = g.gens.len();
    for i in 0..n {
        for j in i + 1..n {
            // Note: Only need to check that the lead term of r isn't in the initial ideal
            let r = g.gens[i].s_poly(&g.gens[j]).reduce(&g);
            if !r.is_zero() {
                return false;
            }
        }
    }
    true
}

impl<'a, P: FPolyRing> Poly<'a, P> {
    fn pop(&mut self) -> Option<Term<P>> {
        self.terms.pop()
    }
    fn s_poly(&self, other: &Self) -> Self {
        let lcm = self.lt().lcm(&other.lt());
        let a_newlead = lcm.div(self.lt()).unwrap();
        let b_newlead = lcm.div(other.lt()).unwrap();
        self.term_scale(&a_newlead)
            .sub(&other.term_scale(&b_newlead))
    }
}

impl<'a, P: FPolyRing> Poly<'a, P> {
    fn reduce(&self, g: &Ideal<'a, P>) -> Self {
        self.divpolys(&g.gens).1
    }

    pub fn divpolys(&self, f: &Vec<Self>) -> (Vec<Self>, Self) {
        let mut p = self.clone();
        let mut r = self.zero(); // This is a zero polynomial in the same ring as f
        let mut q = vec![self.zero(); f.len()];

        'outer: while !p.is_zero() {
            for i in 0..f.len() {
                if let Some(quo) = p.lt().div(f[i].lt()) {
                    p = p.sub(&f[i].term_scale(&quo));
                    q[i] = q[i].add(&Poly::from_terms(vec![quo], &q[i].ring));
                    continue 'outer;
                }
            }
            // Executes if no division occured
            r = r.add(&Poly::from_terms(vec![p.pop().unwrap()], &self.ring));
        }
        (q, r)
    }
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
        let (q, r) = a.divpolys(&vec![b]);
        println!("q = {:?}", q);
        println!("r = {}", r);

        // Multivariate
        let ring = PRDomain::<RR, MultiIndex<U2>, Lex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "3.0x^2y^6 + 5.0x^98y^2").unwrap();
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        println!("a = {}", a);
        println!("b = {}", b);
        let (q, r) = a.divpolys(&vec![b]);
        println!("q = {:?}", q);
        println!("r = {}", r);
    }

    #[test]
    // fn bb_alg_test() {
    //     let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
    //     let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
    //     let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
    //     let r = bb_algorithm(&Ideal::new(vec![a, b]));
    //     println!("r = {:?}", r);
    //     println!("is GB? = {:?}", is_groebner_basis(&r));
    // }
    #[test]
    fn is_gb_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let x = Poly::from_str(&ring, "1.0x^1").unwrap();
        let y = Poly::from_str(&ring, "1.0y^1").unwrap();
        let x_y_ideal = Ideal::new(vec![x, y]);
        assert!(is_groebner_basis(&x_y_ideal));

        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let not_gb = Ideal::new(vec![a, b]);
        assert!(!is_groebner_basis(&not_gb));
    }
}
