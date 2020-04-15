use crate::algebras::*;
use crate::algebras::polyring::*;
use crate::polyu::*;
use crate::polym::*;

pub struct Ideal<'a, P: PolyRing> {
    gens: Vec<Poly<'a, P>>,
}

pub fn BB_algorithm<P: PolyRing>(basis: Ideal<P>) -> Ideal<P> {
    unimplemented!()
}

fn quo_rem<P: PolyRing>(ring: &P, a: Term<P>, b: Term<P>) -> (Term<P>, Term<P>) {
    if a.divides(b) {
        
    }
}

pub fn divpoly<'a, P>(f: Poly<'a, P>, F: Vec<Poly<'a, P>>) -> (Vec<Poly<'a, P>>, Poly<'a, P>) 
    where P: PolyRing {

    let f_vec = f.terms.clone();
    let r = vec![]; // This is a zero polynomial in the same ring as f
    let q: Vec<Vec<Term<P>>> = vec![Vec::new(); F.len()];

    while !f.is_zero() {
        let i = 0;
        let divisionoccured = false;
        while i < F.len() && !divisionoccured {

            let (quo, rem) = quo_rem(f.ring, f.lt(), F[i].lt());
            if rem.coeff == <P::Coeff>::zero() {
                q[i].push(quo);
                f_vec.pop();
                divisionoccured = true;
            } else {
                i += 1;
            }
        }
        if !divisionoccured {
            r.push(f_vec.pop().unwrap()); // Transfer the lead term from f to r
        }
    }
    let q_poly = q.into_iter().map(|v| Poly::from_terms_unchecked(v, f.ring)).collect();
    let r_poly = Poly::from_terms_unchecked(r, f.ring);
    (q_poly, r_poly)
}