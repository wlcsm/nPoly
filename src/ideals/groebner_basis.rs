use crate::algebras::polyring::*;
use crate::algebras::*;

#[derive(Clone)]
pub struct Ideal<'a, P: PolyRing> {
    gens: Vec<Poly<'a, P>>,
}

impl<'a, P: PolyRing> Ideal<'a, P> {
    fn push(&mut self, item: Poly<'a, P>) {
        self.gens.push(item)
    }
}

pub fn bb_algorithm<F: Field, P: PolyRing<Coeff = F>>(g: Ideal<P>) -> Ideal<P> {
    // "New" are the newly added polynomials that we need to check
    // "Acc" is the accumulator for the new Groebner basis
    // "Next" are the next values that we need to transfer into the "New" vector
    let mut new = g.gens.clone();
    let mut acc = g.clone();

    while !new.is_empty() {
        let mut next = Vec::new();
        while let Some(f) = new.pop() {
            for g in &acc.gens {
                let r = f.s_poly(&g).reduce(&acc);
                if r.is_zero() {
                    next.push(r)
                }
            }
            acc.push(f);
        }
        new = next;
    }
    g
}

impl<'a, F: Field, P: PolyRing<Coeff = F>> Poly<'a, P> {
    fn pop(&mut self) -> Option<Term<P>> {
        self.terms.pop()
    }
    fn push(&mut self, item: Term<P>) {
        self.terms.push(item)
    }
    fn s_poly(&self, other: &Self) -> Self {
        let gcd = self.lt().gcd(&other.lt());
        let a_newlead = gcd.div(self.lt()).unwrap();
        let b_newlead = gcd.div(other.lt()).unwrap();
        self.term_scale(&a_newlead)
            .sub(&other.term_scale(&b_newlead))
    }
}

impl<'a, F: Field, P: PolyRing<Coeff = F>> Poly<'a, P> {
    fn reduce(&self, g: &Ideal<'a, P>) -> Self {
        self.divpolys(&g.gens).1
    }

    pub fn divpolys(&self, f: &Vec<Self>) -> (Vec<Self>, Self) {
        let mut p = self.clone();
        let mut r = self.zero(); // This is a zero polynomial in the same ring as f
        let mut q = vec![self.zero(); f.len()];

        'outer: while !p.is_zero() {
            // println!("p = {}", p);
            // println!("r = {}", r);
            // println!("q = {:?}", q);
            for i in 0..f.len() {
                // println!("f[i] = {}", f[i]);
                // println!("p LT = {:?}", p.lt());
                // println!("LT(p) / LT(f[i]) = {:?}", p.lt().div(f[i].lt()));
                if let Some(quo) = p.lt().div(f[i].lt()) {
                    // println!("Inner");
                    // println!("f[i] = {}", f[i]);
                    // println!("quo = {:?}", quo);
                    p = p.sub(&f[i].term_scale(&quo));
                    q[i].push(quo);
                    continue 'outer;
                }
            }
            r.push(p.pop().unwrap()); // Executes if no division occured
        }
        (q, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::RR;
    use crate::parse::*;
    use crate::polyu::*;

    #[test]
    fn division_test() {
        let ring = PRDomain::<RR, UniIndex, UnivarOrder>::new(vec!['x']);
        let a = Poly::from_str(&ring, "3.0x^2 + 5.0x^98").unwrap();
        println!("{:?}", a);
        println!("{}", a);
        let b = Poly::from_str(&ring, "5.0x^6").unwrap();
        let (q, r) = a.divpolys(&vec![b]);
        println!("q = {:?}", q);
        println!("r = {}", r);
    }
}
