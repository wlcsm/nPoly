use crate::algebras::polyring::*;
use crate::polyu::*;
use generic_array::*;
use num_traits::Zero;

/// This will be used to hold the data to perform Kronecker substitution
/// The i^th element corresponds to the value that was substituted into the i^th variable
///
/// This is calculated by taking the maximum degree of each indeterminate in each of the
/// polynomials, adding them togather to obtain the largest possible degree of the indeterminates
/// in the result, and adding one to each value.
///
/// Let d(i) by the result for the ith indeterminate.
/// Then we do the substitution
/// x^{a(1), a(2), \ldots, a(n)}
/// ->
/// x^{a(1) + a(2)*d(1) + a(3)*d(1)d(2) + \ldots + a(n)*d(1)\ldots d(n-1)}
///
/// Therefore we will store (d(1), d(2), \ldots, d(n)) in the 'sub_values' field in
/// the KroneckerData struct.
///
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KroneckerData<'a, P: PolyRing> {
    sub_values: GenericArray<usize, P::NumVar>,
    ring: &'a P,
}

impl<'a, P: PolyRing> KroneckerData<'a, P> {
    fn undo_sub(&self, term: &Term<PRDomain<P::Coeff, UniVarOrder>>) -> Term<P> {
        let mut deg_usize = term.mon.0;
        let new_mon = self
            .sub_values
            .iter()
            .map(|d| {
                let rem = deg_usize.rem_euclid(*d);
                deg_usize = deg_usize.div_euclid(*d);
                rem
            })
            .collect();
        Term::new(term.coeff, new_mon)
    }

    fn sub(&self, term: &Term<P>) -> Term<PRDomain<P::Coeff, UniVarOrder>> {
        let mut acc = 1;
        let new_deg = term.mon.clone().iter().zip(self.sub_values.iter())
            .map(|(a, b)| {
                let prev = acc;
                acc *= b;
                a * prev
            })
            .sum();

        Term {
            coeff: term.coeff,
            mon: UniIndex(new_deg),
        }
    }
}


/// Uses Kronecker substitution to convert the multivariate polynomial to a univariate one
///
/// Need to find the maximum of the individual degree of the monomials, to take the
/// elementwise maximum of two monomial, we can use the lcm function.
pub fn create_kron_data<'a, P>(poly_a: &Poly<'a, P>, poly_b: &Poly<'a, P>) -> KroneckerData<'a, P>
where
    P: PolyRing,
{
    let get_max_degrees = |poly: &Poly<'a, P>| -> P::Mon {
        poly.terms
            .iter()
            .fold(P::Mon::zero(), |acc, term| acc.lcm(&term.mon))
    };

    let max_degs_a = get_max_degrees(poly_a);
    let max_degs_b = get_max_degrees(poly_b);

    // Add them together plus 1. We need a plus one because we need to substitute the term that is
    // strictly larger than the others
    let sub_data = max_degs_a.iter().zip(max_degs_b.iter())
        .map(|(a, b)| a + b + 1)
        .collect();

    KroneckerData {
        sub_values: sub_data,
        ring: poly_a.ring.unwrap(),
    }
}

/// Creates a Univariate polynomial out of a multivariate one and instructions for its
/// Kronecker substitution
///
/// FIXME We clone the monomial each time. This is very inefficient! If we use the
/// MultiIndex monomial struct, then we would easily avoid this because we could implement
/// the Iterator trait, but because of more exotic monomial representations we can't
/// implement Iterator on all of them.
pub fn to_univar<'a, P: PolyRing>(
    poly_a: &Poly<'a, P>,
    kron_data: &KroneckerData<P>,
) -> PolyU<'a, P::Coeff> {
    let new_terms = poly_a
        .terms
        .iter()
        .map(|term| kron_data.sub(term))
        .collect();
    Poly::from_terms(new_terms, None)
}

/// Uses Kronecker substitution to convert the multivariate polynomial to a univariate one
pub fn to_multivar<'a, P: PolyRing>(
    poly_a: PolyU<'a, P::Coeff>,
    kron_data: &KroneckerData<'a, P>,
) -> Poly<'a, P> {
    let terms = poly_a
        .terms
        .iter()
        .map(|term| kron_data.undo_sub(term))
        .collect();
    Poly::from_terms(terms, Some(kron_data.ring))
}

pub fn kronecker<'a, P: PolyRing>(
    poly_a: &Poly<'a, P>,
    poly_b: &Poly<'a, P>,
) -> (
    PolyU<'a, P::Coeff>,
    PolyU<'a, P::Coeff>,
    KroneckerData<'a, P>,
) {
    let kron_data = create_kron_data(&poly_a, &poly_b);
    let poly_a_uni = to_univar(&poly_a, &kron_data);
    let poly_b_uni = to_univar(&poly_b, &kron_data);

    (poly_a_uni, poly_b_uni, kron_data)
}

pub fn kronecker_mult<'a, P: PolyRing>(poly_a: &Poly<'a, P>, poly_b: &Poly<'a, P>) -> Poly<'a, P> {
    let (poly_a_uni, poly_b_uni, kron_data) = kronecker(poly_a, poly_b);
    to_multivar(poly_a_uni * poly_b_uni, &kron_data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::*;
    use crate::parse::*;
    use crate::polym::*;
    use typenum::U2;

    #[test]
    fn create_kron_data_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();

        let kron_data = create_kron_data(&a, &b);

        assert_eq!(
            KroneckerData {
                sub_values: arr![usize; 6, 4],
                ring: &ring
            },
            kron_data
        );
    }

    #[test]
    fn to_univar_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();

        let kron_data = create_kron_data(&a, &b);
        let new_poly_a = to_univar(&a, &kron_data);
        let new_poly_b = to_univar(&b, &kron_data);

        println!("new_poly_a = {}", new_poly_a);
        println!("new_poly_b = {}", new_poly_b);
        // assert_eq!(Poly::from_str(&None, "1.0x^3 - 2x^7").unwrap().terms, new_poly_a.terms);
    }

    #[test]
    fn to_multivar_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();

        let (poly_a_uni, poly_b_uni, kron_data) = kronecker(&a, &b);

        println!("New Poly_a = {}", poly_a_uni);
        println!("New Poly_b = {}", poly_b_uni);

        let poly_a_multi = to_multivar(poly_a_uni.clone(), &kron_data);
        let poly_b_multi = to_multivar(poly_b_uni.clone(), &kron_data);

        assert_eq!(to_univar(&poly_a_multi, &kron_data), poly_a_uni);
        assert_eq!(to_univar(&poly_b_multi, &kron_data), poly_b_uni);
    }

    #[test]
    fn mult_kron_test() {
        let ring = PRDomain::<RR, GLex<MultiIndex<U2>>>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();

        let (poly_a_uni, poly_b_uni, kron_data) = kronecker(&a, &b);

        let res_uni = poly_a_uni * poly_b_uni;
        let res_multi = to_multivar(res_uni, &kron_data);

        assert_eq!(a * b, res_multi);
    }
}
