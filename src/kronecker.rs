use crate::algebras::ScalarRing;

/// This will be used to hold the data to perform Kronecker substitution
/// The i^th element corresponds to the value that was substituted into the i^th variable
///
/// This is calculated by taking the maximum degree of each indeterminate in each of the
/// polynomials, adding them togather to obtain the largest possible degree of the indeterminates
/// in the result, and adding one to each value.
///
/// NOTE: There is a potentially quicker way, it is probably more efficient to just reinterpret the
/// multi-index of the monomial as the mixed-radix representation, rather than actually convert it. 
/// Computationally, nothing would happen, but we would need to make the univariate algorithms
/// accept the interpretation (which is actually not too hard)
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KroneckerData<const VARS: usize>([usize; VARS]);

impl<const VARS: usize> KroneckerData<VARS> {
    pub fn new<R: ScalarRing, const VARS: usize>(
        poly_a: &Poly<R, VARS>,
        poly_b: &Poly<R, VARS>,
    ) -> KroneckerData<VARS> {
        // Find the maximum degree of each indeterminate in the polynomials
        let max_deg = |poly: &Poly<R, VARS>| -> [R; VARS] {
            poly.terms
                .iter()
                .fold([0; VARS], |acc, term| acc.lcm(&term.mon))
        };

        let max_degs_a = max_deg(poly_a);
        let max_degs_b = max_deg(poly_b);

        let max_degrees = izip!(max_degs_a.iter(), max_degs_b.iter())
            .map(|(a, b)| a + b + 1)
            .collect();

        KroneckerData(max_degrees)
    }

    /// Converts a univariate term into a multivariate term
    #[inline(always)]
    fn uni_to_multi(&self, term: &Term<R, 1>) -> Term<R, VARS> {
        let mut deg = term.mon[0];
        let new_mon = self
            .0
            .iter()
            .map(|d| {
                let rem = deg.rem_euclid(*d);
                deg = deg.div_euclid(*d);
                rem
            })
            .collect();
        Term::new(term.coeff, new_mon)
    }

    /// Converts a multivariate term into a univariate term
    #[inline(always)]
    fn multi_to_uni(&self, term: &Term<R, VARS>) -> Term<R, 1> {
        let mut acc = 1;
        let mut sum = 0;

        for i in 0..VARS {
            sum += deg * acc;
            acc *= kron;
        }

        Term::new(term.coeff, [sum])
    }

    /// Creates a Univariate polynomial out of a multivariate one and instructions for its
    /// Kronecker substitution
    pub fn to_univar(&self, poly: &Poly<R, VARS>) -> Poly<R, 1> {
        let new_terms = poly.terms.iter().map(|term| kron_data.multi_to_uni(term)).collect();

        Poly::from_terms(new_terms)
    }

    /// Uses Kronecker substitution to convert the multivariate polynomial to a univariate one
    pub fn to_multivar(&self, poly_a: Poly<R, 1>) -> Poly<R, VARS> {
        let terms = poly_a
            .terms
            .iter()
            .map(|term| kron_data.uni_to_multi(term))
            .collect();

        Poly::from_terms(terms, Some(kron_data.ring))
    }
}

pub fn kronecker_mult<'a, P: PolyRing>(poly_a: &Poly<'a, P>, poly_b: &Poly<'a, P>) -> Poly<'a, P> {
    let (poly_a_uni, poly_b_uni, kron_data) = kronecker(poly_a, poly_b);
    to_multivar(poly_a_uni * poly_b_uni, &kron_data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_kron_data_test() {
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();

        let kron_data = create_kron_data(&a, &b);

        assert_eq!( KroneckerData(arr![usize; 6, 4]), kron_data);
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
