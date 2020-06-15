pub mod complex;
pub mod finite_field;
pub mod integers;
pub mod polyring;
pub mod real;

/// The original motivation for making my own copies of the sub and add traits is so that
/// I can borrow self rather than taking ownership

// The group trait is used in the MonomialIndex trait
pub trait Group: Zero + Sized + Eq + Clone {
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
}

// For the moment I'm also assuming that the rings are integral domains
pub trait Ring: Group + One {
    // Ring operations
    fn mul(&self, other: &Self) -> Self;
}

pub trait ScalarRing:
    Ring + Copy + std::fmt::Debug + std::str::FromStr + std::fmt::Display
{
    const REGEX: &'static str;
    fn add_ass(&mut self, other: &Self);
    fn sub_ass(&mut self, other: &Self);
    fn mul_ass(&mut self, other: &Self);
}

// Don't know if I should have a function for division by multiple arguments.
// It would just be annoying if
pub trait EuclideanDomain: Ring {
    fn gcd(&self, other: &Self) -> Self;
    fn lcm(&self, other: &Self) -> Self;
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)>;

    // Iteratively calls the standard "euclid_div" method on the remainder
    // of the previous division.
    // Returns None if one of the quotients is zero
    fn euclid_div_multi(&self, others: &Vec<&Self>) -> Option<(Vec<Self>, Self)> {
        assert!(!others.is_empty());
        // Guards against the case that any of the quotients are zero
        if others.iter().any(|x| x.is_zero()) {
            None
        } else {
            // The first calculation is done manually because the rem variables
            // needs to be an owned type, and we have only borrowed self.
            let (q_init, r_init) = self.euclid_div(&others[0]).unwrap();
            Some(
                others
                    .iter()
                    .skip(1)
                    .fold((vec![q_init], r_init), |(mut quo, rem), poly| {
                        let (q, r) = rem.euclid_div(poly).unwrap();
                        quo.push(q);
                        (quo, r)
                    }),
            )
        }
    }
    // self | other. Returns None if self is zero
    fn divides(&self, other: &Self) -> Option<bool> {
        if self.is_zero() {
            None
        } else {
            Some(self.euclid_div(other).unwrap().1.is_zero())
        }
    }
}

// Fields are also Euclidean domains and the euclid_div can be easily derived
// from the "div" function for Field
impl<F: Field> EuclideanDomain for F {
    fn gcd(&self, _other: &Self) -> Self {
        self.clone()
    }
    fn lcm(&self, _other: &Self) -> Self {
        self.clone()
    }
    fn euclid_div(&self, other: &Self) -> Option<(Self, Self)> {
        // Remainder is always zero
        self.div(other).map(|s| (s, <F>::zero()))
    }
}

pub trait Field: ScalarRing {
    fn div(&self, other: &Self) -> Option<Self>;
}

pub trait Zero: Eq + PartialEq + Sized {
    fn zero() -> Self;
    fn is_zero(&self) -> bool {
        *self == <Self>::zero()
    }
}

pub trait One: Eq + PartialEq + Sized {
    fn one() -> Self;
}
