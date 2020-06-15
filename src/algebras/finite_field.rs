use crate::algebras::*;
use std::marker::PhantomData;
use typenum::marker_traits::Unsigned;
use typenum::U7;

trait Prime: Unsigned + Copy + Eq {}
impl Prime for U7 {}

#[derive(Clone, Copy)]
struct FF<P: Prime>(i64, PhantomData<P>);

// I need to reimplment this because I'm allowing the numbers in the FF struct to be
// stored positively or negatively, so I need to account for these two forms
impl<P: Prime> PartialEq for FF<P> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0 || <FF<P>>::red(self.0 - other.0) == 0
    }
}

impl<P: Prime> std::cmp::Eq for FF<P> {}

impl<P: Prime> FF<P> {
    // Automatically reduces the number
    pub fn new(q: i64) -> Self {
        FF(q % P::to_i64(), PhantomData)
    }

    #[allow(dead_code)]
    pub fn new_unchecked(q: i64) -> Self {
        FF(q, PhantomData)
    }

    pub fn red(q: i64) -> i64 {
        q % P::to_i64()
    }
}

use std::fmt;

impl<P: Prime> fmt::Debug for FF<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "F({})", self.0)
    }
}

impl<P: Prime> fmt::Display for FF<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<P: Prime> Zero for FF<P> {
    fn zero() -> Self {
        FF::new(0)
    }
}
impl<P: Prime> One for FF<P> {
    fn one() -> Self {
        FF::new(1)
    }
}

impl<P: Prime> Group for FF<P> {
    fn add(&self, other: &Self) -> Self {
        FF::new(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        FF::new(self.0 - other.0)
    }
    fn neg(&self) -> Self {
        FF::new(-self.0)
    }
}

impl<P: Prime> Ring for FF<P> {
    fn mul(&self, other: &Self) -> Self {
        FF::new(self.0 * other.0)
    }
}

impl<P: Prime> std::str::FromStr for FF<P> {
    type Err = std::num::ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.parse::<i64>() {
            Ok(n) => Ok(FF::new(n)),
            Err(e) => Err(e),
        }
    }
}

impl<P: Prime> ScalarRing for FF<P> {
    const REGEX: &'static str = r"-?\d+";
    fn add_ass(&mut self, other: &Self) {
        self.0 = <FF<P>>::red(self.0 + other.0)
    }
    fn sub_ass(&mut self, other: &Self) {
        self.0 = <FF<P>>::red(self.0 - other.0)
    }
    fn mul_ass(&mut self, other: &Self) {
        self.0 = <FF<P>>::red(self.0 * other.0)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use typenum::U7;

    #[test]
    fn basics() {
        let a = FF::<U7>::new(5);
        let b = FF::<U7>::new(-2);
        let c = FF::<U7>::new(10);

        println!("a = {}, b = {}, c = {}", a, b, c);
        assert_eq!(FF::<U7>::new(3), a.add(&b));
        assert_eq!(FF::<U7>::new(0), a.sub(&b));
        assert_eq!(FF::<U7>::new(-2), c.sub(&a));

        // Since one is in a negative representation and the other is positive
        assert_eq!(a, b);
    }
}
