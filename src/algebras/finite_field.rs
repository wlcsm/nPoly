use typenum::marker_traits::Unsigned;
use std::marker::PhantomData;
use crate::algebras::*;
use typenum::{U7};

trait Prime: Unsigned + Copy + Eq {}
impl Prime for U7 {}

#[derive(Eq, PartialEq, Clone, Copy)]
struct FF<U: Prime> (i64, PhantomData<U>);

impl<U: Prime> FF<U> {
    // Automatically reduces the number
    fn new(q: i64) -> Self {
        FF(q % U::to_i64() , PhantomData)
    }

    #[warn(dead_code)]
    fn new_unchecked(q: i64) -> Self {
        FF(q, PhantomData)
    }

    fn red(q: i64) -> i64 {
        q % U::to_i64()
    }
}

use std::fmt;

impl<U: Prime> fmt::Debug for FF<U> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "F({})", self.0)
    }
}

impl<U: Prime> fmt::Display for FF<U> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<U: Prime> Ring for FF<U> {
    type BaseRing = Self;

    fn add(&self, other: &Self) -> Self {
        FF::new(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        FF::new(self.0 - other.0)
    }
    fn neg(&self) -> Self {
        FF::new(-self.0)
    }
    fn zero() -> Self {
        FF::new(0)
    }
    fn mul(&self, other: &Self) -> Self {
        FF::new(self.0 * other.0)
    }
    fn one() -> Self {
        FF::new(1)
    }
}

impl<U: Prime> ScalarRing for FF<U> {
    fn add_ass(&mut self, other: &Self) {
        self.0 = <FF<U>>::red(self.0 + other.0)
    }
    fn sub_ass(&mut self, other: &Self) {
        self.0 = <FF<U>>::red(self.0 - other.0)
    }
    fn mul_ass(&mut self, other: &Self) {
        self.0 = <FF<U>>::red(self.0 * other.0)
    }
}


#[cfg(test)]
mod tests {

    use super::*;
    use typenum::{U7};

    #[test]
    fn basics() {
        let a = FF::<U7>::new(5);
        let b = FF::<U7>::new(-2);
        let c = FF::<U7>::new(10);

        println!("a = {}, b = {}, c = {}", a, b, c);
    }
}