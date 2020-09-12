// TODO Should setup something like Criterion.rs to do benchmarks.

extern crate chrono;
extern crate rand;
use crate::algebras::real::*;
use rand::distributions::uniform::{SampleBorrow, SampleUniform, UniformFloat, UniformSampler};
use rand::prelude::Rng;

// /// Test parameters for RR
// /// I haven't implemented the "dense" option yet, but my idea is that we will be able to specify
// /// the density of the polynomials. By giving a percentage we will then use another random
// /// distribution to choose whether the next number should be zero or not. This is a problem if we
// /// use something like ZZ because if the min is -1 and max is 1, then we already have a 1/3 chance
// /// of choosing a zero anyway, so we need to take this into account.
// struct TestParam<R: ScalarRing> {
//     min: R,
//     max: R,
//     dense: f64,
// }

#[derive(Clone, Copy, Debug)]
pub struct MyUniformDistRR(UniformFloat<f64>);

impl UniformSampler for MyUniformDistRR {
    type X = RR;
    fn new<B1, B2>(low: B1, high: B2) -> Self
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        MyUniformDistRR(UniformFloat::<f64>::new(low.borrow().0, high.borrow().0))
    }
    fn new_inclusive<B1, B2>(low: B1, high: B2) -> Self
    where
        B1: SampleBorrow<Self::X> + Sized,
        B2: SampleBorrow<Self::X> + Sized,
    {
        UniformSampler::new(low, high)
    }
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        RR(self.0.sample(rng))
    }
}

impl SampleUniform for RR {
    type Sampler = MyUniformDistRR;
}

// #[cfg(test)]
// mod test {

//     use super::*;
//     use rand::distributions::uniform::UniformSampler;
//     use crate::algebras::polyring::*;
//     use chrono::Duration;
//     use crate::polyu::*;

//     #[test]
//     fn bench_karatsuba() {
//         let res = bench_dense(vec![2]);
//         println!("{:?}", res);
//     }

//     fn bench_dense(sizes: Vec<usize>) -> Vec<(usize, Duration)> {
//         let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);

//         let mut rng = rand::thread_rng();

//         let distribution = MyUniformDistRR::new(RR(-10.0), RR(10.0));

//         let mut res: Vec<(usize, Duration)> = Vec::new();

//         for size in sizes.into_iter() {
//             let a = PolyU::from_coeff(&ring,
//                           (0..size).map(|_| distribution.sample(&mut rng)).collect());
//             let b = PolyU::from_coeff(&ring,
//                           (0..size).map(|_| distribution.sample(&mut rng)).collect());

//             // res.push((size,
//             //  Duration::span(|| {
//             //     karatsuba::karatsuba(&a, &b);
//             // })))
//         }
//         res
//     }
// }
