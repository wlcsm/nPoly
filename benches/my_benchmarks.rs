// use criterion::{black_box, criterion_group, criterion_main, Criterion};
// use mycrate::fibonacci;

// pub fn criterion_benchmark(c: &mut Criterion) {
//     c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
// }

// criterion_group!(benches, criterion_benchmark);
// criterion_main!(benches);

// extern crate chrono;
// extern crate rand;
// use chrono::*;
// use rand::distributions::{Distribution, Uniform};

// fn bench_dense() {
//     let ring = PRDomain::<RR, UniVarOrder>::new(vec!['x']);

//     // Note all coefficients are nonzero
//     let dist = Uniform::from(1..100);
//     let mut rng = rand::thread_rng();
//     // A function to randomly generate a polynomial with n coefficients
//     let mut make_poly = |n: usize| -> PolyU<RR> {
//         let res_vec = (0..n).map(|_| RR::from_int(dist.sample(&mut rng))).collect();
//         Poly::from_coeff(&ring, res_vec)
//     };

//     // Benches the time required to multiply two arbitrary polynomials of deg = n
//     let mut time_mult = |n: usize| {
//         let a = make_poly(n);
//         let b = make_poly(n);

//         println!("-------------------------------------------");
//         println!("Number of elements = {}", n);
//         println!("-------------------------------------------");
//         println!(
//             "Karatsuba: {:?}",
//             Duration::span(|| {
//                 karatsuba(&a, &b);
//             })
//         );
//         println!("-------------------------------------------");
//     };

//     for i in 5..16 {
//         time_mult(1 << i);
//     }
// }
