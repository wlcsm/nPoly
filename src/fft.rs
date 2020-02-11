
extern crate num_complex;


// use num_complex::Complex;
// use std::f64::consts::PI;


pub mod dft {
    use num_complex::Complex;
    use crate::*;
    use std::f64::consts::PI;

    pub fn dft_mult(poly_a: &PolyU, poly_b: &PolyU) -> PolyU {

        // Maximum number of coefficients in a * b
        let n = poly_a.deg() + poly_b.deg() + 1;

        // Take the DFT of the coefficients of the polynomials
        let a_dft = dft(&poly_a.terms[..], n);
        let b_dft = dft(&poly_b.terms[..], n);

        // Multiply elementwise
        let temp: Vec<Complex<f64>> = a_dft.into_iter().zip(b_dft.into_iter())
                                           .map(|(a,b)| a * b).collect();
        let c = inv_dft(&temp[..], n);
        let c_parsed = c.into_iter().map(|x| x.re.round() as i32).collect();
        PolyU::from_coeff("x".to_string(), c_parsed).unwrap()
    }

    // Performs the DFT on the samples
    pub fn dft(sample: &[Monomial], n: usize) -> Vec<Complex<f64>> {

        let mut result = Vec::new();

        // Generate the nth roots of unity
        let rou = gen_rou(n, false);

        // Evaluates: p * q mod n.
        // Rust's % does remainder (i.e. allows negative numbers) so this
        // corrects that for positive remainders.
        let modx = |p, q| (((p * q) % n) + n) % n;

        // F(k) = \sum^n_{j=0} x_j e^{-2\pi i jk / n}
        for k in 0..n {
            let term = sample.iter()
                             .map(|x| rou[modx(k, x.deg)] * x.coeff as f64)
                             .sum();
            result.push(term);
        }

        result
    }

    pub fn inv_dft(sample: &[Complex<f64>], n: usize) -> Vec<Complex<f64>> {

        let mut result = Vec::new();

        // Generate the nth roots of unity
        let rou = gen_rou(n, true);

        // Evaluates: p * q mod n.
        // Rust's % does remainder (i.e. allows negative numbers) so this
        // corrects that for positive remainders.
        let modx = |p, q| (((p * q) % n) + n) % n;

        // F(k) = \sum^n_{j=0} x_j e^{2\pi i jk / n}
        for k in 0..n {
            // Notice that we iterate over the roots in reverse order
            // because its the inverse DFT
            let term: Complex<f64> = sample.iter().enumerate()
                             .map(|(i, c)| rou[modx(k, i)] * c)
                             .sum();
            result.push(term.scale( 1.0 / n as f64));
        }

        result
    }

    fn gen_rou(n: usize, inverse: bool) -> Vec<Complex<f64>> {
        // Generates all the nth roots of unity 
        // Changes it depending on whether computing the dft or the inverse
        let sign = if inverse {1.0} else {-1.0};
        let base = Complex::new(0.0, sign * 2.0 * PI / n as f64);
        // Iterator magic to make it more idiomatic
        (0..n).map(|k| base.scale(k as f64).exp()).collect()
    }

}

#[cfg(test)]
mod tests {
    use super::dft::*;
    use crate::*;
    use polyu::*;

    #[test]
    fn dft_mult_test() {

        let a = PolyU::from_coeff("x".to_string(), vec![1,1]).unwrap();
        let b = PolyU::from_coeff("x".to_string(), vec![1,3]).unwrap();
        let c = PolyU::from_coeff("x".to_string(), vec![1,2,1]).unwrap();

        assert_eq!(a.mul(&b), dft_mult(&a, &b));
        assert_eq!(b.mul(&c), dft_mult(&b, &c));
        assert_eq!(c.mul(&a), dft_mult(&c, &a));


        let d = PolyU::from_coeff("x".to_string(), vec![-1,3]).unwrap();
        let e = PolyU::from_coeff("x".to_string(), vec![-1,3,4,6]).unwrap();
        assert_eq!(d.mul(&a), dft_mult(&d, &a));
        assert_eq!(d.mul(&e), dft_mult(&d, &e));

        let f = PolyU::from_coeff("x".to_string(), vec![0]).unwrap();
        assert_eq!(a.mul(&f), dft_mult(&a, &f));
        
    }
}

// pub fn fft(sample: &[i32]) -> Vec<Complex<f64>> {
//     let rou = generate_rou(8);
//     let n = sample.len();
//     let mut result: Vec<Complex<f64>> = Vec::new();

//     // Clone and convert everything into Complex<f64>
//     for i in 0..n {
//         result.push(Complex::new(sample[i].into(), 0.0));
//     };

//     // TODO Need to handle the option type better
//     for i in 1..log2(n).unwrap() {
//         let flutter_width = n/i;
//         for j in 0..i {
//             eval_k_flutter(result[flutter_width*j ..  flutter_width * (j+1)]);
//         }
//     }
//     result
// }

// fn eval_k_flutter(&mut flutter: [Complex<f64>]) -> () {
//     let n = flutter.len();
//     let rou = gen_rou_single(log2(n).unwrap());

//     for i in 0..n/2 {
//         let a = flutter[i];
//         let b = flutter[i + n/2];
//         flutter[i] = a + rou * b;
//         flutter[i+n/2] = a - rou * b;
//     }
// }


// fn gen_rou_single(num: usize) -> Complex<f64> {
//     let x = 2.0 * PI / (1 << num) as f64;
//     let cos = x.cos();
//     let sin = x.sin();
//     Complex::new(cos, sin)
// }

// // Does log2 for integers, if parameter is not a power of two
// // then it returns None
// fn log2(num: usize) -> Option<usize> {
//     if (num - 1) & num != 0 {
//         None
//     } else {
//         let result = 0;
//         while num != 0 {
//             num = num >> 1;
//             result += 1;
//         };
//         Some(result)
//     }
// }

// fn rec(sample: &[i32], dep: usize, rou: &Vec<Complex<f64>>)
//                                     -> Vec<Complex<f64>> {
//     if sample.len() == 1 {
//         vec![Complex::new(sample[0] as f64, 0.0)]
//     } else {
//         let a = rec(&sample[..sample.len()/2], dep + 1, rou);
//         let b = rec(&sample[sample.len()/2..], dep + 1, rou);
//         a.iter().zip(b.iter())
//                 .map(|(aa, bb)| aa + rou[dep] * bb)
//                 .collect()
//     }
// }

// pub fn add(slice_a: Vec<f64>, slice_b: Vec<f64>) -> Vec<f64> {
//     slice_a.iter().zip(slice_b.iter())
//                   .map(|(a, b)| a + b).collect()
// }

// pub fn generate_rou(num: usize) -> Vec<Complex<f64>> {
//     let mut result: Vec<Complex<f64>> = Vec::new();

//     for i in 0..num {
//         let x = 2.0 * PI / (1 << i) as f64;
//         let cos = x.cos();
//         let sin = x.sin();
//         let acc = Complex::new(cos, sin);
//         result.push(acc);
//     }
//     result
// }
