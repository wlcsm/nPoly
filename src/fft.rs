
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



use num_complex::Complex;
use crate::*;
use std::f64::consts::PI;

pub fn fft_mult(poly_a: &PolyU, poly_b: &PolyU) -> PolyU {

    // Maximum number of coefficients in a * b
    let n = round_2pow(poly_a.deg() + poly_b.deg() + 1);


    // Take the DFT of the coefficients of the polynomials
    // I feel okay calling unwrap here because there are only two ways an error
    // occurs. 1. The signal is empty which is impossible because the terms of a
    // polynomial are always nonempty, or the length isn't a power of two, but that
    // is guarenteeds by line 7
    let a_fft = perform_fft(&to_coeffs_complex(&poly_a.terms[..], n)[..], false).unwrap();
    let b_fft = perform_fft(&to_coeffs_complex(&poly_b.terms[..], n)[..], false).unwrap();

    // Multiply elementwise
    let temp: Vec<Complex<f64>> = a_fft.into_iter().zip(b_fft.into_iter())
                                        .map(|(a,b)| a * b).collect();
                                        
    // Calculate inverse fft
    let c = perform_fft(&temp[..], true).unwrap();
    // Convert back into i32
    let c_parsed = c.into_iter().map(|x| x.re.round() as i32).collect();

    // Convert back into polynomial type
    PolyU::from_coeff("x".to_string(), c_parsed).unwrap()
}

pub fn to_coeffs_complex(input: &[Monomial], n: usize) -> Vec<Complex<f64>> {
    // Expands the input into a list with only coefficients to a desired 
    // size n. If n is larger than it then it will pad with zeros
    let mut result = Vec::with_capacity(n);
    let mut i = 0;
    
    for mono in input.into_iter() {
        // Between two monomials we fill the vector with zeros
        if mono.deg > i {
            (i .. mono.deg).into_iter().for_each(|_| result.push(Complex::from(0.0)));
        } 
        // Push the monomial and increase the counter
        result.push(Complex::from(mono.coeff as f64));
        i = mono.deg + 1;
    }
    // Pad the rest
    (i .. n).into_iter().for_each(|_| result.push(Complex::from(0.0)));

    result
}


fn perform_fft(sample: &[Complex<f64>], inverse: bool) 
                        -> Result<Vec<Complex<f64>>, &'static str> {
    // Sample is the signal. Note that is fully expanded
    // Note: Assumes that the length of 'sample' is a power of two and greater
    //       than 4
    let not2pow = |n| (n - 1) & n != 0;

    match sample.len() {
        0 => Err("Signal cannot be empty"),
        1 => Ok(vec![sample[0]]),
        2 => Ok(vec![sample[0] + sample[1], sample[0] - sample[1]]),
        n if not2pow(n) => 
             Err("Signal needs to be a power of two"),
        _ => Ok(go_fast(sample, inverse)),
    }
}

fn go_fast(sample: &[Complex<f64>], inverse: bool) -> Vec<Complex<f64>> {

    let n = sample.len();
    let rou = gen_rou(n, inverse);
    let mut result: Vec<Complex<f64>> = Vec::with_capacity(n);
    
    // Preprocess: The first layer is particularly easy because we don't need
    // complex numbers. We then also take this opportunity to put it in the 
    // reverse bit order.
    let mut aux   : Vec<Complex<f64>> = Vec::with_capacity(n/2);

    let mut i = 0;
    while i < n / 2 {
        let x_0 = sample[i];
        let x_1 = sample[i + 1];
        let x_0_n2 = sample[i + (n / 2)];
        let x_1_n2 = sample[i + (n / 2) + 1];
        result.push(x_0 + x_0_n2);
        result.push(x_0 - x_0_n2);
        aux.push(x_1 + x_1_n2);
        aux.push(x_1 - x_1_n2);
        i += 2;
    }
    result.append(&mut aux);

    // From here on we assume the first layer of the bottom up approach has been 
    // completed and it is in the reverse bit order

    // TODO Need to handle the option type better
    // Starts at index two because we already handled the first one
    println!("flutter before = {:?}\n", result);
    // let modx = |p, q| (((p * q) % n) + n) % n;
    for i in 2..=log2(n).unwrap() {
        // Iterate over all the flutters
        for j in 0..(n >> i) {
            // println!("i = {}, j = {}", i, j);
            // println!("flutter before = {:?}", result);
            // Only need to iterate over half of the flutter
            let half_flut = 1 << (i - 1);
            println!("i = {}, j = {}, half_flut = {}", i, j, half_flut);
            for k in (j * (1 << i))..(j * (1 << i))+half_flut {
                println!(" rou index {}\n", k % (1 << i));
                let a = result[k];
                let b = result[k + half_flut];
                result[k] = a + rou[(k % (1 << i)) * (n >> i)] * b; // a + w^j * b
                result[k + half_flut] = a - rou[(k % (1 << i)) * (n >> i)] * b; // a - w^j * b (= a + w^{j + n/2} * b)
            }
            println!("flutter after = {:?}\n", result);
        }
    }
    // println!("Result now = {:?}", result);
    if inverse {
        result.into_iter().map(|x| x / (n as f64)).collect()
    } else {
        result
    }
}

fn round_2pow(n: usize) -> usize {
    // Rounds up to the nearest power of two.
    // Basically, it uses bitwise operations to find the largest non-zero bit in
    // the int, then it return the number that has one higher than that and zeros
    // everywhere else. 
    // I think this could be redone if we rewrote the log2 function to round up 
    // instead of returning a None value if it isn't a power of two
    match n {
        0 => n,
        // If its already a power of two, return it
        x if (x - 1) & n == 0 => x,
        mut x => {
            let mut pow = 0;
            while x != 0 {
                x = x >> 1;
                pow += 1
            }
            (1 << pow)
        }
    }
}

fn gen_rou(n: usize, inverse: bool) -> Vec<Complex<f64>> {
    // Generates all the nth roots of unity 
    // Changes it depending on whether computing the dft or the inverse
    let sign = if inverse {1.0} else {-1.0};
    let base = Complex::new(0.0, sign * 2.0 * PI / n as f64);
    // Iterator magic to make it more idiomatic
    // let scaling = if inverse { 1.0 / n as f64} else {1.0};
    // (0..n).map(|k| base.scale(k as f64).exp().scale(scaling)).collect()
    (0..n).map(|k| base.scale(k as f64).exp()).collect()
}

// Does log2 for integers, if parameter is not a power of two
// then it returns None
// TODO Make this nicer
pub fn log2(n: usize) -> Option<usize> {

    // Very quick test to see if its a power of two
    // The test will actually fail if num = 1, but thats okay because
    // we still want to return Some(0) if num = 1
    match n {
        0 => None,
        1 => Some(0),
        _ if (n - 1) & n != 0 => None,
        mut num => {
            let mut result = 0;
            while num != 1 {
                num = num >> 1;
                result += 1;
            };
            Some(result)
        },
    }
}



#[cfg(test)]
mod tests {
    use super::dft::*;
    use super::fft::*;
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
    #[test]
    fn log() {
        assert_eq!(Some(3), log2(8));
        assert_eq!(Some(6), log2(64));
        assert_eq!(Some(0), log2(1));
        assert_eq!(None, log2(9));
        assert_eq!(None, log2(0));
    }

    #[test]
    fn conversion() {
        let a = PolyU::from_coeff("x".to_string(), vec![1,1]).unwrap();
        let b = PolyU::from_coeff("x".to_string(), vec![1,3]).unwrap();

        let mut a_coeff = vec![Complex::from(1.0); 2];
        assert_eq!(a_coeff, to_coeffs_complex(&a.terms[..], 2));
        a_coeff.append(&mut vec![Complex::from(0.0); 3]);
        assert_eq!(a_coeff, to_coeffs_complex(&a.terms[..], 5));

        let mut b_coeff = vec![Complex::from(1.0), Complex::from(3.0)];
        assert_eq!(b_coeff, to_coeffs_complex(&b.terms[..], 2));
        b_coeff.append(&mut vec![Complex::from(0.0); 5]);
        assert_eq!(b_coeff, to_coeffs_complex(&b.terms[..], 7));
    }

    #[test]
    fn fft_mult_test() {

        let a = PolyU::from_coeff("x".to_string(), vec![1,1]).unwrap();
        let b = PolyU::from_coeff("x".to_string(), vec![1,3]).unwrap();
        let c = PolyU::from_coeff("x".to_string(), vec![1,2,1]).unwrap();

        assert_eq!(a.mul(&b), fft_mult(&a, &b));
        assert_eq!(b.mul(&c), fft_mult(&b, &c));
        assert_eq!(c.mul(&a), fft_mult(&c, &a));


        let d = PolyU::from_coeff("x".to_string(), vec![-1,3]).unwrap();
        let e = PolyU::from_coeff("x".to_string(), vec![-1,3,4,6]).unwrap();
        assert_eq!(d.mul(&a), fft_mult(&d, &a));
        assert_eq!(d.mul(&e), fft_mult(&d, &e));

        let f = PolyU::from_coeff("x".to_string(), vec![0]).unwrap();
        assert_eq!(a.mul(&f), fft_mult(&a, &f));
        
    }
}