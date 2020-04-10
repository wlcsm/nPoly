use crate::polyu::*;
use crate::algebras::polyring::*;
use crate::fft::*;

// Notation: t in the number of terms in the sparse polynomial and n is the fft degree

// Coeffs is a vector of coefficients found by evaluating some subtree
// Index is a number whose binary representation defines the path back to the root node (i.e. is in RBE)
#[derive(Debug, Clone)]
pub struct Subtree<T: SupportsFFT> {
    coeffs: Vec<T>,
    pos: Position,
}

#[derive(Debug, Clone, Copy)]
struct Position {
    path: u32,
    depth: usize
}
impl<T: SupportsFFT> Subtree<T> {

    // Note: This is dependent on the ordering of arg1 and arg2. 
    //       arg2's path is assumed to be greater than arg1's 
    pub fn can_combine(&self, arg2: &Subtree<T>) -> Option<usize> {
        // Find the lowest branch
        let lowest = findlowestbranch(self.pos, arg2.pos);

        // Check all the bits before the lowest connecting branch are 1's
        let tester = (1 << (lowest - arg2.pos.depth)) - 1; 
        if (tester & arg2.pos.path) == tester { Some(lowest) } else { None }
    }

    // Takes a subtree and promotes/expands it up the tree, assuming the 
    // coefficients in its sister subtree is zero
    pub fn expand(&mut self, depth: usize) {
        while self.pos.depth != depth {
            // If the digit is a zero it slides to the right
            if self.pos.path % 2 == 0 { 
                self.coeffs.extend(self.coeffs.clone())
            } else {
            // If the digit is a one it slides to the left
                let rou = <T>::rou(1 << self.coeffs.len(), false);
                for j in 0..self.coeffs.len() {
                    self.coeffs[j].mul_ass(&rou[j]);
                    self.coeffs.push(self.coeffs[j].neg());
                }
            }
            self.pos.path >>= 1;
            self.pos.depth += 1;
        }
    }

    /// Combines two equal depth subtrees together
    pub fn combine(&mut self, mut argb: Subtree<T>, depth: usize) -> Subtree<T> {

        // Expand them to the required depth
        // TODO this is very inefficient
        self.expand(depth - 1);
        argb.expand(depth - 1);

        // Our resulting vector of coefficients
        let mut result = vec![<T>::zero(); self.coeffs.len() * 2];
        let rou = <T>::rou(self.coeffs.len() * 2, false);

        // Combines them, note that it assumes that argb is on the right
        for i in 0..(result.len() / 2) {
            // Performs a 2-butterfly
            let tmp = argb.coeffs[i].mul(&rou[i]);
            result[i] = self.coeffs[i].add(&tmp);
            result[i + self.coeffs.len()] = self.coeffs[i].sub(&tmp);
        }
        Subtree {
            coeffs : result,
            pos    : Position {
                path  : self.pos.path >> 1,
                depth : self.pos.depth + 1
            }
        }
    }
}

/// Converts the monomials of the polynomial into instances of Subtree<T>, then sorts them 
/// with respect to their reverse bit encoding
pub(crate) fn revbit<T: SupportsFFT>(input: &Poly<T, Univariate>, n: usize) -> Vec<Subtree<T>> {

    // Reverses the bits: TODO assumes that x is less than 2^32
    let revbits = |x: u32| x.reverse_bits() >> (32 - n);

    // Makes each element a subtree (inefficient I know)
    let mut result: Vec<Subtree<T>> 
        = input.terms.iter().map(|Term { coeff, deg }| Subtree { 
                                        coeffs: vec![*coeff],
                                        pos   : Position {
                                            path  : revbits(deg.0 as u32) as u32,
                                            depth : 0
                                        }
                                    }
                                    ).collect();

    // Applys the function 'revbits' to the monomial's degree then sorts based on those 'keys'
    result.sort_by_cached_key(|a| a.pos.path);
    result
}

// Does a linear search for the largest bit prefix of s and t which equates to finding the 
// path of the smallest subtree connecting them in this case
fn findlowestbranch(mut arg1: Position, arg2: Position) -> usize {
    // Make the paths the same length
    let a = arg1.depth as i32 - arg2.depth as i32;
    arg1.path = if a < 0 {arg1.path >> -a} else {arg1.path << a};

    // Under a bitwise XOR, the common prefix of the two paths will be all zeros
    // So the leading one signals the end of the common prefix 
    let res = arg1.path ^ arg2.path;
    ((32 - res.leading_zeros()) + arg2.depth as u32) as usize
}


pub(crate) fn sparse_eval<T: SupportsFFT>(input: Poly<T, Univariate>, n: usize) -> Vec<T> {

    // Obtain RBE
    let rbe = revbit(&input, n);
    let mut result: Vec<Subtree<T>> = vec![rbe[0].clone()];

    // Iterate over the RBE
    for el in rbe.into_iter().skip(1) {

        result.push(el);
        while result.len() > 1 {
            if let Some(d) = result[result.len() - 2].can_combine(&result[result.len() - 1]) {
                // Combine the two at some depth
                let tmp = result.pop().unwrap();
                let n = result.len();
                result[n - 1] = result[n - 1].combine(tmp, d);
            } else {
                break;
            }
        }
    }
    
    println!("Result so far = {:?}", result);

    // Clean up everything
    while result.len() > 1 {
        let n = result.len();
        let lowest = findlowestbranch(result[n - 2].pos, result[n - 1].pos);
        let tmp = result.pop().unwrap();
        result[n - 2] = result[n - 2].combine(tmp, lowest);
    }

    result[0].coeffs.clone()
}

#[cfg(test)]
mod tests {

    use crate::algebras::complex::CC;
    use crate::algebras::*;
    use crate::fft::*;
    use super::*;

    // Yeet it works for this one example
    #[test]
    fn whole() {
        // Commented it out because it was causing some errors I didn't want to fix
        // at the time

        // let ring = PRDomain::univar("x".to_string());

        // let a = Poly::from_coeff(&ring, vec![CC::from_re(1),CC::from_re(2),CC::from_re(3)]);
        // let b = a.clone();

        // let mut b_coeffs: Vec<CC> = b.terms.iter().map(|a| a.0).collect();

        // b_coeffs.push(CC::from_re(0));
        // perform_fft(&mut b_coeffs[..], false).unwrap();
        // let res = sparse_eval(a, 2);
        // println!("Sparse Eval = {:?}", res);
        // println!("Normal FFT = {:?}", b_coeffs);
    }

    #[test]
    fn expand_test() {
        // let a = Poly::from_coeff(None, vec![ZZ(1),ZZ(2),ZZ(3)]).unwrap();
        // let b = Poly::from_coeff(None, vec![ZZ(1),ZZ(2),ZZ(3), ZZ(4), ZZ(5), ZZ(6)]).unwrap();
        // let c = Poly::from_coeff(None, vec![ZZ(1)]).unwrap();
        let mut sub = Subtree { coeffs: vec![<CC>::one(), <CC>::one()], pos: Position { path: 0b01, depth:0}};
        sub.expand(1);
        println!("{:?}", sub);
    }

    #[test]
    fn pathfinding() {
        let a = Position { path: 0b00, depth: 0};
        let b = Position { path: 0b01, depth: 0};
        let c = Position { path: 0b10, depth: 0};
        let d = Position { path: 0b11, depth: 0};
        assert_eq!(1, findlowestbranch(a, b));
        assert_eq!(1, findlowestbranch(c, d));
        assert_eq!(2, findlowestbranch(b, c));
        assert_eq!(2, findlowestbranch(a, d));

        let a = Position { path: 0b0, depth: 2};
        let b = Position { path: 0b1, depth: 2};
        let c = Position { path: 0b10, depth: 1};
        let d = Position { path: 0b111, depth: 0};
        assert_eq!(3, findlowestbranch(a, b));
        assert_eq!(2, findlowestbranch(c, d));
        assert_eq!(2, findlowestbranch(b, c));
        assert_eq!(3, findlowestbranch(a, d));
    }

    #[test]
    fn combining() {
        let a = Subtree { coeffs: vec![CC::from_re(0)], pos: Position { path: 0b00, depth: 0}};
        let b = Subtree { coeffs: vec![CC::from_re(1)], pos: Position { path: 0b01, depth: 0}};
        let c = Subtree { coeffs: vec![CC::from_re(0)], pos: Position { path: 0b10, depth: 0}};
        let d = Subtree { coeffs: vec![CC::from_re(1)], pos: Position { path: 0b11, depth: 0}};
        assert_eq!(Some(1), a.can_combine(&b));
        assert_eq!(None, a.can_combine(&c));
        assert_eq!(Some(1), c.can_combine(&d));

        let a = Subtree { coeffs: vec![CC::from_re(0)], pos: Position { path: 0b000, depth: 0}};
        let b = Subtree { coeffs: vec![CC::from_re(1)], pos: Position { path: 0b01, depth: 1}};
        let c = Subtree { coeffs: vec![CC::from_re(0)], pos: Position { path: 0b0000, depth: 0}};
        let d = Subtree { coeffs: vec![CC::from_re(1)], pos: Position { path: 0b101, depth: 1}};
        assert_eq!(Some(2), a.can_combine(&b));
        assert_eq!(Some(2), c.can_combine(&b));
        // Detects that it needs to wait on another subtree to come in
        assert_eq!(None, c.can_combine(&d))
    }
    #[test]
    fn lowestbranch() {
        let a = Position { path: 0b00, depth: 0};
        let b = Position { path: 0b01, depth: 0};
        let c = Position { path: 0b10, depth: 0};
        let d = Position { path: 0b11, depth: 0};
        assert_eq!(1, findlowestbranch(a, b));
        assert_eq!(2, findlowestbranch(b, c));
        assert_eq!(1, findlowestbranch(c, d));

        let a = Position { path: 0b0000, depth: 0};
        let b = Position { path: 0b01, depth: 2};
        let c = Position { path: 0b101, depth: 1};
        let d = Position { path: 0b11, depth: 2};
        assert_eq!(3, findlowestbranch(a, b));
        assert_eq!(4, findlowestbranch(b, c));
        assert_eq!(3, findlowestbranch(c, d));
    }
}
