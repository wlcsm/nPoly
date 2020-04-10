
// Some basic math utilites. This is only really recessary when the base is not
// a power of two, otherwise one should just use fast bit operations

// Rounds down
pub fn logn(num: usize, n: usize) -> usize { 
    if num < n { 0 } else {logn(num / n, n) + 1}
}

pub fn is_n_pow(num: usize, n: usize) -> bool {
    n.pow(logn(num, n) as u32) == num 
}

pub fn next_npow(num: usize, n: usize) -> usize {
    if is_n_pow(num, n) { num } else { n.pow((logn(num, n) + 1) as u32) }
}

// Exponentiation function: n^exp
// pub fn pow(n: usize, exp: usize) -> usize {
//     match exp {
//         0 => 1,
//         1 => n,
//         e => pow(n, e / 2) * pow(n, e / 2) * (if e % 2 == 0 {1} else {n})
//     }
// }