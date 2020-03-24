
// Don't think this abstraction is really necessary but I like the organisation
// I didn't idiot proof these so it is easy to break them with edge cases

// Rounds down
pub fn logn(num: usize, n: usize) -> usize { 
    if num < n { 0 } else {logn(num / n, n) + 1}
}

pub fn is_n_pow(num: usize, n: usize) -> bool {
    if n == 2 {
        (num - 1) & num == 0
    } else {
        pow(n, logn(num, n)) == num 
    }
}

pub fn next_npow(num: usize, n: usize) -> usize {
    if is_n_pow(num, n) { num } else { pow(n, logn(num, n) + 1) }
}

// Exponentiation function: n^exp
pub fn pow(n: usize, exp: usize) -> usize {
    if n == 2 {
        1 << exp
    } else {
        match exp {
            0 => 1,
            1 => n,
            e => pow(n, e / 2) * pow(n, e / 2) * (if e % 2 == 0 {1} else {n})
        }
    }
}