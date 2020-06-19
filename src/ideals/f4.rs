use crate::algebras::polyring::*;
use crate::algebras::*;
use crate::ideals::*;

fn vec_poly_str<'a, P: PolyRing>(poly: &Vec<Poly<'a, P>>) -> String {
    poly.iter()
        .map(|p| format!("{}", p))
        .collect::<Vec<String>>()
        .join("\n")
}

// The F4 algorithm
pub fn f4<'a, P: FPolyRing>(f_gens: Vec<Poly<'a, P>>) -> Ideal<'a, P> {
    let mut g_basis = Ideal::new(f_gens.clone());
    let mut t = f_gens.len();
    // println!("f+gens = {}", vec_poly_str(&f_gens));
    // println!("t = {}", t);
    let mut subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).filter(|(a, b)| a != b).collect();

    while !subsets.is_empty() {
        let selection: Vec<(usize, usize)> = choose_subset(&mut subsets);

        // println!("selection = {:?}", selection);
        // println!("g_basis = \n{}", g_basis);
        let mut l_mat: Vec<Poly<'a, P>> = selection
            .into_iter()
            .map(|(i, j)| left_hand_s_poly(&g_basis.gens[i], &g_basis.gens[j]))
            .collect();

        compute_m(&mut l_mat, &g_basis);
        let m_init = MonomialIdeal::new(l_mat.iter().map(|m| m.lt()).collect());

        // println!("compute_m = \n{}", vec_poly_str(&l_mat));

        // Row reduction and find the new elements
        row_reduce(&mut l_mat);
        let n_plus: Vec<Poly<'a, P>> = l_mat
            .into_iter()
            .filter(|n| !n.is_zero() && !m_init.is_in(&n.lt()))
            .collect();

        // println!("n_plus = \n{}", vec_poly_str(&n_plus));

        // println!("subsets = \n{:?}", subsets);
        for n in n_plus {
            // println!("next to be added to n = {}", n);
            g_basis.add(n);
            // TODO for some reason this failed, i.e. the .add() method didn't actually add
            // anything to the list of generators, this is weird because by the construction
            // of n_plus it definitely should have
            // t += 1; 
            t = g_basis.gens.len();
            let mut new_pairs_l = (0..t-1).map(|new| (new, t-1)).collect();
            let mut new_pairs_r = (0..t-1).map(|new| (t-1, new)).collect();
            subsets.append(&mut new_pairs_l);
            subsets.append(&mut new_pairs_r);
        }
        // println!("len = {}", g_basis.num_gens());
        // println!("t = {}", t);
        // println!("subsets = \n{:?}", subsets);
    }
    println!("Groebner basis = {}", vec_poly_str(&g_basis.gens));
    println!("Is Grobner = {}", g_basis.is_groebner_basis());

    row_reduce(&mut g_basis.gens);
    g_basis.gens = g_basis.gens.into_iter().filter(|x| !x.is_zero()).collect();
    println!("Row reduced = {}", vec_poly_str(&g_basis.gens));
    println!("Is Grobner again = {}", g_basis.is_groebner_basis());
    g_basis
}

// Error log
// mat[i] = -2xy + -1x^2
// mat[j] = -2xy + -1x^2
// Here
// thread 'ideals::f4::tests::bench_f4' panicked at 'attempt to add with overflow', /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/ops/arith.rs:94:45
// stack backtrace:
//    0: backtrace::backtrace::libunwind::trace
//              at /Users/runner/.cargo/registry/src/github.com-1ecc6299db9ec823/backtrace-0.3.46/src/backtrace/libunwind.rs:86
//    1: backtrace::backtrace::trace_unsynchronized
//              at /Users/runner/.cargo/registry/src/github.com-1ecc6299db9ec823/backtrace-0.3.46/src/backtrace/mod.rs:66
//    2: std::sys_common::backtrace::_print_fmt
//              at src/libstd/sys_common/backtrace.rs:78
//    3: <std::sys_common::backtrace::_print::DisplayBacktrace as core::fmt::Display>::fmt
//              at src/libstd/sys_common/backtrace.rs:59
//    4: core::fmt::write
//              at src/libcore/fmt/mod.rs:1069
//    5: std::io::Write::write_fmt
//              at src/libstd/io/mod.rs:1532
//    6: std::sys_common::backtrace::_print
//              at src/libstd/sys_common/backtrace.rs:62
//    7: std::sys_common::backtrace::print
//              at src/libstd/sys_common/backtrace.rs:49
//    8: std::panicking::default_hook::{{closure}}
//              at src/libstd/panicking.rs:198
//    9: std::panicking::default_hook
//              at src/libstd/panicking.rs:218
//   10: std::panicking::rust_panic_with_hook
//              at src/libstd/panicking.rs:477
//   11: rust_begin_unwind
//              at src/libstd/panicking.rs:385
//   12: core::panicking::panic_fmt
//              at src/libcore/panicking.rs:89
//   13: core::panicking::panic
//              at src/libcore/panicking.rs:52
//   14: <usize as core::ops::arith::Add>::add
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/ops/arith.rs:94
//   15: <usize as core::ops::arith::Add<&usize>>::add
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/internal_macros.rs:45
//   16: core::ops::function::FnMut::call_mut
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/ops/function.rs:154
//   17: core::iter::traits::iterator::Iterator::fold::ok::{{closure}}
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2002
//   18: core::iter::traits::iterator::Iterator::try_fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:1878
//   19: core::iter::traits::iterator::Iterator::fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2005
//   20: <usize as core::iter::traits::accum::Sum<&usize>>::sum
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/accum.rs:62
//   21: core::iter::traits::iterator::Iterator::sum
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2805
//   22: n_poly::polym::MultiIndex<N>::new
//              at src/polym.rs:19
//   23: <n_poly::polym::MultiIndex<N> as n_poly::algebras::polyring::Variate>::add
//              at src/polym.rs:136
//   24: n_poly::algebras::polyring::Poly<P>::term_scale::{{closure}}
//              at src/algebras/polyring.rs:343
//   25: core::iter::adapters::map_fold::{{closure}}
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/adapters/mod.rs:785
//   26: core::iter::traits::iterator::Iterator::fold::ok::{{closure}}
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2002
//   27: core::iter::traits::iterator::Iterator::try_fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:1878
//   28: core::iter::traits::iterator::Iterator::fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2005
//   29: <core::iter::adapters::Map<I,F> as core::iter::traits::iterator::Iterator>::fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/adapters/mod.rs:825
//   30: core::iter::traits::iterator::Iterator::for_each
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:658
//   31: <alloc::vec::Vec<T> as alloc::vec::SpecExtend<T,I>>::spec_extend
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/liballoc/vec.rs:2116
//   32: <alloc::vec::Vec<T> as alloc::vec::SpecExtend<T,I>>::from_iter
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/liballoc/vec.rs:2096
//   33: <alloc::vec::Vec<T> as core::iter::traits::collect::FromIterator<T>>::from_iter
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/liballoc/vec.rs:1981
//   34: core::iter::traits::iterator::Iterator::collect
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:1660
//   35: n_poly::algebras::polyring::Poly<P>::term_scale
//              at src/algebras/polyring.rs:340
//   36: n_poly::ideals::f4::left_hand_s_poly
//              at src/ideals/f4.rs:170
//   37: n_poly::ideals::f4::f4::{{closure}}
//              at src/ideals/f4.rs:27
//   38: core::iter::adapters::map_fold::{{closure}}
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/adapters/mod.rs:785
//   39: core::iter::traits::iterator::Iterator::fold::ok::{{closure}}
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2002
//   40: core::iter::traits::iterator::Iterator::try_fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:1878
//   41: core::iter::traits::iterator::Iterator::fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:2005
//   42: <core::iter::adapters::Map<I,F> as core::iter::traits::iterator::Iterator>::fold
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/adapters/mod.rs:825
//   43: core::iter::traits::iterator::Iterator::for_each
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:658
//   44: <alloc::vec::Vec<T> as alloc::vec::SpecExtend<T,I>>::spec_extend
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/liballoc/vec.rs:2116
//   45: <alloc::vec::Vec<T> as alloc::vec::SpecExtend<T,I>>::from_iter
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/liballoc/vec.rs:2096
//   46: <alloc::vec::Vec<T> as core::iter::traits::collect::FromIterator<T>>::from_iter
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/liballoc/vec.rs:1981
//   47: core::iter::traits::iterator::Iterator::collect
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/iter/traits/iterator.rs:1660
//   48: n_poly::ideals::f4::f4
//              at src/ideals/f4.rs:25
//   49: n_poly::ideals::f4::tests::bench_f4::{{closure}}
//              at src/ideals/f4.rs:232
//   50: time::duration::Duration::span
//              at /Users/williamcashman/.cargo/registry/src/github.com-1ecc6299db9ec823/time-0.1.42/src/duration.rs:138
//   51: n_poly::ideals::f4::tests::bench_f4
//              at src/ideals/f4.rs:231
//   52: n_poly::ideals::f4::tests::bench_f4::{{closure}}
//              at src/ideals/f4.rs:213
//   53: core::ops::function::FnOnce::call_once
//              at /Users/williamcashman/.rustup/toolchains/nightly-x86_64-apple-darwin/lib/rustlib/src/rust/src/libcore/ops/function.rs:232
//   54: <alloc::boxed::Box<F> as core::ops::function::FnOnce<A>>::call_once
//              at /rustc/7ebd87a7a1e0e21767422e115c9455ef6e6d4bee/src/liballoc/boxed.rs:1034
//   55: <std::panic::AssertUnwindSafe<F> as core::ops::function::FnOnce<()>>::call_once
//              at /rustc/7ebd87a7a1e0e21767422e115c9455ef6e6d4bee/src/libstd/panic.rs:318
//   56: std::panicking::try::do_call
//              at /rustc/7ebd87a7a1e0e21767422e115c9455ef6e6d4bee/src/libstd/panicking.rs:297
//   57: std::panicking::try
//              at /rustc/7ebd87a7a1e0e21767422e115c9455ef6e6d4bee/src/libstd/panicking.rs:274
//   58: std::panic::catch_unwind
//              at /rustc/7ebd87a7a1e0e21767422e115c9455ef6e6d4bee/src/libstd/panic.rs:394
//   59: test::run_test_in_process
//              at src/libtest/lib.rs:541
//   60: test::run_test::run_test_inner::{{closure}}
//              at src/libtest/lib.rs:450
// note: Some details are omitted, run with `RUST_BACKTRACE=full` for a verbose backtrace.
// test ideals::f4::tests::bench_f4 ... FAILED

// failures:

// failures:
//     ideals::f4::tests::bench_f4

// test result: FAILED. 0 passed; 1 failed; 0 ignored; 0 measured; 14 filtered out

// error: test failed, to rerun pass '--lib'

pub fn row_reduce<'a, P: FPolyRing>(mat: &mut Vec<Poly<'a, P>>) {
    // Standard row reduction, but operates on a vector of polynomials.

    // mat.sort_by(|a, b| <P::Ord>::cmp(&b.lt().deg, &a.lt().deg));
    let mut i = 0;

    while i < mat.len() {
        // Get one row. We then use this is cancel other rows
        let lm = mat[i].lm();
        if lm == <P::Var as Zero>::zero() {
            i += 1;
            continue
        }
        let mut j = i + 1;

        while j < mat.len() && mat[j].lm() == lm {

            let scalar = mat[j].lc().div(&mat[i].lc()).unwrap();
            mat[i].scale_ass(scalar);

            // Cancel the first term
            mat[j] = mat[j].sub(&mat[i]);

            j += 1
        }
        // Preserve the order

        // println!("Matrix so far = [\n{}\n]", vec_poly_str(mat));

        // Goes back up the matrix in order to put it into reduces row echelon
        for j in (0..i).rev() {
            if let Some(c) = mat[j].has(&lm) {
                // Cancel the lead terms
                let scalar = c.div(&mat[i].lc()).unwrap();
                mat[i].scale_ass(scalar);

                mat[j] = mat[j].sub(&mat[i]);
            }
        }
        // println!("Matrix after RREF = [\n{}\n]", vec_poly_str(mat));
        i += 1;
    }

    mat.sort_by(|a, b| <P::Ord>::cmp(&b.lt().deg, &a.lt().deg));
}

pub fn choose_subset(pairs: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    // Returns a subset of B based on G, and also removes that subset from B
    // Currently just processing all of it
    let mut a = Vec::new();
    a.append(pairs);
    a
}

use std::collections::BTreeMap;

pub fn compute_m<'a, P: FPolyRing>(l_mat: &mut Vec<Poly<'a, P>>, g: &Ideal<'a, P>) {
    // This line is a bit of a hack
    let g_init = MonomialIdeal::from(&g);

    let mut i = 0;

    let mut seen_monomials: BTreeMap<Term<P>, Vec<usize>> = BTreeMap::new();

    // Iterates through the monomials of the polynomials in the matrix and adds the extra
    // scaled basis elements if necessary
    while i < l_mat.len() {
        let mut aux = Vec::new();

        // Go through the terms and add any that aren't already in there
        for term in l_mat[i].terms.iter() {
            // If term is in G_init, then give back the scaled f
            if let Some(poly) = g_init.get_poly(&term) {
                // Record the row we found this monomial in
                match seen_monomials.get_mut(&term) {
                    Some(x) => x.push(i),
                    None => {
                        // Evaluate x^\alpha f_ell and pushes it the aux buffer
                        aux.push(poly.term_scale(&term.div(poly.lt()).unwrap()));
                        // Record that we have seen it
                        seen_monomials.insert(term.clone(), vec![i]);
                    }
                }
            }
        }
        l_mat.append(&mut aux);
        i += 1;
    }
}

pub fn left_hand_s_poly<'a, P: FPolyRing>(a: &Poly<'a, P>, b: &Poly<'a, P>) -> Poly<'a, P> {
    // Returns the "Left hand side" of the S-Polynomial of a and b
    let lcm = a.lt().lcm(&b.lt());
    println!("a = {}", a);
    let scale = lcm.div(a.lt()).unwrap();
    a.term_scale(&scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebras::real::*;
    use crate::parse::MyFromStr;
    use crate::polym::*;
    use generic_array::typenum::{U2, U3};

    #[test]
    fn f4_test() {
        let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
                Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
                Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
            ];

        let gb = f4(f_vec);
        assert!(gb.is_groebner_basis());

        let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
                Poly::from_str(&ring, "1.0x^1y^2 + 1.0x^1y^1 - 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2 + 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2 + 1.0x^1y^1 + 1.0").unwrap(),
            ];

        let gb = f4(f_vec);
        assert!(gb.is_groebner_basis());
    }

    use chrono::*;

    #[test]
    fn bench_f4() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);
        let a = Poly::from_str(&ring, "1.0x^3 - 2.0x^1y^1").unwrap();
        let b = Poly::from_str(&ring, "1.0x^2y^1 - 2.0y^2 + 1.0x^1").unwrap();
        let r = &Ideal::new(vec![a, b]);
        println!("BB Alg time = {:?}", 
            Duration::span(|| {
                f4(r.gens.clone());
            }));

        let ring = PRDomain::<RR, MultiIndex<U3>, GLex>::new(vec!['x', 'y', 'z']);
        let f_vec = vec![
                Poly::from_str(&ring, "1.0x^2 + 1.0x^1y^1 - 1.0").unwrap(),
                Poly::from_str(&ring, "1.0x^2 - 1.0z^2").unwrap(),
                Poly::from_str(&ring, "1.0x^1y^1 + 1.0").unwrap(),
            ];

        println!("BB Alg time = {:?}", 
            Duration::span(|| {
                f4(f_vec);
            }));

        // let f_vec = vec![
        //         Poly::from_str(&ring, "1.0x^3 + 1.0x^1y^1 - 1.0").unwrap(),
        //         Poly::from_str(&ring, "1.0x^2y^1 - 1.0z^2 + 1.0").unwrap(),
        //         Poly::from_str(&ring, "1.0x^1y^1 + 1.0 + 2.0z^1").unwrap(),
        //     ];

        // println!("BB Alg time = {:?}", 
        //     Duration::span(|| {
        //         f4(f_vec);
        //     }));
    }

    #[test]
    fn row_reduce_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, GLex>::new(vec!['x', 'y']);

        let mut f_vec = vec![
            Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2 + 2.0x^1y^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2").unwrap(),
            Poly::from_str(&ring, "3.0x^2 + 5.0y^2 + 1.0y^1").unwrap(),
        ];

        row_reduce(&mut f_vec);
        println!("======================");
        println!("{}", vec_poly_str(&f_vec));

    }

    #[test]
    fn compute_m_test() {
        let ring = PRDomain::<RR, MultiIndex<U2>, Lex>::new(vec!['x', 'y']);

        let f_vec = vec![
            Poly::from_str(&ring, "1.0x^1y^1 + 2.0x^1").unwrap(),
            Poly::from_str(&ring, "1.0y^2").unwrap(),
        ];

        let g_basis = Ideal::new(f_vec.clone());
        let t = f_vec.len();
        let subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).filter(|(a, b)| a != b).collect();

        let mut l_mat: Vec<Poly<PRDomain<RR, MultiIndex<U2>, Lex>>> = subsets
            .into_iter()
            .map(|(i, j)| left_hand_s_poly(&g_basis.gens[i], &g_basis.gens[j]))
            .collect();

        println!("{}", vec_poly_str(&l_mat));

        compute_m(&mut l_mat, &g_basis);
        println!("===========================");
        println!("{}", vec_poly_str(&l_mat));
    }
}


// % I have actually made my own way for this. It leverages the fact that we only care about collecting the polynomials whose lead terms were not originally in $G$. Notice that whenever we add a $x^\alpha f_\ell$ into the matrix, this polynomial will be excluded in the end. The only way this wouldn't happen is if while we are reducing the matrix, we find another intermediate polynomial has the same lead monomial as $x^\alpha f_\ell$ and then it can cancel its lead term. But actually, in the process of matrix reduction, we could choose $x^\alpha f_\ell$ to cancel the other polynomials lead term instead. Hence the trick here is to actually avoid unnecessarily making the matrix larger

// % \begin{enumerate}[(1)]
// %     \item Start with the largest polynomial in terms of lead term size in the sorted list
// %     \item Cancel out the lead terms. Remember the lead terms of the resulting polynomials, they are now in the unsorted pile
// %     \item Put this polynomial in the Done list
// %     \item Then add the polynomials with the highest lead term from the sorted list into the unsorted one, also adding its lead term.
// %     \item Select the next biggest lead term, and check if it is in the initial ideal of G. If this was the one from the sorted list, add the next one from the sorted list also into it.
// %     \item  - If it is, then use the polynomial from $G$ to subtract the terms and also from in the done pile. We don't actually include this polynomial in any of the piles because we know that it won't end up in $N^+$.
// %     \item  - If not, then subtract any other polynomials and then add it into the Done pile
// %     \item Remembering to add the new lead terms into the monomial list as well. If the unsorted pile is empty, then add the polynomials will the highest lead terms from the sorted list
// % \end{enumerate}

