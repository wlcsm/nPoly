use crate::ideals::groebner_basis::*;
use crate::algebras::*;
use crate::algebras::polyring::*;


pub fn f4<'a, P: FPolyRing>(F: Vec<Poly<'a, P>>) -> Ideal<'a, P> {
    let mut G = Ideal::new(F.clone());
    let mut t = F.len();
    let mut subsets: Vec<(usize, usize)> = iproduct!(0..t, 0..t).collect();

    while !subsets.is_empty() {
        let selection: Vec<(usize, usize)> = choose_subset(&mut subsets);

        let mut L: Vec<Poly<'a, P>> = selection.into_iter().map(|(i, j)| left_hand_s_poly(&G.gens[i], &G.gens[j]) ).collect();
        ComputeM(&mut L, &G);
        let m_init = MonomialIdeal::new(L.iter().map(|m| m.lt()).collect());

        // Row reduction and find the new elements
        row_reduce(&mut L);
        let n_plus: Vec<Poly<'a, P>> = L.into_iter().filter(|n| !M_init.is_in(&n.lt())).collect();

        for n in n_plus {
            t += 1;
            G.add(n);
            let mut new_pairs = (1..t).map(|n| (n, t)).collect();
            subsets.append(&mut new_pairs);
        }
    }
    G
}

pub fn row_reduce<'a, P: FPolyRing>(mat: &mut Vec<Poly<'a, P>>) {
    // Standard row reduction. Assumes that the polynomial are sorted in descending 
    // order based on their lead coefficient.

    for i in 0..mat.len() {

        let lt = mat[i].lm();
        let mut j = i;

        while mat[j].lm() == lt {
            let scalar = mat[j].lc().div(&mat[i].lc()).unwrap();
            mat[i].scale_ass(scalar);

            // Cancel the first term
            mat[j] = mat[j].sub(&mat[i]);

            j += 1
        }
        // Preserve the order
        mat[i..j+1].sort_by(|a, b| <P::Ord>::cmp(&a.lt().deg, &b.lt().deg));

        // Put it into reduces row echelon
        for j in i..mat.len() {
            if let Some(c) = mat[j].has(&lt) {
                // Cancel the lead terms
                let scalar = mat[j].lc().div(&mat[i].lc()).unwrap();
                mat[i].scale_ass(scalar);

                mat[j] = mat[j].sub(&mat[i]);
            }
        }
    }
}

pub fn choose_subset(B: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    // Returns a subset of B based on G, and also removes that subset from B
    unimplemented!()
}


pub fn ComputeM<'a, P: FPolyRing>(L: &mut Vec<Poly<'a, P>>, G: &Ideal<'a, P>) {

    // This line is a bit of a hack
    let g_init = MonomialIdeal::from(&G);

    L.sort_by(|a, b| <P::Ord>::cmp(&a.lt().deg, &b.lt().deg));
    row_reduce(L);

    let mut i = 0;

    while i < L.len() {
        let mut aux = Vec::new();

        // TODO This isn't removing duplicates
        for term in L[i].terms.iter() {
            // If term is in G_init, then give back the scaled f
            if let Some(poly) = g_init.in_and_get(&term) {
                let p = poly.clone();
                aux.push(p.term_scale(&term.div(poly.lt()).unwrap()))
            }
        }
        L.append(&mut aux);
        i += 1;
    }
}



pub fn left_hand_s_poly<'a, P: FPolyRing>(a: &Poly<'a, P>, b: &Poly<'a, P>) -> Poly<'a, P> {
    // Returns the "Left hand side" of the S-Polynomial of a and b
    let lcm = a.lt().lcm(&b.lt());
    let scale = lcm.div(a.lt()).unwrap();
    a.term_scale(&scale)
}


pub fn pairs_ord(n: usize) -> Vec<(usize, usize)> {
    let mut acc = Vec::new();
    for i in 0..n {
        for j in i+1..n {
            acc.push((i, j))
        }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computeM_test() {
    }
}


