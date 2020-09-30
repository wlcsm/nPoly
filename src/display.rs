use crate::algebras::polyring::*;

pub fn show_vec_poly<'a, P: PolyRing>(poly: &Vec<Poly<'a, P>>) -> String {
    format!(
        "[\n {} \n ]",
        poly.iter()
            .map(|p| format!("{}", p))
            .collect::<Vec<String>>()
            .join("\n")
    )
}

pub fn show_term<'a, P: PolyRing>(term: &Term<P>, ring: &Option<&P>) -> String {
    let mut acc = term.coeff.to_string();

    // Add the extra variables if it is in a defined ring
    match ring {
        Some(r) => {
            for (i, symb) in r.symb().iter().enumerate() {
                match term.mon.get(i).unwrap() {
                    0 => {}
                    1 => acc.push(*symb),
                    d => acc.push_str(&format!("{}^{}", symb, d)),
                }
            }
        }
        // Can only be a univariate polynomial
        None => match term.mon.get(0).unwrap() {
            0 => {}
            1 => acc.push('?'),
            d => acc.push_str(&format!("?^{}", d)),
        },
    }
    acc
}
