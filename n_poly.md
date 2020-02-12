# Items
## Algorithms 
### Karatsuba's
### [ ] DFT Multiplication 
### FFT Multiplication 
So the idea for this is a bottom up approach since its not very nice in the
recursive approach in Rust. For one, Rust does easily allow Python's slice 
with step operation e.g. [a:b:2]. This makes it more difficult, and upon 
thought this way would be less eff

In the prepocessing stage, it would be quicker if we could use the "with_capacity" 
command rather than having to fill it with zeros. However it should be evident that
with the current algorithm that won't work because we can only work on elements in
sequence in the "with_capacity" pushing a value each time. However, what we can do
is manage two vectors, one for the first half and the other for the second, and then
append them in the end. You can get around the two vectors thing using a recursive
function but I suspect that would be less efficient

I could scale all the roots of unity but that would be less accurate than just
scaling it at the end

Instead of adding `n / 2` I could do a bitwise set operation
#### Recursive implementation
A recursive implementation, not working though. 
```rust
fn fft_rec(sample: &[Complex<f64>], step: usize, offset: usize) -> Vec<Complex<f64>> {
    let n = sample.len();
    if step = n {
        // Base case. There is only one element left
        // This should be optimised to see what the best base case should be
        vec![sample[step]]
    } else {
        let mut result = Vec::new();

        for i in 0..n {
            fft_rec(sample, step: step << 2, offset) + 
            fft_rec(sample, step: step << 2, offset + 1)
        }
    }
}
```
### DONE Log2
A very quick test to see if a number is a power of two is
```rust
((num - 1) & num) != 0
```
Where & is the bitwise AND operation.
Note: This only works when num > 1