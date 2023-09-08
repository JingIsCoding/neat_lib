pub fn sigmod(value: f64) -> f64 {
    let z = (5.0 * value).min(60.0).max(-60.0);
    1.0 / (1.0 + (-z).exp())
} 

pub fn sin(value: f64) -> f64 {
    let z = (5.0 * value).min(60.0).max(-60.0);
    z.sin()
}

#[cfg(test)]
mod tests {
}
