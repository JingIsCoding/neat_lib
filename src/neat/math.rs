pub fn clamp(value: f64, min: f64, max: f64) -> f64 {
    value.min(max).max(min)
}
