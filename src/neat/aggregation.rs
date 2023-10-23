pub type AggregationFunc = fn(Vec<f64>) -> f64;

pub fn sum(values: Vec<f64>) -> f64 {
    values.iter().sum()
}

pub fn max(values: Vec<f64>) -> f64 { 
    *values.iter().max_by(|a, b| a.total_cmp(b)).unwrap()
}

pub fn min(values: Vec<f64>) -> f64 {
    *values.iter().min_by(|a, b| a.total_cmp(b)).unwrap()
}

pub fn mean(values: Vec<f64>) -> f64 {
    let sum: f64 = values.iter().sum();
    sum / values.len() as f64
}

pub fn get_aggregation(name: &str) -> AggregationFunc {
    match name {
        "sum" => sum,
        "max" => max,
        "min" => min,
        "mean" => mean,
        _ => sum,
    }
}

