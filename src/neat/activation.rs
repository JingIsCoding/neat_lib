pub type ActivationFunc = fn(f64) -> f64;

pub fn binary_step(value: f64) -> f64 {
    if value < 0.0 {
        0.0
    } else {
        1.0
    }
}

pub fn sigmod(value: f64) -> f64 {
    let z = (5.0 * value).min(60.0).max(-60.0);
    1.0 / (1.0 + (-z).exp())
} 

pub fn sin(value: f64) -> f64 {
    (5.0 * value).min(60.0).max(-60.0).sin()
}

pub fn tanh(value: f64) -> f64 {
    (2.5 * value).min(60.0).max(-60.0).tanh()
}

pub fn gauss(value: f64) -> f64 {
    (value.max(-3.4).min(3.4).powi(2) * -5.0).exp()
}

pub fn relu(value: f64) -> f64 {
    value.max(0.0)
}

pub fn get_activation(name: &str) -> ActivationFunc {
    match name {
        "sigmod" => sigmod,
        "binary_step" => binary_step,
        "sin" => sin,
        "tanh" => tanh,
        "gauss" => gauss,
        "relu" => relu,
        _ => sin,
    }
}

#[cfg(test)]
mod tests {
}
