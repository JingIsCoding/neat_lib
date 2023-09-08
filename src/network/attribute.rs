use serde::{Deserialize, Serialize};
use crate::network::math::clamp;
use rand_distr::{Normal, Distribution};
use std::ops::{Add, Sub, Mul, Div, AddAssign};
use rand::Rng;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FloatAttribute(f64);

impl FloatAttribute {
    pub fn new(value: f64) -> Self {
        FloatAttribute(value)
    }

    pub fn new_norm(mean: f64, std_dev: f64) -> Self {
        let normal = Normal::new(mean, std_dev).unwrap();
        FloatAttribute(normal.sample(&mut rand::thread_rng()))
    }

    pub fn abs(&self) -> f64 {
        self.0.abs()
    }

    pub fn mutate(&mut self, replace_rate: f64, mutate_rate: f64, mutate_power: f64, min: f64, max: f64, mean: f64, std_dev: f64 ) {
        let mut rng = rand::thread_rng();
        let rand = rng.gen_range(0.0..1.0);
        if rand < mutate_rate {
            let normal = Normal::new(0.0, mutate_power).unwrap();
            self.0 = clamp(self.0 + normal.sample(&mut rng), min, max);
            return
        }

        if rand < replace_rate + mutate_rate {
            self.0 = FloatAttribute::init_value(&mut rng, mean, std_dev, min, max)
        }
    }

    pub fn set_value(&mut self, value: f64) {
        self. 0 = value
    }

    pub fn value(&self) -> f64 {
        self.0
    }

    fn init_value(rng: &mut impl Rng, mean: f64, std_dev: f64, min: f64, max: f64) -> f64 {
        let normal = Normal::new(mean, std_dev).unwrap();
        clamp(normal.sample(rng), min, max)
    }
}

impl Add for FloatAttribute {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        return FloatAttribute(self.value() + rhs.value())
    }
}

impl Sub for FloatAttribute {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        return FloatAttribute(self.value() - rhs.value())
    }
}

impl Mul for FloatAttribute {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        return FloatAttribute(self.value() * rhs.value())
    }
}


impl Div for FloatAttribute {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        return FloatAttribute(self.value() / rhs.value())
    }
}

impl AddAssign for FloatAttribute {
    fn add_assign(&mut self, rhs: Self) {
        self.set_value(self.value() + rhs.value())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_float_attribute_mutate() {
        let mut float = FloatAttribute::new(0.0);
        float.mutate(0.0, 1.0, 1.0, -1.0, 1.0, 0.0, 0.0);
        assert!(float.value() < 1.0, "float should be smaller than 1.0")
    }

    #[test]
    fn test_float_attribute_ops() {
        let float1 = FloatAttribute::new(1.0);
        let float2 = FloatAttribute::new(1.0);
        assert!((float1 + float2).value() == 2.0, "+");
        assert!((float1 - float2).value() == 0.0, "-");
        assert!((float1 * float2).value() == 1.0, "*");
        assert!((float1 / float2).value() == 1.0, "/");
    }
}
