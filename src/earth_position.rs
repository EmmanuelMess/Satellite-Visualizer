use std::f64::consts::PI;

#[derive(Debug, Copy, Clone, Default)]
pub(crate) struct EcefPosition {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

#[derive(Debug, Copy, Clone, Default)]
pub(crate) struct LlhPosition {
    pub(crate) latitude: f64,
    pub(crate) longitude: f64,
    pub(crate) height: f64
}

impl EcefPosition {
    pub(crate) fn new(x: f64, y: f64, z: f64) -> EcefPosition {
        EcefPosition {x, y, z}
    }
}

impl LlhPosition {
    pub(crate) fn new(latitude: f64, longitude: f64, height: f64) -> LlhPosition {
        LlhPosition {latitude, longitude, height}
    }
    
    pub(crate) fn on_surface(latitude: f64, longitude: f64) -> LlhPosition {
        LlhPosition {latitude, longitude, height: 0.0}
    }
}

impl Into<LlhPosition> for EcefPosition {
    fn into(self) -> LlhPosition {
        // See https://gssc.esa.int/navipedia/index.php?title=Ellipsoidal_and_Cartesian_Coordinates_Conversion

        // WGS84 constants
        let a = 6_378_137.0;
        let f = 1.0 / 298.257_223_563;
        let b = a * (1.0 - f);
        let e_2 = 2.0 * f - f * f;

        let l = f64::atan2(self.y, self.x);
        let p = f64::hypot(self.x, self.y);

        if p < 1e-20 {
            return if self.z >= 0.0 {
                LlhPosition::new(PI / 2.0, 0.0, self.z - b)
            } else {
                LlhPosition::new(-PI / 2.0, 0.0, -self.z - b)
            };
        }

        let mut theta = (self.z / ((1.0 - e_2) * p)).atan();
        let mut h = 0.0;

        for _ in 0..100 {
            let n = a / (1.0 - e_2 * theta.sin().powi(2)).sqrt();
            h = p / theta.cos() - n;
            let new_theta = (self.z / ((1.0 - e_2 * (n / (n + h))) * p)).atan();

            // See https://wiki.openstreetmap.org/wiki/Precision_of_coordinates
            if (theta - new_theta).abs() < 1e-9 {
                break;
            }
            theta = new_theta;
        }

        LlhPosition::new(theta, l, h)
    }
}

impl Into<EcefPosition> for LlhPosition {
    fn into(self) -> EcefPosition {
        // See https://gssc.esa.int/navipedia/index.php?title=Ellipsoidal_and_Cartesian_Coordinates_Conversion
        
        // WGS84 constants
        let a: f64 = 6_378_137.0;
        let f: f64 = 1.0 / 298.257_223_563;
        //let b: f64 = a * (1.0 - f);
        let e_2 = 2.0 * f - f.powi(2);

        let n = a / (1.0 - e_2 * self.latitude.sin().powi(2)).sqrt();

        let x = (n + self.height) * self.latitude.cos() * self.longitude.cos();
        let y = (n + self.height) * self.latitude.cos() * self.longitude.sin();
        let z = ((1.0 - e_2) * n + self.height) * self.latitude.sin();
        
        EcefPosition::new(x, y, z)
    }
}

#[cfg(test)]
mod tests {
    use assert_float_eq::assert_float_absolute_eq;
    use super::*;

    #[test]
    fn test_cordoba() {
        let positionEcef = EcefPosition::new(2345503.9452, -4910842.9601, -3316365.5474);
        let positionLlh: LlhPosition = positionEcef.into();

        // See https://gssc.esa.int/navipedia/index.php?title=Ellipsoidal_and_Cartesian_Coordinates_Conversion
        assert_float_absolute_eq!(positionLlh.latitude, f64::to_radians(-31.5284356), f64::to_radians(1e-6));
        assert_float_absolute_eq!(positionLlh.longitude, f64::to_radians(-64.4700483), f64::to_radians(1e-6));
        assert_float_absolute_eq!(positionLlh.height, 747.07, 0.1);
    }
}