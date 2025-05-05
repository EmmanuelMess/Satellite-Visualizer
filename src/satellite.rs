use std::fmt::Debug;
use std::num::ParseFloatError;
use hifitime::{Duration, Epoch};
use hifitime::TimeScale::GPST;
use regex::Regex;
use crate::string_format::fmt_f64;
use crate::time::{get_gpst_seconds_of_week, get_gpst_week};

pub(crate) struct Satellite {
    epoch: Epoch,
    // TODO actually pass the semi major axis
    inclination: f64,
    longitude_of_ascending_node: f64,
    eccentricity: f64,
    argument_of_perigee: f64,
    mean_anomaly: f64,
    sqrt_semi_major_axis: f64,
    delta_mean_motion: f64,
    correction_latitude_cos: f64,
    correction_latitude_sin: f64,
    correction_radius_cos: f64,
    correction_radius_sin: f64,
    correction_inclination_cos: f64,
    correction_inclination_sin: f64,
    rate_of_inclination: f64,
    rate_of_right_ascension: f64,
}

impl Satellite {
    pub(crate) fn crude_model(epoch: Epoch, inclination: f64, longitude_of_ascending_node: f64,
                              eccentricity: f64, argument_of_perigee: f64, mean_anomaly: f64,
                              sqrt_semi_major_axis: f64) -> Satellite {
        Satellite {
            epoch: epoch.to_time_scale(GPST),
            inclination,
            longitude_of_ascending_node,
            eccentricity,
            argument_of_perigee,
            mean_anomaly,
            sqrt_semi_major_axis,
            delta_mean_motion: 0.0,
            correction_latitude_cos: 0.0,
            correction_latitude_sin: 0.0,
            correction_radius_cos: 0.0,
            correction_radius_sin: 0.0,
            correction_inclination_cos: 0.0,
            correction_inclination_sin: 0.0,
            rate_of_inclination: 0.0,
            rate_of_right_ascension: 0.0,
        }
    }

    pub(crate) fn full_model(epoch: Epoch, inclination: f64, longitude_of_ascending_node: f64,
                             eccentricity: f64, argument_of_perigee: f64, mean_anomaly: f64,
                             sqrt_semi_major_axis: f64, delta_mean_motion: f64,
                             correction_latitude_cos: f64, correction_latitude_sin: f64,
                             correction_radius_cos: f64, correction_radius_sin: f64,
                             correction_inclination_cos: f64, correction_inclination_sin: f64,
                             rate_of_inclination: f64, rate_of_right_ascension: f64) -> Satellite {
        Satellite {
            epoch: epoch.to_time_scale(GPST), inclination, longitude_of_ascending_node, 
            eccentricity, argument_of_perigee, mean_anomaly, sqrt_semi_major_axis, 
            delta_mean_motion, correction_latitude_cos, correction_latitude_sin, 
            correction_radius_cos, correction_radius_sin, correction_inclination_cos,
            correction_inclination_sin, rate_of_inclination, rate_of_right_ascension
        }
    }

    pub(crate) fn from_rinex(rinex: &str) -> Satellite {
        let regex: Regex = Regex::new(r#"(?m)^G(?<sv>\d{2}) (?<year>\d{4}) (?<month>\d{2}) (?<day>\d{2}) (?<hour>\d{2}) (?<minute>\d{2}) (?<second>\d{2})( ?(?<clock_bias>-?\d.\d{12}E-?\+?\d{2}))( ?(?<clock_drift>-?\d.\d{12}E-?\+?\d{2}))( ?(?<clock_drift_rate>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<iode>-?\d.\d{12}E-?\+?\d{2}))( ?(?<crs>-?\d.\d{12}E-?\+?\d{2}))( ?(?<delta_n>-?\d.\d{12}E-?\+?\d{2}))( ?(?<m0>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<cuc>-?\d.\d{12}E-?\+?\d{2}))( ?(?<e>-?\d.\d{12}E-?\+?\d{2}))( ?(?<cus>-?\d.\d{12}E-?\+?\d{2}))( ?(?<sqrt_a>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<toe>-?\d.\d{12}E-?\+?\d{2}))( ?(?<cic>-?\d.\d{12}E-?\+?\d{2}))( ?(?<omega0>-?\d.\d{12}E-?\+?\d{2}))( ?(?<cis>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<i0>-?\d.\d{12}E-?\+?\d{2}))( ?(?<crc>-?\d.\d{12}E-?\+?\d{2}))( ?(?<omega>-?\d.\d{12}E-?\+?\d{2}))( ?(?<omega_dot>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<idot>-?\d.\d{12}E-?\+?\d{2}))( ?(?<codes_l2>-?\d.\d{12}E-?\+?\d{2}))( ?(?<gps_week>-?\d.\d{12}E-?\+?\d{2}))( ?(?<l2_data>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<accuracy>-?\d.\d{12}E-?\+?\d{2}))( ?(?<health>-?\d.\d{12}E-?\+?\d{2}))( ?(?<tgd>-?\d.\d{12}E-?\+?\d{2}))( ?(?<iodc>-?\d.\d{12}E-?\+?\d{2}))([\n ]*)( ?(?<transmission_time>-?\d.\d{12}E-?\+?\d{2}))( ?(?<fit_interval>-?\d.\d{12}E-?\+?\d{2}))"#).unwrap();

        let captures = &regex.captures(rinex).unwrap();

        
        let seconds_of_week = captures.name("toe").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let week = captures.name("gps_week").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let epoch = Epoch::from_gpst_seconds(seconds_of_week + week * 7.0 * 24.0 * 60.0 * 60.0).to_time_scale(GPST);
        let correction_radius_sin: f64 = captures.name("crs").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let delta_mean_motion: f64 = captures.name("delta_n").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let mean_anomaly: f64 = captures.name("m0").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let correction_latitude_cos: f64 = captures.name("cuc").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let eccentricity: f64 = captures.name("e").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let correction_latitude_sin: f64 = captures.name("cus").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let sqrt_semi_major_axis: f64 = captures.name("sqrt_a").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let correction_inclination_cos: f64 = captures.name("cic").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let longitude_of_ascending_node: f64 = captures.name("omega0").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let correction_inclination_sin: f64 = captures.name("cis").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let inclination: f64 = captures.name("i0").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let correction_radius_cos: f64 = captures.name("crc").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let argument_of_perigee: f64 = captures.name("omega").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let rate_of_right_ascension: f64 = captures.name("omega_dot").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();
        let rate_of_inclination: f64 = captures.name("idot").map(|x| { x.as_str().parse::<f64>().unwrap() }).unwrap();

        Satellite {
            epoch,
            inclination,
            longitude_of_ascending_node,
            eccentricity,
            argument_of_perigee,
            mean_anomaly,
            sqrt_semi_major_axis,
            delta_mean_motion,
            correction_latitude_cos,
            correction_latitude_sin,
            correction_radius_cos,
            correction_radius_sin,
            correction_inclination_cos,
            correction_inclination_sin,
            rate_of_inclination,
            rate_of_right_ascension,
        }
    }
    
    pub(crate) fn to_rinex(self) -> String {
        let (year, month, day, hour, minute, second, _) = self.epoch.to_gregorian_tai();
        let seconds_of_week = get_gpst_seconds_of_week(self.epoch);
        let week = get_gpst_week(self.epoch);
        let satellite_accuracy = 2.0; // URA index see section 20.3.3.3.1 of IS-GPS-200M
        
        format!("\
G00 {year:04} {month:02} {day:02} {hour:02} {minute:02} {second:02}{}{}{}
    {}{}{}{}
    {}{}{}{}
    {}{}{}{}
    {}{}{}{}
    {}{}{}{}
    {}{}{}{}
    {}{}\
            ",
            fmt_f64(0.0,  19, 12, 2),
            fmt_f64(0.0, 19, 12, 2),
            fmt_f64(0.0, 19, 12, 2),
            fmt_f64(0.0, 19, 12, 2), // TODO
            fmt_f64(self.correction_radius_sin, 19, 12, 2),
            fmt_f64(self.delta_mean_motion, 19, 12, 2),
            fmt_f64(self.mean_anomaly, 19, 12, 2),
            fmt_f64(self.correction_latitude_cos, 19, 12, 2),
            fmt_f64(self.eccentricity, 19, 12, 2),
            fmt_f64(self.correction_latitude_sin, 19, 12, 2),
            fmt_f64(self.sqrt_semi_major_axis, 19, 12, 2),
            fmt_f64(seconds_of_week as f64, 19, 12, 2),
            fmt_f64(self.correction_inclination_cos, 19, 12, 2),
            fmt_f64(self.longitude_of_ascending_node, 19, 12, 2),
            fmt_f64(self.correction_inclination_sin, 19, 12, 2),
            fmt_f64(self.inclination, 19, 12, 2),
            fmt_f64(self.correction_radius_cos, 19, 12, 2),
            fmt_f64(self.argument_of_perigee, 19, 12, 2),
            fmt_f64(self.rate_of_right_ascension, 19, 12, 2),
            fmt_f64(self.rate_of_inclination, 19, 12, 2),
            fmt_f64(0.0, 19, 12, 2),
            fmt_f64(week  as f64, 19, 12, 2),
            fmt_f64(0.0, 19, 12, 2),
            fmt_f64(2.0, 19, 12, 2),
            fmt_f64(0.0, 19, 12, 2),
            fmt_f64(0.0, 19, 12, 2), // TODO
            fmt_f64(0.0, 19, 12, 2), // TODO
            fmt_f64(0.0, 19, 12, 2), // TODO
            fmt_f64(4.0, 19, 12, 2),
        )
    }

    ///
    /// Get the satellite position in ECEF in meters
    ///
    pub(crate) fn position_velocity(&self, time: Epoch) -> ((f64, f64, f64), (f64, f64, f64)) {
        // Gravitational effect of the Earth
        let mu = 3.986005e14;
        // Rate of rotation of the Earth
        let omega_dot_e = 7.2921151467e-5;

        let e = self.eccentricity;
        let sqrt_a =  self.sqrt_semi_major_axis;
        let omega_0 =  self.longitude_of_ascending_node;
        let delta_n = self.delta_mean_motion;
        let m_0 = self.mean_anomaly;
        let i_0 = self.inclination;
        let omega = self.argument_of_perigee;
        let omega_dot = self.rate_of_right_ascension;
        let idot = self.rate_of_inclination;
        let c_uc = self.correction_latitude_cos;
        let c_us = self.correction_latitude_sin;
        let c_rc = self.correction_radius_cos;
        let c_rs = self.correction_radius_sin;
        let c_ic = self.correction_inclination_cos;
        let c_is = self.correction_inclination_sin;

        // Semi-major axis
        let a = sqrt_a.powi(2);
        // Mean motion
        let n_0 = (mu / a.powi(3)).sqrt();
        let t_oe = get_gpst_seconds_of_week(self.epoch) as i64;
        let t = get_gpst_seconds_of_week(time) as i64;
        let tk: f64 = Satellite::gpst_seconds_wrap(t - t_oe) as f64;
        // Corrected mean motion
        let n = n_0 + delta_n;
        // Mean anomaly
        let m_k = m_0 + n * tk;
        // Kepler to solve for Eccentric Anomaly
        let mut ea: [f64; 4] = [m_k, 0.0, 0.0, 0.0];
        for i in 1..=3 {
            ea[i] = ea[i-1] + (m_k - ea[i-1] + e * ea[i-1].sin()) / (1.0 - e * ea[i-1].cos());
        }
        let e_k = ea[3];
        // True anomaly
        let vk = 2.0 * (((1.0+e) / (1.0-e)).sqrt() * (e_k /2.0).tan()).atan();
        // Argument of latitude
        let phi_k = vk + omega;
        let phi_k_2 = phi_k * 2.0;

        // Second harmonic perturbations
        // argument_of_latitude_correction
        let delta_u_k = c_us * phi_k_2.sin() + c_uc * phi_k_2.cos();
        // radius_correction
        let delta_r_k = c_rs * phi_k_2.sin() + c_rc * phi_k_2.cos();
        // inclination_correction
        let delta_i_k = c_is * phi_k_2.sin() + c_ic * phi_k_2.cos();

        // corrected_argument_latitude
        let u_k = phi_k + delta_u_k;
        // corrected_radius
        let r_k = a * (1.0 - e * e_k.cos()) + delta_r_k;
        // corrected_inclination
        let i_k = i_0 + delta_i_k + idot * tk;

        let plane_x = r_k * u_k.cos();
        let plane_y = r_k * u_k.sin();

        // Corrected longitude of the ascending node
        let omega_k = omega_0 + (omega_dot - omega_dot_e) * tk - omega_dot_e * t_oe as f64;

        let x = plane_x * omega_k.cos() - plane_y * i_k.cos() * omega_k.sin();
        let y = plane_x * omega_k.sin() + plane_y * i_k.cos() * omega_k.cos();
        let z = plane_y * i_k.sin();

        let position = (x, y, z);

        // eccentric_anomaly_rate
        let e_dot_k = n / (1.0 - e * e_k.cos());
        // true_anomaly_rate
        let v_dot_k = e_dot_k * ((1.0 - e.powi(2)).sqrt() / (1.0 - e * e_dot_k.cos()));

        // corrected_inclination_angle_rate
        let di_k_dt = idot + 2.0 * v_dot_k * (c_is * phi_k_2.cos() - c_ic * phi_k_2.sin());
        // corrected_argument_latitude_rate
        let u_dot_k = v_dot_k + 2.0 * v_dot_k * (c_us * phi_k_2.cos() - c_uc * phi_k_2.sin());
        // corrected_argument_latitude_rate
        let r_dot_k = e * a * e_dot_k * e_k.sin() 
            + 2.0 * v_dot_k * (c_rs * phi_k_2.cos() - c_rc * phi_k_2.sin());

        // longitude_ascending_node_rate
        let omega_dot_k = omega_dot - omega_dot_e;

        // Plane velocity
        let x_dot_prime = r_dot_k * u_k.cos()
            - r_k * u_dot_k * u_k.sin();
        let y_dot_prime = r_dot_k * u_k.sin()
            + r_k * u_dot_k * u_k.cos();

        let x_velocity = 
            - plane_x * omega_dot_k * omega_k.sin()
            + x_dot_prime * omega_k.cos()
            - y_dot_prime * omega_k.sin() * delta_i_k.cos()
            - plane_y * (omega_dot_k * omega_dot_k.cos() * delta_i_k.cos()
                        - di_k_dt * omega_dot_k.sin() * delta_i_k.sin());
        let y_velocity = 
            plane_x * omega_dot_k * omega_k.cos()
            + x_dot_prime * omega_k.sin()
            + y_dot_prime * omega_k.cos() * delta_i_k.cos()
            - plane_y * (omega_dot_k * omega_dot_k.sin() * delta_i_k.cos()
                        + di_k_dt * omega_dot_k.cos() * delta_i_k.sin());
        let z_velocity = y_dot_prime * delta_i_k.sin() + plane_y * di_k_dt * delta_i_k.cos();

        let velocity = (x_velocity, y_velocity, z_velocity);

        (position, velocity)
    }

    pub(crate) fn get_orbit(&self, point_number: u64) -> Vec<(f64, f64, f64)>  {
        let time_start = (-Duration::from_hours(2.0)).to_seconds().floor() as i64;
        let time_end = (Duration::from_hours(2.0)).to_seconds().floor() as i64;

        let time_difference = (time_end - time_start) / (point_number-1) as i64;
        
        (0..point_number)
            .map(|i| time_start + i as i64 * time_difference)
            .map(|t| self.position_velocity(self.epoch + Duration::from_seconds(t as f64)))
            .map(|(p, _)| p)
            .collect()
    }

    fn gpst_seconds_wrap(t: i64) -> i64 {
        let seconds_per_week = 7 * 24 * 60 * 60;
        
        if t > seconds_per_week / 2 {
            t - seconds_per_week
        } else if t < -seconds_per_week / 2 {
            t + seconds_per_week
        } else {
            t
        }
    }
}



#[cfg(test)]
mod tests {
    use assert_float_eq::assert_float_absolute_eq;
    use super::*;

    #[test]
    fn test_rinex() {
        let rinex_line = "\
G06 2025 04 09 23 59 44-2.959421835840E-04-2.148681232939E-11 0.000000000000E+00
     5.000000000000E+00-9.875000000000E+00 3.730869691412E-09-1.899624399301E+00
    -5.550682544708E-07 3.463134868070E-03 1.232884824276E-05 5.153554180145E+03
     3.455840000000E+05-9.499490261078E-08 2.686281956689E+00-3.539025783539E-08
     9.890951254507E-01 1.600000000000E+02-6.486031939322E-01-7.595316375414E-09
    -1.135761594744E-10 1.000000000000E+00 2.361000000000E+03 0.000000000000E+00
     2.000000000000E+00 0.000000000000E+00 3.725290298462E-09 5.000000000000E+00
     3.398820000000E+05 4.000000000000E+00\
        ";
        let satellite = Satellite::from_rinex(rinex_line);

        let epoch = Epoch::from_gpst_seconds(3.455840000000E+05 + 2.361000000000E+03 * 7.0 * 24.0 * 60.0 * 60.0);
        assert_eq!(satellite.epoch, epoch);
        assert_float_absolute_eq!(satellite.correction_radius_sin, -9.875000000000E+00, 0.1);
        assert_float_absolute_eq!(satellite.delta_mean_motion, 3.730869691412E-09, 0.1);
        assert_float_absolute_eq!(satellite.mean_anomaly, -1.899624399301E+00, 0.1);
        assert_float_absolute_eq!(satellite.correction_latitude_cos, -5.550682544708E-07, 0.1);
        assert_float_absolute_eq!(satellite.eccentricity, 3.463134868070E-03, 0.1);
        assert_float_absolute_eq!(satellite.correction_latitude_sin, 1.232884824276E-05, 0.1);
        assert_float_absolute_eq!(satellite.sqrt_semi_major_axis, 5.153554180145E+03, 0.1);
        assert_float_absolute_eq!(satellite.correction_inclination_cos, 9.499490261078E-08, 0.1);
        assert_float_absolute_eq!(satellite.longitude_of_ascending_node, 2.686281956689E+00, 0.1);
        assert_float_absolute_eq!(satellite.correction_inclination_sin, -3.539025783539E-08, 0.1);
        assert_float_absolute_eq!(satellite.inclination, 9.890951254507E-01, 0.1);
        assert_float_absolute_eq!(satellite.correction_radius_cos, 1.600000000000E+02, 0.1);
        assert_float_absolute_eq!(satellite.argument_of_perigee, -6.486031939322E-01, 0.1);
        assert_float_absolute_eq!(satellite.rate_of_right_ascension, -7.595316375414E-09, 0.1);
        assert_float_absolute_eq!(satellite.rate_of_inclination, -1.135761594744E-10, 0.1);
    }

    #[test]
    fn test_crude() {
        let inclination: f64 = 0.958_138_635_455;
        let longitude_of_ascending_node: f64 = 0.635_682_074_832;
        let eccentricity: f64 = 0.003_710_3;
        let argument_of_perigee: f64 = 4.057_441_961_27;
        let mean_anomaly: f64 = 2.367_757_296_49;
        let sqrt_semi_major_axis: f64 = 6_499.315_173_94;
        let satellite = Satellite::crude_model(
            Epoch::from_gregorian_utc(2025, 4, 19,  9,  46,  50, 0),
            inclination,
            longitude_of_ascending_node,
            eccentricity,
            argument_of_perigee,
            mean_anomaly,
            sqrt_semi_major_axis,
        );
        let calculated_position_time = Epoch::from_gregorian_utc(2025, 4, 19,  9, 0,  0, 0);
        let (position, velocity) = satellite.position_velocity(calculated_position_time);

        // Calculated by hand
        assert_float_absolute_eq!(position.0, -12_214_024.99, 0.1);
        assert_float_absolute_eq!(position.1, -40_481_865.45, 0.1);
        assert_float_absolute_eq!(position.2, -1_945_051.56, 0.1);
        // Velocity makes no sense for this example
    }

    #[test]
    fn test_precise() {
        // Data from COD0MGXFIN_20251000000_01D_05M_ORB.SP3 and CORD00ARG_R_20251000000_01D_GN.rnx
        // Courtesy of NASA's CDDIS database
        
        // 2025 04 09 23 59 44 GPS
        let satellite = Satellite::from_rinex("\
G06 2025 04 09 23 59 44-2.959421835840E-04-2.148681232939E-11 0.000000000000E+00
     5.000000000000E+00-9.875000000000E+00 3.730869691412E-09-1.899624399301E+00
    -5.550682544708E-07 3.463134868070E-03 1.232884824276E-05 5.153554180145E+03
     3.455840000000E+05-9.499490261078E-08 2.686281956689E+00-3.539025783539E-08
     9.890951254507E-01 1.600000000000E+02-6.486031939322E-01-7.595316375414E-09
    -1.135761594744E-10 1.000000000000E+00 2.361000000000E+03 0.000000000000E+00
     2.000000000000E+00 0.000000000000E+00 3.725290298462E-09 5.000000000000E+00
     3.398820000000E+05 4.000000000000E+00\
        ");
        // 2025  4 10  0  0  0.00000000 GPS
        let calculated_position_time = Epoch::from_gpst_seconds(1428278400.0);
        let (position, velocity) = satellite.position_velocity(calculated_position_time);
        // Here error is larger because the estimated satellite position is an amalgamation of data
        // PG06  23201.392605  -4035.022068 -12344.501318   -295.941970
        assert_float_absolute_eq!(position.0, 23_201_392.605, 5.0);
        assert_float_absolute_eq!(position.1, -4_035_022.068, 5.0);
        assert_float_absolute_eq!(position.2, -12_344_501.318 , 5.0);
    }
}