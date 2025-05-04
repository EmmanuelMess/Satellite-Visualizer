mod satellite;
mod earth_position;
mod ui;
mod time;

use std::cmp::min;
use std::collections::HashMap;
use std::f32::consts::TAU;
use hifitime::{Duration, Epoch};
use hifitime::TimeScale::GPST;
use raylib::prelude::*;
use raylib::prelude::MouseButton::MOUSE_BUTTON_LEFT;
use crate::earth_position::{EcefPosition, LlhPosition};
use crate::satellite::Satellite;
use crate::ui::{grid_anchor, grid_anchor_from, grid_rectangle, grid_value, gui_scroll_panel_helper, gui_slider_helper};

// TODO delete this function
fn vec(values: (f64, f64, f64)) -> Vector3 {
    Vector3::new(values.0 as f32, values.1 as f32, values.2 as f32)
}

fn add(rectangle: Rectangle, vector: Vector2) -> Rectangle{
    Rectangle::new(vector.x + rectangle.x, vector.y + rectangle.y, rectangle.width, rectangle.height)
}


//TODO move ui to ui module
fn main() {
    let width = 1200;
    let height = 800;

    let scale = 1.0/1_000_000.0;
    
    let initial_camera_distance = 70.0;
    let initial_tilt_angle = -TAU/4.0;
    let initial_yaw_angle = -TAU/4.0;
    
    let orbit_point_number = 100;

    // Initialization
    let (mut rl, thread) = raylib::init()
        .size(width, height)
        .title("Satellite Visualizer")
        .build();

    // System sizes
    let earth_radius: f32 =  6_378_137.0;
    let satellite_radius = 100000.0;

    // Camera setup
    let initial_camera_position = Vector3::new(0.0, 0.0, initial_camera_distance)
        .transform_with(Matrix::rotate_xyz(Vector3::new(initial_tilt_angle, 0.0, initial_yaw_angle)));
    let mut camera = Camera3D::perspective(
        initial_camera_position,
        Vector3::new(0.0,  0.0,  0.0),
        Vector3::new(0.0,  0.0,  1.0),
        45.0,
    );

    // Camera state
    let mut camera_distance = initial_camera_distance;
    let mut tilt_angle = initial_tilt_angle;
    let mut orbit_angle = initial_yaw_angle;

    // Data out
    let top_text = grid_anchor(1.0, 1.0);
    let satellite_text = grid_anchor_from(top_text, 0.0, 2.0);
    let orbit_text = grid_anchor_from(top_text, 0.0, 8.0);
    let reference_text = grid_anchor_from(top_text, 0.0, 44.0);

    let top_text_layout = HashMap::from([
        ("fps_text",                         grid_rectangle(top_text, 0.0, 0.0, 64.0, 3.0)),
        ("satellite_text",                   grid_rectangle(satellite_text, 0.0, 0.0,64.0, 3.0)),
        ("satellite_position",               grid_rectangle(satellite_text, 1.0, 2.0, 64.0, 3.0)),
        ("satellite_velocity",               grid_rectangle(satellite_text, 1.0, 4.0, 64.0, 3.0)),

        ("orbit_text",                       grid_rectangle(orbit_text, 0.0, 0.0, 64.0, 3.0)),
        ("epoch_date_text",                  grid_rectangle(orbit_text, 1.0, 2.0, 64.0, 3.0)),
        ("epoch_seconds_text",               grid_rectangle(orbit_text, 1.0, 4.0, 64.0, 3.0)),

        ("correction_radius_sin_text",       grid_rectangle(orbit_text, 1.0, 6.0, 64.0, 3.0)),
        ("delta_mean_motion_text",           grid_rectangle(orbit_text, 1.0, 8.0, 64.0, 3.0)),
        ("mean_anomaly_text",                grid_rectangle(orbit_text, 1.0, 10.0, 64.0, 3.0)),
        ("correction_latitude_cos_text",     grid_rectangle(orbit_text, 1.0, 12.0, 64.0, 3.0)),
        ("eccentricity_text",                grid_rectangle(orbit_text, 1.0, 14.0, 64.0, 3.0)),
        ("correction_latitude_sin_text",     grid_rectangle(orbit_text, 1.0, 16.0, 64.0, 3.0)),
        ("sqrt_semi_major_axis_text",        grid_rectangle(orbit_text, 1.0, 18.0, 64.0, 3.0)),
        ("correction_inclination_cos_text",  grid_rectangle(orbit_text, 1.0, 20.0, 64.0, 3.0)),
        ("longitude_of_ascending_node_text", grid_rectangle(orbit_text, 1.0, 22.0, 64.0, 3.0)),
        ("correction_inclination_sin_text",  grid_rectangle(orbit_text, 1.0, 24.0, 64.0, 3.0)),
        ("inclination_text",                 grid_rectangle(orbit_text, 1.0, 26.0, 64.0, 3.0)),
        ("correction_radius_cos_text",       grid_rectangle(orbit_text, 1.0, 28.0, 64.0, 3.0)),
        ("argument_of_perigee_text",         grid_rectangle(orbit_text, 1.0, 30.0, 64.0, 3.0)),
        ("rate_of_right_ascension_text",     grid_rectangle(orbit_text, 1.0, 32.0, 64.0, 3.0)),
        ("rate_of_inclination_text",         grid_rectangle(orbit_text, 1.0, 34.0, 64.0, 3.0)),

        ("reference_text",                   grid_rectangle(reference_text, 0.0, 0.0, 64.0, 3.0)),
        ("reference_position_ecef",          grid_rectangle(reference_text, 1.0, 2.0, 64.0, 3.0)),
        ("reference_position_llh",           grid_rectangle(reference_text, 1.0, 4.0, 64.0, 3.0)),
    ]);

    // UI state
    let sidebar_box = grid_anchor(115.0, 0.0);
    let sidebar_content = grid_anchor_from(sidebar_box, 5.0, 2.0);
    let reference = grid_anchor_from(sidebar_content, 0.0, 1.0);
    let satellite = grid_anchor_from(sidebar_content, 0.0, 22.0);

    let mut sidebar_scroll = Vector2::default();
    let mut sidebar_visible = Rectangle::default();
    let mut vb1_edit = false;  let mut year_value: i32 = 2025;
    let mut vb2_edit = false;  let mut month_value: i32 = 4;
    let mut vb3_edit = false;  let mut day_value: i32 = 9;
    let mut vb4_edit = false;  let mut hour_value: i32 = 23;
    let mut vb5_edit = false;  let mut minutes_value: i32 = 59;
    let mut vb6_edit = false;  let mut seconds_value: i32 = 42;

    let mut position_along_orbit: f32 = 0.0;
    let mut correction_radius_sin: f32 = -9.875000000000E+00;
    let mut delta_mean_motion: f32 = 3.730869691412E-09;
    let mut mean_anomaly: f32 = -1.899624399301E+00;
    let mut correction_latitude_cos: f32 = -5.550682544708E-07;
    let mut eccentricity: f32 = 3.463134868070E-03;
    let mut correction_latitude_sin: f32 = 1.232884824276E-05;
    let mut sqrt_semi_major_axis: f32 = 5.153554180145E+03;
    let mut correction_inclination_cos: f32 = 9.499490261078E-08;
    let mut longitude_of_ascending_node: f32 = 2.686281956689E+00;
    let mut correction_inclination_sin: f32 = -3.539025783539E-08;
    let mut inclination: f32 = 9.890951254507E-01;
    let mut correction_radius_cos: f32 = 1.600000000000E+02;
    let mut argument_of_perigee: f32 = -6.486031939322E-01;
    let mut rate_of_right_ascension: f32 = -7.595316375414E-09;
    let mut rate_of_inclination: f32 = -1.135761594744E-10;

    let mut reference_point_latitude_deg: f32 = -31.5284356;
    let mut reference_point_longitude_deg: f32 = -64.4700483;

    let sidebar_layout = HashMap::from([
        ("sidebar", grid_rectangle(sidebar_box, 0.0, 0.0, grid_value(width.as_f32() - sidebar_box.x), grid_value(height.as_f32()))),
        ("sidebar_content", grid_rectangle(sidebar_content, 0.0, 0.0, grid_value(width.as_f32() - sidebar_box.x) - 2.0, 200.0)),

        ("reference_line", grid_rectangle(reference, 0.0, 0.0, 24.0, 3.0)),

        ("reference_latitude_slider", grid_rectangle(reference, 0.0, 3.0, 24.0, 9.0)),

        ("reference_longitude_slider", grid_rectangle(reference, 0.0, 12.0, 24.0, 9.0)),

        // Epoch
        ("satellite_line",                     grid_rectangle(satellite, 0.0, 0.0, 24.0, 3.0)),
        ("epoch_text",                         grid_rectangle(satellite, 0.0, 3.0, 24.0, 3.0)),
        ( "2",                                 grid_rectangle(satellite, 0.0, 6.0, 6.0, 3.0)),
        ( "3",                                 grid_rectangle(satellite, 9.0, 6.0, 6.0, 3.0)),
        ( "4",                                 grid_rectangle(satellite, 18.0, 6.0, 6.0, 3.0)),
        ( "5",                                 grid_rectangle(satellite, 0.0, 10.0, 6.0, 3.0)),
        ( "6",                                 grid_rectangle(satellite, 9.0, 10.0, 6.0, 3.0)),
        ( "7",                                 grid_rectangle(satellite, 18.0, 10.0, 6.0, 3.0)),
        // Satellite position
        ("position_slider",                    grid_rectangle(satellite, 0.0, 13.0, 24.0, 9.0)),

        ("correction_radius_sin_slider",       grid_rectangle(satellite, 0.0,  22.0, 24.0, 9.0)),
        ("delta_mean_motion_slider",           grid_rectangle(satellite, 0.0,  31.0, 24.0, 9.0)),
        ("mean_anomaly_slider",                grid_rectangle(satellite, 0.0,  40.0, 24.0, 9.0)),
        ("correction_latitude_cos_slider",     grid_rectangle(satellite, 0.0,  49.0, 24.0, 9.0)),
        ("eccentricity_slider",                grid_rectangle(satellite, 0.0,  58.0, 24.0, 9.0)),
        ("correction_latitude_sin_slider",     grid_rectangle(satellite, 0.0,  67.0, 24.0, 9.0)),
        ("sqrt_semi_major_axis_slider",        grid_rectangle(satellite, 0.0,  76.0, 24.0, 9.0)),
        ("correction_inclination_cos_slider",  grid_rectangle(satellite, 0.0,  85.0, 24.0, 9.0)),
        ("longitude_of_ascending_node_slider", grid_rectangle(satellite, 0.0,  94.0, 24.0, 9.0)),
        ("correction_inclination_sin_slider",  grid_rectangle(satellite, 0.0, 103.0, 24.0, 9.0)),
        ("inclination_slider",                 grid_rectangle(satellite, 0.0, 112.0, 24.0, 9.0)),
        ("correction_radius_cos_slider",       grid_rectangle(satellite, 0.0, 121.0, 24.0, 9.0)),
        ("argument_of_perigee_slider",         grid_rectangle(satellite, 0.0, 130.0, 24.0, 9.0)),
        ("rate_of_right_ascension_slider",     grid_rectangle(satellite, 0.0, 139.0, 24.0, 9.0)),
        ("rate_of_inclination_slider",         grid_rectangle(satellite, 0.0, 148.0, 24.0, 9.0)),
    ]);

    rl.set_target_fps(60);

    while !rl.window_should_close() {
        // Update
        // Toggle camera controls
        let moving_camera = rl.is_mouse_button_down(MOUSE_BUTTON_LEFT) || rl.get_mouse_wheel_move().abs() > 0.1;
        let in_gui = sidebar_layout["sidebar"].check_collision_point_rec(rl.get_mouse_position());
        let in_window = Rectangle::new(0.0, 0.0, width.as_f32(), height.as_f32()).check_collision_point_rec(rl.get_mouse_position());

        if moving_camera && in_window && !in_gui {
            let mouse_delta = rl.get_mouse_delta();
            tilt_angle += mouse_delta.y * 0.01; // TODO fix the rotation math to not have weird singularities
            orbit_angle += -mouse_delta.x * 0.01;

            let mouse_wheel_delta = rl.get_mouse_wheel_move();
            camera_distance -= mouse_wheel_delta;

            let camera_position = Vector3::new(0.0, 0.0, camera_distance);
            let matrix = Matrix::rotate_xyz(Vector3::new(tilt_angle, 0.0, orbit_angle));

            camera.position = camera_position.transform_with(matrix);
        }

        let epoch = Epoch::from_gregorian_utc(
            year_value, month_value as u8, day_value as u8, hour_value as u8, minutes_value as u8,
            seconds_value as u8, 0);
        // TODO allow for satellite to be empty, if position fails
        let satellite = Satellite::full_model(
            epoch,
            inclination.into(),
            longitude_of_ascending_node.into(),
            eccentricity.into(),
            argument_of_perigee.into(),
            mean_anomaly.into(),
            sqrt_semi_major_axis.into(),
            delta_mean_motion.into(),
            correction_latitude_cos.into(),
            correction_latitude_sin.into(),
            correction_radius_cos.into(),
            correction_radius_sin.into(),
            correction_inclination_cos.into(),
            correction_inclination_sin.into(),
            rate_of_inclination.into(),
            rate_of_right_ascension.into(),
        );
        let (satellite_position, satellite_velocity) =
            satellite.position_velocity(epoch + Duration::from_hours(position_along_orbit.into()));

        let satellite_position_vector =
            Vector3::new(satellite_position.0 as f32, satellite_position.1 as f32,
                         satellite_position.2 as f32);

        let satellite_orbit = satellite.get_orbit(orbit_point_number);

        let reference_point_llh =
            LlhPosition::on_surface(reference_point_latitude_deg.to_radians().into(),
                                    reference_point_longitude_deg.to_radians().into());
        let reference_point_ecef: EcefPosition = reference_point_llh.into();

        // Begin Drawing
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::RAYWHITE);

        // 3D mode
        {
            // Enter 3D mode (on drop, end_mode3d is called automatically)
            let mut d3d = d.begin_mode3D(&camera);

            unsafe {
                ffi::rlPushMatrix();
                ffi::rlScalef(scale, scale, scale);
            }

            // Draw Earth at origin
            d3d.draw_sphere_ex(Vector3::zero(), earth_radius, 256, 256,
                Color::new(0x00 ,0x00, 0xA0, 0x88));

            // Satellite position along its orbit
            d3d.draw_sphere_ex(satellite_position_vector, satellite_radius, 16, 16, Color::RED);

            // Reference point at Earth
            d3d.draw_sphere_ex(
                Vector3::new(reference_point_ecef.x as f32, reference_point_ecef.y as f32,
                reference_point_ecef.z as f32),
                100_000.0,
                16,
                16,
                Color::LIMEGREEN
            );


            // Draw Satellite's orbit circle
            for i in 0..((orbit_point_number-1) as usize) {
                let start = satellite_orbit[i];
                let end = satellite_orbit[i+1];
                d3d.draw_line_3D(
                    vec(start),
                    vec(end),
                    Color::LIGHTGRAY
                );
            }

            // Line to Earth
            d3d.draw_cylinder_ex(Vector3::zero(), satellite_position_vector, satellite_radius / 2.0,
                satellite_radius / 2.0, 16, Color::RED.alpha(0.75));

            // Grid and axis
            d3d.draw_cylinder_ex(Vector3::zero(), Vector3::new(1_000_000.0 * 10.0, 0.0, 0.0),
                                 50_000.0, 50_000.0, 32, Color::RED);
            d3d.draw_cylinder_ex(Vector3::zero(), Vector3::new(0.0, 1_000_000.0 * 10.0, 0.0),
                                 50_000.0, 50_000.0, 32, Color::GREEN);
            d3d.draw_cylinder_ex(Vector3::zero(), Vector3::new(0.0, 0.0, 1_000_000.0 * 10.0),
                                 50_000.0, 50_000.0, 32, Color::BLUE);
            
            unsafe {
                ffi::rlPopMatrix();
            }
        }

        // Gui
        {
            // Top Text
            {
                d.gui_label(top_text_layout["fps_text"], format!("FPS {}", d.get_fps()).as_str());
                d.gui_label(top_text_layout["satellite_text"], "Satellite:");
                d.gui_label(top_text_layout["satellite_position"],
                            format!("Position: x={:8.1}m y={:8.1}m z={:8.1}m", satellite_position.0,
                                    satellite_position.1, satellite_position.2).as_str());
                d.gui_label(top_text_layout["satellite_velocity"],
                            format!("Velocity: x={:8.1}m/s y={:8.1}m/s z={:8.1}m/s", satellite_velocity.0,
                                    satellite_velocity.1, satellite_velocity.2).as_str());
                d.gui_label(top_text_layout["orbit_text"], "Orbit:");
                d.gui_label(top_text_layout["epoch_date_text"],
                            format!("Epoch: {}", epoch.to_time_scale(GPST).to_string()).as_str());
                d.gui_label(top_text_layout["epoch_seconds_text"],
                            format!("Epoch seconds: {}s", epoch.to_gpst_seconds()).as_str());
                d.gui_label(top_text_layout["inclination_text"],
                            format!("Inclination: {:e}", inclination).as_str());
                d.gui_label(top_text_layout["longitude_of_ascending_node_text"],
                            format!("Longitude of the ascending node: {:e}", longitude_of_ascending_node).as_str());
                d.gui_label(top_text_layout["eccentricity_text"],
                            format!("Eccentricity: {:e}", eccentricity).as_str());
                d.gui_label(top_text_layout["argument_of_perigee_text"],
                            format!("Argument of perigee: {:e}", argument_of_perigee).as_str());
                d.gui_label(top_text_layout["mean_anomaly_text"],
                            format!("Mean anomaly: {:e}", mean_anomaly).as_str());
                d.gui_label(top_text_layout["sqrt_semi_major_axis_text"],
                            format!("Square root of the semi-major axis: {:e}", sqrt_semi_major_axis).as_str());
                d.gui_label(top_text_layout["delta_mean_motion_text"],
                            format!("Delta mean motion: {:e}", delta_mean_motion).as_str());
                d.gui_label(top_text_layout["correction_latitude_cos_text"],
                            format!("Correction of latitude (cosine): {:e}", correction_latitude_cos).as_str());
                d.gui_label(top_text_layout["correction_latitude_sin_text"],
                            format!("Correction latitude (sine): {:e}", correction_latitude_sin).as_str());
                d.gui_label(top_text_layout["correction_radius_cos_text"],
                            format!("Correction radius (cosine): {:e}", correction_radius_cos).as_str());
                d.gui_label(top_text_layout["correction_radius_sin_text"],
                            format!("Correction radius (sine): {:e}", correction_radius_sin).as_str());
                d.gui_label(top_text_layout["correction_inclination_cos_text"],
                            format!("Correction inclination (cosine): {:e}", correction_inclination_cos).as_str());
                d.gui_label(top_text_layout["correction_inclination_sin_text"],
                            format!("Correction inclination (sine): {:e}", correction_inclination_sin).as_str());
                d.gui_label(top_text_layout["rate_of_inclination_text"],
                            format!("Rate of inclination: {:e}", rate_of_inclination).as_str());
                d.gui_label(top_text_layout["rate_of_right_ascension_text"],
                            format!("Rate of right ascension: {:e}", rate_of_right_ascension).as_str());
                d.gui_label(top_text_layout["reference_text"], "Reference:");
                d.gui_label(top_text_layout["reference_position_ecef"],
                            format!("x={:8.1}m y={:8.1}m z={:8.1}m", reference_point_ecef.x,
                                    reference_point_ecef.y, reference_point_ecef.z).as_str());
                d.gui_label(top_text_layout["reference_position_llh"],
                            format!("lat={:1.3}rad lon={:1.3}rad alt={:8.1}m",
                                    reference_point_llh.latitude, reference_point_llh.longitude,
                                    reference_point_llh.height).as_str());
            }

            //Sidebar
            {
                d.draw_rectangle_rec(sidebar_layout["sidebar"], Color::LIGHTGRAY.alpha(0.3));

                gui_scroll_panel_helper(sidebar_layout["sidebar"], "Parameters",
                                        sidebar_layout["sidebar_content"],
                                        &mut sidebar_scroll, &mut sidebar_visible);

                d.draw_scissor_mode(sidebar_visible.x as i32, sidebar_visible.y as i32,
                                    sidebar_visible.width as i32, sidebar_visible.height as i32,
                                    |mut s| {
                    s.gui_line(add(sidebar_layout["reference_line"], sidebar_scroll), "Reference position");

                    gui_slider_helper(&mut s, add(sidebar_layout["reference_latitude_slider"], sidebar_scroll),
                                      "Latitude [deg]",
                                      format!("{:3.2}°", reference_point_latitude_deg).as_str(),
                                      &mut reference_point_latitude_deg, -90.0, 90.0, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["reference_longitude_slider"], sidebar_scroll),
                                      "Longitude [deg]",
                                      format!("{:3.2}°", reference_point_longitude_deg).as_str(),
                                      &mut reference_point_longitude_deg, -180.0, 180.0, false);

                    s.gui_line(add(sidebar_layout["satellite_line"], sidebar_scroll), "Satellite");

                    s.gui_label(add(sidebar_layout["epoch_text"], sidebar_scroll), "Epoch (UTC)");
                    if s.gui_value_box(add(sidebar_layout["2"], sidebar_scroll), "Y", &mut year_value, 2000, 2050, vb1_edit) { vb1_edit = !vb1_edit; }
                    if s.gui_value_box(add(sidebar_layout["3"], sidebar_scroll), "M", &mut month_value, 1, 12, vb2_edit) { vb2_edit = !vb2_edit; }
                    if s.gui_value_box(add(sidebar_layout["4"], sidebar_scroll), "D", &mut day_value, 1, 31, vb3_edit) { vb3_edit = !vb3_edit; }
                    if s.gui_value_box(add(sidebar_layout["5"], sidebar_scroll), "H", &mut hour_value, 0, 23, vb4_edit) { vb4_edit = !vb4_edit; }
                    if s.gui_value_box(add(sidebar_layout["6"], sidebar_scroll), "m", &mut minutes_value, 0, 59, vb5_edit) { vb5_edit = !vb5_edit; }
                    if s.gui_value_box(add(sidebar_layout["7"], sidebar_scroll), "s", &mut seconds_value, 0, 59, vb6_edit) { vb6_edit = !vb6_edit; }

                    // TODO use the correct min and max values
                    gui_slider_helper(&mut s, add(sidebar_layout["position_slider"], sidebar_scroll),
                                      "Position at time [hours]",
                                      format!("{} hs", position_along_orbit).as_str(),
                                      &mut position_along_orbit, -2.0, 2.0, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["correction_radius_sin_slider"], sidebar_scroll),
                                      "Correction of radius (sine) [m]",
                                      format!("{} m", correction_radius_sin).as_str(),
                                      &mut correction_radius_sin, -(1 << 10).as_f32(), (1 << 10).as_f32(), true);

                    gui_slider_helper(&mut s, add(sidebar_layout["delta_mean_motion_slider"], sidebar_scroll),
                                      "Delta mean motion [rad]",
                                      format!("{} rad", delta_mean_motion).as_str(),
                                      &mut delta_mean_motion, 0.0, TAU, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["mean_anomaly_slider"], sidebar_scroll),
                                      "Mean anomaly [rad]",
                                      format!("{} rad", mean_anomaly).as_str(),
                                      &mut mean_anomaly, -TAU, TAU, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["correction_latitude_cos_slider"], sidebar_scroll),
                                      "Correction of latitude (cosine) [rad]",
                                      format!("{} rad", correction_latitude_cos).as_str(),
                                      &mut correction_latitude_cos, - (1 << 10).as_f32(), (1 << 10).as_f32(), true);

                    gui_slider_helper(&mut s, add(sidebar_layout["eccentricity_slider"], sidebar_scroll),
                                      "Eccentricity [1]",
                                      format!("{}", eccentricity).as_str(),
                                      &mut eccentricity, 0.0, 0.03, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["correction_latitude_sin_slider"], sidebar_scroll),
                                      "Correction of latitude (sine) [rad]",
                                      format!("{} rad", correction_latitude_sin).as_str(),
                                      &mut correction_latitude_sin, - (1 << 10).as_f32(), (1 << 10).as_f32(), true);

                    gui_slider_helper(&mut s, add(sidebar_layout["sqrt_semi_major_axis_slider"], sidebar_scroll),
                                      "Square root of semi-major axis [m^½]",
                                      format!("{} m^½", sqrt_semi_major_axis).as_str(),
                                      &mut sqrt_semi_major_axis, 0.0, 10_000.0, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["correction_inclination_cos_slider"], sidebar_scroll),
                                      "Correction of inclination (cosine) [rad]",
                                      format!("{} rad", correction_inclination_cos).as_str(),
                                      &mut correction_inclination_cos, - (1 << 10).as_f32(), (1 << 10).as_f32(), true);

                    gui_slider_helper(&mut s, add(sidebar_layout["longitude_of_ascending_node_slider"], sidebar_scroll),
                                      "Longitude of ascending node [rad]",
                                      format!("{} rad", longitude_of_ascending_node).as_str(),
                                      &mut longitude_of_ascending_node, 0.0, TAU, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["correction_inclination_sin_slider"], sidebar_scroll),
                                      "Correction of inclination (sine) [rad]",
                                      format!("{} rad", correction_inclination_sin).as_str(),
                                      &mut correction_inclination_sin, -TAU, TAU, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["inclination_slider"], sidebar_scroll),
                                      "Inclination [rad]",
                                      format!("{} rad", inclination).as_str(),
                                      &mut inclination, 0.0, TAU, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["correction_radius_cos_slider"], sidebar_scroll),
                                      "Correction of radius (cosine) [rad]",
                                      format!("{} rad", correction_radius_cos).as_str(),
                                      &mut correction_radius_cos, - (1 << 10).as_f32(), (1 << 10).as_f32(), true);

                    gui_slider_helper(&mut s, add(sidebar_layout["argument_of_perigee_slider"], sidebar_scroll),
                                      "Argument of perigee [rad]",
                                      format!("{} rad", argument_of_perigee).as_str(),
                                      &mut argument_of_perigee, -TAU, TAU, false);

                    gui_slider_helper(&mut s, add(sidebar_layout["rate_of_right_ascension_slider"], sidebar_scroll),
                                      "Rate of right ascension [rad/s]",
                                      format!("{:+e} rad/s", rate_of_right_ascension).as_str(),
                                      &mut rate_of_right_ascension, -2e-6, 0.0, true);

                    gui_slider_helper(&mut s, add(sidebar_layout["rate_of_inclination_slider"], sidebar_scroll),
                                      "Rate of inclination [rad/s]",
                                      format!("{:+e} rad/s", rate_of_inclination).as_str(),
                                      &mut rate_of_inclination, -1.0, 1.0, true);
                });
            }
        }
    }
}
