use std::collections::HashMap;
use std::f32::consts::TAU;
use std::ffi::CString;
use arboard::Clipboard;
use hifitime::{Duration, Epoch};
use hifitime::TimeScale::GPST;
use raylib::camera::Camera3D;
use raylib::color::Color;
use raylib::consts::MouseButton::MOUSE_BUTTON_LEFT;
use raylib::drawing::{RaylibDraw, RaylibDraw3D, RaylibDrawHandle, RaylibMode3DExt, RaylibScissorModeExt};
use raylib::{ffi, RaylibHandle, RaylibThread};
use raylib::math::{Matrix, Rectangle, Vector2, Vector3};
use raylib::misc::AsF32;
use raylib::prelude::RaylibDrawGui;
use crate::earth_position::{EcefPosition, LlhPosition};
use crate::satellite::Satellite;

const UI_GRID_SIZE: f32 = 8.0;

fn gui_slider_bar_helper(d: &mut RaylibDrawHandle, rect: impl Into<ffi::Rectangle>,
                                    value: &mut f32, min: f32, max: f32) {
    d.gui_slider_bar(rect, format!("{min:2.1}").as_str(), format!("{max:2.1}").as_str(), value, min, max);
}

fn gui_slider_bar_helper_sci(d: &mut RaylibDrawHandle, rect: impl Into<ffi::Rectangle>,
                                         value: &mut f32, min: f32, max: f32) {
    d.gui_slider_bar(rect, format!("{min:+e}").as_str(), format!("{max:+e}").as_str(), value, min, max);
}

fn gui_slider_helper(d: &mut RaylibDrawHandle, rectangle: Rectangle, title: &str,
                                text: &str, value: &mut f32, min: f32, max: f32, scientific: bool) {
    let anchor = Vector2::new(rectangle.x, rectangle.y);
    let title_rectangle = grid_rectangle(anchor, 2.0, 0.0, grid_value(rectangle.width), 3.0);
    let slider_rectangle = grid_rectangle(anchor, 2.0, 3.0, grid_value(rectangle.width) - 4.0, 3.0);
    let text_rectangle = grid_rectangle(anchor, 2.0, 6.0, grid_value(rectangle.width) - 4.0, 3.0);

    d.gui_label(title_rectangle, title);
    if !scientific {
        gui_slider_bar_helper(d, slider_rectangle, value, min, max);
    } else {
        gui_slider_bar_helper_sci(d, slider_rectangle, value, min, max);
    }
    d.gui_label(text_rectangle, text);
}

fn gui_scroll_panel_helper(bounds: impl Into<ffi::Rectangle>,
                                      text: &str, content: impl Into<ffi::Rectangle>,
                                      scroll: &mut (impl Into<ffi::Vector2> + From<ffi::Vector2> + Copy),
                                      view: &mut (impl Into<ffi::Rectangle> + From<ffi::Rectangle> + Copy)) -> bool {
    let c_text = CString::new(text).unwrap();
    // TODO fix the clone
    let mut scroll_out = scroll.clone().into();
    let mut view_out= view.clone().into();
    let result = unsafe {
        ffi::GuiScrollPanel(
            bounds.into(),
            c_text.as_ptr(),
            content.into(),
            &mut scroll_out,
            &mut view_out,
        )
    };
    *scroll = scroll_out.into();
    *view = view_out.into();
    result > 0
}

fn grid_value(value: f32) -> f32 {
    value / UI_GRID_SIZE
}

fn grid_anchor(x: f32, y: f32) -> Vector2 {
    Vector2::new(x * UI_GRID_SIZE, y * UI_GRID_SIZE)
}

fn grid_anchor_from(reference: Vector2, x: f32, y: f32) -> Vector2 {
    reference + grid_anchor(x, y)
}

fn grid_rectangle(anchor: Vector2, x: f32, y: f32, w: f32, h: f32) -> Rectangle {
    Rectangle::new(
        anchor.x + x * UI_GRID_SIZE,
        anchor.y + y * UI_GRID_SIZE,
        w * UI_GRID_SIZE,
        h * UI_GRID_SIZE
    )
}

// TODO delete this function
fn vec(values: (f64, f64, f64)) -> Vector3 {
    Vector3::new(values.0 as f32, values.1 as f32, values.2 as f32)
}

fn add(rectangle: Rectangle, vector: Vector2) -> Rectangle{
    Rectangle::new(vector.x + rectangle.x, vector.y + rectangle.y, rectangle.width, rectangle.height)
}


pub(crate) struct Window<'a> {
    rl: RaylibHandle,
    thread: RaylibThread,

    camera: Camera3D,
    camera_distance: f32,
    tilt_angle: f32,
    orbit_angle: f32,

    top_text_layout: HashMap<&'a str, Rectangle>,

    sidebar_scroll: Vector2,
    sidebar_visible: Rectangle,
    vb1_edit: bool, year_value: i32,
    vb2_edit: bool, month_value: i32,
    vb3_edit: bool, day_value: i32,
    vb4_edit: bool, hour_value: i32,
    vb5_edit: bool, minutes_value: i32,
    vb6_edit: bool, seconds_value: i32,

    simplified_mode: bool,
    simplified_mode_edited: bool,

    position_along_orbit: f32,
    correction_radius_sin: f32,
    delta_mean_motion: f32,
    mean_anomaly: f32,
    correction_latitude_cos: f32,
    eccentricity: f32,
    correction_latitude_sin: f32,
    sqrt_semi_major_axis: f32,
    correction_inclination_cos: f32,
    longitude_of_ascending_node: f32,
    correction_inclination_sin: f32,
    inclination: f32,
    correction_radius_cos: f32,
    argument_of_perigee: f32,
    rate_of_right_ascension: f32,
    rate_of_inclination: f32,

    reference_point_latitude_deg: f32,
    reference_point_longitude_deg: f32,

    sidebar_layout: HashMap<&'a str, Rectangle>,

    simple_parameters_layout: HashMap<&'a str, Rectangle>,

    complex_parameters_layout: HashMap<&'a str, Rectangle>,

    update_clipboard: bool,

    copy_rinex_rectangle: Rectangle,

    satellite_state: Option<SatelliteState>,

    reference_point_llh: LlhPosition,
}

struct SatelliteState {
    epoch: Epoch,

    satellite: Satellite,
    satellite_position: (f64, f64, f64),
    satellite_velocity: (f64, f64, f64),
    satellite_position_vector: Vector3,
    satellite_orbit: Vec<(f64, f64, f64)>,
}

impl<'a> Window<'a> {
    const WIDTH: i32 = 1200;
    const HEIGHT: i32 = 800;

    const SCALE: f32 = 1.0/1_000_000.0;

    const INITIAL_CAMERA_DISTANCE: f32 = 70.0;
    const INITIAL_TILT_ANGLE: f32 = -TAU/4.0;
    const INITIAL_YAW_ANGLE: f32 = -TAU/4.0;

    const ORBIT_POINT_NUMBER: i32 = 100;

    // System sizes
    const EARTH_RADIUS: f32 =  6_378_137.0;
    const SATELLITE_RADIUS: f32 = 100000.0;

    const FPS: u32 = 60;

    pub(crate) fn init() -> Window<'a> {
        // Initialization
        let (mut rl, thread) = raylib::init()
            .size(Self::WIDTH, Self::HEIGHT)
            .title("Satellite Visualizer")
            .build();

        // Camera setup
        let initial_camera_position = Vector3::new(0.0, 0.0, Self::INITIAL_CAMERA_DISTANCE)
            .transform_with(Matrix::rotate_xyz(Vector3::new(Self::INITIAL_TILT_ANGLE, 0.0, Self::INITIAL_YAW_ANGLE)));
        let camera = Camera3D::perspective(
            initial_camera_position,
            Vector3::new(0.0,  0.0,  0.0),
            Vector3::new(0.0,  0.0,  1.0),
            45.0,
        );

        // Camera state
        let camera_distance = Self::INITIAL_CAMERA_DISTANCE;
        let tilt_angle = Self::INITIAL_TILT_ANGLE;
        let orbit_angle = Self::INITIAL_YAW_ANGLE;

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

        // Sidebar state
        let sidebar_box = grid_anchor(115.0, 0.0);
        let sidebar_content = grid_anchor_from(sidebar_box, 5.0, 2.0);
        let reference = grid_anchor_from(sidebar_content, 0.0, 1.0);
        let satellite = grid_anchor_from(sidebar_content, 0.0, 22.0);

        let sidebar_scroll = Vector2::default();
        let sidebar_visible = Rectangle::default();
        let vb1_edit = false;  let year_value: i32 = 2025;
        let vb2_edit = false;  let month_value: i32 = 5;
        let vb3_edit = false;  let day_value: i32 = 1;
        let vb4_edit = false;  let hour_value: i32 = 9;
        let vb5_edit = false;  let minutes_value: i32 = 0;
        let vb6_edit = false;  let seconds_value: i32 = 0;

        let simplified_mode = true;

        let position_along_orbit: f32 = 0.0;
        let correction_radius_sin: f32 = 0.0;
        let delta_mean_motion: f32 = 0.0;
        let mean_anomaly: f32 = -1.899624399301E+00;
        let correction_latitude_cos: f32 = 0.0;
        let eccentricity: f32 = 0.0;
        let correction_latitude_sin: f32 = 0.0;
        let sqrt_semi_major_axis: f32 = 5.153554180145E+03;
        let correction_inclination_cos: f32 = 0.0;
        let longitude_of_ascending_node: f32 = 0.0;
        let correction_inclination_sin: f32 = 0.0;
        let inclination: f32 = 9.890951254507E-01;
        let correction_radius_cos: f32 = 0.0;
        let argument_of_perigee: f32 = -6.486031939322E-01;
        let rate_of_right_ascension: f32 = 0.0;
        let rate_of_inclination: f32 = 0.0;

        let reference_point_latitude_deg: f32 = -32.9575;
        let reference_point_longitude_deg: f32 = -60.639444444444;

        let sidebar_layout = HashMap::from([
            ("sidebar", grid_rectangle(sidebar_box, 0.0, 0.0, grid_value(Self::WIDTH.as_f32() - sidebar_box.x), grid_value(Self::HEIGHT.as_f32()))),
            ("sidebar_content", grid_rectangle(sidebar_content, 0.0, 0.0, grid_value(Self::WIDTH.as_f32() - sidebar_box.x) - 2.0, 200.0)),

            ("reference_line", grid_rectangle(reference, 0.0, 0.0, 24.0, 3.0)),

            ("reference_latitude_slider", grid_rectangle(reference, 0.0, 3.0, 24.0, 9.0)),

            ("reference_longitude_slider", grid_rectangle(reference, 0.0, 12.0, 24.0, 9.0)),

            // Epoch
            ("satellite_line",                     grid_rectangle(satellite, 0.0,   0.0, 24.0, 3.0)),
            ("simplified_mode",                    grid_rectangle(satellite, 4.0,   3.0, 3.0, 3.0)),
            ("epoch_text",                         grid_rectangle(satellite, 0.0,   6.0, 24.0, 3.0)),
            ( "2",                                 grid_rectangle(satellite, 0.0,   9.0, 6.0, 3.0)),
            ( "3",                                 grid_rectangle(satellite, 9.0,   9.0, 6.0, 3.0)),
            ( "4",                                 grid_rectangle(satellite, 18.0,  9.0, 6.0, 3.0)),
            ( "5",                                 grid_rectangle(satellite, 0.0,  13.0, 6.0, 3.0)),
            ( "6",                                 grid_rectangle(satellite, 9.0,  13.0, 6.0, 3.0)),
            ( "7",                                 grid_rectangle(satellite, 18.0, 13.0, 6.0, 3.0)),
            // Satellite position
            ("position_slider",                    grid_rectangle(satellite, 0.0,  16.0, 24.0, 9.0)),
        ]);

        let simple_parameters_layout = HashMap::from([
            ("mean_anomaly_slider",                grid_rectangle(satellite, 0.0,  25.0, 24.0, 9.0)),
            ("sqrt_semi_major_axis_slider",        grid_rectangle(satellite, 0.0,  34.0, 24.0, 9.0)),
            ("longitude_of_ascending_node_slider", grid_rectangle(satellite, 0.0,  43.0, 24.0, 9.0)),
            ("inclination_slider",                 grid_rectangle(satellite, 0.0,  52.0, 24.0, 9.0)),
            ("argument_of_perigee_slider",         grid_rectangle(satellite, 0.0,  61.0, 24.0, 9.0)),
        ]);

        let complex_parameters_layout = HashMap::from([
            ("correction_radius_sin_slider",       grid_rectangle(satellite, 0.0,  25.0, 24.0, 9.0)),
            ("delta_mean_motion_slider",           grid_rectangle(satellite, 0.0,  34.0, 24.0, 9.0)),
            ("mean_anomaly_slider",                grid_rectangle(satellite, 0.0,  43.0, 24.0, 9.0)),
            ("correction_latitude_cos_slider",     grid_rectangle(satellite, 0.0,  52.0, 24.0, 9.0)),
            ("eccentricity_slider",                grid_rectangle(satellite, 0.0,  61.0, 24.0, 9.0)),
            ("correction_latitude_sin_slider",     grid_rectangle(satellite, 0.0,  70.0, 24.0, 9.0)),
            ("sqrt_semi_major_axis_slider",        grid_rectangle(satellite, 0.0,  79.0, 24.0, 9.0)),
            ("correction_inclination_cos_slider",  grid_rectangle(satellite, 0.0,  88.0, 24.0, 9.0)),
            ("longitude_of_ascending_node_slider", grid_rectangle(satellite, 0.0,  97.0, 24.0, 9.0)),
            ("correction_inclination_sin_slider",  grid_rectangle(satellite, 0.0, 106.0, 24.0, 9.0)),
            ("inclination_slider",                 grid_rectangle(satellite, 0.0, 115.0, 24.0, 9.0)),
            ("correction_radius_cos_slider",       grid_rectangle(satellite, 0.0, 124.0, 24.0, 9.0)),
            ("argument_of_perigee_slider",         grid_rectangle(satellite, 0.0, 133.0, 24.0, 9.0)),
            ("rate_of_right_ascension_slider",     grid_rectangle(satellite, 0.0, 142.0, 24.0, 9.0)),
            ("rate_of_inclination_slider",         grid_rectangle(satellite, 0.0, 151.0, 24.0, 9.0)),
        ]);

        // Copy rinex UI state
        let update_clipboard= false;

        let copy_button_height = 3.0;
        let copy_rinex_anchor = grid_anchor(0.0, grid_value(Self::HEIGHT.as_f32()) - copy_button_height);
        let copy_rinex_rectangle = grid_rectangle(copy_rinex_anchor, 0.0, 0.0, 10.0, copy_button_height);

        let reference_point_llh =
            LlhPosition::on_surface(reference_point_latitude_deg.to_radians().into(),
                                    reference_point_longitude_deg.to_radians().into());

        rl.set_target_fps(Self::FPS);

        Window {
            rl,
            thread,
            camera,
            camera_distance,
            tilt_angle,
            orbit_angle,
            top_text_layout,
            sidebar_scroll,
            sidebar_visible,
            vb1_edit, year_value,
            vb2_edit, month_value,
            vb3_edit, day_value,
            vb4_edit, hour_value,
            vb5_edit, minutes_value,
            vb6_edit, seconds_value,
            simplified_mode,
            simplified_mode_edited: simplified_mode,
            position_along_orbit,
            correction_radius_sin,
            delta_mean_motion,
            mean_anomaly,
            correction_latitude_cos,
            eccentricity,
            correction_latitude_sin,
            sqrt_semi_major_axis,
            correction_inclination_cos,
            longitude_of_ascending_node,
            correction_inclination_sin,
            inclination,
            correction_radius_cos,
            argument_of_perigee,
            rate_of_right_ascension,
            rate_of_inclination,
            reference_point_latitude_deg,
            reference_point_longitude_deg,
            sidebar_layout,
            simple_parameters_layout,
            complex_parameters_layout,
            update_clipboard,
            copy_rinex_rectangle,
            satellite_state: Option::None,
            reference_point_llh
        }
    }

    pub(crate) fn run(&mut self) {
        while !self.rl.window_should_close() {
            self.update();
            self.draw();
        }
    }

    fn update(&mut self) {
        // Update
        self.simplified_mode = self.simplified_mode_edited;

        // Toggle camera controls
        let moving_camera = self.rl.is_mouse_button_down(MOUSE_BUTTON_LEFT) || self.rl.get_mouse_wheel_move().abs() > 0.1;
        let in_gui = self.sidebar_layout["sidebar"].check_collision_point_rec(self.rl.get_mouse_position());
        let in_window = Rectangle::new(0.0, 0.0, Self::WIDTH.as_f32(), Self::HEIGHT.as_f32()).check_collision_point_rec(self.rl.get_mouse_position());

        if moving_camera && in_window && !in_gui {
            let mouse_delta = self.rl.get_mouse_delta();
            self.tilt_angle += mouse_delta.y * 0.01; // TODO fix the rotation math to not have weird singularities
            self.orbit_angle += -mouse_delta.x * 0.01;

            let mouse_wheel_delta = self.rl.get_mouse_wheel_move();
            self.camera_distance -= mouse_wheel_delta;

            let camera_position = Vector3::new(0.0, 0.0, self.camera_distance);
            let matrix = Matrix::rotate_xyz(Vector3::new(self.tilt_angle, 0.0, self.orbit_angle));

            self.camera.position = camera_position.transform_with(matrix);
        }

        let epoch = Epoch::maybe_from_gregorian_utc(
            self.year_value, self.month_value as u8, self.day_value as u8, self.hour_value as u8,
            self.minutes_value as u8, self.seconds_value as u8, 0);

        self.satellite_state = epoch.map(|epoch| {
            let satellite = Satellite::full_model(
                epoch,
                self.inclination.into(),
                self.longitude_of_ascending_node.into(),
                self.eccentricity.into(),
                self.argument_of_perigee.into(),
                self.mean_anomaly.into(),
                self.sqrt_semi_major_axis.into(),
                self.delta_mean_motion.into(),
                self.correction_latitude_cos.into(),
                self.correction_latitude_sin.into(),
                self.correction_radius_cos.into(),
                self.correction_radius_sin.into(),
                self.correction_inclination_cos.into(),
                self.correction_inclination_sin.into(),
                self.rate_of_inclination.into(),
                self.rate_of_right_ascension.into(),
            );
            let (satellite_position, satellite_velocity) =
                satellite.position_velocity(epoch + Duration::from_hours(self.position_along_orbit.into()));

            let satellite_position_vector =
                Vector3::new(satellite_position.0 as f32, satellite_position.1 as f32,
                             satellite_position.2 as f32);

            let satellite_orbit = satellite.get_orbit(Self::ORBIT_POINT_NUMBER);

            SatelliteState {
                epoch,
                satellite,
                satellite_position,
                satellite_velocity,
                satellite_position_vector,
                satellite_orbit
            }
        }).ok();

        self.reference_point_llh =
            LlhPosition::on_surface(self.reference_point_latitude_deg.to_radians().into(),
                                    self.reference_point_longitude_deg.to_radians().into());
    }

    fn draw(&mut self) {

        let parameters_layout =
            if self.simplified_mode { &self.simple_parameters_layout } else { &self.complex_parameters_layout };

        let reference_point_ecef: EcefPosition = self.reference_point_llh.into();

        let mut d = self.rl.begin_drawing(&self.thread);
        d.clear_background(Color::RAYWHITE);

        // 3D mode
        {
            // Enter 3D mode (on drop, end_mode3d is called automatically)
            let mut d3d = d.begin_mode3D(&self.camera);

            unsafe {
                ffi::rlPushMatrix();
                ffi::rlScalef(Self::SCALE, Self::SCALE, Self::SCALE);
            }

            // Draw Earth at origin
            d3d.draw_sphere_ex(Vector3::zero(), Self::EARTH_RADIUS, 256, 256,
                               Color::new(0x00 ,0x00, 0xA0, 0x88));

            if let Some(ref satellite_state) = self.satellite_state {
                // Satellite position along its orbit
                d3d.draw_sphere_ex(satellite_state.satellite_position_vector, Self::SATELLITE_RADIUS, 16, 16, Color::RED);

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
                for i in 0..((Self::ORBIT_POINT_NUMBER - 1) as usize) {
                    let start = satellite_state.satellite_orbit[i];
                    let end = satellite_state.satellite_orbit[i + 1];
                    d3d.draw_line_3D(
                        vec(start),
                        vec(end),
                        Color::LIGHTGRAY
                    );
                }

                // Line to Earth
                d3d.draw_cylinder_ex(Vector3::zero(), satellite_state.satellite_position_vector,
                                     Self::SATELLITE_RADIUS / 2.0,
                                     Self::SATELLITE_RADIUS / 2.0, 16, Color::RED.alpha(0.75));
            }

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
                d.gui_label(self.top_text_layout["fps_text"], format!("FPS {}", d.get_fps()).as_str());
                d.gui_label(self.top_text_layout["satellite_text"], "Satellite:");
                match self.satellite_state {
                    Some(ref satellite_state) => {
                        d.gui_label(self.top_text_layout["satellite_position"],
                                    format!("Position: x={:8.1}m y={:8.1}m z={:8.1}m", satellite_state.satellite_position.0,
                                            satellite_state.satellite_position.1, satellite_state.satellite_position.2).as_str());
                        d.gui_label(self.top_text_layout["satellite_velocity"],
                                    format!("Velocity: x={:8.1}m/s y={:8.1}m/s z={:8.1}m/s", satellite_state.satellite_velocity.0,
                                            satellite_state.satellite_velocity.1, satellite_state.satellite_velocity.2).as_str());
                    }
                    None => {
                        d.gui_label(self.top_text_layout["satellite_position"], "Position: undefined");
                        d.gui_label(self.top_text_layout["satellite_velocity"], "Velocity: undefined");
                    }
                }
                d.gui_label(self.top_text_layout["orbit_text"], "Orbit:");
                match self.satellite_state {
                    Some(ref satellite_state) => {
                        d.gui_label(self.top_text_layout["epoch_date_text"],
                                    format!("Epoch: {}", satellite_state.epoch.to_time_scale(GPST).to_string()).as_str());
                        d.gui_label(self.top_text_layout["epoch_seconds_text"],
                                    format!("Epoch seconds: {}s", satellite_state.epoch.to_gpst_seconds()).as_str());
                    }
                    None => {
                        d.gui_label(self.top_text_layout["epoch_date_text"], "Epoch: undefined");
                        d.gui_label(self.top_text_layout["epoch_seconds_text"], "Epoch seconds: undefined");
                    }
                }
                d.gui_label(self.top_text_layout["inclination_text"],
                            format!("Inclination: {:e}", self.inclination).as_str());
                d.gui_label(self.top_text_layout["longitude_of_ascending_node_text"],
                            format!("Longitude of the ascending node: {:e}", self.longitude_of_ascending_node).as_str());
                d.gui_label(self.top_text_layout["eccentricity_text"],
                            format!("Eccentricity: {:e}", self.eccentricity).as_str());
                d.gui_label(self.top_text_layout["argument_of_perigee_text"],
                            format!("Argument of perigee: {:e}", self.argument_of_perigee).as_str());
                d.gui_label(self.top_text_layout["mean_anomaly_text"],
                            format!("Mean anomaly: {:e}", self.mean_anomaly).as_str());
                d.gui_label(self.top_text_layout["sqrt_semi_major_axis_text"],
                            format!("Square root of the semi-major axis: {:e}", self.sqrt_semi_major_axis).as_str());
                d.gui_label(self.top_text_layout["delta_mean_motion_text"],
                            format!("Delta mean motion: {:e}", self.delta_mean_motion).as_str());
                d.gui_label(self.top_text_layout["correction_latitude_cos_text"],
                            format!("Correction of latitude (cosine): {:e}", self.correction_latitude_cos).as_str());
                d.gui_label(self.top_text_layout["correction_latitude_sin_text"],
                            format!("Correction latitude (sine): {:e}", self.correction_latitude_sin).as_str());
                d.gui_label(self.top_text_layout["correction_radius_cos_text"],
                            format!("Correction radius (cosine): {:e}", self.correction_radius_cos).as_str());
                d.gui_label(self.top_text_layout["correction_radius_sin_text"],
                            format!("Correction radius (sine): {:e}", self.correction_radius_sin).as_str());
                d.gui_label(self.top_text_layout["correction_inclination_cos_text"],
                            format!("Correction inclination (cosine): {:e}", self.correction_inclination_cos).as_str());
                d.gui_label(self.top_text_layout["correction_inclination_sin_text"],
                            format!("Correction inclination (sine): {:e}", self.correction_inclination_sin).as_str());
                d.gui_label(self.top_text_layout["rate_of_inclination_text"],
                            format!("Rate of inclination: {:e}", self.rate_of_inclination).as_str());
                d.gui_label(self.top_text_layout["rate_of_right_ascension_text"],
                            format!("Rate of right ascension: {:e}", self.rate_of_right_ascension).as_str());
                d.gui_label(self.top_text_layout["reference_text"], "Reference:");
                d.gui_label(self.top_text_layout["reference_position_ecef"],
                            format!("x={:8.1}m y={:8.1}m z={:8.1}m", reference_point_ecef.x,
                                    reference_point_ecef.y, reference_point_ecef.z).as_str());
                d.gui_label(self.top_text_layout["reference_position_llh"],
                            format!("lat={:1.3}rad lon={:1.3}rad alt={:8.1}m",
                                    self.reference_point_llh.latitude, self.reference_point_llh.longitude,
                                    self.reference_point_llh.height).as_str());
            }

            //Sidebar
            {
                d.draw_rectangle_rec(self.sidebar_layout["sidebar"], Color::LIGHTGRAY.alpha(0.3));

                gui_scroll_panel_helper(self.sidebar_layout["sidebar"], "Parameters",
                                        self.sidebar_layout["sidebar_content"],
                                        &mut self.sidebar_scroll, &mut self.sidebar_visible);

                d.draw_scissor_mode(self.sidebar_visible.x as i32, self.sidebar_visible.y as i32,
                                    self.sidebar_visible.width as i32, self.sidebar_visible.height as i32,
                    |mut s| {
                        s.gui_line(add(self.sidebar_layout["reference_line"], self.sidebar_scroll), "Reference position");

                        gui_slider_helper(&mut s, add(self.sidebar_layout["reference_latitude_slider"], self.sidebar_scroll),
                                          "Latitude [deg]",
                                          format!("{:3.2}°", self.reference_point_latitude_deg).as_str(),
                                          &mut self.reference_point_latitude_deg, -90.0, 90.0, false);

                        gui_slider_helper(&mut s, add(self.sidebar_layout["reference_longitude_slider"], self.sidebar_scroll),
                                          "Longitude [deg]",
                                          format!("{:3.2}°", self.reference_point_longitude_deg).as_str(),
                                          &mut self.reference_point_longitude_deg, -180.0, 180.0, false);

                        s.gui_line(add(self.sidebar_layout["satellite_line"], self.sidebar_scroll), "Satellite");

                        s.gui_check_box(add(self.sidebar_layout["simplified_mode"], self.sidebar_scroll), "Simple orbit mode", &mut self.simplified_mode_edited);

                        s.gui_label(add(self.sidebar_layout["epoch_text"], self.sidebar_scroll), "Epoch (UTC)");
                        if s.gui_value_box(add(self.sidebar_layout["2"], self.sidebar_scroll), "Y", &mut self.year_value, 2000, 2050, self.vb1_edit) { self.vb1_edit = !self.vb1_edit; }
                        if s.gui_value_box(add(self.sidebar_layout["3"], self.sidebar_scroll), "M", &mut self.month_value, 1, 12, self.vb2_edit) { self.vb2_edit = !self.vb2_edit; }
                        if s.gui_value_box(add(self.sidebar_layout["4"], self.sidebar_scroll), "D", &mut self.day_value, 1, 31, self.vb3_edit) { self.vb3_edit = !self.vb3_edit; }
                        if s.gui_value_box(add(self.sidebar_layout["5"], self.sidebar_scroll), "H", &mut self.hour_value, 0, 23, self.vb4_edit) { self.vb4_edit = !self.vb4_edit; }
                        if s.gui_value_box(add(self.sidebar_layout["6"], self.sidebar_scroll), "m", &mut self.minutes_value, 0, 59, self.vb5_edit) { self.vb5_edit = !self.vb5_edit; }
                        if s.gui_value_box(add(self.sidebar_layout["7"], self.sidebar_scroll), "s", &mut self.seconds_value, 0, 59, self.vb6_edit) { self.vb6_edit = !self.vb6_edit; }

                        // TODO use the correct min and max values
                        gui_slider_helper(&mut s, add(self.sidebar_layout["position_slider"], self.sidebar_scroll),
                                          "Position at time [hours]",
                                          format!("{} hs", self.position_along_orbit).as_str(),
                                          &mut self.position_along_orbit, -2.0, 2.0, false);

                        if !self.simplified_mode {
                            gui_slider_helper(&mut s, add(parameters_layout["correction_radius_sin_slider"], self.sidebar_scroll),
                                              "Correction of radius (sine) [m]",
                                              format!("{} m", self.correction_radius_sin).as_str(),
                                              &mut self.correction_radius_sin, -(1 << 10).as_f32(), (1 << 10).as_f32(), true);


                            gui_slider_helper(&mut s, add(parameters_layout["delta_mean_motion_slider"], self.sidebar_scroll),
                                              "Delta mean motion [rad]",
                                              format!("{} rad", self.delta_mean_motion).as_str(),
                                              &mut self.delta_mean_motion, 0.0, TAU, false);
                        }

                        gui_slider_helper(&mut s, add(parameters_layout["mean_anomaly_slider"], self.sidebar_scroll),
                                          "Mean anomaly [rad]",
                                          format!("{} rad", self.mean_anomaly).as_str(),
                                          &mut self.mean_anomaly, -TAU, TAU, false);

                        if !self.simplified_mode {
                            gui_slider_helper(&mut s, add(parameters_layout["correction_latitude_cos_slider"], self.sidebar_scroll),
                                              "Correction of latitude (cosine) [rad]",
                                              format!("{} rad", self.correction_latitude_cos).as_str(),
                                              &mut self.correction_latitude_cos, -(1 << 10).as_f32(), (1 << 10).as_f32(), true);

                            gui_slider_helper(&mut s, add(parameters_layout["eccentricity_slider"], self.sidebar_scroll),
                                              "Eccentricity [1]",
                                              format!("{}", self.eccentricity).as_str(),
                                              &mut self.eccentricity, 0.0, 0.03, false);

                            gui_slider_helper(&mut s, add(parameters_layout["correction_latitude_sin_slider"], self.sidebar_scroll),
                                              "Correction of latitude (sine) [rad]",
                                              format!("{} rad", self.correction_latitude_sin).as_str(),
                                              &mut self.correction_latitude_sin, -(1 << 10).as_f32(), (1 << 10).as_f32(), true);
                        }

                        gui_slider_helper(&mut s, add(parameters_layout["sqrt_semi_major_axis_slider"], self.sidebar_scroll),
                                          "Square root of semi-major axis [m^½]",
                                          format!("{} m^½", self.sqrt_semi_major_axis).as_str(),
                                          &mut self.sqrt_semi_major_axis, 0.0, 10_000.0, false);

                        if !self.simplified_mode {
                            gui_slider_helper(&mut s, add(parameters_layout["correction_inclination_cos_slider"], self.sidebar_scroll),
                                              "Correction of inclination (cosine) [rad]",
                                              format!("{} rad", self.correction_inclination_cos).as_str(),
                                              &mut self.correction_inclination_cos, -(1 << 10).as_f32(), (1 << 10).as_f32(), true);
                        }

                        gui_slider_helper(&mut s, add(parameters_layout["longitude_of_ascending_node_slider"], self.sidebar_scroll),
                                          "Longitude of ascending node [rad]",
                                          format!("{} rad", self.longitude_of_ascending_node).as_str(),
                                          &mut self.longitude_of_ascending_node, 0.0, TAU, false);

                        if !self.simplified_mode {
                            gui_slider_helper(&mut s, add(parameters_layout["correction_inclination_sin_slider"], self.sidebar_scroll),
                                              "Correction of inclination (sine) [rad]",
                                              format!("{} rad", self.correction_inclination_sin).as_str(),
                                              &mut self.correction_inclination_sin, -TAU, TAU, false);
                        }

                        gui_slider_helper(&mut s, add(parameters_layout["inclination_slider"], self.sidebar_scroll),
                                          "Inclination [rad]",
                                          format!("{} rad", self.inclination).as_str(),
                                          &mut self.inclination, 0.0, TAU, false);

                        if !self.simplified_mode {
                            gui_slider_helper(&mut s, add(parameters_layout["correction_radius_cos_slider"], self.sidebar_scroll),
                                              "Correction of radius (cosine) [rad]",
                                              format!("{} rad", self.correction_radius_cos).as_str(),
                                              &mut self.correction_radius_cos, - (1 << 10).as_f32(), (1 << 10).as_f32(), true);
                        }

                        gui_slider_helper(&mut s, add(parameters_layout["argument_of_perigee_slider"], self.sidebar_scroll),
                                          "Argument of perigee [rad]",
                                          format!("{} rad", self.argument_of_perigee).as_str(),
                                          &mut self.argument_of_perigee, -TAU, TAU, false);

                        if !self.simplified_mode {
                            gui_slider_helper(&mut s, add(parameters_layout["rate_of_right_ascension_slider"], self.sidebar_scroll),
                                              "Rate of right ascension [rad/s]",
                                              format!("{:+e} rad/s", self.rate_of_right_ascension).as_str(),
                                              &mut self.rate_of_right_ascension, -2e-6, 0.0, true);

                            gui_slider_helper(&mut s, add(parameters_layout["rate_of_inclination_slider"], self.sidebar_scroll),
                                              "Rate of inclination [rad/s]",
                                              format!("{:+e} rad/s", self.rate_of_inclination).as_str(),
                                              &mut self.rate_of_inclination, -1.0, 1.0, true);
                        }
                    });
            }

            if d.gui_button(self.copy_rinex_rectangle, "Copy RINEX") { self.update_clipboard = true; };
        }

        {
            if self.update_clipboard {
                match self.satellite_state {
                    Some(ref satellite_state) => {
                        let clipboard_result =
                            Clipboard::new().and_then(|mut clipboard| {
                                let rinex = satellite_state.satellite.to_rinex();
                                println!("{}", rinex);
                                clipboard.set_text(rinex)
                            });

                        if let Err(error) = clipboard_result {
                            println!("Error copying to clipboard: {}", error);
                        }
                    }
                    None => {
                        println!("Error copying to clipboard: Missing satellite");
                    }
                }

                self.update_clipboard = false;
            }
        }
    }
}