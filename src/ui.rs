use std::ffi::CString;
use raylib::drawing::RaylibDrawHandle;
use raylib::ffi;
use raylib::math::{Rectangle, Vector2, Vector3};
use raylib::prelude::RaylibDrawGui;
use crate::add;

const UI_GRID_SIZE: f32 = 8.0;

pub(crate) fn gui_slider_bar_helper(d: &mut RaylibDrawHandle, rect: impl Into<ffi::Rectangle>,
                                    value: &mut f32, min: f32, max: f32) {
    d.gui_slider_bar(rect, format!("{min:2.1}").as_str(), format!("{max:2.1}").as_str(), value, min, max);
}

pub(crate) fn gui_slider_bar_helper_sci(d: &mut RaylibDrawHandle, rect: impl Into<ffi::Rectangle>,
                                         value: &mut f32, min: f32, max: f32) {
    d.gui_slider_bar(rect, format!("{min:+e}").as_str(), format!("{max:+e}").as_str(), value, min, max);
}

pub(crate) fn gui_slider_helper(d: &mut RaylibDrawHandle, rectangle: Rectangle, title: &str,
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

pub(crate) fn gui_scroll_panel_helper(bounds: impl Into<ffi::Rectangle>,
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

pub(crate) fn grid_value(value: f32) -> f32 {
    value / UI_GRID_SIZE
}

pub(crate) fn grid_anchor(x: f32, y: f32) -> Vector2 {
    Vector2::new(x * UI_GRID_SIZE, y * UI_GRID_SIZE)
}

pub(crate) fn grid_anchor_from(reference: Vector2, x: f32, y: f32) -> Vector2 {
    reference + grid_anchor(x, y)
}

pub(crate) fn grid_rectangle(anchor: Vector2, x: f32, y: f32, w: f32, h: f32) -> Rectangle {
    Rectangle::new(
        anchor.x + x * UI_GRID_SIZE,
        anchor.y + y * UI_GRID_SIZE,
        w * UI_GRID_SIZE,
        h * UI_GRID_SIZE
    )
}
