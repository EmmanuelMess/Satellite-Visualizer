mod satellite;
mod earth_position;
mod ui;
mod time;

mod string_format;

fn main() {
    let mut window = ui::Window::init();

    window.run();
}
