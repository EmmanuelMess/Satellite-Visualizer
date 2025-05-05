use hifitime::{Epoch, TimeScale};

pub(crate) fn get_gpst_seconds_of_week(time: Epoch) -> u64 {
    let (_, nanoseconds) = time.to_time_scale(TimeScale::GPST).to_time_of_week();
    nanoseconds / 1_000_000_000
}

pub(crate) fn get_gpst_week(time: Epoch) -> u64 {
    let (week, _) = time.to_time_scale(TimeScale::GPST).to_time_of_week();
    week.into()
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seconds_of_week() {
        let epoch = Epoch::from_gregorian_utc(2025, 4, 19,  9,  46,  50, 0);
        let seconds = get_gpst_seconds_of_week(epoch);
        assert_eq!(seconds, 553628);
    }
}