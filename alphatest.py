import math
from bisect import bisect_left

time_altitude_dict = {}

MOMENT_OF_INERTIA = 2.2527   # kg*m^2
MAX_TORQUE = 0.8             # N*m
GROUND_STATION_FROM_POLE = 1609.34  # horizontal distance in meters

def build_dict():
    with open("timeoveraltitude.csv", 'r') as filename:
        for line in filename:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            t = float(parts[0])
            altitude = float(parts[1])
            time_altitude_dict[t] = altitude
    return time_altitude_dict


def get_altitude_at_time_interpolated(time_keys, time_altitude_dict, t: float):
    if t <= time_keys[0]:
        return time_altitude_dict[time_keys[0]]
    if t >= time_keys[-1]:
        return time_altitude_dict[time_keys[-1]]
    idx = bisect_left(time_keys, t)
    t_before = time_keys[idx - 1]
    alt_before = time_altitude_dict[t_before]
    t_after = time_keys[idx]
    alt_after = time_altitude_dict[t_after]
    slope = (alt_after - alt_before) / (t_after - t_before)
    time_progress = t - t_before
    return alt_before + (slope * time_progress)


def main():
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())
    # Resample between CSV entries using linear slope (20 ms default)
    resample_interval = 0.02  # seconds

    angle_deg_dict = {}
    altitude_dict = {}
    resampled_times = []

    for i in range(len(time_keys) - 1):
        t0 = time_keys[i]
        t1 = time_keys[i + 1]
        alt0 = time_altitude_dict[t0]
        alt1 = time_altitude_dict[t1]

        # include the start point
        if t0 not in altitude_dict:
            altitude_dict[t0] = alt0
            angle_deg_dict[t0] = round(math.degrees(math.atan(alt0 / GROUND_STATION_FROM_POLE)), 2)
            resampled_times.append(t0)

        # linear slope between the two CSV points
        dt_segment = t1 - t0
        if dt_segment <= 0:
            continue
        slope = (alt1 - alt0) / dt_segment

        # add intermediate points at resample_interval
        t = t0 + resample_interval
        while t < t1:
            alt = alt0 + slope * (t - t0)
            altitude_dict[t] = alt
            angle_deg_dict[t] = round(math.degrees(math.atan(alt / GROUND_STATION_FROM_POLE)), 2)
            resampled_times.append(t)
            t += resample_interval

    # include final CSV time
    last_t = time_keys[-1]
    altitude_dict[last_t] = time_altitude_dict[last_t]
    angle_deg_dict[last_t] = round(math.degrees(math.atan(altitude_dict[last_t] / GROUND_STATION_FROM_POLE)), 2)
    resampled_times.append(last_t)

    # sort resampled times
    resampled_times = sorted(set(resampled_times))

    # Convert rounded/displayed angles to radians for dynamics
    theta_rad = {t: math.radians(angle_deg_dict[t]) for t in resampled_times}

    # compute angular velocity at each resampled time
    omega_dict = {}
    for i in range(1, len(resampled_times)):
        t0 = resampled_times[i - 1]
        t1 = resampled_times[i]
        dt = t1 - t0
        if dt == 0:
            omega = 0.0
        else:
            omega = (theta_rad[t1] - theta_rad[t0]) / dt
        omega_dict[t1] = omega

    # compute angular acceleration and torque using resampled times
    for i in range(2, len(resampled_times)):
        t_prev = resampled_times[i - 1]
        t_curr = resampled_times[i]
        dt = t_curr - t_prev
        if dt == 0:
            continue

        omega_prev = omega_dict.get(t_prev, 0.0)
        omega_curr = omega_dict.get(t_curr, 0.0)

        alpha = (omega_curr - omega_prev) / dt
        torque = MOMENT_OF_INERTIA * alpha

        altitude = altitude_dict[t_curr]
        angle_deg = angle_deg_dict[t_curr]

        status = "OK" if abs(torque) <= MAX_TORQUE else "EXCEEDS"

        # single consolidated print (starred summary + dynamics)
        print(f"T: {t_curr:05.2f} | Altitude: {int(altitude):05d} | Angle: {angle_deg:05.2f} *** | Alpha: {alpha: .6f} | Torque: {torque: .6f} | {status}")


if __name__ == "__main__":
    main()