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

    # Step 1: compute angle at each time
    angle_dict = {}
    altitude_dict = {}
    for t in sampled_times:
        alt = get_altitude_at_time_interpolated(time_keys, time_altitude_dict, t)
        angle_rad = math.atan(alt / GROUND_STATION_FROM_POLE)
        altitude_dict[t] = alt
        angle_dict[t] = angle_rad

    # Angular velocity (omega) on sampled times
    omega_dict = {}
    for i in range(1, len(sampled_times)):
        t0 = sampled_times[i - 1]
        t1 = sampled_times[i]
        dt = t1 - t0
        theta0 = angle_dict[t0]
        theta1 = angle_dict[t1]
        omega = (theta1 - theta0) / dt
        omega_dict[t1] = omega

    # Compute angular acceleration and torque, print results
    print("Time(s) | Alt(m) | Angle(deg) | Alpha(rad/s^2) | Torque(Nm) | Status")
    print("-" * 75)

    for i in range(2, len(sampled_times)):
        t_prev = sampled_times[i - 1]
        t_curr = sampled_times[i]

        dt = t_curr - t_prev
        omega_prev = omega_dict[t_prev]
        omega_curr = omega_dict[t_curr]

        alpha = (omega_curr - omega_prev) / dt
        torque = MOMENT_OF_INERTIA * alpha

        altitude = altitude_dict[t_curr]
        angle_deg = math.degrees(angle_dict[t_curr])

        status = "OK"
        if abs(torque) > MAX_TORQUE:
            status = "EXCEEDS"

        print(f"{t_curr:6.3f} | "
              f"{altitude:7.3f} | "
              f"{angle_deg:10.4f} | "
              f"{alpha:14.6f} | "
              f"{torque:10.6f} | "
              f"{status}")


if __name__ == "__main__":
    main()