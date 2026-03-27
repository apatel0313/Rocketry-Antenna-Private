import math
from bisect import bisect_left

MOMENT_OF_INERTIA = 2.2527   # kg*m^2
MAX_TORQUE = 0.8             # N*m
GROUND_STATION_FROM_POLE = 1609.34  # horizontal distance in meters


def build_dict():
    data = {}
    with open("timeoveraltitude.csv", "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            t = float(parts[0])
            altitude = float(parts[1])
            data[t] = altitude
    return data


def main():
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())

    # Using CSV points directly for dynamics
    theta_rad = {t: math.atan(time_altitude_dict[t]/GROUND_STATION_FROM_POLE) for t in time_keys}
    angle_deg_rounded = {t: round(math.degrees(theta_rad[t]), 2) for t in time_keys}

    # Compute omega via first-difference at CSV points (one-sided at ends, central inside)
    omega = {}
    for i in range(len(time_keys)):
        if i == 0:
            t0, t1 = time_keys[i], time_keys[i+1]
            omega[time_keys[i]] = (theta_rad[t1] - theta_rad[t0]) / (t1 - t0)
        elif i == len(time_keys)-1:
            t0, t1 = time_keys[i-1], time_keys[i]
            omega[time_keys[i]] = (theta_rad[t1] - theta_rad[t0]) / (t1 - t0)
        else:
            t0, t1, t2 = time_keys[i-1], time_keys[i], time_keys[i+1]
            o0 = (theta_rad[t1] - theta_rad[t0]) / (t1 - t0)
            o1 = (theta_rad[t2] - theta_rad[t1]) / (t2 - t1)
            omega[time_keys[i]] = (o0 + o1) / 2

    # Compute alpha via first differences of omega at CSV points
    alpha = {}
    for i in range(1, len(time_keys)):
        t0, t1 = time_keys[i-1], time_keys[i]
        alpha[t1] = (omega[t1] - omega[t0]) / (t1 - t0)

    # Print only for points where accel is defined (skip first point alpha not available)
    for t in time_keys[1:]:
        alt = time_altitude_dict[t]
        angle = angle_deg_rounded[t]
        a = alpha[t]
        torque = MOMENT_OF_INERTIA * a
        status = "OK" if abs(torque) <= MAX_TORQUE else "EXCEEDS"
        print(f"T: {t:05.2f} | Altitude: {int(alt):05d} | Angle: {angle:05.2f} *** | Alpha: {a: .4f} | Torque: {torque: 4.3f} | {status}")


if __name__ == "__main__":
    main()
