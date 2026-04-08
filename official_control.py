import sys
import select
import time
import math
from bisect import bisect_left
import pigpio

pi = pigpio.pi()

try:
    import tty
    import termios
    has_tty = True
except ImportError:
    has_tty = False

# Constants from alphatest_csv.py
MOMENT_OF_INERTIA = 0.2   # kg*m^2
MAX_TORQUE = 0.4          # N*m
GROUND_STATION_FROM_POLE = 1609.34  # horizontal distance in meters

# Constants from servo_control.py
PIN = 18
max_pulse = 1250 # Set to pulsewidth corresponding to 90 deg - should be about 45 deg from min
rest_pulse = 1750 # Set to pulsewidth corresponding to 0 deg - should be about 45 deg max
deg_step_pulse = 5.56 # set to the pulsewidth corresponding to a 1 degree increase

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

def get_altitude_at_time_interpolated(time_keys, time_altitude_dict, t:float):
    if t <= time_keys[0]:
        return time_altitude_dict[time_keys[0]]
    if t >= time_keys[-1]:
        return time_altitude_dict[time_keys[-1]]
    idx = bisect_left(time_keys, t)
    t_before = time_keys[idx - 1]
    alt_before = time_altitude_dict[t_before]
    t_after = time_keys[idx]
    alt_after = time_altitude_dict[t_after]
    slope = (alt_after - alt_before)/(t_after - t_before)
    time_progress = t - t_before
    alt_interpolated = alt_before + (slope * time_progress)
    return alt_interpolated

def _solve_linear_system(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        piv = i
        for j in range(i + 1, n):
            if abs(M[j][i]) > abs(M[piv][i]):
                piv = j
        if abs(M[piv][i]) < 1e-12:
            raise ValueError("Singular system")
        M[i], M[piv] = M[piv], M[i]
        pivot_val = M[i][i]
        M[i] = [v / pivot_val for v in M[i]]
        for j in range(n):
            if j == i:
                continue
            factor = M[j][i]
            M[j] = [M[j][k] - factor * M[i][k] for k in range(n + 1)]
    return [M[i][-1] for i in range(n)]

def fit_polynomial(t_values, theta_values, center, degree):
    x = [t - center for t in t_values]
    n = len(x)
    m = degree + 1
    X = [[xi ** k for k in range(m)] for xi in x]
    AtA = [[0.0]*m for _ in range(m)]
    Atb = [0.0]*m
    for i in range(n):
        for p in range(m):
            for q in range(m):
                AtA[p][q] += X[i][p] * X[i][q]
            Atb[p] += X[i][p] * theta_values[i]
    for i in range(m):
        AtA[i][i] += 1e-6
    return _solve_linear_system(AtA, Atb)

def main():
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())
    
    # Polynomial settings
    window_s = 2
    degree = 4
    
    last_print_str = ""
    fail_count = 0
    last_altitude = -1
    last_servo_time = 0.0
    last_sent_pulse = None
    servo_interval = 0.02
    
    if has_tty:
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        tty.setcbreak(fd)
    else:
        old_settings = None

    try:
        last_monotonic = time.monotonic()
        sim_time = 0.0
        direction = 1.0

        while True:
            now = time.monotonic()
            dt = now - last_monotonic
            last_monotonic = now

            sim_time += direction * dt
            if sim_time < 0.0:
                sim_time = 0.0
            
            # Check for 's' keypress to enter reverse playback
            if has_tty:
                if select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], []):
                    c = sys.stdin.read(1)
                    if c and c.lower() == 's':
                        direction = -1.0
            
            # Use sim_time plus the first time key to get actual CSV time
            csv_time = time_keys[0] + sim_time
            lookup_time = min(csv_time, time_keys[-1])

            try:
                # 1. Altitude & Angle logic (from servo_control.py)
                altitude = get_altitude_at_time_interpolated(time_keys, time_altitude_dict, lookup_time)
                angle_deg = round(math.degrees(math.atan(altitude/GROUND_STATION_FROM_POLE)), 2)
                pulse = rest_pulse + deg_step_pulse * angle_deg
                
                # 2. Omega, Alpha, Torque logic (from alphatest_csv.py)
                # Compute polynomial dynamics based on real-time rolling window
                window_start = lookup_time - window_s / 2
                window_end = lookup_time + window_s / 2
                window_keys = [tk for tk in time_keys if window_start <= tk <= window_end]

                if len(window_keys) < degree + 1:
                    # Fallback to nearest points if window too small
                    idx = bisect_left(time_keys, lookup_time)
                    num_needed = degree + 1
                    start = max(0, idx - num_needed // 2)
                    end = min(len(time_keys), start + num_needed)
                    window_keys = time_keys[start:end]

                if len(window_keys) >= degree + 1:
                    theta_window = [math.atan(time_altitude_dict[tk] / GROUND_STATION_FROM_POLE) for tk in window_keys]
                    coeffs = fit_polynomial(window_keys, theta_window, lookup_time, degree)
                    omega = coeffs[1] if len(coeffs) > 1 else 0.0
                    alpha = 2 * coeffs[2] if len(coeffs) > 2 else 0.0
                else:
                    omega = 0.0
                    alpha = 0.0

                torque = MOMENT_OF_INERTIA * alpha
                status = "OK" if abs(torque) <= MAX_TORQUE else "EXCEEDS"

            except Exception as e:
                print(f"Error at t={lookup_time:.3f}: {e}\r")
                fail_count += 1
                if fail_count > 5:
                    break
                continue

            # Update servo and print if altitude moves or we are reversing continuously
            if abs(altitude - last_altitude) >= 0.1 or direction == -1.0:
                if now - last_servo_time >= servo_interval and pulse != last_sent_pulse:
                    # pi.set_servo_pulsewidth(PIN, pulse) # Uncomment in hardware
                    last_servo_time = now
                    last_sent_pulse = pulse
                    indicate_servo_update_str = "***"
                else:
                    indicate_servo_update_str = ""

                last_altitude = altitude
                # Show reversing status
                rev_str = "[REVERSING] " if direction == -1.0 else ""
                pi.set_servo_pulsewidth(PIN, pulse)
                print_str = f"{rev_str}T:{sim_time:05.2f} |Alt:{int(altitude):05d} |Ang:{angle_deg:05.2f} |Om:{math.degrees(omega):.2f} |Al:{alpha:.3f} |Trq:{torque:.4f} |{status} {indicate_servo_update_str}"
                
                if print_str != last_print_str:
                    print(print_str + "\r")
                    last_print_str = print_str
            
            time.sleep(0.003)
            
            if sim_time <= 0.0 and direction == -1.0:
                time.sleep(0.01) # Yield CPU if fully reversed back to zero

            if direction == 1.0 and sim_time >= (time_keys[-1] - time_keys[0]):
                time.sleep(0.01) # Yield CPU if finished standard run without quitting

    finally:
        if old_settings and has_tty:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        print("\nProcess Completed.")

if __name__ == "__main__":
    main()
