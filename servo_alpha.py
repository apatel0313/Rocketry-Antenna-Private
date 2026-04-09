import math
import time
from bisect import bisect_left
import lgpio

# Constants
MOMENT_OF_INERTIA = 0.2   # kg*m^2
MAX_TORQUE = 0.4          # N*m
GROUND_STATION_FROM_POLE = 1609.34  # horizontal distance in meters

# Servo constants
PIN = 18
max_pulse = 1250  # pulsewidth corresponding to 90 deg
rest_pulse = 1750  # pulsewidth corresponding to 0 deg
deg_step_pulse = 5.56  # pulsewidth per degree
servo_interval = 0.02  # 20 ms between servo updates

# Polynomial fitting parameters
WINDOW_S = 5  # 2 second window for smoothing
POLYNOMIAL_DEGREE = 2
DT_OUT = 0.1  # output/analysis interval

# Diagnostic/debug mode - set to False to disable torque, omega calculations for reduced CPU usage
COMPUTE_DIAGNOSTICS = False


def build_dict():
    """Load altitude vs time data from CSV."""
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


def get_altitude_at_time_interpolated(time_keys, time_altitude_dict, t):
    """Get altitude at time t via linear interpolation."""
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
    alt_interpolated = alt_before + (slope * time_progress)
    
    return alt_interpolated


def _solve_linear_system(A, b):
    """Solve linear system Ax = b using Gaussian elimination."""
    n = len(A)
    # augmented matrix
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    
    for i in range(n):
        # pivot
        piv = i
        for j in range(i + 1, n):
            if abs(M[j][i]) > abs(M[piv][i]):
                piv = j
        if abs(M[piv][i]) < 1e-12:
            raise ValueError("Singular system")
        M[i], M[piv] = M[piv], M[i]
        
        # normalize row i
        pivot_val = M[i][i]
        M[i] = [v / pivot_val for v in M[i]]
        
        # eliminate other rows
        for j in range(n):
            if j == i:
                continue
            factor = M[j][i]
            M[j] = [M[j][k] - factor * M[i][k] for k in range(n + 1)]
    
    return [M[i][-1] for i in range(n)]


def fit_polynomial(t_values, theta_values, center, degree):
    """Fit polynomial of given degree to theta values centered at given time."""
    x = [t - center for t in t_values]
    n = len(x)
    m = degree + 1  # number of terms
    X = [[xi ** k for k in range(m)] for xi in x]
    
    # normal equations: AtA coeff matrix and Atb RHS
    AtA = [[0.0]*m for _ in range(m)]
    Atb = [0.0]*m
    for i in range(n):
        for p in range(m):
            for q in range(m):
                AtA[p][q] += X[i][p] * X[i][q]
            Atb[p] += X[i][p] * theta_values[i]
    
    # Add small regularization to prevent singularity
    for i in range(m):
        AtA[i][i] += 1e-6
    
    return _solve_linear_system(AtA, Atb)


def compute_angle_and_derivatives(t_eval, t_window, theta_window, coeffs, center, degree):
    """
    Evaluate polynomial at t_eval and compute angle, omega (angular velocity), 
    and alpha (angular acceleration).
    """
    tau = t_eval - center
    
    # Position: theta
    theta = sum(coeffs[k] * (tau ** k) for k in range(degree + 1))
    
    # Velocity: omega = d(theta)/dt
    omega = sum(k * coeffs[k] * (tau ** (k - 1)) for k in range(1, degree + 1))
    
    # Acceleration: alpha = d(omega)/dt = d²(theta)/dt²
    alpha = sum(k * (k - 1) * coeffs[k] * (tau ** (k - 2)) for k in range(2, degree + 1))
    
    return theta, omega, alpha


def main():
    """Real-time servo control with polynomial smoothing."""
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())
    
    launch_time = time.monotonic()
    fail_count = 0
    last_servo_time = 0.0
    last_sent_pulse = None
    last_angle = None
    last_print_str = ""
    
    # Track data for polynomial fitting
    data_buffer = []  # list of (t, altitude, theta) tuples
    last_fit_time = -1
    current_coeffs = None
    current_center = None
    max_torque_session = 0.0  # Track maximum torque observed during session
    
    # Initialize GPIO chip for servo control
    h = lgpio.gpiochip_open(0)
    lgpio.gpio_claim_output(h, PIN)
    
    try:
        while True:
            now = time.monotonic()
            time_from_launch = now - launch_time
            
            try:
                altitude = get_altitude_at_time_interpolated(time_keys, time_altitude_dict, time_from_launch)
            except Exception as e:
                print(f"Lookup error at t={time_from_launch:.3f}: {e}")
                fail_count += 1
                if fail_count > 5:
                    print("Too many lookup errors. Exiting.")
                    break
                time.sleep(0.003)
                continue
            
            # Convert altitude to angle
            theta = math.atan(altitude / GROUND_STATION_FROM_POLE)
            
            # Add to buffer
            data_buffer.append((time_from_launch, altitude, theta))
            
            # Fit polynomial approximately every DT_OUT seconds
            if time_from_launch - last_fit_time >= DT_OUT:
                # Remove points older than window
                window_start = time_from_launch - WINDOW_S / 2
                window_end = time_from_launch + WINDOW_S / 2
                data_buffer = [(t, alt, th) for t, alt, th in data_buffer 
                               if window_start <= t <= window_end]
                
                # Try to fit polynomial if we have enough points
                if len(data_buffer) >= POLYNOMIAL_DEGREE + 1:
                    try:
                        t_window = [t for t, _, _ in data_buffer]
                        theta_window = [th for _, _, th in data_buffer]
                        current_center = time_from_launch
                        current_coeffs = fit_polynomial(t_window, theta_window, current_center, POLYNOMIAL_DEGREE)
                        last_fit_time = time_from_launch
                    except ValueError as e:
                        print(f"Polynomial fitting error at t={time_from_launch:.3f}: {e}")
                        current_coeffs = None
            
            # Calculate angle and derivatives using fitted polynomial
            angle_deg = None
            if current_coeffs is not None:
                try:
                    # Always compute angle (needed for servo control)
                    if COMPUTE_DIAGNOSTICS:
                        # Compute angle and derivatives for diagnostics
                        theta_pred, omega, alpha = compute_angle_and_derivatives(
                            time_from_launch, 
                            [t for t, _, _ in data_buffer],
                            [th for _, _, th in data_buffer],
                            current_coeffs,
                            current_center,
                            POLYNOMIAL_DEGREE
                        )
                        angle_deg = round(math.degrees(theta_pred), 2)
                        torque = MOMENT_OF_INERTIA * alpha
                        torque_status = "OK" if abs(torque) <= MAX_TORQUE else "EXCEEDS"
                        # Track maximum torque
                        if abs(torque) > max_torque_session:
                            max_torque_session = abs(torque)
                    else:
                        # Only compute angle, skip expensive derivatives
                        tau = time_from_launch - current_center
                        theta_pred = sum(current_coeffs[k] * (tau ** k) for k in range(POLYNOMIAL_DEGREE + 1))
                        angle_deg = round(math.degrees(theta_pred), 2)
                    
                    # Command servo if angle changed
                    if last_angle is None or abs(angle_deg - last_angle) >= 0.3:
                        if now - last_servo_time >= servo_interval:
                            pulse = rest_pulse - deg_step_pulse * angle_deg
                            duty_cycle = (pulse / 20000) * 100
                            lgpio.tx_pwm(h, PIN, 50, duty_cycle)
                            last_servo_time = now
                            last_sent_pulse = pulse
                            
                            # Print only when servo is actually commanded
                            if COMPUTE_DIAGNOSTICS:
                                print_str = (f"T:{time_from_launch:05.2f} | Alt:{altitude:5.1f}m | "
                                           f"Angle:{angle_deg:05.2f}° | Omega:{math.degrees(omega):.2f}°/s | "
                                           f"Torque:{torque:+.4f}Nm [{torque_status}] ***")
                            else:
                                print_str = (f"T:{time_from_launch:05.2f} | Alt:{altitude:5.1f}m | "
                                           f"Angle:{angle_deg:05.2f}° ***")
                            if print_str != last_print_str:
                                print(print_str)
                                last_print_str = print_str
                        
                        last_angle = angle_deg
                except Exception as e:
                    print(f"Angle calculation error at t={time_from_launch:.3f}: {e}")
            
            # Check if we've run out of data
            if time_from_launch > time_keys[-1]:
                print("Reached end of trajectory data.")
                break
            
            time.sleep(0.003)  # Small sleep to avoid excessive CPU usage
    
    except KeyboardInterrupt:
        print("\nShutdown requested by user.")
    finally:
        lgpio.gpiochip_close(h)
        if COMPUTE_DIAGNOSTICS:
            print(f"Max torque during session: {max_torque_session:.4f} Nm")
        print("Servo control finished.")


if __name__ == "__main__":
    main()
