import math
import time
from bisect import bisect_left
import pigpio

# Constants
MOMENT_OF_INERTIA = 0.2   # kg*m^2
MAX_TORQUE = 0.4          # N*m
GROUND_STATION_FROM_POLE = 1609.34  # horizontal distance in meters

# Servo constants
PIN = 18
min_pulse = 500     # pulsewidth corresponding to -90 deg
max_pulse = 2500    # pulsewidth corresponding to 90 deg
rest_pulse = 1500   # pulsewidth corresponding to 0 deg
deg_step_pulse = 11.11  # pulsewidth per degree (2000μs / 180°)
servo_interval = 0.10  # 100 ms between servo updates (was 20ms)
max_angle_rate = 0.3  # max degrees per servo update to limit torque

# Antenna angle limits (degrees)
MIN_ANGLE = 0
MAX_ANGLE = 90

# Button constants
BUTTON_UP_PIN = 17      # GPIO pin for antenna up button
BUTTON_DOWN_PIN = 22    # GPIO pin for antenna down button
BUTTON_START_PIN = 23   # GPIO pin for start trajectory tracking button
BUTTON_DEBOUNCE = 0.05  # Debounce time in seconds
MANUAL_DEG_PER_SEC = 10  # Manual hold speed (degrees/second)
MANUAL_ACCEL_DEG_PER_S2 = 60  # Acceleration while holding (deg/s^2)
MANUAL_DECEL_DEG_PER_S2 = 90  # Deceleration when released or reversing (deg/s^2)

# Polynomial fitting parameters
WINDOW_S = 5  # 5-second window for smoothing
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


def update_servo(pi, angle_deg, last_sent_pulse=None):
    """Send servo command for given angle in degrees."""
    pulse = rest_pulse - deg_step_pulse * angle_deg
    # Clamp to valid range
    pulse = max(min_pulse, min(max_pulse, pulse))
    pi.set_servo_pulsewidth(PIN, int(pulse))
    return pulse


def setup_buttons(pi):
    """Setup GPIO pins for buttons as inputs."""
    pi.set_mode(BUTTON_UP_PIN, pigpio.INPUT)
    pi.set_mode(BUTTON_DOWN_PIN, pigpio.INPUT)
    pi.set_mode(BUTTON_START_PIN, pigpio.INPUT)
    # Enable internal pull-ups for active-low buttons
    pi.set_pull_up_down(BUTTON_UP_PIN, pigpio.PUD_UP)
    pi.set_pull_up_down(BUTTON_DOWN_PIN, pigpio.PUD_UP)
    pi.set_pull_up_down(BUTTON_START_PIN, pigpio.PUD_UP)


class ButtonState:
    """Track button states and debounce."""
    def __init__(self):
        self.last_press_time = {
            'up': -float('inf'),
            'down': -float('inf'),
            'start': -float('inf')
        }
        self.button_held = {
            'up': False,
            'down': False,
            'start': False
        }
        self.mode = 'manual'  # 'manual' or 'auto'
        self.manual_angle = 0  # Current manual angle in degrees
        self.tracking_started = False
        self.manual_vel = 0.0  # deg/s, ramps with accel/decel
    
    def check_button_held(self, pi, button_pin, button_name):
        """Check if button is currently held down (continuous movement)."""
        try:
            state = pi.read(button_pin)
            return state == 0  # Button pressed (active low)
        except Exception as e:
            print(f"Error reading button {button_name}: {e}")
        return False
    
    def check_button_press(self, pi, button_pin, button_name, current_time):
        """Check if button was just pressed (with debounce) - for toggle buttons."""
        if current_time - self.last_press_time[button_name] < BUTTON_DEBOUNCE:
            return False
        
        try:
            state = pi.read(button_pin)
            if state == 0:  # Button pressed (active low)
                self.last_press_time[button_name] = current_time
                return True
        except Exception as e:
            print(f"Error reading button {button_name}: {e}")
        
        return False
    
    def handle_buttons(self, pi, current_time, dt):
        """Handle button presses and continuous holds. Time-based movement using dt (s)."""
        movement = 0.0  # degrees to move this tick
        should_print = False
        
        if self.mode == 'manual':
            # UP button - continuous movement while held
            up_held = self.check_button_held(pi, BUTTON_UP_PIN, 'up')
            if up_held and not self.button_held['up']:
                print(f"[MANUAL] Antenna UP (holding...)")
                self.button_held['up'] = True
            elif not up_held and self.button_held['up']:
                print(f"[MANUAL] Antenna UP released → {self.manual_angle}°")
                self.button_held['up'] = False
            
            # DOWN button - continuous movement while held
            # DOWN button - continuous movement while held
            down_held = self.check_button_held(pi, BUTTON_DOWN_PIN, 'down')
            if down_held and not self.button_held['down']:
                print(f"[MANUAL] Antenna DOWN (holding...)")
                self.button_held['down'] = True
            elif not down_held and self.button_held['down']:
                print(f"[MANUAL] Antenna DOWN released → {self.manual_angle}°")
                self.button_held['down'] = False

            # Desired velocity from held buttons
            if up_held and not down_held:
                desired_vel = MANUAL_DEG_PER_SEC
            elif down_held and not up_held:
                desired_vel = -MANUAL_DEG_PER_SEC
            else:
                desired_vel = 0.0

            # Choose accel or decel limit
            same_dir = (self.manual_vel == 0.0 and desired_vel != 0.0) or (self.manual_vel * desired_vel > 0)
            speeding_up = same_dir and abs(desired_vel) > abs(self.manual_vel)
            accel_limit = MANUAL_ACCEL_DEG_PER_S2 if speeding_up else MANUAL_DECEL_DEG_PER_S2

            # Ramp velocity toward desired with limits
            dv = desired_vel - self.manual_vel
            max_step = accel_limit * dt
            if dv > max_step:
                dv = max_step
            elif dv < -max_step:
                dv = -max_step
            self.manual_vel += dv

            # Movement this tick
            movement += self.manual_vel * dt
        
        # START button - toggle between manual and auto (debounced single press)
        if self.check_button_press(pi, BUTTON_START_PIN, 'start', current_time):
            self.tracking_started = not self.tracking_started
            if self.tracking_started:
                self.mode = 'auto'
                print("[AUTO MODE] Trajectory tracking started...")
            else:
                self.mode = 'manual'
                self.manual_angle = 0
                print("[MANUAL MODE] Trajectory tracking stopped. Manual control active.")
            return movement, True
        
        # Apply continuous movement
        if movement != 0:
            old_angle = self.manual_angle
            self.manual_angle = self.manual_angle + movement
            # Clamp to allowed [MIN_ANGLE, MAX_ANGLE]
            if self.manual_angle < MIN_ANGLE:
                self.manual_angle = MIN_ANGLE
                if self.manual_vel < 0:
                    self.manual_vel = 0.0
            elif self.manual_angle > MAX_ANGLE:
                self.manual_angle = MAX_ANGLE
                if self.manual_vel > 0:
                    self.manual_vel = 0.0
            return movement, True
        
        return movement, False


def main():
    """Real-time servo control with button input and manual/auto modes."""
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())
    
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
    
    # Initialize pigpio for servo and button control
    pi = pigpio.pi()
    if not pi.connected:
        raise RuntimeError("pigpio daemon not running or Pi not reachable. Start it with: sudo pigpiod")
    pi.set_mode(PIN, pigpio.OUTPUT)
    
    # Setup button inputs
    setup_buttons(pi)
    button_state = ButtonState()
    
    # Set servo to rest position
    pi.set_servo_pulsewidth(PIN, int(rest_pulse))
    
    print("Setup complete. Press START button to begin trajectory tracking, or use UP/DOWN for manual control...")
    print(f"Using GPIO pins: UP={BUTTON_UP_PIN}, DOWN={BUTTON_DOWN_PIN}, START={BUTTON_START_PIN}")
    
    launch_time = None
    last_loop_time = time.monotonic()
    
    try:
        while True:
            now = time.monotonic()
            dt = now - last_loop_time
            if dt < 0:
                dt = 0.0
            elif dt > 0.2:
                dt = 0.2  # cap to avoid large jumps if paused
            last_loop_time = now
            
            # Handle button presses
            button_state.handle_buttons(pi, now, dt)
            
            # Start/reset timer when tracking begins
            if button_state.tracking_started and launch_time is None:
                launch_time = now
                print(f"Launch time recorded. Tracking begins now.")
            
            # MANUAL MODE: Just accept up/down button input
            if not button_state.tracking_started:
                if now - last_servo_time >= servo_interval:
                    pulse = update_servo(pi, button_state.manual_angle)
                    last_servo_time = now
                    last_angle = button_state.manual_angle
                time.sleep(0.01)
                continue
            
            # AUTOMATIC MODE: Track trajectory
            if launch_time is None:
                time.sleep(0.01)
                continue
            
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
                    
                    # Command servo if angle changed significantly
                    if last_angle is None or abs(angle_deg - last_angle) >= 0.15:
                        if now - last_servo_time >= servo_interval:
                            # Apply slew rate limiting to avoid high torque
                            if last_angle is None:
                                commanded_angle = angle_deg
                            else:
                                angle_change = angle_deg - last_angle
                                angle_change = max(-max_angle_rate, min(max_angle_rate, angle_change))
                                commanded_angle = last_angle + angle_change
                            # Clamp commanded angle to allowed [MIN_ANGLE, MAX_ANGLE]
                            if commanded_angle < MIN_ANGLE:
                                commanded_angle = MIN_ANGLE
                            elif commanded_angle > MAX_ANGLE:
                                commanded_angle = MAX_ANGLE
                            
                            pulse = update_servo(pi, commanded_angle)
                            last_servo_time = now
                            last_sent_pulse = pulse
                            
                            # Print only when servo is actually commanded
                            if COMPUTE_DIAGNOSTICS:
                                print_str = (f"T:{time_from_launch:05.2f} | Alt:{altitude:5.1f}m | "
                                           f"Angle:{angle_deg:05.2f}° | Omega:{math.degrees(omega):.2f}°/s | "
                                           f"Torque:{torque:+.4f}Nm [{torque_status}] ***")
                            else:
                                print_str = (f"T:{time_from_launch:05.2f} | Alt:{altitude:5.1f}m | "
                                           f"Angle:{angle_deg:05.2f}° → {commanded_angle:05.2f}° Pulse:{pulse} ***")
                            if print_str != last_print_str:
                                print(print_str)
                                last_print_str = print_str
                        
                        last_angle = commanded_angle
                except Exception as e:
                    print(f"Angle calculation error at t={time_from_launch:.3f}: {e}")
            
            # Check if we've run out of data
            if time_from_launch > time_keys[-1]:
                print("Reached end of trajectory data.")
                button_state.tracking_started = False
                launch_time = None
                button_state.mode = 'manual'
                button_state.manual_angle = 0
                print("Reverting to manual mode.")
            
            time.sleep(0.003)  # Small sleep to avoid excessive CPU usage
    
    except KeyboardInterrupt:
        print("\nShutdown requested by user.")
    finally:
        # Stop servo pulses and close pigpio connection
        try:
            pi.set_servo_pulsewidth(PIN, 0)
        except Exception:
            pass
        pi.stop()
        if COMPUTE_DIAGNOSTICS:
            print(f"Max torque during session: {max_torque_session:.4f} Nm")
        print("Servo control finished.")


if __name__ == "__main__":
    main()
