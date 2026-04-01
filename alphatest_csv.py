import math
from bisect import bisect_left
import time

MOMENT_OF_INERTIA = 0.2   # kg*m^2
MAX_TORQUE = 0.4             # N*m
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


def get_altitude_at_time_interpolated(time, time_altitude_dict, time_keys):
    if time in time_altitude_dict:
        return time_altitude_dict[time]
    # find the nearest
    idx = bisect_left(time_keys, time)
    if idx == 0:
        return time_altitude_dict[time_keys[0]]
    if idx == len(time_keys):
        return time_altitude_dict[time_keys[-1]]
    t1 = time_keys[idx-1]
    t2 = time_keys[idx]
    if abs(time - t1) < abs(time - t2):
        return time_altitude_dict[t1]
    else:
        return time_altitude_dict[t2]


def _solve_linear_system(A, b):
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


def main():
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())

    window_s = 2
    dt_out = 0.02
    degree = 4  # adjustable polynomial degree (2 = quadratic, 3 = cubic, 4 = quartic, etc.)
    top_torque = 0.0
    output = []
    t = time_keys[0]
    max_t = time_keys[-1]

    while t <= max_t - window_s:
        window_start = t - window_s / 2
        window_end = t + window_s / 2
        window_keys = [tk for tk in time_keys if window_start <= tk <= window_end]

        if len(window_keys) < degree + 1:
            # fallback: pick the closest points around current t
            idx = bisect_left(time_keys, t)
            num_needed = degree + 1
            start = max(0, idx - num_needed // 2)
            end = min(len(time_keys), start + num_needed)
            window_keys = time_keys[start:end]
            if len(window_keys) < num_needed:
                # If still not enough, skip this t
                t += dt_out
                continue

        theta_window = [math.atan(time_altitude_dict[tk] / GROUND_STATION_FROM_POLE) for tk in window_keys]
        center = t  # center at t for symmetric window
        try:
            coeffs = fit_polynomial(window_keys, theta_window, center, degree)
        except ValueError:
            t += dt_out
            continue

        tau = 0  # at t
        theta_pred = sum(coeffs[k] * (tau ** k) for k in range(degree + 1))
        angle_deg = round(math.degrees(theta_pred), 2)

        alt = get_altitude_at_time_interpolated(t, time_altitude_dict, time_keys)
        output.append({'t': t, 'theta': theta_pred, 'angle_deg': angle_deg, 'alt': alt, 'coeffs': coeffs, 'center': center, 'window_keys': window_keys, 'theta_window': theta_window})
        t += dt_out
        if t > 5: break

    # Compute derivatives analytically
    for item in output:
        tau = item['t'] - item['center']
        coeffs = item['coeffs']
        omega = sum(k * coeffs[k] * (tau ** (k - 1)) for k in range(1, degree + 1))
        alpha = sum(k * (k - 1) * coeffs[k] * (tau ** (k - 2)) for k in range(2, degree + 1))
        item['omega'] = omega
        item['alpha'] = alpha

    # Now compute torque and print
    top_torque = 0.0
    for item in output:
        torque = MOMENT_OF_INERTIA * item['alpha']
        status = "OK" if abs(torque) <= MAX_TORQUE else "EXCEEDS"
        print(f"T:{item['t']:05.2f} |Altitude:{float(item.get('alt') or 0.0):5.2f} |Angle:{item['angle_deg']:04.2f} |Omega:{math.degrees(item['omega']):.2f} |Alpha:{item['alpha']:.3f} |Torque:{torque:.4f} |{status}")
        if abs(torque) > top_torque:
            top_torque = abs(torque)
            if abs(torque) > MAX_TORQUE:
                break
    print(f"Max torque observed: {top_torque:.4f} N*m")
if __name__ == "__main__":
    main()
    