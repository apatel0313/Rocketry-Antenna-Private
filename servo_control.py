#import pigpio
import time
import math
from bisect import bisect_left

#pi = pigpio.pi()
PIN = 18

max_pulse = 1250 # Set to pulsewidth corresponding to 90 deg - should be about 45 deg from min
rest_pulse = 1750 # Set to pulsewidth corresponding to 0 deg - should be about 45 deg max
deg_step_pulse = 5.56 # set to the pulsewidth corresponding to a 1 degree increase
pulse = rest_pulse
time_altitude_dict={}

def build_dict():
    with open("timeoveraltitude.csv", 'r') as filename:
        for line in filename:
            line = line.split(',')
            time_from_launch_predicted = float(line[0])
            altitude_predicted = float(line[1])
            time_altitude_dict[time_from_launch_predicted] = altitude_predicted
        return time_altitude_dict

# def get_altitude_at_time(time_keys, t:float):
#     #   Return altitude corresponding to the closest time fig in CSV
#     idx = bisect_left(time_keys, t)
#     if idx == 0:
#         nearest_key = time_keys[0]
#     elif idx == len(time_keys):
#         nearest_key = time_keys[-1]
#     else:
#         before = time_keys[idx - 1]
#         after = time_keys[idx]
#         #   pick whichever is closer to t
#         if t - before <= after - t:
#             nearest_key = before
#         else:
#             nearest_key = after
#     return time_altitude_dict[nearest_key]

def get_altitude_at_time_interpolated(time_keys, time_altitude_dict, t:float):
    if t <= time_keys[0]:
        return time_altitude_dict[time_keys[0]]
    if t >= time_keys[-1]:
        return time_altitude_dict[time_keys[-1]]
    idx = bisect_left(time_keys, t) # Find the negibborign point
    t_before = time_keys[idx - 1] #first point
    alt_before = time_altitude_dict[t_before]
    t_after = time_keys[idx] #second point
    alt_after = time_altitude_dict[t_after]
    #math
    slope = (alt_after - alt_before)/(t_after - t_before)
    time_progress = t - t_before
    alt_interpolated = alt_before + (slope * time_progress)

    return alt_interpolated



    



def main():
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())
    ground_station_from_pole = 1609.34
    last_print_str = ""
    launch_time = time.monotonic()
    fail_count = 0
    last_altitude = -1
    last_servo_time = 0.0
    last_sent_pulse = None
    servo_interval = 0.02  # 20 ms
    while True:
        now = time.monotonic()
        time_from_launch = now - launch_time

        try:
            altitude = get_altitude_at_time_interpolated(time_keys, time_altitude_dict, time_from_launch)
        except Exception as e:
            print(f"Lookup error at t={time_from_launch:.3f}: {e}")
            fail_count += 1
            if fail_count > 5:
                break
            continue

        angle = round(math.degrees(math.atan(altitude/ground_station_from_pole)), 2)
        pulse = rest_pulse + deg_step_pulse * angle

        if abs(altitude - last_altitude) >= 0.1:
            # throttle servo updates to once every `servo_interval` seconds
            if now - last_servo_time >= servo_interval and pulse != last_sent_pulse:
                #pi.set_servo_pulsewidth(PIN, pulse)
                last_servo_time = now
                last_sent_pulse = pulse
                indicate_servo_update_str = "***"
            else: indicate_servo_update_str = ""

            last_altitude = altitude
            print_str = f"T: {time_from_launch:05.2f} | Altitude: {int(altitude):05d} | Angle: {angle:05.2f} {indicate_servo_update_str}"
            if print_str != last_print_str:
                print(print_str)
                last_print_str = print_str
        time.sleep(0.003) # sleep 3ms to avoid eating up CPU power

if __name__ == "__main__":
    main()