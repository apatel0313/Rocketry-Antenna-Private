import pigpio
import time
import math
from bisect import bisect_left

pi = pigpio.pi()
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

def get_altitude_at_time(time_keys, t:float):
    #   Return altitude corresponding to the closest time fig in CSV
    idx = bisect_left(time_keys, t)
    if idx == 0:
        nearest_key = time_keys[0]
    elif idx == len(time_keys):
        nearest_key = time_keys[-1]
    else:
        before = time_keys[idx - 1]
        after = time_keys[idx]
        #   pick whichever is closer to t
        if t - before <= after - t:
            nearest_key = before
        else:
            nearest_key = after
    return time_altitude_dict[nearest_key]



def main():
    time_altitude_dict = build_dict()
    time_keys = sorted(time_altitude_dict.keys())
    ground_station_from_pole = 1609.34
    launch_time = time.time() 
    time_from_last = 1
    fail_count=0
    last_altitude=-1
    while True:
        if time.time() - time_from_last > 0.03: # ensures at least 3 ms between loops
            time_from_last = time.time()
            time_from_launch = time.time() - launch_time
            
        time_from_launch = time.time() - launch_time 

        try:
            altitude = get_altitude_at_time(time_keys, time_from_launch)
        except Exception as e:
            print(f"Lookup error at t={time_from_launch:.3f}: {e}")
            fail_count += 1
            if fail_count > 5:
                break
            continue

        angle = round(math.degrees(math.atan(altitude/ground_station_from_pole)),2)
        if altitude != last_altitude:
            pulse = rest_pulse + deg_step_pulse * angle
            pi.set_servo_pulsewidth(PIN, pulse)                
            last_altitude = altitude
            print_str = f"T: {time_from_launch:05.2f} | Altitude: {int(altitude):04d} | Angle: {angle:05.2f}"
            print(print_str)  
        else:
            pass

if __name__ == "__main__":
    main()