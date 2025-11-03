This MATLAB project implements a complete strapdown Inertial Navigation System (INS) mechanization, translating raw IMU sensor data (specific force and angular rate) into a full 9-state navigation solution. 
Based on the fundamental equations of motion detailed in E. H. Shin's research, the script integrates sensor measurements using a quaternion-based attitude update to continuously track a vehicle's position
(latitude, longitude, height), velocity (north, east, down), and attitude (roll, pitch, yaw), demonstrating a core algorithm used in obotics and autonomous vehicle guidance.
