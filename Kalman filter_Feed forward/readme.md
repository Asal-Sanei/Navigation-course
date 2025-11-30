This MATLAB code simulates an Inertial Navigation System (INS) and improves its accuracy using an Extended Kalman Filter (EKF) with GPS corrections:

1.Computes position, velocity, and orientation from IMU data.

2.Applies EKF to correct INS drift using GPS measurements.

3.Plots comparisons between truth, INS-only, and EKF-corrected results.

4.Includes utility functions for gravity, coordinate conversions, and quaternions.

Purpose: Demonstrates how combining IMU and GPS data with EKF improves navigation accuracy.
