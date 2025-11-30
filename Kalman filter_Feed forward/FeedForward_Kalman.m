clc; clear; close all;

% Load data
load("path2sensor_output.mat") 

% INS Mechanization with EKF
[pos_history, vel_history, euler_history, pos_kf, vel_kf, euler_kf] = ...
    INS_EKF(t, f_b_IMU, w_b_ib_IMU, r_n, v_n, roll, pitch, yaw, GPS, var_GPS, Q);

% Position plots with Kalman
figure('Name','Position Comparison','NumberTitle','off');
subplot(3,1,1);
plot(t, r_n(1,:), 'r--', t, pos_history(:,1), 'b', t, pos_kf(:,1), 'g-.','LineWidth',1.5);
title('Latitude'); xlabel('Time (s)'); ylabel('rad'); legend('Truth','INS','EKF'); grid on;
subplot(3,1,2);
plot(t, r_n(2,:), 'r--', t, pos_history(:,2), 'b', t, pos_kf(:,2), 'g-.','LineWidth',1.5);
title('Longitude'); xlabel('Time (s)'); ylabel('rad'); legend('Truth','INS','EKF'); grid on;
subplot(3,1,3);
plot(t, r_n(3,:), 'r--','LineWidth',2)
hold on
plot(t, pos_history(:,3), 'b', t, pos_kf(:,3), 'g-.','LineWidth',1.5);
title('Height'); xlabel('Time (s)'); ylabel('m'); legend('Truth','INS','EKF'); grid on;
% Velocity plots with Kalman
figure('Name','Velocity Comparison','NumberTitle','off');
subplot(3,1,1);
plot(t, vel_history(:,1), 'b', t, v_n(1,:), 'r--', t, vel_kf(:,1), 'g-.','LineWidth',1.5);
title('v_N'); xlabel('Time (s)'); ylabel('m/s'); legend('INS','Truth','EKF'); grid on;
subplot(3,1,2);
plot(t, vel_history(:,2), 'b', t, v_n(2,:), 'r--', t, vel_kf(:,2), 'g-.','LineWidth',1.5);
title('v_E'); xlabel('Time (s)'); ylabel('m/s'); legend('INS','Truth','EKF'); grid on;
subplot(3,1,3);
plot(t, vel_history(:,3), 'b', t, v_n(3,:), 'r--', t, vel_kf(:,3), 'g-.','LineWidth',1.5);
title('v_D'); xlabel('Time (s)'); ylabel('m/s'); legend('INS','Truth','EKF'); grid on;

% Attitude plots with Kalman
figure('Name','Attitude Comparison','NumberTitle','off');
subplot(3,1,1);
plot(t, euler_history(:,1), 'b', t, roll, 'r--', t, euler_kf(:,1), 'g-.','LineWidth',1.5);
title('Roll (\phi)'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth','EKF'); grid on;
subplot(3,1,2);
plot(t, euler_history(:,2), 'b', t, pitch, 'r--', t, euler_kf(:,2), 'g-.','LineWidth',1.5);
title('Pitch (\theta)'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth','EKF'); grid on;
subplot(3,1,3);
plot(t, euler_history(:,3), 'b', t, yaw, 'r--', t, euler_kf(:,3), 'g-.','LineWidth',1.5);
title('Yaw (\psi)'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth','EKF'); grid on;


% INS Mechanization with EKF Function 
function [pos_history, vel_history, euler_history, pos_kf, vel_kf, euler_kf] = ...
    INS_EKF(t, f_b, w_b_ib, r_n, v_n, roll, pitch, yaw, GPS, var_GPS, Q)

omega_e = 7.292115e-5;

% Initial states
pos = r_n(:,1);
vel = v_n(:,1);
att = [roll(1); pitch(1); yaw(1)];
q = euler2quat(att(1), att(2), att(3));

N = length(t);
pos_history = zeros(N,3);
vel_history = zeros(N,3);
euler_history = zeros(N,3);
pos_kf = zeros(N,3);
vel_kf = zeros(N,3);
euler_kf = zeros(N,3);

pos_history(1,:) = pos';
vel_history(1,:) = vel';
euler_history(1,:) = att';
pos_kf(1,:) = pos';
vel_kf(1,:) = vel';
euler_kf(1,:) = att';

% EKF states
pos_c = pos;
vel_c = vel;
q_c = q;

% Kalman filter initialization
P = eye(9) * 1e-6;
x_hat = zeros(9,1);  % State vector
% Process noise 
Q_cont = Q;
% Measurement noise
R = diag(var_GPS);
% Measurement matrix
H = [eye(3), zeros(3), zeros(3)];

for k = 2:N
    dt = t(10) - t(9);
    
    fb_k = f_b(:,k-1);
    wib_b_k = w_b_ib(:,k-1);
    
    % INS mechanization 
    [pos, vel, q] = dynamic(pos, vel, q, fb_k, wib_b_k, dt, omega_e);
    pos_history(k,:) = pos';
    vel_history(k,:) = vel';
    euler_history(k,:) = dcm2euler(quat2dcm(q))';
    
    %  transition matrix PHI and process noise Q_k
    Cbn = quat2dcm(q);
    f_n = Cbn * fb_k;
    [PHI, G] = calculate_phi(pos, vel, f_n, dt, omega_e, Cbn);
    Q_k = PHI * G * Q_cont * G' * PHI' * dt;
    
    % Prediction 
    x_hat = PHI * x_hat;
    P = PHI * P * PHI' + Q_k;
    
    % Update when GPS available
    if ~isnan(GPS(1,k))
        % GPS measurement (position)
        z_gps = GPS(1:3,k);
        
        % Innovation: z = INS_position - GPS_position 
        phi_cur = pos(1); h_cur = pos(3);
        [M_r, N_r] = calM_N(phi_cur);
        
        % Scaled measurement equation 
        z = [(M_r+h_cur)*(pos(1) - z_gps(1));
             (N_r+h_cur)*cos(phi_cur)*(pos(2) - z_gps(2));
             pos(3) - z_gps(3)];
        
        H_scaled = [M_r+h_cur, 0, 0, zeros(1,6);
                    0, (N_r+h_cur)*cos(phi_cur), 0, zeros(1,6);
                    0, 0, 1, zeros(1,6)];
        
        R_scaled = diag([R(1,1)*(M_r+h_cur)^2, R(2,2)*((N_r+h_cur)*cos(phi_cur))^2, R(3,3)]);
        
        % Kalman gain
        S = H_scaled * P * H_scaled' + R_scaled;
        K = P * H_scaled' / S;
        
        % State update
        x_hat = x_hat + K * (z - H_scaled * x_hat);
        
        % Covariance update
        P = (eye(9) - K * H_scaled) * P;
    end
    
    % Feedforward
    pos_corrected = pos - x_hat(1:3);
    vel_corrected = vel - x_hat(4:6);
    
    % Attitude correction
    eps = x_hat(7:9);
    E_n = [0, -eps(3), eps(2);
           eps(3), 0, -eps(1);
           -eps(2), eps(1), 0];
    Cbn_corrected = (eye(3) + E_n) * Cbn;
    euler_corrected = dcm2euler(Cbn_corrected);
    
    pos_kf(k,:) = pos_corrected';
    vel_kf(k,:) = vel_corrected';
    euler_kf(k,:) = euler_corrected';
end
end

% Calculate PHI and G matrices
function [PHI, G] = calculate_phi(pos, vel, f_n, dt, omega_e, Cbn)
phi = pos(1); h = pos(3);
vn = vel(1); ve = vel(2); vd = vel(3);
[M, N] = calM_N(phi);

% F_rr 
F_rr = [0, 0, -vn/(M+h)^2;
        ve*sin(phi)/((N+h)*cos(phi)^2), 0, -ve/((N+h)^2*cos(phi));
        0, 0, 0];

% F_rv
F_rv = [1/(M+h), 0, 0;
        0, 1/((N+h)*cos(phi)), 0;
        0, 0, -1];

% F_vr 
a1 = 9.7803267715; a2 = 0.0052790414; a3 = 0.0000232718;
a4 = -0.0000030876910891; a5 = 0.0000000043977311; a6 = 0.0000000000007211;
gamma = a1*(1 + a2*(sin(phi))^2 + a3*(sin(phi))^4) + (a4 + a5*(sin(phi))^2)*h + a6*h^2;
%gamma = 9.8;
R_e = sqrt(M*N);
F_vr = [-2*ve*omega_e*cos(phi) - ve^2/((N+h)*cos(phi)^2), 0, ...
        -vn*vd/(M+h)^2 + ve^2*tan(phi)/(N+h)^2;
        2*omega_e*(vn*cos(phi)-vd*sin(phi)) + ve*vn/((N+h)*cos(phi)^2), 0, ...
        -ve*vd/(N+h)^2 - vn*ve*tan(phi)/(N+h)^2;
        2*ve*omega_e*sin(phi), 0, ve^2/(N+h)^2 + vn^2/(M+h)^2 - 2*gamma/(R_e+h)];

% F_vv 
F_vv = [vd/(M+h), -2*omega_e*sin(phi)-2*ve*tan(phi)/(N+h), vn/(M+h);
        2*omega_e*sin(phi)+ve*tan(phi)/(N+h), (vd+vn*tan(phi))/(N+h), 2*omega_e*cos(phi)+ve/(N+h);
        -2*vn/(M+h), -2*omega_e*cos(phi)-2*ve/(N+h), 0];

% F_ve (f_n x)
F_ve = [0, -f_n(3), f_n(2);
        f_n(3), 0, -f_n(1);
        -f_n(2), f_n(1), 0];

% F_er 
F_er = [-omega_e*sin(phi), 0, -ve/(N+h)^2;
        0, 0, vn/(M+h)^2;
        -omega_e*cos(phi)-ve/((N+h)*cos(phi)^2), 0, ve*tan(phi)/(N+h)^2];

% F_ev 
F_ev = [0, 1/(N+h), 0;
        -1/(M+h), 0, 0;
        0, -tan(phi)/(N+h), 0];

% omega_in_n
omega_in_n = [omega_e*cos(phi) + ve/(N+h);
              -vn/(M+h);
              -omega_e*sin(phi) - ve*tan(phi)/(N+h)];

% F_ee = -(omega_in_n x)
F_ee = -[0, -omega_in_n(3), omega_in_n(2);
         omega_in_n(3), 0, -omega_in_n(1);
         -omega_in_n(2), omega_in_n(1), 0];

% Full F matrix
F = [F_rr, F_rv, zeros(3);
     F_vr, F_vv, F_ve;
     F_er, F_ev, F_ee];

% Transition matrix 
PHI = eye(9) + F * dt;

% G matrix for process noise
G = [zeros(3), zeros(3);
     Cbn, zeros(3);
     zeros(3), -Cbn];
end

% Sub functions
function [M,N] = calM_N(phi)
a = 6378136.6; b = 6356751.9;
e = sqrt(1-(b/a)^2);
N = a / sqrt(1 - (e*sin(phi))^2);
M = a*(1 - e^2) / (1 - (e*sin(phi))^2)^(1.5);
end

function g_n = gravity(phi,h)
a1 = 9.7803267715; a2 = 0.0052790414; a3 = 0.0000232718;
a4 = -0.0000030876910891; a5 = 0.0000000043977311; a6 = 0.0000000000007211;
gamma = a1*(1 + a2*(sin(phi))^2 + a3*(sin(phi))^4) + (a4 + a5*(sin(phi))^2)*h + a6*h^2;
g_n = [0;0;gamma];
end

function [pos, vel, q] = dynamic(pos, vel, q, f_b, w_b_ib, dt, omega_e)
phi = pos(1); h = pos(3);
vn = vel(1); ve = vel(2);

[M, N] = calM_N(phi);
g_n = gravity(phi,h);

omega_ie_n = [omega_e*cos(phi); 0; -omega_e*sin(phi)];
omega_en_n = [ve/(N+h); -vn/(M+h); -ve*tan(phi)/(N+h)];
omega_in_n = omega_ie_n + omega_en_n;

Cbn = quat2dcm(q);
Cold = Cbn;

omega_nb_b = w_b_ib - Cbn' * omega_in_n;
delta_theta = omega_nb_b * dt;
delta_theta_mag = norm(delta_theta);

if delta_theta_mag > 1e-8
    s = (2/delta_theta_mag)*sin(delta_theta_mag/2);
    c = 2*(cos(delta_theta_mag/2)-1);
else
    s = 1 - delta_theta_mag^2/24;
    c = -delta_theta_mag^2/4;
end

Omega = [c, s*delta_theta(3), -s*delta_theta(2), s*delta_theta(1);
        -s*delta_theta(3), c, s*delta_theta(1), s*delta_theta(2);
        s*delta_theta(2), -s*delta_theta(1), c, s*delta_theta(3);
        -s*delta_theta(1), -s*delta_theta(2), -s*delta_theta(3), c];

q = q + 0.5*Omega*q;
q = q / norm(q);

Cbn = quat2dcm(q);
matrix = eye(3) + 0.5*[0, delta_theta(3), -delta_theta(2);
                       -delta_theta(3), 0, delta_theta(1);
                       delta_theta(2), -delta_theta(1), 0];

delta_vf_n = Cold * matrix * (f_b * dt);
delta_v_n = delta_vf_n - cross((2*omega_ie_n + omega_en_n), vel) * dt + g_n * dt;

vel_new = vel + delta_v_n;

D_inv = [1/(M+h), 0, 0;
         0, 1/((N+h)*cos(phi)), 0;
         0, 0, -1];
pos = pos + 0.5 * D_inv * (vel + vel_new) * dt;
vel = vel_new;
end

function Cbn = quat2dcm(q)
q1=q(1); q2=q(2); q3=q(3); q4=q(4);
Cbn = [1-2*q2^2-2*q3^2, 2*(q1*q2-q3*q4), 2*(q1*q3+q2*q4);
       2*(q1*q2+q3*q4), 1-2*q1^2-2*q3^2, 2*(q2*q3-q1*q4);
       2*(q1*q3-q2*q4), 2*(q2*q3+q1*q4), 1-2*q1^2-2*q2^2];
end

function euler = dcm2euler(Cbn)
c11=Cbn(1,1); c21=Cbn(2,1); c31=Cbn(3,1);
c32=Cbn(3,2); c33=Cbn(3,3);
pitch=-asin(c31); roll=atan2(c32,c33); yaw=atan2(c21,c11);
euler=[roll; pitch; yaw];
end

function q = euler2quat(roll,pitch,yaw)
cy=cos(yaw/2); sy=sin(yaw/2);
cp=cos(pitch/2); sp=sin(pitch/2);
cr=cos(roll/2); sr=sin(roll/2);
q=zeros(4,1);
q(4)=cr*cp*cy + sr*sp*sy;
q(1)=sr*cp*cy - cr*sp*sy;
q(2)=cr*sp*cy + sr*cp*sy;
q(3)=cr*cp*sy - sr*sp*cy;
q=q/norm(q);
end