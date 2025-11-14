clc; clear; close all;

% Load data
load("path2sensor_output.mat") 

% INS Mechanization
[pos_history, vel_history, euler_history] = INS_mechanization(t, f_b, w_b_ib, r_n, v_n, roll, pitch, yaw);

% Position plots
figure('Name','Position','NumberTitle','off');
subplot(3,1,1);
plot(t, pos_history(:,1), 'b', t, r_n(1,:), 'r--','LineWidth',1.5);
title('Latitude'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth'); grid on;

subplot(3,1,2);
plot(t, pos_history(:,2), 'b', t, r_n(2,:), 'r--','LineWidth',1.5);
title('Longitude'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth'); grid on;

subplot(3,1,3);
plot(t, pos_history(:,3), 'b', t, r_n(3,:), 'r--','LineWidth',1.5);
title('Height'); xlabel('Time (s)'); ylabel('m'); legend('INS','Truth'); grid on;

% Velocity plots
figure('Name','Velocity','NumberTitle','off');
subplot(3,1,1);
plot(t, vel_history(:,1), 'b', t, v_n(1,:), 'r--','LineWidth',1.5);
title('v_N'); xlabel('Time (s)'); ylabel('m/s'); legend('INS','Truth'); grid on;

subplot(3,1,2);
plot(t, vel_history(:,2), 'b', t, v_n(2,:), 'r--','LineWidth',1.5);
title('v_E'); xlabel('Time (s)'); ylabel('m/s'); legend('INS','Truth'); grid on;

subplot(3,1,3);
plot(t, vel_history(:,3), 'b', t, v_n(3,:), 'r--','LineWidth',1.5);
title('v_D'); xlabel('Time (s)'); ylabel('m/s'); legend('INS','Truth'); grid on;

% Attitude plots
figure('Name','Attitude','NumberTitle','off');
subplot(3,1,1);
plot(t, euler_history(:,1), 'b', t, roll, 'r--','LineWidth',1.5);
title('Roll (\phi)'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth'); grid on;

subplot(3,1,2);
plot(t, euler_history(:,2), 'b', t, pitch, 'r--','LineWidth',1.5);
title('Pitch (\theta)'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth'); grid on;

subplot(3,1,3);
plot(t, euler_history(:,3), 'b', t, yaw, 'r--','LineWidth',1.5);
title('Yaw (\psi)'); xlabel('Time (s)'); ylabel('rad'); legend('INS','Truth'); grid on;

%Error
pos_error   = pos_history - r_n';
vel_error   = vel_history - v_n';
euler_error = euler_history - [roll; pitch; yaw]';

% Position Error
figure('Name','Position Error','NumberTitle','off');
subplot(3,1,1); plot(t, pos_error(:,1), 'b','LineWidth',1.5); title('Latitude Error'); xlabel('Time (s)'); ylabel('rad'); grid on;
subplot(3,1,2); plot(t, pos_error(:,2), 'b','LineWidth',1.5); title('Longitude Error'); xlabel('Time (s)'); ylabel('rad'); grid on;
subplot(3,1,3); plot(t, pos_error(:,3), 'b','LineWidth',1.5); title('Height Error'); xlabel('Time (s)'); ylabel('m'); grid on;

% Velocity Error
figure('Name','Velocity Error','NumberTitle','off');
subplot(3,1,1); plot(t, vel_error(:,1), 'b','LineWidth',1.5); title('v_N Error'); xlabel('Time (s)'); ylabel('m/s'); grid on;
subplot(3,1,2); plot(t, vel_error(:,2), 'b','LineWidth',1.5); title('v_E Error'); xlabel('Time (s)'); ylabel('m/s'); grid on;
subplot(3,1,3); plot(t, vel_error(:,3), 'b','LineWidth',1.5); title('v_D Error'); xlabel('Time (s)'); ylabel('m/s'); grid on;

% Attitude Error
figure('Name','Attitude Error','NumberTitle','off');
subplot(3,1,1); plot(t, euler_error(:,1), 'b','LineWidth',1.5); title('Roll Error'); xlabel('Time (s)'); ylabel('rad'); grid on;
subplot(3,1,2); plot(t, euler_error(:,2), 'b','LineWidth',1.5); title('Pitch Error'); xlabel('Time (s)'); ylabel('rad'); grid on;
subplot(3,1,3); plot(t, euler_error(:,3), 'b','LineWidth',1.5); title('Yaw Error'); xlabel('Time (s)'); ylabel('rad'); grid on;

%INS Mechanization Function 
function [pos_history, vel_history, euler_history] = INS_mechanization(t, f_b, w_b_ib, r_n, v_n, roll, pitch, yaw)
omega_e = 7.292115e-5; % Earth rotation rate

% Initial states
pos = r_n(:,1);
vel = v_n(:,1);
att = [roll(1); pitch(1); yaw(1)];
q = euler2quat(att(1), att(2), att(3));   % scalar-last [qx;qy;qz;qw]

N = length(t);
pos_history = zeros(N,3);
vel_history = zeros(N,3);
euler_history = zeros(N,3);

pos_history(1,:) = pos';
vel_history(1,:) = vel';
euler_history(1,:) = att';

for k = 2:N
    dt = t(10) - t(9);   

    fb_k    = f_b(:,k-1);
    wib_b_k = w_b_ib(:,k-1);

    [pos, vel, q] = dynamic(pos, vel, q, fb_k, wib_b_k, dt, omega_e);

    pos_history(k,:)   = pos';
    vel_history(k,:)   = vel';
    euler_history(k,:) = dcm2euler(quat2dcm(q))';
end
end

%Dynamics Function 
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

Omega = [  c,                 s*delta_theta(3), -s*delta_theta(2),  s*delta_theta(1);
          -s*delta_theta(3),  c,                 s*delta_theta(1),  s*delta_theta(2);
           s*delta_theta(2), -s*delta_theta(1),  c,                 s*delta_theta(3);
          -s*delta_theta(1), -s*delta_theta(2), -s*delta_theta(3),  c ];

q = q + 0.5*Omega*q;
q = q / norm(q);

Cbn = quat2dcm(q);
matrix = eye(3) + 0.5*[  0,               delta_theta(3), -delta_theta(2);
                        -delta_theta(3),   0,               delta_theta(1);
                         delta_theta(2),  -delta_theta(1),  0 ];

delta_vf_n = Cold * matrix * (f_b * dt);
delta_v_n  = delta_vf_n - cross((2*omega_ie_n + omega_en_n), vel) * dt + g_n * dt;

vel_new = vel + delta_v_n;

D_inv = [1/(M+h), 0, 0;
         0, 1/((N+h)*cos(phi)), 0;
         0, 0, -1];
pos = pos + 0.5 * D_inv * (vel + vel_new) * dt;

vel = vel_new;
end

%Sub functions
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

function Cbn = quat2dcm(q)
q1=q(1); q2=q(2); q3=q(3); q4=q(4);
Cbn = [ 1-2*q2^2-2*q3^2, 2*(q1*q2-q3*q4), 2*(q1*q3+q2*q4);
        2*(q1*q2+q3*q4), 1-2*q1^2-2*q3^2, 2*(q2*q3-q1*q4);
        2*(q1*q3-q2*q4), 2*(q2*q3+q1*q4), 1-2*q1^2-2*q2^2 ];
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
