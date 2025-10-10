% Rotation Conversion inputs
clc; clear; close all;

disp('Select input type:');
disp('1: Euler angles (phi, theta, psi)');
disp('2: Quaternion [q1; q2; q3; q4] (scalar last)');
disp('3: DCM (3x3)');
disp('4: Axis-Angle [angle, axis]'); 
inputType = input('Enter choice (1–4): ');

switch inputType
    case 1
        unitChoice = input('Are the angles in degrees or radians? Enter 1 for degrees, 2 for radians: ');
        phi = input('Enter roll φ: ');
        theta = input('Enter pitch θ: ');
        psi = input('Enter yaw ψ: ');
        
        if unitChoice == 1
            phi = deg2rad(phi);
            theta = deg2rad(theta);
            psi = deg2rad(psi);
        end

        disp('Select output type:');
        disp('1: Quaternion');
        disp('2: DCM');
        choice = input('Enter choice (1–2): ');
        if choice == 1
            q = euler_to_quat(phi, theta, psi);
            disp('Quaternion [q1; q2; q3; q4]:');
            disp(q);
        elseif choice == 2
            Cbn = euler_to_dcm(phi, theta, psi);
            disp('Direction Cosine Matrix (Cbn):');
            disp(Cbn);
        end

    case 2
        q = input('Enter quaternion [q1; q2; q3; q4] (column or row): ');
        q = q(:);
        disp('Select output type:');
        disp('1: Euler angles');
        disp('2: DCM');
        choice = input('Enter choice (1–2): ');
        if choice == 1
            [phi, theta, psi] = quat_to_euler(q);
            outputUnit = input('Output in degrees or radians? Enter 1 for degrees, 2 for radians: ');
            if outputUnit == 1
                disp('Euler angles [deg]:');
                disp(rad2deg([phi, theta, psi]));
            else
                disp('Euler angles [rad]:');
                disp([phi, theta, psi]);
            end
        elseif choice == 2
            Cbn = quat_to_dcm(q);
            disp('Direction Cosine Matrix (Cbn):');
            disp(Cbn);
        end

    case 3
        disp('Enter DCM row by row:');
        Cbn = zeros(3);
        for i = 1:3
            Cbn(i,:) = input(['Row ' num2str(i) ': ']);
        end
        disp('Select output type:');
        disp('1: Euler angles');
        disp('2: Quaternion');
        choice = input('Enter choice (1–2): ');
        if choice == 1
            [phi, theta, psi] = dcm_to_euler(Cbn);
            outputUnit = input('Output in degrees or radians? Enter 1 for degrees, 2 for radians: ');
            if outputUnit == 1
                disp('Euler angles [deg]:');
                disp(rad2deg([phi, theta, psi]));
            else
                disp('Euler angles [rad]:');
                disp([phi, theta, psi]);
            end
        elseif choice == 2
            q = dcm_to_quat(Cbn);
            disp('Quaternion [q1; q2; q3; q4]:');
            disp(q);
        end

    case 4
        unitChoice = input('Is the angle in degrees or radians? Enter 1 for degrees, 2 for radians: ');
        angle = input('Enter rotation angle: ');
        axis = input('Enter rotation axis [kx; ky; kz]: ');
        axis = axis(:) / norm(axis); % normalize axis
        
        if unitChoice == 1
            angle = deg2rad(angle); % convert to radians
        end

        Cbn = axis_angle_to_dcm(angle, axis);
        disp('Direction Cosine Matrix (Cbn):');
        disp(Cbn);
        
        disp('Select output type:');
        disp('1: Euler angles');
        disp('2: Quaternion');
        choice = input('Enter choice (1–2): ');
        if choice == 1
            [phi, theta, psi] = dcm_to_euler(Cbn);
            outputUnit = input('Output in degrees or radians? Enter 1 for degrees, 2 for radians: ');
            if outputUnit == 1
                disp('Euler angles [deg]:');
                disp(rad2deg([phi, theta, psi]));
            else
                disp('Euler angles [rad]:');
                disp([phi, theta, psi]);
            end
        elseif choice == 2
            q = dcm_to_quat(Cbn);
            disp('Quaternion [q1; q2; q3; q4]:');
            disp(q);
        end
end


% Conversion Functions

function Cbn = euler_to_dcm(phi, theta, psi)
Cbn = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
       sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), ...
       sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), sin(phi)*cos(theta);
       cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), ...
       cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), cos(phi)*cos(theta)];
end

function [phi, theta, psi] = dcm_to_euler(Cbn)
theta = -asin(Cbn(1,3));
phi   = atan2(Cbn(2,3), Cbn(3,3));
psi   = atan2(Cbn(1,2), Cbn(1,1));
end

function Cbn = quat_to_dcm(q)
qv = q(1:3);
qs = q(4);
Cbn = (qs^2 - qv.'*qv)*eye(3) + 2*(qv*qv.') - 2*qs*skew_symmetric(qv);
end

function q = dcm_to_quat(Cbn)
tr = trace(Cbn);
if tr > 0
    s = 0.5 / sqrt(tr + 1);
    qs = 0.25 / s;
    qv = s * [Cbn(3,2) - Cbn(2,3);
              Cbn(1,3) - Cbn(3,1);
              Cbn(2,1) - Cbn(1,2)];
else
    [~, i] = max(diag(Cbn));
    switch i
        case 1
            s = sqrt(1 + Cbn(1,1) - Cbn(2,2) - Cbn(3,3)) * 2;
            qv = [0.25*s;
                  (Cbn(1,2)+Cbn(2,1))/s;
                  (Cbn(1,3)+Cbn(3,1))/s];
            qs = (Cbn(3,2)-Cbn(2,3))/s;
        case 2
            s = sqrt(1 - Cbn(1,1) + Cbn(2,2) - Cbn(3,3)) * 2;
            qv = [(Cbn(1,2)+Cbn(2,1))/s;
                  0.25*s;
                  (Cbn(2,3)+Cbn(3,2))/s];
            qs = (Cbn(1,3)-Cbn(3,1))/s;
        otherwise
            s = sqrt(1 - Cbn(1,1) - Cbn(2,2) + Cbn(3,3)) * 2;
            qv = [(Cbn(1,3)+Cbn(3,1))/s;
                  (Cbn(2,3)+Cbn(3,2))/s;
                  0.25*s];
            qs = (Cbn(2,1)-Cbn(1,2))/s;
    end
end
q = [qv; qs];
q = q / norm(q);
end

function q = euler_to_quat(phi, theta, psi)
c1 = cos(phi/2); s1 = sin(phi/2);
c2 = cos(theta/2); s2 = sin(theta/2);
c3 = cos(psi/2); s3 = sin(psi/2);
q = [s1*c2*c3 - c1*s2*s3;
     c1*s2*c3 + s1*c2*s3;
     c1*c2*s3 - s1*s2*c3;
     c1*c2*c3 + s1*s2*s3];
q = q / norm(q);
end

function [phi, theta, psi] = quat_to_euler(q)
q1=q(1); q2=q(2); q3=q(3); q4=q(4);
phi   = atan2(2*(q4*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));
theta = asin( 2*(q4*q2 - q3*q1) );
psi   = atan2(2*(q4*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));
end  

function Cbn = axis_angle_to_dcm(angle, k)
    K = skew_symmetric(k);
    Cbn = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;
end

function S = skew_symmetric(v)
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end
