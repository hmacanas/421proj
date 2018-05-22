% Henry Macanas
% Eric Ekstrom
% Michael Johnston
% 421 Assignment 6
clc, close all, clear all

addpath('overhead_functions/')

%% orbital parameters
h = 53335.2;
ecc = 0;
RAAN = 0;
inc = deg2rad(98.43);
W = 0;
nu = 0;

% initial r and v vectors
[r0,v0 ] = coes2rv(ecc,inc,RAAN,W,h,nu );

% retrieves period
[ ~,~,~,~,~,~,P,~ ] = coes( r0,v0 );

% angular velcoity of lvlh frame
w_lvlh_eci_eci = cross(r0,v0)/norm(r0)^2; % this is wrong, should be a 3x1 not 1x1

%% rigid body prop

% lvlh to eci tranformation 
C_lvlh_eci= eci2lvlh(r0,v0);
C_eci_lvlh = C_lvlh_eci';
C_body_lvlh = eye(3);
C_body_eci = C_body_lvlh*C_lvlh_eci;

% initial states in reference to eci        
w0_body_eci = [0;-2*pi/P;0];

r0_eci_eci = r0;
v0_eci_eci = v0;
euler_angles0_eci = euler_angs(C_body_eci);
q0_eci_eci = quaternion(C_body_eci);

% states relative to lvlh
euler_angs0_lvlh = [0;0;0];
w_body_eci = w0_body_eci;
w_body_lvlh0 = w_body_eci - C_body_eci*w_lvlh_eci_eci;
q0_body_lvlh = [0;0;0;1];

state = [euler_angles0_eci;w0_body_eci;q0_eci_eci;r0_eci_eci;v0_eci_eci;...
    euler_angs0_lvlh;w_body_lvlh0;q0_body_lvlh];

% Distances to individual centers of mass
Rbb = zeros(3,1); % Distance from bus COM to bus COM
Rbsp = [0; 2.5; 0]; % Distance from bus COM to solar panel COM
Rbsens = [0; 0; 1.5]; % Distance from bus COM to sensor COM
masses = [500,20,20,100]; % Component masses
dims = [2 2 2 0.25; 2 3 3 0.25; 2 0.05 0.05 1]; % Component dimensions
[consts.COM, consts.I] = getCOM(masses,...
	[Rbb Rbsp -Rbsp Rbsens],...
	dims); % Spacecraft center of mass

% Normal vectors
consts.n   = [1 0 -1  0 0  0    1 1 -1 -1 0 0  0  0    1 -1 0  0;
       0 1  0 -1 0  0    0 0  0  0 0 0  0  0    0  0 1 -1;
       0 0  0  0 1 -1    0 0  0  0 1 1 -1 -1    0  0 0  0]; % Bus    Solar Panels    Sensor
% Rho vectors
consts.rho = consts.COM - [2 0 -2  0 0  0    2  2 -2 -2    0    0     0     0    .125 -.125    0     0;
             0 2  0 -2 0  0    4 -4  4 -4    4   -4     4    -4       0     0 .125 -.125;
             0 0  0  0 2 -2    0  0  0  0 .025 .025 -.025 -.025     1.5   1.5  1.5   1.5]; % Bus    Solar Panels    Sensor
consts.A   = [4 4 4 4 4 4    .15 .15 .15 .15 6 6 6 6    .25 .25 .25 .25]; % Bus    Solar Panels    Sensor

% -- ODE call
Torque = 'yes';
tspan = [0 3*P];
options = odeset('RelTol',1e-8,'AbsTol',1e-8, 'OutputFcn',@(t,y,flag,varargin) odeOutFunc(t,y,flag));
[tnew, statenew] = ode45(@day_func,tspan,state,options,Torque,consts);
% Save and load solutions for speed
% save('soln','tnew','statenew')
% load('soln')

%% Body rel to ECI Plots
figure
subplot(2,3,1)
set(groot,'DefaultAxesXGrid','on', 'DefaultAxesYGrid','on')
plot(tnew,statenew(:,4:6),'LineWidth',2)
title('Absolute Angular Velocity of Spacecraft: F_b relative to F_{ECI}')
xlabel('Time (s)')
ylabel('Angular Velocity (rads/s)')
legend('\omega_x','\omega_y','\omega_z')

subplot(2,3,2)
plot(tnew,rad2deg(statenew(:,1:3)),'LineWidth',2)
title('Euler Angles from F_b to F_{ECI}')
xlabel('Time (s)')
ylabel('Angle (degrees)')
legend('\phi','\theta','\psi')

subplot(2,3,3)
plot(tnew,statenew(:,7:10),'LineWidth',2)
title('Quaternion Components from F_b to F_{ECI}')
xlabel('Time (s)')
ylabel('Magnitude (None)')
legend('\epsilon_x','\epsilon_y','\epsilon_z','\eta')

%% Body rel to LVLH Plots
subplot(2,3,4)
plot(tnew,statenew(:,20:22),'LineWidth',2)
title('Absolute Angular Velocity of Spacecraft: F_b relative to F_{LVLH}')
xlabel('Time (s)')
ylabel('Angular Velocity (rads/s)')
legend('\omega_x','\omega_y','\omega_z')

subplot(2,3,5)
plot(tnew,rad2deg(statenew(:,17:19)),'LineWidth',2)
title('Euler Angles from F_b to F_{LVLH}')
xlabel('Time (s)')
ylabel('Angle (degrees)')
legend('\phi','\theta','\psi')

subplot(2,3,6)
plot(tnew,statenew(:,23:26),'LineWidth',2)
title('Quaternion Components from F_b to F_{LVLH}')
xlabel('Time (s)')
ylabel('Magnitude (None)')
legend('\epsilon_x','\epsilon_y','\epsilon_z','\eta')

%% Orbit plot
figure
hold on
plot3(statenew(:,11),statenew(:,12),statenew(:,13))
plot3(statenew(1,11),statenew(1,12),statenew(1,13),'*')

%% Total Torque
figure
plot(Torques.tot(:,1),Torques.tot(:,2:4), 'lineWidth', 2)
grid on
title('Disturbance Torques')
xlabel('Time (seconds)'), ylabel('Disturbance Torque [Nm]')
legend('Tx', 'Ty', 'Tz')

%% Individual Torques
figure
subplot(2,2,1)
plot(Torques.grav(:,1),Torques.grav(:,2:4), 'lineWidth', 2)
grid on
title('Gravity Torques')
xlabel('Time (seconds)'), ylabel('Disturbance Torque [Nm]')
legend('Tx', 'Ty', 'Tz')

subplot(2,2,2)
plot(Torques.srp(:,1),Torques.srp(:,2:4), 'lineWidth', 2)
grid on
title('Solar Radiation Pressure Torques')
xlabel('Time (seconds)'), ylabel('Disturbance Torque [Nm]')
legend('Tx', 'Ty', 'Tz')

subplot(2,2,3)
plot(Torques.drag(:,1),Torques.drag(:,2:4), 'lineWidth', 2)
grid on
title('Drag Torques')
xlabel('Time (seconds)'), ylabel('Disturbance Torque [Nm]')
legend('Tx', 'Ty', 'Tz')

subplot(2,2,4)
plot(Torques.mag(:,1),Torques.mag(:,2:4), 'lineWidth', 2)
grid on
title('Magnetic Field Torques')
xlabel('Time (seconds)'), ylabel('Disturbance Torque [Nm]')
legend('Tx', 'Ty', 'Tz')

