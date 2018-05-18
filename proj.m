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

% initial states inn reference to eci        
w0_body_eci = [0;-2*pi/P;0];

w_body_eci = w0_body_eci;
w_body_lvlh0 = w_body_eci - C_lvlh_eci*w_lvlh_eci_eci;

r0_eci_eci = r0;
v0_eci_eci = v0;
euler_angles0_eci = euler_angs(C_lvlh_eci);
q0_eci_eci = quaternion(C_lvlh_eci);
state = [euler_angles0_eci; w0_body_eci;q0_eci_eci;r0_eci_eci;v0_eci_eci;...
	C_lvlh_eci*euler_angles0_eci; C_lvlh_eci*w0_body_eci];

% ode call
Torque = 'no';
tspan = [0 3*P];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tnew0, statenew0] = ode45(@day_func,tspan,state,options,Torque, C_lvlh_eci);

%% no torque plots
figure
plot(tnew0,statenew0(:,4:6),'LineWidth',2)
title('Absolute Angular Velocity of Spacecraft with no Torque: F_b relative to F_eci')
xlabel('Time (s)')
ylabel('Angular Velocity (rads/s)')
legend('\omega_x','\omega_y','\omega_z')

figure
plot(tnew0,rad2deg(statenew0(:,1:3)),'LineWidth',2)
title('Euler Angles from F_b to F_eci')
xlabel('Time (s)')
ylabel('Angular Velocity (rads/s)')
legend('\phi','\theta','\psi')

figure
plot(tnew0,statenew0(:,7:10),'LineWidth',2)
title('Quaternion Components from F_b to F_eci')
xlabel('Time (s)')
ylabel('Magnitude (None)')
legend('\epsilon_x','\epsilon_y','\epsilon_z','\eta')

figure
plot(tnew0,rad2deg(statenew0(:,17:19)),'LineWidth',2)
title('Euler Angles from F_b to F_(LVLH)')
xlabel('Time (s)')
ylabel('Angular Velocity (rads/s)')
legend('\phi','\theta','\psi')

%% Functions
function [y]=day_func(t,state,Torque,C_lvlh_eci)

eulerAngs_eci = state(1:3);
w_b_eci = state(4:6);
q_eci_eci = state(7:10);
r_eci_eci = state(11:13);
v_eci_eci = state(14:16);

% Constants
ns_eci = [-1;0;0]; % Sun vector (constant in ECI)
I = diag([857.091666666667 590.425 626.666666666667]); % SC inertia matrix

% Transformation matrix from ECI to body
C_b_ECI = cx(eulerAngs_eci(1))*cy(eulerAngs_eci(2))*cz(eulerAngs_eci(3));

% Sun vector in body
ns_b = C_b_ECI*ns_eci;

% Velocity vector in body
v_b =  C_b_ECI*r_eci_eci;

% if no torque set torques to zero
if strcmp(Torque,'no')
    T = [0;0;0];
end

% Equations for attitude motion in ECI 
wdot = I\(T - cross(w_b_eci,I*w_b_eci));
eulrates = euler_rates(w_b_eci,eulerAngs_eci(1),eulerAngs_eci(2));
quaternion_rates = quatrates(w_b_eci,q_eci_eci);

% Equations for attitude motion in LVLH
eulerAngs_lvlh = C_lvlh_eci*eulerAngs_eci;
wdot_lvlh = I\(T - cross(C_lvlh_eci*w_b_eci,I*C_lvlh_eci*w_b_eci));
eulrates_lvlh = euler_rates(C_lvlh_eci*w_b_eci,eulerAngs_eci(1),eulerAngs_eci(2));
% quaternion_rates_lvlh = quatrates(w_b_eci,q_eci_eci);

% orbital motion equations
muearth = 398600;
r = norm(r_eci_eci);
acc = -muearth*r_eci_eci./r^3;

% outputs that will be intergrated
% Euler angles, euler rates, quaternions, position, velocity
y = [eulrates;wdot;quaternion_rates;v_eci_eci;acc;eulrates_lvlh;wdot_lvlh];

end
