%% Henry Macanas
% 421 Assignment 6
clc, close all, clear all

addpath('overhead_functions/')

%% eom symbolically
 
syms Ixx Iyy Izz wx wy wz dwx dwy dwz Tx Ty Tz

I = [Ixx 0 0;0 Iyy 0;0 0 Izz];
T = [Tx;Ty;Tz];
dw = [dwx;dwy;dwz];
w = [wx;wy;wz];

eq1 = T == I*dw+cross(w,I*w);
eq2 = eq1(1);
eq3 = eq1(2);
eq4 = eq1(3);

[dwx,dwy,dwz] = solve([eq2 eq3 eq4],[dwx dwy dwz]);

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
w_lvlh_eci = 2*pi/P; % this is wrong, should be a 3x1 not 1x1

%% rigid body prop

% lvlh to eci tranformation 
C_lvlh_eci= eci2lvlh(r0,v0);
C_eci_lvlh = C_lvlh_eci';

% initial states inn reference to eci        
w0_body_eci = [0;-2*pi/P;0];

w_body_eci = w0_body_eci;
w_body_lvlh = C_lvlh_eci*w_body_eci - C_lvlh_eci*w_lvlh_eci;

r0_eci_eci = r0;
v0_eci_eci = v0;
euler_angles0_eci = euler_angs(C_eci_lvlh);
q0_eci_eci = quaternion(C_eci_lvlh);
state = [euler_angles0_eci;w0_body_eci;q0_eci_eci;r0_eci_eci;v0_eci_eci];

% ode call
Torque = 'no';
tspan = [0 3*P];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tnew0, statenew0] = ode45(@day_func,tspan,state,options,Torque);

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

%% Functions
function [y]=day_func(t,state,Torque)

% Constants
I = diag([644.5917 377.9250 626.6667]);

% if no torque set torques to zero
if strcmp(Torque,'no')
    T = [0;0;0];
end
    
% attitude motion equatiuons 
wdot = I\(T - cross(state(4:6),I*state(4:6)));
eulrates = euler_rates(state(4:6),state(1),state(2));
quaternion_rates = quatrates(state(4:6),state(7:10));

% orbital motion equations
muearth = 398600;
r = norm(state(11:13));
acc = -muearth*state(11:13)./r^3;

% outputs that will be intergrated 
y = [eulrates;wdot;quaternion_rates;state(14:16); acc];

end