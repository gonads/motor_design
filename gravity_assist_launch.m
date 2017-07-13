close all; clear all; clc;

global G0 M0 ME R0 T A Cd Ht tf mdot 

G0 = 6.67408E-11; % Gravitational constant, (m^3/kg.s^2).
ME = 5.972E24; % Mass of the Earth, (kg).     
R0 = 6371E3; % Radius of the Earth, (m).
H = 0; % Initial height above mean sea level, (m).    

Mfrac1 = 0.8; % Stage 1 mass fraction.
M0 = 2000; % Stage 1 take-off mass, (kg).
Mp1 = M0 * Mfrac1; % Stage 1 propellant mass, (kg).
Mf1 = M0 * (1-Mfrac1); % Stage 1 hardware mass, (kg).
  
Ht = 200; % Stage 1 gravity turn height, (m).

ISP = 200 ; % Stage 1 motor specific impulse, (s). TO BE MADE DYNAMIC USING MOTOR SIMULATION.       
tb = 60; % Stage 1 motor specific burn time, (s). TO BE MADE DYNAMIC USING MOTOR SIMULATION.       
mdot = Mp1 / tb; % Stage 1 motor mass flow rate, (kg/s). TO BE MADE DYNAMIC USING MOTOR SIMULATION. 
T = mdot * ISP * 9.807; % Stage 1 motor thrust, (N). TO BE MADE DYNAMIC USING MOTOR SIMULATION.    

r = 0.25; % Rocket radius, (m).           
A = pi*r^2; % Rocket frontal area, (m^2).
Cd = 0.8 ; % Rocket drag coefficient.

t0 = 0; % Set launch at time equal to zero, (s).             
tf = t0 + tb; % Set time when propellant tank empty, (s).
tr = [t0, tf]; % Integration interval, (s).

gamma0 = 89.5/180*pi; % Set initial flight angle, (rad).     
v0 = 0; % Rocket velocity taking Earth's rotation into account, (m/s).
x0 = 0; % Down range distance, (km).
h0 = H; % Set initial altitude to launch site altitide, (km).
vD0 = 0; % Loss due to drag, (N).
vG0 = 0; % Loss due to gravity, (N).
state0 = [v0, gamma0, x0, h0, vD0, vG0];

% Now solve problem using ODE approach.
opts = odeset('Reltol',1e-12,'AbsTol',1e-13); % Set error tolerances.
[t,state] = ode45(@gravity_assist_eq, tr, state0, opts); 
v = state(:,1)/1000;    
gamma = state(:,2)*180/pi;    
x = state(:,3)/1000;      
h = state(:,4)/1000;     
vD = -state(:,5)/1000;   
vG = -state(:,6)/1000;    

% Plotting
figure;
yyaxis right;
plot(t, h, t, x, '-b');
title('Position and Angular Profile');
xlabel('Time, (s)');
ylabel('Distance, (km)');
yyaxis left;
plot(t, gamma, 'k');
ylabel('Angle, (deg)');
set(gca,'TickDir','out');
legend('Flight Angle', 'Altitude','Downrange');
grid minor;

figure;
plot(t, v, t, vD, t, vG);
title('Velocity Profile');
xlabel('Time, (s)');
ylabel('Velocity, (km/s)');
set(gca,'TickDir','out');
legend('Velocity', 'Drag Loss','Gravity Loss');
grid minor;

