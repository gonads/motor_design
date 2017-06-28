% Written by Nadeem Gabbani - 28/06/2017

clear all; close all; clc;

G0 = 6.67408E-11; % Gravitational constant, (m^3/kg.s^2).
M0 = 5.972E24; % Mass of Earth, (kg).
R0 = 6371E3; % Radius of Earth, (m).

T = []; % Set initial temperature to zero.
p = []; % Set initial pressure to zero.
a = []; % Set initial gravitation acceleration to zero.
k = 0; % Start at step zero.

for H = 0:100000
   k = k + 1;
  
if H < 11000
    T(k) = 273.1+(15.04-(0.00649*H));
    p(k) = 101.29*(((T(k))/288.08).^5.256);
    rho(k) = p(k)/(0.2869*(T(k)));  
    a(k) = (G0*M0)/(R0+H)^2;
    
elseif H < 25000 
    T(k) = 273.1-56.46;
    p(k) = 22.65*exp(1.73-0.000157*H);
    rho(k) = p(k)/(0.2869*(T(k)));  
    a(k) = (G0*M0)/(R0+H)^2;
    
else
    T(k) = 273.1-131.21+(0.00299*H);
    p(k) = 2.488*(((T(k))/216.6)^-11.388); 
    rho(k) = p(k)/(0.2869*(T(k)));
    a(k) = (G0*M0)/(R0+H)^2;
    
    end
end

H = 0:100000;

figure;
plot(H, T); 
title('Temperature vs. Altitude');
xlabel('Altitude (m)');
ylabel('Temperature (K)');
set(gca,'TickDir','out');
grid minor;
  
figure; 
plot(H, p);   
title('Pressure vs. Altitude');
xlabel('Altitude (m)');
ylabel('Pressure (mbar)');
set(gca,'TickDir','out');
grid minor;

figure; 
plot(H, rho); 
title('Temperature vs. Altitude');
xlabel('Altitude (m)');
ylabel('Density (kg/m^3)');
set(gca,'TickDir','out');
grid minor;

figure; 
plot(H, a); 
title('Gravity vs. Altitude');
xlabel('Altitude (m)');
ylabel('Gravitational Acceleration (m/s^2)');
set(gca,'TickDir','out');
grid minor;

