% Written by Nadeem Gabbani - 28/06/2017

clear all; close all; clc;

T = []; % Set initial temperature to zero.
p = []; % Set initial pressure to zero.
k = 0; % Start at step zero.

for H = 0:100000
   k = k + 1;
  
if H < 11000
    T(k) = 273.1+(15.04-(0.00649*H));
    p(k) = 101.29*(((T(k))/288.08).^5.256);
    rho(k) = p(k)/(0.2869*(T(k)));  
    
elseif H < 25000 
    T(k) = 273.1-56.46;
    p(k) = 22.65*exp(1.73-0.000157*H);
    rho(k) = p(k)/(0.2869*(T(k)));  
    
else
    T(k) = 273.1-131.21+(0.00299*H);
    p(k) = 2.488*(((T(k))/216.6)^-11.388); 
    rho(k) = p(k)/(0.2869*(T(k)));
    
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



