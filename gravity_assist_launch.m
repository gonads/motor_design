close all; clear all; clc;

global M0 G0 T A Cd R0 Ht tf mdot

G0 = 9.81;        
R0 = 6371E3;         

H = 0;             

Mfrac1 = 0.8;
M0 = 1000; 
Mf1 = M0 * (1-Mfrac1);
Mp1 = M0 * Mfrac1;       

ISP = 200 ;        
tb = 100;       
Ht = 200;      

r = 0.25;           
A = pi*r^2;      
Cd = 0.8 ;            

mdot = Mp1/tb; 
T = mdot*ISP*G0;   
Mf1 = M0 - Mp1;     

t0 = 0;             
tf = t0 + tb;      
tr = [t0,tf];  

gamma0 = 89.5/180*pi;      
v0 = 0;   
x0 = 0;   
h0 = H;
vD0 = 0;  
vG0 = 0;
state0 = [v0, gamma0, x0, h0, vD0, vG0];

[t,state] = ode45(@gravity_assist_eq, tr, state0);
v = state(:,1)/1000;    
gamma = state(:,2)*180/pi;    
x = state(:,3)/1000;      
h = state(:,4)/1000;     
vD = -state(:,5)/1000;   
vG = -state(:,6)/1000;    

fprintf('\n Final speed               = %4.2f, (km/s)',v(end))
fprintf('\n Final flight path angle   = %4.2f,(deg)',gamma(end))
fprintf('\n Drag loss                 = %4.2f, (km/s)',vD(end))
fprintf('\n Gravity loss              = %4.2f, (km/s)',vG(end))
fprintf('\n');

figure;
plot(t, h, t, x);
title('Flight Profile');
xlabel('Time, (s)');
ylabel('Distance, (km)');
set(gca,'TickDir','out');
legend('Altitude','Range');
grid minor;

