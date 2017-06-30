function [] = plot_atmosphere_model()

    close all; clear all; clc;

    Ta = [];
    p = [];
    a = [];

    k = 0;
    for H = 0:100000
        k = k + 1;
        [Ta(k), p(k), rho(k), a(k)] = atmosphere_model(H);
    end

    H = 0:100000;

    figure;
    plot(H, Ta); 
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

end
