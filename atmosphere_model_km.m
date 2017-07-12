% Written by Nadeem Gabbani - 12/07/2017

function [Ta, p, rho, g] = atmosphere_model_km(h)

global G0 R0 M0

    G0 = 6.67408E-11; % Gravitational constant, (m^3/kg.s^2).
    M0 = 5.972E24; % Mass of Earth, (kg).
    R0 = 6371E3; % Radius of Earth, (m).

    % Can these constants be named meanifully (i.e., define them in variables?)
    if h < 11
        Ta = 273.1 + (15.04 - (0.00649 * (h*1000)));
        p = 101.29 * ((Ta / 288.08)^5.256);
        rho = p / (0.2869 * Ta);  
        g = (G0*M0) / (R0+(h*1000))^2;

    elseif h < 25
        Ta = 273.1 - 56.46;
        p = 22.65 * exp(1.73 - 0.000157 * (h*1000));
        rho = p / (0.2869 * Ta);  
        g = (G0 * M0) / (R0 + (h*1000))^2;

    else
        Ta = 273.1 - 131.21 + (0.00299 * (h*1000));
        p = 2.488 * ((Ta / 216.6)^(-11.388)); 
        rho = p / (0.2869 * Ta);
        g = (G0 * M0) / (R0 + (h*1000))^2;
    end
end