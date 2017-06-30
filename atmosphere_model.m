% Written by Nadeem Gabbani - 28/06/2017

function [Ta, p, rho, a] = atmosphere_model(H)

    G0 = 6.67408E-11; % Gravitational constant, (m^3/kg.s^2).
    M0 = 5.972E24; % Mass of Earth, (kg).
    R0 = 6371E3; % Radius of Earth, (m).

    % Can these constants be named meanifully (i.e., define them in variables?)
    if H < 11000
        Ta = 273.1 + (15.04 - (0.00649 * H));
        p = 101.29 * ((Ta / 288.08)^5.256);
        rho = p / (0.2869 * Ta);  
        a = (G0*M0) / (R0+H)^2;

    elseif H < 25000 
        Ta = 273.1 - 56.46;
        p = 22.65 * exp(1.73 - 0.000157 * H);
        rho = p / (0.2869 * Ta);  
        a = (G0 * M0) / (R0 + H)^2;

    else
        Ta = 273.1 - 131.21 + (0.00299 * H);
        p = 2.488 * ((Ta / 216.6)^(-11.388)); 
        rho = p / (0.2869 * Ta);
        a = (G0 * M0) / (R0 + H)^2;
    end
end
