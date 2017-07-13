% Written by Nadeem Gabbani - 12/07/2017

function rho = atmosphere_model_km(h)
    
    if h < 11
        Ta = 273.1 + (15.04 - (0.00649 * (h)));
        p = 101.29 * ((Ta / 288.08)^5.256);
        rho = p / (0.2869 * Ta);  
        
    elseif h < 25
        Ta = 273.1 - 56.46;
        p = 22.65 * exp(1.73 - 0.000157 * (h));
        rho = p / (0.2869 * Ta);  
   
    else
        Ta = 273.1 - 131.21 + (0.00299 * (h));
        p = 2.488 * ((Ta / 216.6)^(-11.388)); 
        rho = p / (0.2869 * Ta);
    end
end