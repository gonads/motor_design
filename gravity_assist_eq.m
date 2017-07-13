function dfdt = gravity_assist_eq(t,y)

global G0 M0 ME R0 T A Cd Ht mdot 

G0 = 6.67408E-11; % Gravitational constant, (m^3/kg.s^2).
ME = 5.972E24; % Mass of Earth, (kg).     
R0 = 6371E3; % Radius of the Earth, (m).

v  =  y(1); % Rocket Velocity, (m/s).     
gm =  y(2); % Flight angle, (rad). 
x  =  y(3); % Down range distance, (km).   
h  =  y(4); % Altitude, (km).   
vD =  y(5); % Loss due to drag, (N).    
vG =  y(6); % Loss due to gravity, (N).   

m = M0 - mdot*t; % Change of rocket mass during flight, (kg).
rho  = atmosphere_model_km(h); % Change of atmospheric density with altitude, (kg/m^3).
g = (G0*ME) / (R0+(h))^2; % Change of gravitational acceleration with altitude, (m/s^2).
D = 1/2 * rho*v^2 * A * Cd; % Change of drag during flight, (N).

% Behaviour before gravity turn.
if h <= Ht 
    dv_dt = T/m - D/m - g;
    dgm_dt = 0;
    dx_dt = 0;
    dh_dt = v;
    dvG_dt = -g;

% Behaviour during gravity turn.
else                               
    dv_dt = T/m - D/m - g*sin(gm);
    dgm_dt = -1/v*(g - v^2/(R0 + h))*cos(gm);
    dx_dt = R0/(R0 + h)*v*cos(gm) + 463*sin(gm);
    dh_dt = v*sin(gm) + 463*cos(gm); 
    dvG_dt = -g*sin(gm);             
end

% Rate of drag loss.
dvD_dt = -D/m;           
dfdt = [dv_dt,dgm_dt, dx_dt,dh_dt, dvD_dt, dvG_dt]';

return

