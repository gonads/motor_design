function dfdt = gravity_assist_eq(t,y)

global M0 T A Cd R0 Ht mdot 

v  =  y(1);     
gm =  y(2);    
x  =  y(3);     
h  =  y(4);    
vD =  y(5);     
vG =  y(6);    

m = M0 - mdot*t;  

g  = atmosphere_model_km(h);       
rho = atmosphere_model_km(h);


D = 1/2 * rho*v^2 * A * Cd;  
if h <= Ht 
    dv_dt = T/m - D/m - g;
    dgm_dt = 0;
    dx_dt = 0;
    dh_dt = v;
    dvG_dt = -g;
    
else                               
    dv_dt = T/m - D/m - g*sin(gm);
    dgm_dt = -1/v*(g - v^2/(R0 + h))*cos(gm);
    dx_dt = R0/(R0 + h)*v*cos(gm) + 463*sin(gm);
    dh_dt = v*sin(gm) + 463*cos(gm); 
    dvG_dt = -g*sin(gm);             
end

dvD_dt = -D/m;           
dfdt = [dv_dt,dgm_dt, dx_dt,dh_dt, dvD_dt, dvG_dt]';

return

