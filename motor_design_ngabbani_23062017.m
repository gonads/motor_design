clear all; close all; clc; format long;

%Set constants and conditions.
G0 = 9.807; %Gravitation acceleration MSL, (m/s^2).
T0 = 293; %Initial Temperature, (K).
P0 = 101325;% Initial Pressure, (Pa).

%Set fuel properties
G_F = 1.25; %Fuel specific heat ratio.
M_N = 0.0382; %Fuel molecular weight, (kg/mol).
T_F = 1770; %Fuel steady combustion temperature, (K).
R_U = 8.315; %Universal gas constant, (J/mol.K).
R_F = R_U/M_N; %Fuel gas constant, (J/kg.K).
Rho_F = 1800; %Fuel density, (kg/m^3).
a = 8.763; %Fuel burn rate coefficient.
n = -0.314; %Fuel pressure exponent.

%Set motor properties.
Ro = 0.0348/2; %Fuel core outer radius, (m).
Ri = 0.009/2; %Fuel core inner radius, (m).
L = 0.0705; %Fuel core length, (m).
A = 2*pi*Ri*L; %Fuel burn area, (m^2).
V = (pi*(Ro^2)*L)-(pi*(Ri^2)*L); %Fuel core volume
M_F = Rho_F*V; %Fuel mass, (kg).
Rt = 0.0045/2; %Throat radius, (m).
At = pi*(Rt^2); %Throat area, (m^2).
Re = 0.0045/2; %Exit radius, (m).
Ae = pi*(Re^2); %Exit area, (m^2).

%Set initial conditions.
P(1) = P0; %Chamber pressure to atmospheric, (Pa).
Pe(1) = P0; %Exit pressure to atmospheric, (Pa).
T(1) = T0; %Chamber temperature to atmospheric, (K).
M_F(1) = M_F; %Initial mass of fuel, (kg).
L(1) = L; %Initial length of fuel cure, (m).
Ri(1) = Ri; %Initial burn radius, (m).
V(1) = pi*(Ri(1)^2)*L(1); %Initial burn volume, (m^3).
N(1) = (P(1)*V(1))/(R_F*T(1)); %Initial number of moles.
M_G(1) = N(1)*M_N;
Rho_G(1) = N(1)/V(1);

t(1) = 0; %Time starts at zero, (s).
dt = 0.0001; %Time steps, (s).

%Combustion Process.
i = 0;

if Ae == At
 Me = 1; %Choked flow condition
else
 
 
 
 for Me = 0.001:0.00001:20
     eq1 = (Ae/At)^2;
     eq2 = (1/Me^2)*((2/(G_F+1))*(1+((G_F-1)/2)*Me^2))^((G_F+1)/(G_F-1));
 
    if eq1 > eq2
        var1 = eq1/eq2;
    else
        var1 = eq2/eq1;
    end
    
    if var1 > 1.0001
    else
        Mexit = Me;
        break;
    end
 end
 
 Me = Mexit;
end

Te = 2*T_F/(2+Me^2*G_F-Me^2); %exit temperature during combustion, (K).


while Ri(i+1) <= Ro;
 i = i+1;
    
      if i == 1
     V(i) = pi*(Ri(1)^2)*L(1);
      else
     V(i) = pi*(Ri(1)^2)*L(1);
      end
      
        
        P(i) = (N(i)*R_F*T_F)/V(i);
        r(i) = a*P(i)^n;
        Ab(i) = 2*pi*Ri(i)*L(1);

         mdot1(i) = Ab(i)*Rho_F*r(i); %Mass flow rate of fuel, (kg/s).
         mdot3(i) = ((At*P(i))/(T_F))*sqrt(G_F/R_F)*((G_F+1)/2)^(-(G_F+1)/(2*(G_F-1))); %Mass flow rate through nozzle, (kg/s).
         mdot2(i) = mdot1(i) - mdot3(i); %Mass flow rate through chamber, (kg/s).

             N(i) = P(i)*V(i)/(R_F*T_F); %Number of moles during combustion.
             %M_G(i) = N(i)/M_N; %Mass of fuel during combustion, (kg).
             Rho_G(i) = N(i)/V(i); %Density of gas during combustion, (kg/m^3).
 
  
 Pe(i) = P(i)/((1+(1/(2*Me^2*G_F-1)/2)*Me^2)^(G_F/(G_F-1))); %Exit presurre during combustion, (Pa).
 
    
   
 Ve(i)= sqrt((2*G_F*R_F*T_F)/(G_F-1)*(1-(Pe(i)/P(i))^((G_F-1)/G_F)))/12; %Gas exit velocity, (m/s^2).
 F(i) = mdot3(i)*Ve(i)+(Pe(i) - P(1))*Ae; %Thrust produced, (N).
 ISP(i) = F(i)/(mdot3(i)*G0); %Resulant specigic impulse, (s).
 
 if i >= 2
    dP(i)=(P(i)-P(i-1))/dt; %Pressure change, (Pa).
 else
 end
 
 N(i+1) = N(i) + mdot2(i)*dt;
 Ri(i+1) = Ri(i) + r(i)*dt;
 L(i+1) = L(i) - r(i)*dt;
 t(i+1) = t(i) + dt
 
end

t_plot1 = i;

while P(i) >= P(1)
 i = i+1;
 
 Ab(i) = 0;
 V(i) = V(i-1);
 P(i) = (N(i)*R_F*T_F)/(V(i));
 mdot1(i) = 0;
 mdot3(i) = ((At*P(i))/(T_F))*sqrt(G_F/R_F)*((G_F+1)/2)^(-(G_F+1)/(2*(G_F-1))); 
 mdot2(i) = mdot1(i) - mdot3(i);
 
 Rho_G(i) = N(i)/V(i);
 Rho_G(i)
 Pe(i) = P(i)/((1+(1/(2*Me^2*G_F-1)/2)*Me^2)^(G_F/(G_F-1))); 
 Ve(i)= sqrt((2*G_F*R_F*T_F)/(G_F-1)*(1-(Pe(i)/P(i))^((G_F-1)/G_F)))/12;

 F(i) = mdot3(i)*Ve(i)+(Pe(i)-P(1))*Ae; %[lbf] Sea Level
 ISP(i) = (F(i))/(mdot3(i)*G0); %[s]
 N(i+1) = N(i) + mdot2(i)*dt;
 t(i+1) = t(i) + dt;
 
end

[var1 k] = size(r); 
j = 1;
while j <= k
 dr(j) = r(j)*dt;
 j = j + 1;
end
netdr = sum(dr);
clear var1 k

[var1 k] = size(F);
j = 1;
while j <= k
 I(j) = F(j)*dt;
 j = j + 1;
end
Impulse = sum(I);



figure('Name','Combustion Pressure vs. Time');
plot(t(1:i),P);
title('Combustion Pressure vs. Time');
xlabel('T [s]');
ylabel('Combustion Presure [psi]');
figure('Name','Regression Rate vs. Time');
plot(t(1:t_plot1),r(1:t_plot1));
title('Regression Rate vs. Time');
xlabel('T [s]');
ylabel('Regression Rate [in/s]');
figure('Name','Burning Area vs. Time');
plot(t(1:t_plot1),Ab(1:t_plot1));
title('Burning Area vs. Time');
xlabel('T [s]');
ylabel('Burning Area [in^2]');
figure('Name','Combustion Volume vs. Time');
plot(t(1:t_plot1),V(1:t_plot1));
title('Combustion Volume vs. Time');
xlabel('T [s]');
ylabel('Combustion Volume [in^3]');
figure('Name','Mass in Chamber vs. Time');
plot(t(1:t_plot1),N(1:t_plot1));
title('Mass in Chamber vs. Time');
xlabel('T [s]');
ylabel('Mass in Chamber [lbm]');
figure('Name','Nozzle Mass Flow Rate vs. Time');
plot(t(1:i),mdot3);
title('Nozzle Mass Flow Rate vs. Time');
xlabel('T [s]');
ylabel('Nozzle Mass Flow Rate [lbm\s]'); 
figure('Name','Combustion Chamber Volume vs. Time');
plot(t(1:i-1),V(1:i-1));
title('Combustion Chamber Volume vs. Time');
xlabel('T [s]');
ylabel('Combustion Chamber Volume [in^3]');
figure('Name','Gas Density in Chamber vs. Time');
plot(t(1:i-1),Rho_G(1:i-1));
title('Gas Density in Chamber vs. Time');
xlabel('T [s]');
ylabel('Gas Density in Chamber [lbm/in^3]');
figure('Name','Nozzle Exit Pressure vs. Time');
plot(t(1:i),Pe);
title('Nozzle Exit Pressure vs. Time');
xlabel('T [s]');
ylabel('Nozzle Exit Presure [psi]');
ylim([0 1]);
figure('Name','Nozzle Exit Velocity vs. Time');
plot(t(1:i),Ve);
title('Nozzle Exit Velocity vs. Time');
xlabel('T [s]');
ylabel('Nozzle Exit Velocity [ft/s]');
figure('Name','Thrust vs. Time');
plot(t(1:i-1),F(1:i-1));
title('Thrust vs. Time');
xlabel('T [s]');
ylabel('Thrust [lbf]');
figure('Name','Specific Impulse vs. Time');
plot(t(1:i-1),ISP(1:i-1));
title('Specific Impulse vs. Time');
xlabel('T [s]');
ylabel('Specific Impulse [s]');
ylim([0 350]); 
