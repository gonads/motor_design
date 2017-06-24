%Written by Nadeem Gabbani - 24/06/2017
clear all; close all; clc;

%Set constants and conditions.
G0 = 9.807; %Gravitation acceleration MSL, (m/s^2).
T0 = 293; %Initial Temperature, (K).
P0 = 101325;% Initial Pressure, (Pa).

%Set fuel properties
G_F = 1.24; %Fuel specific heat ratio (https://www.nakka-rocketry.net/techs2.html).
M_N = 0.0382; %Fuel molecular weight, (kg/mol).
T_F = 1600; %Fuel steady-state combustion temperature, (K).
R_U = 8.315; %Universal gas constant, (J/mol.K).
R_F = R_U/M_N; %Fuel gas constant, (J/kg.K).
Rho_F = 1800; %Fuel density, (kg/m^3).
a = 7.77; %Fuel burn rate coefficient.
n = -0.013; %Fuel pressure exponent.

%Set motor properties.
Ro = 0.0342/2; %Fuel core outer radius, (m).
Ri = 0.015/2; %Fuel core inner radius, (m).
L = 0.0705; %Fuel core length, (m).
A = 2*pi*Ri*L; %Fuel burn area, (m^2).
V = (pi*(Ro^2)*L)-(pi*(Ri^2)*L); %Fuel core volume
M_F = Rho_F*V; %Fuel mass, (kg).
Rt = 0.0045/2; %Throat radius, (m).
At = pi*(Rt^2); %Throat area, (m^2).
Re = 0.0045/2; %Exit radius, (m).
Ae = pi*(Re^2); %Exit area, (m^2).
ARact = Ae/At; %Exit area to throat area ratio.

%Set initial conditions.
P(1) = P0; %Chamber pressure to atmospheric, (Pa).
Pt(1) = P0; %Exit pressure to atmospheric, (Pa).
Pe(1) = P0; %Exit pressure to atmospheric, (Pa).
T(1) = T0; %Chamber temperature to atmospheric, (K).
Tt(1) = T0; %Chamber temperature to atmospheric, (K).
Te(1) = T0; %Chamber temperature to atmospheric, (K).
L(1) = L; %Initial length of fuel cure, (m).
Ri(1) = Ri; %Initial burn radius, (m).
V(1) = pi*(Ri(1)^2)*L(1); %Initial burn volume, (m^3).
N(1) = (P(1)*V(1))/(R_F*T(1)); %Initial number of moles.
Ve(1) = 0; %Initial exhaust velocity, (m/s^2).
Rho_G(1) = N(1)/V(1); %Initial gas density, (kg/m^3).

%Set time and loop counter.
t(1) = 0; %Time starts at zero, (s).
dt = 1E-4; %Time steps, (s).
i = 0; %Start at step zero.

%Identify nozzle behaviour.
if Ae==At
 Me = 1; %Choked flow conditions at the nozzle best describes this situation.
    else
     fMe = 0.001:0.001:10; 
     ARopt = (1./fMe).*(((1+(G_F-1./2).*fMe.^2)./(1+(G_F-1/2))).^((G_F+1/(2.*(G_F-1))))); %Find the relationship between Mach number and area ratio.
     AM = vertcat(fMe, ARopt);      
     Me = interp1(fMe, ARopt, ARact); %Find the Mach number the selected area ratio satisfies. Remember that Pe needs to match P0 for optimum performance, use Pratio to help.
end
      
while Ri(i+1) <= Ro; %Combustion can occur if fuel is remaining.      
 i = i+1;
     
     V(i) = pi*(Ri(i)^2)*L(i); %Dynamic chamber volume.
     P(i) = (N(i)*R_F*T_F)/V(i); %Dynamic chamber pressure.
     Pt(i) = P(i)*(1+(G_F-1/2)^(-G_F/(G_F-1))); %Dynamic throat pressure.
     Pe(i) = Pt(i)*((1+(G_F-1/2)*Me^2).^(-G_F/(G_F-1))); %Dynamic exit pressure.
     
         Rb(i) = (a/1000)*(P(i))^n; %Dynamic burn rate.
         Ab(i) = 2*pi*Ri(i)*L(i); %Dynamic burn area.
     
             mdot1(i) = Ab(i)*Rho_F*Rb(i); %Mass flow rate through chamber.
             mdot3(i) = ((At*P(i))/(sqrt(T_F)))*sqrt(G_F/R_F)*((G_F+1)/2)^(-(G_F+1)/(2*(G_F-1))); %Mass flow rate through nozzle.
             mdot2(i) = mdot1(i) - mdot3(i); %Mass flow rate causing back pressure.
    
                 T(i) = (P(i)*V(i))/(N(i)*R_F); %Combustion temperature.
                 Tt(i) = T(i)/(1+(G_F-1/2)); %Throat temperature.
                 Te(i) = Tt(i)*((1+(G_F-1/2)*Me^2)^-1); %Exit temperature.
     
                    Rho_G(i) = mdot1(i)/V(i); %Density of gas during combustion, (kg/m^3).          
   
                         Ve(i) = sqrt(((2*G_F)/(G_F-1))*R_F*T(i)*((1-((Pe(i)/P(i))^(G_F-1/G_F))))); %Exhaust velocity.
                         F(i) = (mdot3(i)*Ve(i))+((Pe(i)-P0)*Ae); %Dynamic thrust produced.
                         ISP(i) = F(i)/(mdot3(i)*G0); %Dynamic specific impulse,.
                         Cs (i) = (At*P(i))/mdot3(i); %Characteristic exhaust velocity.
 
    if i >= 2
         dP(i)=(P(i)-P(i-1))/dt; %Pressure change.
        else
    end
 
             N(i+1) = N(i)+(mdot2(i)*dt);
             Ri(i+1) = Ri(i)+(Rb(i)*dt);
             L(i+1) = L(i)-(Rb(i)*dt);
             t(i+1) = t(i)+dt;
 
end


while P(i) >= P(1) %Tail down occurs when fuel runs out.
 i = i+1;
     
     Ab(i) = 0;
     Rb(i) = 0;
     V(i) = V(i-1);
     P(i) = (N(i)*R_F*T_F)/V(i); %Chamber Pressure, (Pa).
     Pt(i) = P(i)*(1+(G_F-1/2)^(-G_F/(G_F-1)));
     Pe(i) = Pt(i)*((1+(G_F-1/2)*Me^2)^(-G_F/(G_F-1))); %Exit pressure, (Pa).
     
         mdot1(i) = 0; %Mass flow rate of fuel, (kg/s).
         mdot3(i) = ((At*P(i))/(sqrt(T_F)))*sqrt(G_F/R_F)*((G_F+1)/2)^(-(G_F+1)/(2*(G_F-1))); %Mass flow rate through nozzle, (kg/s).
         mdot2(i) = mdot1(i) - mdot3(i); %Mass flow rate through chamber, (kg/s).
    
             T(i) = (P(i)*V(i))/(N(i)*R_F);
             Tt(i) = T(i)/(1+(G_F-1/2));
             Te(i) = Tt(i)*((1+(G_F-1/2)*Me^2)^-1);
     
                Rho_G(i) = mdot1(i)/V(i); %Density of gas during combustion, (kg/m^3).          
   
                     Ve(i) = sqrt(((2*G_F)/(G_F-1))*R_F*T(i)*((1-((Pe(i)/P(i))^(G_F-1/G_F)))));
                     F(i) = (mdot3(i)*Ve(i))+((Pe(i)-P0)*Ae); %Thrust produced, (N).
                     ISP(i) = F(i)/(mdot3(i)*G0); %Resulant specific impulse, (s).
                     Cs (i) = (At*P(i))/mdot3(i);

                        N(i+1) = N(i)+(mdot2(i)*dt);
                        t(i+1) = t(i)+dt;
  
end

if Me ~= 1
    plot(fMe, ARopt);
    title('Mach Number vs. Exit/Throat Area Ratio'); xlabel('Mach Number'); ylabel('Area Ratio'); set(gca,'TickDir','out'); grid minor;
end
    
figure; plot(t(1:end-1), repelem(P0, length(t(1:end-1))), t(1:end-1), P, t(1:end-1), Pt, t(1:end-1), Pe);
title('Pressures vs. Time'); xlabel('Time (s)'); ylabel('Pressure (Pa)'); set(gca,'TickDir','out'); grid minor; legend('P0','PC', 'PT', 'PE');

figure; plot(t(1:end-1), repelem(T0, length(t(1:end-1))), t(1:end-1), T, t(1:end-1), Tt, t(1:end-1), Te); 
title('Temperatures vs. Time'); xlabel('Time (s)'); ylabel('Temperature (K)'); set(gca,'TickDir','out'); grid minor; legend('T0','TC', 'TT', 'TE');

figure; plot(t(1:end-1), mdot1, t(1:end-1), abs(mdot2), t(1:end-1), mdot3);
title('Mass-Flow Rates vs. Time'); xlabel('Time (s)'); ylabel('Mass Flow Rate (kg/s)'); set(gca,'TickDir','out'); grid minor; legend('MdotC','MdotE', 'MdotT');

figure; plot(t(1:end-1), Rho_G);
title('Gas Density vs. Time'); xlabel('Time (s)'); ylabel('Gas Density (kg/m^3)'); set(gca,'TickDir','out'); grid minor;

figure; plot(t(1:end-1), Rb);
title('Burn Rate vs. Time'); xlabel('Time (s)'); ylabel('Burn Rate (m/s)'); set(gca,'TickDir','out'); grid minor;

figure; plot(t(1:end-1), V);
title('Chamber Volume vs. Time'); xlabel('Time (s)'); ylabel('Volume (m^3)'); set(gca,'TickDir','out'); grid minor;

figure; plot(t(1:end-(length(t(1:end))-length(Ri))), Ri);
title('Chamber Radius vs. Time'); xlabel('Time (s)'); ylabel('Radius (m)'); set(gca,'TickDir','out'); grid minor;

figure; plot(t(1:end-1), F);
title('Thrust vs. Time'); xlabel('Time (s)'); ylabel('Thrust (N)'); set(gca,'TickDir','out'); grid minor;

figure; plot(t(1:end-1), ISP);
title('Specific Impukse vs. Time'); xlabel('Time (s)'); ylabel('Specific Impulse (s)'); set(gca,'TickDir','out'); grid minor;

figure; plot(t(1:end-1), Ve, t(1:end-1), Cs);
title('Velocities vs. Time'); xlabel('Time (s)'); ylabel('Velocity (m/s)'); set(gca,'TickDir','out'); grid minor; legend('Ve','C*');

PR = mean(P0./Pe); fprintf('Pratio = %f', PR); fprintf('\n');






