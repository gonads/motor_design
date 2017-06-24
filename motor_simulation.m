% Written by Nadeem Gabbani - 24/06/2017
clear all; close all; clc;

% Set constants and conditions
G0 = 9.807; %Gravitation acceleration MSL, (m/s^2).
T0 = 293; %Initial Temperature, (K).
P0 = 101325;% Initial Pressure, (Pa).

% Fuel properties
G_F = 1.24; % Fuel specific heat ratio (https://www.nakka-rocketry.net/techs2.html).
M_N = 0.0382; % Fuel molecular weight, (kg/mol).
T_F = 1600; % Fuel steady-state combustion temperature, (K).
R_U = 8.315; % Universal gas constant, (J/mol.K).
R_F = R_U/M_N; % Fuel gas constant, (J/kg.K).
Rho_F = 1800; % Fuel density, (kg/m^3).
a = 7.77; % Fuel burn rate coefficient.
n = -0.013; % Fuel pressure exponent.

% Motor properties
Ro = 0.0342/2; % Fuel core outer radius, (m).
Ri = 0.015/2; % Fuel core inner radius, (m).
L = 0.0705; % Fuel core length, (m).
A = 2*pi*Ri*L; % Fuel burn area, (m^2).
V = (pi*(Ro^2)*L)-(pi*(Ri^2)*L); % Fuel core volume
M_F = Rho_F*V; % Fuel mass, (kg).
Rt = 0.0045/2; % Throat radius, (m).
At = pi*(Rt^2); % Throat area, (m^2).
Re = 0.0045/2; % Exit radius, (m).
Ae = pi*(Re^2); % Exit area, (m^2).
ARact = Ae/At; % Exit area to throat area ratio.

% Set initial conditions for arrays
P(1) = P0; % Chamber pressure to atmospheric, (Pa).
Pt(1) = P0; % Exit pressure to atmospheric, (Pa).
Pe(1) = P0; % Exit pressure to atmospheric, (Pa).
T(1) = T0; % Chamber temperature to atmospheric, (K).
Tt(1) = T0; % Chamber temperature to atmospheric, (K).
Te(1) = T0; % Chamber temperature to atmospheric, (K).
L(1) = L; % Initial length of fuel cure, (m).
Ri(1) = Ri; % Initial burn radius, (m).
V(1) = pi*(Ri(1)^2)*L(1); % Initial burn volume, (m^3).
N(1) = (P(1)*V(1))/(R_F*T(1)); % Initial number of moles.
Ve(1) = 0; % Initial exhaust velocity, (m/s^2).
Rho_G(1) = N(1)/V(1); % Initial gas density, (kg/m^3).

% Set time and loop counter.
t(1) = 0; % Time starts at zero, (s).
dt = 1E-4; % Time steps, (s).
i = 1; % Index

% Identify nozzle behaviour.
if Ae==At
    % Choked flow conditions at the nozzle best describes this situation.
    Me = 1;
else

    % Find the relationship between Mach number and area ratio.
    fMe = 0.001:0.001:10; 
    ARopt = (1./fMe).*(((1+(G_F-1./2).*fMe.^2)./(1+(G_F-1/2))).^((G_F+1/(2.*(G_F-1)))));
    AM = vertcat(fMe, ARopt);      
    %{ Find the Mach number the selected area ratio satisfies. Remember that Pe
       needs to match P0 for optimum performance, use Pratio to help. %}
    Me = interp1(fMe, ARopt, ARact);
end

% Integrate from start of combustion to when chamber pressure back to external
while true
    i = i + 1;
    t(i) = t(i-1) + dt;

    % Update chamber conditions and fuel properties
    N(i) = N(i - 1) + (mdot2(i - 1)*dt); % Number of moles

    % Combustion phase
    if Ri(i + 1) <= Ro
        Ri(i) = Ri(i-1) + (Rb(i-1)*dt);
        L(i) = L(i-1) - (Rb(i-1)*dt);
        V(i) = pi * (Ri(i)^2) * L(i); % Dynamic chamber volume.
        P(i) = (N(i)*R_F*T_F) / V(i); % Dynamic chamber pressure.
        Rb(i) = (a/1000) * (P(i))^n; % Dynamic burn rate.
        Ab(i) = 2 * pi * Ri(i) * L(i); % Dynamic burn area.
        mdot1(i) = Ab(i) * Rho_F * Rb(i); % Mass flow rate through chamber.

    % Tail phase
    else
        V(i) = V(i-1); % So is this just 0 now?
        P(i) = (N(i)*R_F*T_F) / V(i); % Chamber Pressure, (Pa).
        Rb(i) = 0;
        Ab(i) = 0;
        mdot1(i) = 0; % Mass flow rate of fuel, (kg/s).
    end
    Rho_G(i) = mdot1(i)/V(i); % Density of gas during combustion, (kg/m^3).          

    % Update nozzle conditions
    Pt(i) = P(i)*(1+(G_F-1/2)^(-G_F/(G_F-1))); % Dynamic throat pressure.
    Pe(i) = Pt(i)*((1+(G_F-1/2)*Me^2).^(-G_F/(G_F-1))); % Dynamic exit pressure.
    mdot3(i) = ((At*P(i))/(sqrt(T_F)))*sqrt(G_F/R_F)*((G_F+1)/2)^...
            (-(G_F+1)/(2*(G_F-1))); % Mass flow rate through nozzle.
    mdot2(i) = mdot1(i) - mdot3(i); % Mass flow rate causing back pressure.

    % Update temperatures
    T(i) = (P(i)*V(i))/(N(i)*R_F); % Combustion temperature.
    Tt(i) = T(i)/(1+(G_F-1/2)); % Throat temperature.
    Te(i) = Tt(i)*((1+(G_F-1/2)*Me^2)^-1); % Exit temperature.

    % Update exhaust conditions and properties
    Ve(i) = sqrt(((2*G_F)/(G_F-1))*R_F*T(i)*((1-((Pe(i)/P(i))^(G_F-1/G_F))))); % Exhaust velocity.
    F(i) = (mdot3(i)*Ve(i))+((Pe(i)-P0)*Ae); % Dynamic thrust produced.
    ISP(i) = F(i)/(mdot3(i)*G0); % Dynamic specific impulse,.
    Cs (i) = (At*P(i))/mdot3(i); % Characteristic exhaust velocity.

    % Check if pressure returned to initial conditions
    if P(i) >= P(1):
        break
    end
end

% Calculate average ratio of initial pressure to exit pressure
PR = mean(P0./Pe);
fprintf('Pratio = %f', PR);
fprintf('\n');

% Only plot if not choked flow
if Me ~= 1
    title('Mach Number vs. Exit/Throat Area Ratio');
    xlabel('Mach Number');
    ylabel('Area Ratio');
    set(gca,'TickDir','out');
    grid minor;
    plot(fMe, ARopt);
end
    
% Plotting
figure;
title('Pressures vs. Time');
xlabel('Time (s)');
ylabel('Pressure (Pa)');
set(gca,'TickDir','out');
grid minor;
legend('P0','PC', 'PT', 'PE');
plot(t, repelem(P0, length(t)), t), P, t), Pt, t, Pe);

figure;
title('Temperatures vs. Time');
xlabel('Time (s)');
ylabel('Temperature (K)');
set(gca,'TickDir','out');
grid minor;
legend('T0','TC', 'TT', 'TE');
plot(t, repelem(T0, length(t)), t, T, t, Tt, t, Te); 

figure;
title('Mass-Flow Rates vs. Time');
xlabel('Time (s)');
ylabel('Mass Flow Rate (kg/s)');
set(gca,'TickDir','out');
grid minor;
legend('MdotC','MdotE', 'MdotT');
plot(t, mdot1, t, abs(mdot2), t, mdot3);

figure;
title('Gas Density vs. Time');
xlabel('Time (s)');
ylabel('Gas Density (kg/m^3)');
set(gca,'TickDir','out');
grid minor;
plot(t, Rho_G);

figure;
plot(t, Rb);
title('Burn Rate vs. Time');
xlabel('Time (s)');
ylabel('Burn Rate (m/s)');
set(gca,'TickDir','out');
grid minor;

figure;
title('Chamber Volume vs. Time');
xlabel('Time (s)');
ylabel('Volume (m^3)');
set(gca,'TickDir','out');
grid minor;
plot(t, V);

figure;
title('Chamber Radius vs. Time');
xlabel('Time (s)');
ylabel('Radius (m)');
set(gca,'TickDir','out');
grid minor;
plot(t(1:length(Ri)), Ri);

figure;
title('Thrust vs. Time');
xlabel('Time (s)');
ylabel('Thrust (N)');
set(gca,'TickDir','out');
grid minor;
plot(t, F);

figure;
title('Specific Impukse vs. Time');
xlabel('Time (s)');
ylabel('Specific Impulse (s)');
set(gca,'TickDir','out');
grid minor;
plot(t, ISP);

figure;
title('Velocities vs. Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
set(gca,'TickDir','out');
grid minor;
legend('Ve','C*');
plot(t, Ve, t, Cs);
