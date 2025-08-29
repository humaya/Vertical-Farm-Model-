%%
% Copyright (c) 2025 Homa Esmaeli, Wageningen University Foundation All rights reserved
% Matlab code for Model developed for Vertical Farming systems expained at 
% paper: {A Dynamic Control-Oriented Model of Climate and Crop Growth in Vertical FarmsA Dynamic Control-Oriented Model of Climate and Crop Growth in Vertical Farms}

clear all;
close all; clc;
addpath(genpath('bin'));

% model parameters 
p = DefineParametersShare;  % {adjust for your case, only relevant parameters}

ops.nx  = 4;    % number of states
ops.ny  = 5;    % number of measurements
ops.nu  = 7;    % number of controllable inputs
ops.nd  = 0;    % number of disturbance inputs

% simulation parameters
c         = 86400;              % One day om seconds (24*60*60)
nDays     = 21;                 % days in simulation (max = 46 weeks) Time: from Transplant to Harvest {adjust for your case} 
h_int     = 1*3;                % time intervals in calculations  [min]                                 {adjust for your case} 
ops.h     = h_int*60;           % sample period in seconds: #minutes*60 (so if you want h (time step) to be 1 hr (60 min), then ops.h= 60*60 means each time step is 3600 s or an hour, for h=15 min then ops.h=15*60, and so on)
ops.L     = nDays*c;            % final time simulation (second)
ops.t     = 0:ops.h:ops.L;      % initial time vector
ops.N     = length(ops.t);      % number of samples in initial time vector
ops.Np    = 0;                  % simulation horizen
ops.U     = 3;                  % number of resources 

% signals
d          = zeros(ops.nd,ops.N);     % d = LoadDisturbances(ops);  
x          = zeros(ops.nx,ops.N+1);
y          = zeros(ops.ny,ops.N);
U          = zeros(ops.U,ops.N);
R          = zeros(ops.U,1);
x(:,1)     = [0.0035 2.34e-3 24 0.015];  % initial condition [dry_weight CO2 temp AH]  {adjust for your case} 

%% define climate set points here  {adjust for your case} 
PPFD         = 200;             % set point             [umol.m-2.s-1] 
photo_period = 16;              % set point             [hr]
air_speed    = 0.30;            % set point             [m.s-1]
SP_CO2_day     = 2.16e-3;       % set point 1200 ppm    [kg.m-3]  
SP_CO2_night   = 0* 2.16e-3;       % set point 1200 ppm    [kg.m-3]  
SP_CO2ppm_day     = 1200;       % set point
SP_CO2ppm_night   = 0*1200;       % set point 
SP_t_day   = 24;                % set point           [C]
SP_t_night = 20;                % set point           [C]
SP_h_day   = 0.016;             % RH 75%  @ 24 C                  [kg.m-3]
SP_h_night = 0.013;             % RH 75%  @ 20 C                  [kg.m-3]
SP_rh_day   = 75;               % RH 75%  @ 24 C                  [%]
SP_rh_night = 75;               % RH 75%  @ 20 C                  [%]
xi_red     = 0.88;                      % red fraction in total spectrum
xi_blue    = 0.12;                      % blue fraction in total spectrum
xi_green   = 0;                         % green fraction in total spectrum    
xi_FR      = 0;                         % FR fraction in total spectrum 

CPUTime    = zeros(1,ops.N);

% Initialization u(1)
error_u1      =  zeros(1, ops.N);            %  Initialization of error
integral_u1   =  zeros(1, ops.N);            %  Initialization of integral
derivative_u1 =  zeros(1, ops.N);            %  Initialization of derivative
output_u1     =  zeros(1, ops.N);            %  Initialization of output

% Initialization u(4)
error_u4      =  zeros(1, ops.N);            %  Initialization of error
integral_u4   =  zeros(1, ops.N);            %  Initialization of integral
derivative_u4 =  zeros(1, ops.N);            %  Initialization of derivative
output_u4     =  zeros(1, ops.N);            %  Initialization of output

% Initialization u(5)
error_u5      =  zeros(1, ops.N);            %  Initialization of error
integral_u5   =  zeros(1, ops.N);            %  Initialization of integral
derivative_u5 =  zeros(1, ops.N);            %  Initialization of derivative
output_u5     =  zeros(1, ops.N);            %  Initialization of output
u = zeros(5,ops.N);

% time-loop
disp('Time loop')
for kk=1:ops.N

    % this section contrllers are defined and tuned
    u(6,kk) = photo_period ; 
    u(7,kk) = ((xi_red*p.landa_red)+(xi_blue*p.landa_blue)+(xi_green*p.landa_green))+(xi_FR*p.landa_FR) ;
    p.xi_lampPAR = u(7,kk)/(p.h_plank*p.c_speedlight*6.023e+23*1e-6);     % if we consider the light spectrum portions
    
    % on-off signal PPFD 
    u(2,kk)    = PPFD*(mod(kk,24*(60/h_int))<u(6,kk)*(60/h_int) ) + 0*(mod(kk,24*(60/h_int) )>=u(6,kk)*(60/h_int) );  %on-off signal for light (PPFD)
    
    % air speed signal (p controller)
    Kp_u3      =  1; 
    u(3,:)     =  Kp_u3  * air_speed; 
   
    u(1,1)     = 0;  
    % PID Controller Parameters for u(1)
    Kp_u1 = 15*10;                     % Proportional gain
    Ki_u1 = 1*0.01;                    % Integral gain
    Kd_u1 = 0.01*0.1;                  % Derivative gain
    SP_CO2 = SP_CO2_day*(mod(kk,24*(60/h_int) )<photo_period*(60/h_int) ) + SP_CO2_night*(mod(kk,24*(60/h_int) )>=photo_period*(60/h_int)); % when we have 2 set points for CO2 (day and night)
    error_u1(kk+1)       =  SP_CO2 - x(2,kk);                                                                     % Error calculation
    integral_u1(kk+1)    =  integral_u1(kk) + error_u1(kk+1) * ops.h;                                             % Integral term
    derivative_u1(kk+1)  =  (error_u1(kk+1) - error_u1(kk)) / ops.h;                                              % Derivative term
    output_u1(kk+1)      =  Kp_u1 * error_u1(kk+1) + Ki_u1 * integral_u1(kk+1) + Kd_u1 * derivative_u1(kk+1);     % PID output
    u(1,kk+1)            =  max(0,((output_u1(kk+1))*p.volum)/(p.A_cultivation*3600)) ;   % use when we have 2 set points for CO2 (day and night)

    SP_t  = SP_t_day*(mod(kk,24*(60/h_int) )<photo_period*(60/h_int) ) + SP_t_night*(mod(kk,24*(60/h_int) )>=photo_period*(60/h_int));  % tempeature set point day and night
    SP_h  = SP_h_day*(mod(kk,24*(60/h_int) )<photo_period*(60/h_int) ) + SP_h_night*(mod(kk,24*(60/h_int) )>=photo_period*(60/h_int));  % humidity set point day and night

    % PID Controller Parameters for u(4) temp of AC  Ziegler-Nichols method
    k_uzn = 0.5;
    P_uzn = 0.152*60*60;
    T_i   = P_uzn/2; 
    T_d   = P_uzn/8;
    Kp_u4 = 0.6*k_uzn;           % Proportional gain 
    Ki_u4 = Kp_u4/T_i ;          % Integral gain
    Kd_u4 = Kp_u4*T_d;           % Derivative gain

    error_u4(kk+1)       =  SP_t - x(3,kk); % output_u4(kk);                                                % Error calculation
    integral_u4(kk+1)    =  integral_u4(kk) + error_u4(kk+1) * ops.h;                                       % Integral term
    derivative_u4(kk+1)  =  (error_u4(kk+1) - error_u4(kk)) / ops.h;                                        % Derivative term
    output_u4(kk+1)      =  Kp_u4 * error_u4(kk+1) + Ki_u4 * integral_u4(kk+1) + Kd_u4 * derivative_u4(kk+1);     % PID output
    u(4,1)     = 20; 
    u(4,kk+1) = min (25, max(18, output_u4(kk+1)));       % pid-controller  signal

    % PID Controller Parameters for u(5) humidity of AC Ziegler-Nichols method
    k_uznh = 2;
    P_uznh = 0.15*60*60;
    T_ih   = P_uznh/2; 
    T_dh   = P_uznh/8;
    Kp_u5 =  0.6*k_uznh;           % Proportional gain 
    Ki_u5 =  Kp_u5/T_ih ;          % Integral gain
    Kd_u5 =  Kp_u5*T_dh;           % Derivative gain
 
    error_u5(kk+1)       =  SP_h - x(4,kk); % output_u5(kk);                                                      % Error calculation
    integral_u5(kk+1)    =  integral_u5(kk) + error_u5(kk+1) * ops.h;                                             % Integral term
    derivative_u5(kk+1)  =  (error_u5(kk+1) - error_u5(kk)) / ops.h;                                              % Derivative term
    output_u5(kk+1)      =  Kp_u5 * error_u5(kk+1) + Ki_u5 * integral_u5(kk+1) + Kd_u5 * derivative_u5(kk+1);     % PID output
    
    u(5,kk+1) =  max (0, output_u5(kk+1)); 
    tic

    % propagate model one time step ahead
    x(:,kk+1)  = f(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    y(:,kk)    = g(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
    U(:,kk)    = e(x(:,kk), u(:,kk), d(:,kk), p, ops.h);
 
    CPUTime(kk)  = toc;

end

disp(' ')
disp(['mean CPU time: ',num2str(mean(CPUTime)),' (sec).'])
disp(' ')

%% plotting 
UserColor    = 1/255*[0,0,0];

figure(1);clf;
ops.y2_min_k(1:kk) = 0; 
ops.y2_max_k(1:kk) = 1200; 
ops.y3_min_k(1:kk)  = 20; 
ops.y3_max_k(1:kk)  = 24; 
ops.y4_min_k(1:kk)  = 75; 
ops.y4_max_k(1:kk)  = 75; 
ops.u2_min_k(1:kk)  = 0;
ops.u2_max_k(1:kk)  = 200;
ops.u3_min_k(1:kk)  = 0.2;
ops.u3_max_k(1:kk)  = 0.2;

for ll=1:ops.ny
    subplot(3,ops.ny,ll)
    stairs(ops.t(1:kk)/c,y(ll,1:kk),'linewidth',1.5);grid;hold on
    axis tight
    if  ll==1
        ylabel('$y_1$ DW (g/m$^2$)','fontsize',12,'interpreter','latex')
    elseif ll==2
        stairs(ops.t(1:kk)/c,ops.y2_min_k(1:kk),'--','color',UserColor,'linewidth',1.5);
        stairs(ops.t(1:kk)/c,ops.y2_max_k(1:kk),'r--','linewidth',1.5);        
        ylabel('$y_2$ CO2 (ppm)','fontsize',12,'interpreter','latex')
    elseif ll==3
        stairs(ops.t(1:kk)/c,ops.y3_min_k(1:kk),'--','color',UserColor,'linewidth',1.5);
        stairs(ops.t(1:kk)/c,ops.y3_max_k(1:kk),'r--','linewidth',1.5);        
        ylabel('$y_3$ Temp ($^\circ$C)','fontsize',12,'interpreter','latex')
    elseif ll==4
        stairs(ops.t(1:kk)/c,ops.y4_min_k(1:kk),'--','color',UserColor,'linewidth',1.5);
        stairs(ops.t(1:kk)/c,ops.y4_max_k(1:kk),'r--','linewidth',1.5);        
        ylabel('$y_4$ RH (\%)','fontsize',12,'interpreter','latex')
        ylim([0 100])
    elseif ll==5
        ylabel('$y_5$ FW (kg/m$^2$)','fontsize',12,'interpreter','latex')      
    end
end


for ll=1:ops.ny
    subplot(3,ops.ny,ops.ny+ll)
    stairs(ops.t(1:kk)/c,u(ll,1:kk),'linewidth',1.5);grid;hold on
    axis tight
    if ll==1      
        ylabel('$u_1$ (kg/m$^2$.s)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    elseif ll==2
        stairs(ops.t(1:kk)/c,ops.u2_min_k(1:kk),'--','color',UserColor,'linewidth',1.5);
        stairs(ops.t(1:kk)/c,ops.u2_max_k(1:kk),'r--','linewidth',1.5);  
        ylabel('$u_2$ (micromol/m$^2$.s)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    elseif ll==3
        stairs(ops.t(1:kk)/c,ops.u3_min_k(1:kk),'--','color',UserColor,'linewidth',1.5);
        stairs(ops.t(1:kk)/c,ops.u3_max_k(1:kk),'r--','linewidth',1.5);  
        ylabel('$u_3$ (m/s)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    elseif ll==4
        ylabel('$u_4$ ($^\circ$C)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    elseif ll==5
        ylabel('$u_5$ (kg/m$^2$)','fontsize',13,'interpreter','latex')
        xlabel('Time (days)','fontsize',13,'interpreter','latex')
    end
end


% Calculation of Resource Use 
R(1)    = sum (U(1,:));                  % electricty consumption by LED   [W.m-2]
R_1kwh  = R(1)*(ops.h/3600)*1.0e-3;      % electricty consumption by LED   [sum of *(W.m-2)*(time interval in hours)/1000 = kWh.m-2]
R(2)    = sum (U(2,:));                  % electricty consumption by AC    [W.m-2]
R_2kwh  = R(2)*(ops.h/3600)*1.0e-3;      % electricty consumption by AC    [kWh.m-2]
R(3)    = sum (U(3,:));                  % CO2 consumption                 [kg.m-2] 
R_3ppm  = 24.5*(R(3)*10^6)/44.01;        % CO2 consumption [ppm.m-2]

Resorces = [R(1) " LED [W.m-2]", R_1kwh "LED [kWh.m-2]"; R(2)  " AC [W.m-2]", R_2kwh   "AC [kWh.m-2]"; R(3)   "CO2 [kg.m-2]",  R_3ppm   "CO2 [ppm.m-2]"];

R_eff = {
    R(1)/y(5,end), 'LED [W.kg fw-1]', R_1kwh/y(5,end), 'LED [kWh.kg fw-1]';
    R(2)/y(5,end), 'AC [W.kg fw-1]',  R_2kwh/y(5,end), 'AC [kWh.kg fw-1]';
    R(3)/y(5,end), 'CO2 [kg.kg fw-1]', R_3ppm/y(5,end), 'CO2 [ppm.kg fw-1]'
};

Resorces
R_eff

% solver RK4th order 
function dx = f(x,u,d,p,h)

k1  = F(x,u,d,p,h); 
k2  = F(x + h/2 * k1,u,d,p,h);
k3  = F(x + h/2 * k2,u,d,p,h);
k4  = F(x + h * k3,u,d,p,h);
x   = x + h/6*(k1 + 2*k2 + 2*k3 + k4);

dx  = x;

end

% model equations
function ki = F(x,u,d,p,h)       

v_dot        = ((u(3)*p.A_duct));                                     % volumetric flow rate of the air
Gamma        =  (p.gamma*(p.Q10_gamma^((0.1*x(3))-2)));
epsilon      =  (p.epsi*((x(2)-Gamma)/(x(2)+2*Gamma)));
sigma_CO2    =  (-p.photCO2_1*x(3).^2+p.photCO2_2*x(3)-p.photCO2_3);
phi_phot     =  (1-exp(-p.k*p.lar_d*10^3*(1-p.tao)*x(1)));
phi_phot_max =  (epsilon.*(u(2)/p.xi_lampPAR).*sigma_CO2.*(x(2)-Gamma))./(epsilon.*(u(2)/p.xi_lampPAR) +sigma_CO2.*(x(2)-Gamma));
h_sat=((p.v1/(p.R*(x(3)+p.T)))*exp((p.v2*x(3))/(x(3)+p.v3)));
phi_transp   =  ((1-exp(-p.pl_d*x(1)))*p.mass*(h_sat-x(4)));
phi_ac       =  v_dot*(x(4)-u(5))/p.A_cultivation;
phi_resp     =  ((p.resp_s*(1-p.tao))+(p.resp_r*p.tao))*x(1).*p.Q10_resp.^(0.1*x(3)-2.5);
q_ac_sens    =  v_dot*p.ro*p.c_p*(x(3)-u(4))/p.A_cultivation;
q_ac_latent  =  p.L*phi_ac;
q_ac         =  q_ac_sens + q_ac_latent;
q_trans      =  p.L*phi_transp;
q_led        =  ((u(2)/p.PPE)*p.eta_conv)+((u(2)/p.PPE)*p.eta_NIR)+((u(2)/p.PPE)*p.eta_PAR*abs(1-p.k*(epsilon*p.delta_G/p.MW_CO2)));
phi_phot_total = phi_phot*phi_phot_max;

ki =  [ p.beta * (p.alfa*( phi_phot_total - phi_resp));               % crop growth model, x_1 is state varible X_d : crop dry weight    [kg. m-2.s-1]

1/p.CO2cap * (u(1) - phi_phot_total + phi_resp );                     % CO2 balance,  x_2 is state varible X_c : CO2 concentration       [kg. m-3.s-1]

1/p.aCap * ( -q_ac - q_trans + q_led);                                % energy balance,  x_3 is state varible X_t : VF air temperature   [C.s-1]

1/p.H2Ocap * (phi_transp - phi_ac)];                                  % humidity balance,  x_4 is state varible X_h : absolute humidity  [kg. m-3.s-1]


end


function y = g(x,u,d,p,h)         % unit change for state variables + fresh weight  

y = [1e3*x(1);                    % g/m2 dry weight
     co2dens2ppm(x(3),x(2));      % ppm  CO2
     x(3);                        % temp C
     vaporDens2rh(x(3), x(4));    % RH % 
     p.fw *x(1)*(1-p.tao)];       % kg/m2 fw
end



function U = e (x,u,d,p,h)  % calcuation of resource use 
 
 U(1)         = u(2)/p.PPE;                                    % electricty consumption by LED at each time step  [W.m-2]
 v_dot        =  (u(3)*p.A_duct);
 q_ac_sens    =  v_dot*p.ro*p.c_p*(x(3)-u(4))/p.A_cultivation;
 phi_ac       =  v_dot*(x(4)-u(5))/p.A_cultivation;
 q_ac_latent  =  p.L*phi_ac;
 q_ac         =  q_ac_sens + q_ac_latent;
 U(2)         = abs(q_ac)/p.cop;                              % electricty consumption by AC  at each time step  [W.m-2]
 U(3) = u(1)*h ;                                              % CO2 consumption at each time step [kg.m-2.s-1 * sec (of each step size)]

end  
