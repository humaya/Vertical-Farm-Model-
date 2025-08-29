function p = DefineParametersShare % parameters 

%% physical properties (non changable) 
% parameter 					description 									[unit] 					nominal value
p.satH2O1 = 9348; 				% saturation water vapour parameter 			[J m^{-3}] 				9348
p.satH2O2 = 17.4; 				% saturation water vapour parameter 			[-] 					17.4
p.satH2O3 = 239; 				% saturation water vapour parameter 			[°C] 					239
p.satH2O4 = 10998;  			% saturation water vapour parameter 			[J m^{-3}] 				10998
p.R = 8314; 					% ideal gas constant 							[J K^{-1} kmol^{-1}] 	8314
p.T = 273.15; 					% conversion from C to K 						[K] 					273.15
p.ro=1.225;                     % air density                                   [kg m^{-3}]
p.c_p=1000;                     % air heat capacity                             [J kg^{-1} K^{-1}]
p.L=2.25e6;                     % latent heat of water vaporization             [J kg^{-1}]
p.CO2cap = 4.1; 				% CO2 capacity of the greenhouse 				[m^3{air} m^{-2}{gh}]   4.1
p.H2Ocap = 4.1; 				% Vapor capacity of the greenhouse 				[m^3{air} m^{-2}{gh}]   4.1
p.aCap = 3e4; 					% effective heat capacity of the greenhouse air [J m^{-2}{gh} °C^{-1}]  3e4
p.ventCap = 1290; 				% heat capacity per volume of greenhouse air 	[J m^{-3}{gh} °C^{-1}]  1290

% spectrum of light
p.h_plank = 6.62607015e-34;                   % Planck constant            [m2.kg.s-1]
p.c_speedlight = 299792458;                   % speed of light             [m.s-1]
p.avag = 6.023e+23;                           % Avogadro constant          [photon.mol-1]
p.landa_red = 685e-9;                         % wavelength of red photons  [m]  620-750 nm 
p.landa_blue = 472.5e-9;                      % wavelength of blue photons  [m]  450-495 nm
p.landa_green = 532.5e-9;                     % wavelength of green photons  [m]  495-570 nm 
p.landa_FR    = 750e-9;                       % wavelength of Far Red  photons  [m]  700-800 nm 

%% crop growth model parametres 
% parameter 					description 									[unit] 					nominal value
p.alfaBeta = 0.544; 			% yield factor 	(replaced by p.alfa and p.beta)								[-] 					0.544
p.alfa=0.68;                    % Physical constant                             [-] 
p.beta=0.8;                     % yield factor     0-1                              [-] 
p.gamma=7.32e-5;                % CO2 compensation point at 20?                 [kg m^{-2}]
p.Q10_gamma=2;                  % Q10 factor, for every temperature increase of 10?, the CO2 compensation point increase by a factor cQ10,?      [-]
p.epsi=17e-9;                   % Light use efficiency at very high CO2 concentration in the absence of photorespiration             [kg J^{-1}]          
p.resp_s=3.47e-7;               % Maintenance respiration rates for shoot at 25? expressed in the mass of glucose consumed           [s^{-1}]
p.resp_r=1.16e-7;               % Maintenance respiration rates for root at 25? expressed in the mass of glucose consumed            [s^{-1}]
p.Wc_a = 2.65e-7; 				% respiration rate 								[s^{-1}] 				2.65e-7
p.CO2c_a = 4.87e-7; 			% respiration coefficient 						[s^{-1}]  				4.87e-7
p.k=0.9;                        % Extinction coefficient of canopy              [-] 
p.lar_d=62.5e-3;                % Shoot leaf area ratio                         [m^{2} kg}
p.tao=0.07;                     % Ratio of the root dry weight to the total crop dry weight                 []
p.mass=3.6e-3;                  % mass transfer coefficient (same c_v_pl_ai)    [m s^{-1}
p.v1=9348;                      %                                               [J kg{-1} m^{-3} kmol {-1}]
p.v2=17.4;                      %                                               [-]
p.v3=239;                       %                                               [C]
p.Q10_resp=2;
p.photI0 = 3.55e-9; 			% light use efficiency 							[kg{CO2} J^{-1}]  		3.55e-9
p.photCO2_1 = 5.11e-6;  	    % temperature influence on photosynthesis 		[m s^{-1} °C^{-2}] 		5.11e-6
p.photCO2_2 = 2.3e-4;		    % temperature influence on photosynthesis 		[m s^{-1} °C^{-1}] 		2.3e-4
p.photCO2_3 = 6.29e-4; 			% temperature influence on photosynthesis 		[m s^{-1}] 				6.29e-4
p.photGamma = 5.2e-5; 			% carbon dioxide compensation point 			[kg{CO2} m^{-3}{air}] 	5.2e-5
p.evap_c_a = 3.6e-3; 			% coefficient of leaf-air vapor flow 			[m s^{-1}] 				3.6e-3
p.delta_G=0.476e6;              % Gibbs free energy [J.mol-1]
p.MW_CO2=44.01e-3;              % molecular weight CO2 {kg.mol-1}
p.pl_d=53;                      % effective canopy surface m2.kg-1
p.resp=4.87e-7;
p.fw=22.5;                      % the ration of crop fresh weight to crop dry weight (cultivar related) 

%% LED parametres 
% parameter 					description 									[unit] 					nominal value
p.eta_PAR = 0.31;               % led lamp effecinecy for PAR
p.eta_FIR = 0.0237;             % led lamp waste of FIR
p.eta_NIR = 0.02;               % led lamp waste of NIR
p.eta_convection = 0.0163;    
p.eta_cooling = 0.63;
p.eta_conv = p.eta_cooling + p.eta_convection;                     % LED with no cooling/passive cooling
p.PPE = 3;                       % possible amounts  2; 2.5; 3;                       % micromol (PAR)/J (input)
p.xi_lampPAR = 5.4;                % micromol (PAR)/J(PAR) ---  average of PAR

%% Size and Layout of VF 
% parameter 					description 									[unit] 					nominal value
p.N_compartments=1;                       % number of compartments in VF (only first compartment) 
p.Hieght =0.50;                           % height of each compartment 
p.Weidth =1.5 ;                           % Wedith of each compartment 
p.Length =1.5 ;                           % Lenghth of each compartment 
p.c_d = 0.7 ;                             % discharge coefficent   (perforated wall)
p.OAR = 0.20 ;                            % open air ratio 20-50%  (perforated wall)
p.A_wall = p.Hieght*p.Weidth;             % area of the input air wall - this case is a porforated wall                   [m^{2}]
p.A_duct = p.c_d* p.OAR * p.A_wall;       % area of the input air duct - this case is a porforated wall                   [m^{2}]
p.A_compartment = p.Weidth * p.Length;    % [m^{2}] (2.4*1)=(length*width) 
p.A_cultivation = p.A_compartment * p.N_compartments;              % area under cultivation              
p.volum = p.A_cultivation*p.Hieght;       % air volume  [m^{3}]
p.cop = 3;                                % COP of air conditiong system 
