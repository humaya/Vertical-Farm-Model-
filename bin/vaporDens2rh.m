function rh = vaporDens2rh(temp, vaporDens)
% vaporDens2rh Convert vapor density [kg{H2O} m^{-3}] to relative humidity [%]
%
% Usage:
%   rh = vaporDens2rh(temp, vaporDens)
% Inputs:
%   temp        given temperatures [°C] (numeric vector)
%   vaporDens   absolute humidity [kg{H20} m^{-3}] (numeric vector)
%   Inputs should have identical dimensions
% Outputs:
%   rh          relative humidity [%] between 0 and 100 (numeric vector)
%
% Calculation based on 
%   http://www.conservationphysics.org/atmcalc/atmoclc2.pdf

% David Katzin, Wageningen University
% david.katzin@wur.nl

    % constants
    R = 8.3144598; % molar gas constant [J mol^{-1} K^{-1}]
    C2K = 273.15; % conversion from Celsius to Kelvin [K]
    Mw = 18.01528e-3; % molar mass of water [kg mol^-{1}]
    
    % parameters used in the conversion
    p = [610.78 238.3 17.2694 -6140.4 273 28.916];
        % default value is [610.78 238.3 17.2694 -6140.4 273 28.916]
    
    satP = p(1)*exp(p(3)*temp./(temp+p(2))); 
        % Saturation vapor pressure of air in given temperature [Pa]
       
    % convert to relative humidity using the ideal gas law pV=nRT => n=pV/RT 
    % so n=p/RT is the number of moles in a m^3, and Mw*n=Mw*p/(R*T) is the 
    % number of kg in a m^3, where Mw is the molar mass of water.
        
    rh = min(100,(100*R*(temp+C2K)./(Mw*satP)).*vaporDens);
end