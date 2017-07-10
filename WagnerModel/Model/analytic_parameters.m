% Iceberg Parameters ------------------------------------------------------
%
% This script computes the iceberg parameters of WDE17 and
% analytical expressions of gamma, Lambda, alpha, beta (eqns 7,8,9)
%
% Till Wagner, Oct 2016, tjwagner@ucsd.edu
% -------------------------------------------------------------------------

% set constants -----------------------------------------------------------
R = 6378*1e3;           % earth radius in m
rhow = 1027;            % density of water (kg/m^3)
rhoa = 1.2;             % density of air   (kg/m^3)
rhoi = 850;             % density of ice   (kg/m^3)
drho = rhow-rhoi;
Cw = 0.9;               % bulk coefficient water  (Bigg et al 1997)
Ca = 1.3;               % bulk coefficient air    (Bigg et al 1997)
om = 7.2921*1e-5;       % rotation rate of earth (rad/s)

% define necessary functions ----------------------------------------------
ff = @(lati) 2*om*sin(abs(lati)*pi/180);        % latitude in degrees
g = sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw));          % gamma = sqrt(ca/cw)
S = @(l,w) pi*l.*w./(l+w);                      % harmonic mean length
Ut = @(u,lati,S) sqrt(2)*Cw*g/ff(lati).*u/S;    % Lambda in the papers

% alpha and beta functions
% these are slightly differently in the papers (sqrt(2)) -- not sure why
% b_big is accurate for U greater than 0.1 otherwise b_small is to be used
a_naive = @(U) sqrt(2)./U.^3.*(1-sqrt(1+U.^4));
a_taylor = @(U) sqrt(2).*(-U./2 + U.^5/8 - U.^9./16 + 5.*U.^13./128 - 7.*U.^17./256);
%b_small = @(U) real(U.^3./sqrt(8).*sqrt(1-3.*U.^4./4+9.*U.^8./16-7.*U.^12./16+45.*U.^16./256));
b_naive = @(U) real(1./U.^3.*sqrt((4+U.^4).*sqrt(1+U.^4)-3*U.^4-4));
b_taylor = @(U) U.^3.*(U.^4.*(0.0745776683282687.*U.^4 - 0.132582521472478) ...
                + 0.353553390593274) - 4.93696020934508e-17;


b23 = @(U) real(sqrt(2)*(1/4*U^3 - 3/32*U^7 + 27/512*U^11 -143/4096*U^15 - ...
                    3315/131072*U^19));
b27 = @(U) real(sqrt(2)*(1/4*U^3 - 3/32*U^7 + 27/512*U^11 -143/4096*U^15 + ...
                    3315/131072*U^19 + 17667/1048576*U^23));
b31 = @(U) real(sqrt(2)*(1/4*U^3 - 3/32*U^7 + 27/512*U^11 -143/4096*U^15 + ...
                    3315/131072*U^19 + 17667/1048576*U^23 + 252497/16777216*U^27));
b35 = @(U) real(sqrt(2)*(1/4*U^3 - 3/32*U^7 + 27/512*U^11 -143/4096*U^15 + ...
                    3315/131072*U^19 - 17667/1048576*U^23 + 260015/16777216*U^27 + ...
                    1803513/134217728*U^31));

% Melt parameters (WDE16 Appendix) ----------------------------------------
Ti0 = -4;
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1;
CMv1 = 7.62e-3; CMv2 = 1.29e-3;
CMe1 = 0.5;
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2;
