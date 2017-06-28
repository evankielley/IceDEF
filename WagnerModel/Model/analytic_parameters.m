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
a = @(U) sqrt(2)./U.^3.*(1-sqrt(1+U.^4));
b_big = @(U) real(1./U.^3.*sqrt((4+U.^4).*sqrt(1+U.^4)-3*U.^4-4));
b_small = @(U) real(U.^3./sqrt(8).*sqrt(1-3.*U.^4./4+9.*U.^8./16-7.*U.^12./16+45.*U.^16./256));

% Melt parameters (WDE16 Appendix) ----------------------------------------
Ti0 = -4;
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1;
CMv1 = 7.62e-3; CMv2 = 1.29e-3;
CMe1 = 0.5;
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2;
