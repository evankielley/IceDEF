% Melt Model --------------------------------------------------------------
%
% This script computes the iceberg melting component, and is called upon at
% every timestep of "iceberg_shell.m". A discussion of this melting module
% is given in the Appendix of WDE17. 
%
% Till Wagner, Oct 2016, tjwagner@ucsd.edu
% -------------------------------------------------------------------------
I = i; 
mmelted(I)=false;
% Compute Melt Terms-------------------------------------------------------
Me = CMe1*(Cs1*Ua^Cs2 + Cs3*Ua);      % wind driven erosion
Mv = CMv1*SST + CMv2*SST^2;           % thermal side wall erosion
Mb = CMb1*sqrt((ui - uw).^2 + (vi - vw).^2)^CMb2...
    *(SST-Ti0)/l(I)^CMb3;             % turbulent basal melt
% Apply Melt Rates --------------------------------------------------------
dldt = - Mv - Me; dhdt = - Mb;
l(I+1) = l(I) + dldt*Dt; 
w(I+1) = w(I)+dldt*Dt; 
h(I+1)=h(I)+dhdt*Dt;
% -------------------------------------------------------------------------
% Make sure the berg is not negative size
if l(I+1)<0 || w(I+1)<0 || h(I+1)<0
    l(I+1)=0; w(I+1)=0; h(I+1)=0;
    melted = 1;
    mmelted(I)=true;
    mm = mm+1;
end
% Rollovers:
if w(I+1) < 0.85 * h(I+1)
    hn = w(I+1); w(I+1) = h(I+1); h(I+1) = hn;
end
% Make sure length>width:
if w(I+1)>l(I+1)
    wn = l(I+1); l(I+1)=w(I+1); w(I+1) = wn;
end
% save length, width, and height
ml(I+1)=l(I+1); mw(I+1)=w(I+1);mh(I+1)=h(I+1);

% Compute new volume and dvol
v(I+1) = l(I+1)*w(I+1)*h(I+1);
dv(I+1) = v(I)-v(I+1);
% Check whether iceberg survived-------------------------------------------
if I == lt-1 && v(I+1) > 0
    ss = ss+1;
end
% -------------------------------------------------------------------------
Mev(I) = Me; Mvv(I) = Mv; Mbv(I) = Mb;  %store melt rates