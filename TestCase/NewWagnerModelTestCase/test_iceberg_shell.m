% TEST CASE

% set directories ---------------------------------------------------------
root_dir = '/home/evankielley/IceDEF/TestCase/';  % root dir for project
input_dir = strcat(root_dir,'Inputs/'); % input directory
output_dir = strcat(root_dir,'Outputs/'); % output directory

% load input fields -------------------------------------------------------
tic
load(strcat(input_dir,'mask.mat'));  %landmask
load(strcat(input_dir,'E2_vels_1992.mat'));
load(strcat(input_dir,'E2_sst_1992.mat'));
fprintf('model data loaded \n')
toc

% read in all parameters and analytic expressions for alpha and beta ------
test_analytic_parameters

% specify the space domains -----------------------------------------------
LAT = double(vel.latw); LON = double(vel.lonw);
minLAT = min(LAT(:)); maxLAT = max(LAT(:));
minLON = min(LON(:)); maxLON = max(LON(:));

% set run parameters ------------------------------------------------------
final_t = 122;           % number of input field time steps
startrange = final_t/2;  % input field start range
tres = 1;                % time resoln such that "model Dt"="input DT"/tres
DT = 1;                  % Input fields time step
Dt = DT/tres;            % model timestep in days
dt = Dt*24*3600;         % model timestep in seconds
R = 6378*1e3;            % earth radius in m
dtR = dt/R*180/pi;       % need this ratio for distances in "drifting.m"
t = 1:final_t;                  % how long is the run
nt= length(t)*tres;             % number of model timesteps
tt = linspace(1,length(t),nt);  % model time

T_ALL = linspace(1, final_t, 1);

% these are the circulation fields-----------------------------------------
uwF = vel.uw(:,:,t); vwF = vel.vw(:,:,t);   % water velocities input
uaF = vel.ua(:,:,t); vaF = vel.va(:,:,t);   % air velocities input
sst = double(sst(:,:,t));                   % sst input


% set output arrays -------------------------------------------------------
XIL = nan(1,nt); YIL = XIL;   
mXI = XIL; mYI = XIL;
VOL = XIL; DVOL = VOL;                           
mL = XIL; mW = XIL; mH = XIL;
mUa = XIL; mUT = XIL;
UI = XIL; UA = XIL; UW = XIL;                    
VI = XIL; VA = XIL; VW = XIL;                    
TE = XIL;                                        
mMe = XIL; mMv = XIL; mMb = XIL;
mOB = XIL; mGROUNDED = XIL; mMELTED = XIL;
mALPHA = XIL; mBETA = XIL;

% initialize the iceberg-----------------------------------------------

% run drift and melt---------------------------------------------------
tic
mm = 0; ss = 0; ob = 0;  % melted, survived, out of bounds status'

% initialize output vectors
xil = nan(1,nt); yil = xil;
mxi = xil; myi = xil;
v = xil; dv = xil;
ml = xil; mw = xil; mh = xil;
mua = xil; mut = xil;
uiv = v; uav = v; uwv = v;
viv = v; vav = v; vwv = v;
temp = v;
Mev = v; Mvv = v; Mbv = v;
mob = v; mgrounded = v; mmelted = v;
malpha = v; mbeta = v;


% set initial values for output vectors
xil(1) = 310; 
yil(1) = 50;
L = 600;
W = 500;
H = 400;
l = L*ones(1,nt); w = W*ones(1,nt); h = H*ones(1,nt); 
v(1) = L*W*H; dv(1) = 0;                
ml(1) = L; mw(1) = W; mh(1) = H;

% set initial conditions for iceberg status (boolean)
outofbound = 0; melted = 0;

% initialize iterating variable
i = 0; 


while outofbound == 0 && melted == 0 && i < nt-1
    i = i+1;
    test_drifting
    test_melting
end

ind = 1:i+1;

% store output vectors into output matrices
XIL(ind)=xil(ind); YIL(ind)=yil(ind);
mXI(ind) = mxi(ind); mYI(ind) = myi(ind);
VOL(ind)=v(ind); DVOL(ind)=dv(ind);  
mL(ind) = ml(ind); mW(ind) = mw(ind); mH(ind) = mh(ind);
mUa(ind) = mua(ind); mUT(ind) = mut(ind);
UI(ind) = uiv(ind); VI(ind) = viv(ind);
UA(ind) = uav(ind); VA(ind) = vav(ind);
UW(ind) = uwv(ind); VW(ind) = vwv(ind);
TE(ind) =temp(ind);                       
mMe(ind) = Mev(ind); mMv(ind) = Mvv(ind); mMb(ind) = Mbv(ind);
mGROUNDED(ind)=mgrounded(ind);mOB(ind)=mob(ind);
mMELTED(ind)=mmelted(ind);
mALPHA(ind)=malpha(ind); mBETA(ind)=mbeta(ind);

% save inputs and outputs for analysis and model comparison
mLAT = LAT; mLON = LON; mmsk = msk;
save(strcat(output_dir,'NewWagnerTestCaseOutput'),...
    'XIL','YIL','mXI','mYI','VOL','DVOL','mUa','mUT','UI','VI','UA','VA','UW','VW',...
    'mL','mW','mH','mLAT','mLON','mmsk','TE','mMe','mMv','mMb',...
    'mGROUNDED','mOB','mMELTED','mALPHA','mBETA');
