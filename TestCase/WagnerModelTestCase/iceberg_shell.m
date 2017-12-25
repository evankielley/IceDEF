% TEST CASE

% set directories ---------------------------------------------------------
root_dir = '/home/evankielley/IceDEF/TestCase/';  % root dir for project
input_dir = strcat(root_dir,'Inputs/'); % input directory
output_dir = strcat(root_dir,'Outputs/'); % output directory

% load input fields -------------------------------------------------------
tic
load(strcat(input_dir,'bergdims.mat'));  %landmask
load(strcat(input_dir,'mask.mat'));  %landmask
load(strcat(input_dir,'Laurent_Seed.mat'));  %landmask
load(strcat(input_dir,'E2_vels_1992.mat'));
load(strcat(input_dir,'E2_sst_1992.mat'));
load(strcat(input_dir,'fixed.mat'),'ts_all','randoX_all','randoY_all');
fprintf('model data loaded \n')
toc

% read in all parameters and analytic expressions for alpha and beta ------
test_analytic_parameters

% specify the space domain ------------------------------------------------
LAT = double(vel.latw); LON = double(vel.lonw);
minLAT = min(LAT(:)); maxLAT = max(LAT(:));
minLON = min(LON(:)); maxLON = max(LON(:));

% set run parameters ------------------------------------------------------
trajnum = 25;            % total number of iceberg trajectories to compute
final_t = 122;           % number of input field time steps
startrange = final_t/2;  % input field start range
tres = 3;                % time resoln such that "model Dt"="input DT"/tres
DT = 3;                  % Input fields time step
Dt = DT/tres;            % model timestep in days
dt = Dt*24*3600;         % model timestep in seconds
R = 6378*1e3;            % earth radius in m
dtR = dt/R*180/pi;       % need this ratio for distances in "drifting.m"
t = 1:final_t;                  % how long is the run
nt= length(t)*tres;             % number of model timesteps
tt = linspace(1,length(t),nt);  % model time

% Load Seeding fields -----------------------------------------------------
% From Laurent_Seed.mat
seed_X = repmat(Seed_X(:),[100,1]); % cycle through each location 100x
seed_Y = repmat(Seed_Y(:),[100,1]); % i.e. this can run 3600 icebergs

% these are the circulation fields-----------------------------------------
uwF = vel.uw(:,:,t); vwF = vel.vw(:,:,t);   % water velocities input
uaF = vel.ua(:,:,t); vaF = vel.va(:,:,t);   % air velocities input
sst = double(sst(:,:,t));                   % sst input


% Load berg sizes (classification from Bigg et al 1997) -------------------
bvec = 1:10;   % vector of which size classes to compute - has to be [1,10]

% set output arrays -------------------------------------------------------
XIL = nan(length(bvec),trajnum,nt); YIL = XIL;   
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

% loop over individual initial iceberg size classes -----------------------
for bb = bvec
    bergsize = bb;   % current berg size class
    fprintf('run bergsize B%d \n',bergsize)
    
    % initialize the iceberg-----------------------------------------------
    L = bergdims(bergsize,1);
    W = bergdims(bergsize,2);
    H = bergdims(bergsize,3);
    
    % run drift and melt---------------------------------------------------
    tic
    mm = 0; ss = 0; ob = 0;  % melted, survived, out of bounds status'
    for j = 1:trajnum
        if mod(j,10)==0; toc; fprintf('%d trajectories computed \n',j); end
        % pick a random trajectory start time (of Input field)
        % (note: this is not the best setup for validation runs, where you
        % may want specific trajectories to compare!)
        %ts = randi([0,round(startrange)],1);
        ts = ts_all(bb,j);
        tts= ts*tres;                  %trajectory start time (of model)
        lt = nt-tts;                   %trajectory run length
        
        % initialize output vectors
        xil = nan(1,lt); yil = xil;
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
        
        % pick fixed seeding location (for model comparison)
        randoX = randoX_all(bb,j);  % fixed set of random int inputs
        randoY = randoY_all(bb,j);  % fixed set of random int inputs
        yig = seed_Y(randoY); xig = seed_X(randoX);
                        
        % set initial values for output vectors
        xil(1) = LON(xig); yil(1) = LAT(yig);   %initial lon and lat
        l = L*ones(1,lt); w = l*W/L; h = l*H/L; %initial berg dimensions
        v(1) = L*W*H; dv(1) = 0;                %initial volume and dvol
        ml(1) = L; mw(1) = W; mh(1) = H;
        
        % set initial conditions for iceberg status (boolean)
        outofbound = 0; melted = 0;
        
        % initialize iterating variable
        i = 0; 
        
        % now integrate as long as the iceberg is in the domain and not
        % melted and over the time period specified above
        
        while outofbound == 0 && melted == 0 && i<lt-1
            i = i+1;
            test_drifting
            test_melting
        end
        
        ind = 1:i+1;
        
        % store output vectors into output matrices
        XIL(bb,j,ind)=xil(ind); YIL(bb,j,ind)=yil(ind);
        mXI(bb,j,ind) = mxi(ind); mYI(bb,j,ind) = myi(ind);
        VOL(bb,j,ind)=v(ind); DVOL(bb,j,ind)=dv(ind);  
        mL(bb,j,ind) = ml(ind); mW(bb,j,ind) = mw(ind); mH(bb,j,ind) = mh(ind);
        mUa(bb,j,ind) = mua(ind); mUT(bb,j,ind) = mut(ind);
        UI(bb,j,ind) = uiv(ind); VI(bb,j,ind) = viv(ind);
        UA(bb,j,ind) = uav(ind); VA(bb,j,ind) = vav(ind);
        UW(bb,j,ind) = uwv(ind); VW(bb,j,ind) = vwv(ind);
        TE(bb,j,ind) =temp(ind);                       
        mMe(bb,j,ind) = Mev(ind); mMv(bb,j,ind) = Mvv(ind); mMb(bb,j,ind) = Mbv(ind);
        mGROUNDED(bb,j,ind)=mgrounded(ind);mOB(bb,j,ind)=mob(ind);
        mMELTED(bb,j,ind)=mmelted(ind);
        mALPHA(bb,j,ind)=malpha(ind); mBETA(bb,j,ind)=mbeta(ind);
        
    end
    % ---------------------------------------------------------------------
    fprintf('%d icebergs died, %d lived, %d left the domain \n',mm,ss,ob)
    % ---------------------------------------------------------------------  
end


% save inputs and outputs for analysis and model comparison
mLAT = LAT; mLON = LON; mmsk = msk;
save(strcat(output_dir,'WagnerTestCaseOutput'),...
    'XIL','YIL','mXI','mYI','VOL','DVOL','mUa','mUT','UI','VI','UA','VA','UW','VW',...
    'mL','mW','mH','mLAT','mLON','mmsk','TE','mMe','mMv','mMb',...
    'mGROUNDED','mOB','mMELTED','mALPHA','mBETA');
