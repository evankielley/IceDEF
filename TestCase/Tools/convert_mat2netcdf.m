clear;

root_dir = '/home/evankielley/IceDEF/TestCase/';
input_dir = strcat(root_dir,'Inputs/');
output_dir = strcat(root_dir,'Outputs/');

tic
load(strcat(input_dir,'mask.mat'));
load(strcat(input_dir,'E2_vels_1992.mat'));
load(strcat(input_dir,'E2_sst_1992.mat'));
fprintf('model data loaded \n')
toc

LAT = double(vel.latw); LON = double(vel.lonw);
minLAT = min(LAT(:)); maxLAT = max(LAT(:));
minLON = min(LON(:)); maxLON = max(LON(:));

final_t = 122;
t = 1:final_t;

uwF = vel.uw(:,:,t); vwF = vel.vw(:,:,t);
uaF = vel.ua(:,:,t); vaF = vel.va(:,:,t);
sst = double(sst(:,:,t));
mask = int64(msk);

%uwF = vel.uw(1:120,1:120,1:120); vwF = vel.vw(1:120,1:120,1:120);
%uaF = vel.ua(1:120,1:120,1:120); vaF = vel.va(1:120,1:120,1:120);
%sst = double(sst(1:120,1:120,1:120));

%--------------------------------------------------------------------------

ncid = netcdf.create('../Inputs/test_data.nc','NC_WRITE');

dimidLON = netcdf.defDim(ncid,'LonDim',340);
dimidLAT = netcdf.defDim(ncid,'LatDim',360);
dimidTIME = netcdf.defDim(ncid,'TimeDim',122);

%dimidLON = netcdf.defDim(ncid,'LonDim',120);
%dimidLAT = netcdf.defDim(ncid,'LatDim',120);
%dimidTIME = netcdf.defDim(ncid,'TimeDim',120);

varid0 = netcdf.defVar(ncid,'uaF','NC_DOUBLE',[dimidLON dimidLAT dimidTIME]);
varid1 = netcdf.defVar(ncid,'vaF','NC_DOUBLE',[dimidLON dimidLAT dimidTIME]);
varid2 = netcdf.defVar(ncid,'uwF','NC_DOUBLE',[dimidLON dimidLAT dimidTIME]);
varid3 = netcdf.defVar(ncid,'vwF','NC_DOUBLE',[dimidLON dimidLAT dimidTIME]);
varid4 = netcdf.defVar(ncid,'sst','NC_DOUBLE',[dimidLON dimidLAT dimidTIME]);
varid5 = netcdf.defVar(ncid,'mask','NC_INT',[dimidLON dimidLAT]);

netcdf.endDef(ncid);

netcdf.putVar(ncid, varid0, uaF);
netcdf.putVar(ncid, varid1, vaF);
netcdf.putVar(ncid, varid2, uwF);
netcdf.putVar(ncid, varid3, vwF);
netcdf.putVar(ncid, varid4, sst);
netcdf.putVar(ncid, varid5, mask);


netcdf.close(ncid);

ncid2 = netcdf.open('../Inputs/test_data.nc','NC_NOWRITE');

uaF_copy = netcdf.getVar(ncid2, varid0);
vaF_copy = netcdf.getVar(ncid2, varid1);
uwF_copy = netcdf.getVar(ncid2, varid2);
vwF_copy = netcdf.getVar(ncid2, varid3);
sst_copy = netcdf.getVar(ncid2, varid4);
mask_copy = netcdf.getVar(ncid2, varid5);


if isequal(uaF, uaF_copy)
      disp('uaF Data match');
else
      disp('uaF Data mis-match');
end

if isequal(vaF, vaF_copy)
      disp('vaF Data match');
else
      disp('vaF Data mis-match');
end

if isequal(uwF, uwF_copy)
      disp('uwF Data match');
else
      disp('uwF Data mis-match');
end

if isequal(vwF, vwF_copy)
      disp('vwF Data match');
else
      disp('vwF Data mis-match');
end

if isequal(sst, sst_copy)
      disp('sst Data match');
else
      disp('sst Data mis-match');
end

if isequal(mask, mask_copy)
      disp('mask Data match');
else
      disp('mask Data mis-match');
end