clear;

root_dir = '/home/evankielley/IceDEF/TestCase/';
input_dir = strcat(root_dir,'Inputs/');

load(strcat(input_dir,'mask.mat'));
intmsk = int64(msk);

ncid = netcdf.create('../Inputs/mask.nc','NC_WRITE');

dimidLON = netcdf.defDim(ncid,'LonDim',340);
dimidLAT = netcdf.defDim(ncid,'LatDim',360);

varid0 = netcdf.defVar(ncid,'intmsk','NC_INT',[dimidLON dimidLAT]);

netcdf.endDef(ncid);

netcdf.putVar(ncid, varid0, intmsk);

netcdf.close(ncid);

ncid2 = netcdf.open('../Inputs/mask.nc','NC_NOWRITE');

intmsk_copy = netcdf.getVar(ncid2, varid0);

if isequal(intmsk, intmsk_copy)
      disp('msk Data match');
else
      disp('msk Data mis-match');
end