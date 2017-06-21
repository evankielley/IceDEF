% Construct Seeding -------------------------------------------------------
%
% This script constructs a grid of release locations near the outlet
% glacier in question. See WDE17.
%
% Till Wagner, Oct 2016, tjwagner@ucsd.edu
% -------------------------------------------------------------------------
glacier = 'H';  %pick the glacier you want to construct the seeding for

if glacier == 'H'
    % Helheim Ice Sheet Seeding Locations ---------------------------------
    Glac_Y = 259:1:262; %[261,190] 66.3500° N, 38.2000° W coord, Helheim
    Glac_X = 188:1:192;
    glac_name = 'Helheim';
elseif glacier == 'J'
    % Jakobshaven Ice Sheet Seeding Locations -----------------------------
    Glac_Y = 276:1:280; %[277,120] 69°10 N 49°50 W coord, Jakobshavn
    Glac_X = 119:1:123;
    glac_name = 'Jakobsh';
elseif glacier == 'K'
    % Kangerd Ice Sheet Seeding Locations ---------------------------------
    Glac_Y = 267:1:271; %[277,120] 68°38 N 33°0 W coord, Kangerd
    Glac_X = 213:1:218;
    glac_name = 'Kangerd';
elseif glacier == 'L';
    % Laurentide Ice Sheet Seeding Locations ------------------------------
    Glac_Y = 245:2:255; %Laurentide LAT seed locations (in ECCO 2 grid)
    Glac_X = 86:2:101;   %Laurentide LON seed locations (in ECCO 2 grid)
    glac_name = 'Laurent';
end

seed_LAT{1} =  [LAT(Glac_Y(1))*[1 1]; LAT(Glac_Y(end))*[1 1]];
seed_LON{1} = repmat([LON(Glac_X(1)) LON(Glac_X(end))],[2,1]);
seed_LAT{2} =  repmat([LAT(Glac_Y(1)) LAT(Glac_Y(end))],[2,1]);
seed_LON{2} = [LON(Glac_X(1))*[1 1]; LON(Glac_X(end))*[1 1]];
[Seed_Y, Seed_X] = meshgrid(Glac_Y, Glac_X);
% save(strcat(sprintf('%s',glac_name),'_Seed'),'Seed_Y','Seed_X','seed_LAT','seed_LON')

% figure(4); %clf
% latlim = [50 72]; lonlim = [-70 -25];
% setmap4
% hold on;
% plotm(LAT(Seed_Y), LON(Seed_X),'x')
% plotm(seed_LAT{1},seed_LON{1},'k','linewidth',2)
% plotm(seed_LAT{2},seed_LON{2},'k','linewidth',2)