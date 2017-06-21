% -------------------------------------------------------------------------
% Plot Iceberg Trajectories------------------------------------------------
%
% This loads in output files from "iceberg_shell.m" for different iceberg
% sizes and plots them (colorcoded). 
%
% Requires the matlab mapping toolbox. I made the atl.mat file to quickly
% plot a low-resolution landmask of the north atlantic. This is called on 
% in "setmap4.m" Uncomment the "plotm(berglimit ..." line to plot the 
% approximate normal iceberg range from the International Ice Patrol.
%
% See 
%
% Till Wagner, Oct 2016, tjwagner@ucsd.edu
% -------------------------------------------------------------------------

load atl
load berglimit

figure(1); clf
latlim = [33 68];
lonlim = [-80 0];

modelfull = 'ECCO_20th';
modelshort= 'E2';
root   = '~/WagnerModel';  % root directory for project
outloc = strcat(root,'/output/',modelfull,'/'); % output directory

N = 25;   %number of trajectories to be plotted

bvec = [9 7 5 3 1];  %specify the different iceberg sizes you're interested in

ii = 0;
setmap4
for i  = bvec
    load(strcat(outloc,sprintf('%s_B%d_full',modelshort,i)))
    bb = i; bs = bb;
    XL = XIL;
    YL = YIL;
    %         ind = randi([1,size(XL,1)],min(size(XL,1),N),1);  %pick random index
    ind = 1:N;
    ii = ii+1;
    for tind = 1:length(ind)
        vend = find(VOL(ind(tind),:)<.1*VOL(ind(tind),1),1,'first'); %plot until iceberg is 90% decayed
        if isempty(vend)
            vend = size(VOL,2);
        end
        if tind == 1
            p{ii} = plotm(YL(ind(1),1:vend)',XL(ind(1),1:vend)','-','col',1-[bb/15 1-bb/15 1-bb/15]);
        else
            plotm(YL(ind(tind),1:vend)',XL(ind(tind),1:vend)','-','col',1-[bb/15 1-bb/15 1-bb/15]);
        end
    end
end

for i = 1:ii; pp = p{i}; pl{i}=pp(1); leg{i} = sprintf('%d',bvec(i)); end
pos = get(gca,'position');
legend([pl{:}],{leg{:}},'location','northeastoutside','interpreter','latex');
set(gca,'position',pos)
%         plotm(berglimit(:,2),berglimit(:,1),'k--')
title(sprintf('%s - 1992 - %d icebergs per size class',modelshort,N))
gridm on