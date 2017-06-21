load atl

if ~exist('latlim','var')
    latlim = [25 68];
    lonlim = [-80 0];
end
worldmap(latlim,lonlim)
fra = framem('on');
set(fra,'clipping','on');
mlp = mlabel('on');
set(mlp,'clipping','on'); plp = plabel('on');
set(plp,'clipping','on'); gri = gridm('off');
set(gri,'clipping','on');
geo = geoshow(atl(1:10), 'FaceColor', [.9 .9 .9]);
setm(gca,'FFaceColor',[1 1 1])
hold on

hm = mlabel('on');
set(hm(:),'fontsize',13);
set(hm(:),'verticalalignment','baseline');
pm = plabel('on');
set(pm(:),'fontsize',13);
set(pm(:),'verticalalignment','baseline');