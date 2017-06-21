x = linspace(1,10,1);
y = linspace(1,10,1);
z = linspace(1,10,1);
Vi = interp3(x,y,z,uaF(x,y,z),[4,4.5],[4,4.5],[1,2]);

%[xx,yy,zz]=meshgrid(x,y,z)
%V0=meshgrid(x,y,z);
%V=uaF(1:10,1:10,1);

%V=griddedInterpolant(xx,yy,z);
%surf(V);