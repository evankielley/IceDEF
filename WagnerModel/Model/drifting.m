% Drift Model -------------------------------------------------------------
%
% This script computes the iceberg drift component, and is called upon at
% every timestep of "iceberg_shell.m". The drift integration is discussed
% in WDE17, Section 3.
%
% Till Wagner, Oct 2016, tjwagner@ucsd.edu
% -------------------------------------------------------------------------
I = i;
interpolate = false;

if interpolate
    
    YI = dsearchn(LAT,yil(I));  % just for saving purposes
    XI = dsearchn(LON,xil(I));  % likewise

    timestep = tt(tts + I); t=timestep;                                                      
    t1 = floor(timestep); t2 = t1 + 1;                      
    ti1=(t2-t)/(t2-t1); ti2=(t-t1)/(t2-t1);                                              
    ti=t1*ti1+t2*ti2;
    
    YI1 = find(LAT<=yil(I),1,'last');
    YI2 = find(LAT>yil(I),1,'first');
    XI1 = find(LON<=xil(I),1,'last');
    XI2 = find(LON>xil(I),1,'first');
    
    x=xil(I);                                                      
    xi1=(LON(XI2)-x)/(LON(XI2)-LON(XI1)); xi2=(x-LON(XI1))/(LON(XI2)-LON(XI1));                            
    xi=XI1*xi1+XI2*xi2; 
    
    y=yil(I);                                                      
    yi1=(LAT(YI2)-y)/(LAT(YI2)-LAT(YI1)); yi2=(y-LAT(YI1))/(LAT(YI2)-LAT(YI1));                            
    yi=YI1*yi1+YI2*yi2; 
    
    %points=[[XI1,XI2],[YI1,YI2],[t1,t2]];                          
    xpoints = [XI1,XI2];
    ypoints = [YI1,YI2];
    tpoints = [t1,t2];
    %xyti=[xi,yi,ti];                                               
    
    uaValues = uaF(XI1:XI2,YI1:YI2,t1:t2);                                                                    
    uaGRID = griddedInterpolant({xpoints,ypoints,tpoints},uaValues);
    ua = uaGRID({xi,yi,ti});
    vaValues = vaF(XI1:XI2,YI1:YI2,t1:t2);                                                                    
    vaGRID = griddedInterpolant({xpoints,ypoints,tpoints},vaValues);
    va = vaGRID({xi,yi,ti});
    uwValues = uwF(XI1:XI2,YI1:YI2,t1:t2);                                                                    
    uwGRID = griddedInterpolant({xpoints,ypoints,tpoints},uwValues);
    uw = uwGRID({xi,yi,ti});
    vwValues = vwF(XI1:XI2,YI1:YI2,t1:t2);                                                                    
    vwGRID = griddedInterpolant({xpoints,ypoints,tpoints},vwValues);
    vw = vwGRID({xi,yi,ti});    
    sstValues = sst(XI1:XI2,YI1:YI2,t1:t2);                                                                    
    sstGRID = griddedInterpolant({xpoints,ypoints,tpoints},sstValues);
    SST = sstGRID({xi,yi,ti});    
   
    
else
    % find nearest neighbour using dsearchn,
    % CAUTION: this only works on a regular grid!
    % if strcmp(modelshort,'E2')
    YI = dsearchn(LAT,yil(I));
    XI = dsearchn(LON,xil(I));

    % else  %CCSM is on a irregular grid!
    %     indic = dsearchn([LAT2,LON2],[yil(I),xil(I)]);
    %     [XI,YI]=find(indic==vec);
    % end*
    % now interpolate fields linearly between timesteps------------------------
    timestep = tt(tts+I);
    t1  = floor(timestep); t2 = t1+1;
    dt1 = timestep-t1; dt2 = t2-timestep;

    ua = uaF(XI,YI,t1)*dt1 + uaF(XI,YI,t2)*dt2;
    va = vaF(XI,YI,t1)*dt1 + vaF(XI,YI,t2)*dt2;
    uw = uwF(XI,YI,t1)*dt1 + uwF(XI,YI,t2)*dt2;
    vw = vwF(XI,YI,t1)*dt1 + vwF(XI,YI,t2)*dt2;
    SST= sst(XI,YI,t1)*dt1 + sst(XI,YI,t2)*dt2;
end
    
% compute wind speed and "U tilde" at location (for given iceberg size)----
Ua = sqrt(ua^2+va^2);
UT = Ut(Ua,yil(I),S(l(I),w(I)));   %U tilde is \Lambda in the papers
% now compute analytic iceberg velocity solution---------------------------
ui = uw - g*a(UT)*va + g*b(UT)*ua;
vi = vw + g*a(UT)*ua + g*b(UT)*va;

% iceberg translation (note that I'm converting from m to deg lat/lon)-----
dlon = ui*dtR;
dlat = vi*dtR;

mxi(I) = XI; myi(I) = YI;
uiv(I) = ui; viv(I) = vi;
uav(I) = ua; vav(I) = va;
uwv(I) = uw; vwv(I) = vw;
temp(I) = SST;

yil(I+1) = yil(I) + dlat;
xil(I+1) = xil(I) + dlon/cos((yil(I+1)+yil(I))/2*pi/180);

% check you haven't gone out of bounds-------------------------------------
if xil(I+1)>maxLON || xil(I+1)<minLON || yil(I+1)>maxLAT || yil(I+1)<minLAT
    outofbound = 1;
    ob = ob+1;
    fprintf('iceberg %d left boundary at timestep %d \n',j,I);
else % now check you didn't send the iceberg on land-----------------------
    %     if strcmp(modelshort,'E2')
    yi2(1) = find(LAT<=yil(I+1),1,'last');
    yi2(2) = find(LAT>yil(I+1),1,'first');
    xi2(1) = find(LON<=xil(I+1),1,'last');
    xi2(2) = find(LON>xil(I+1),1,'first');
    if any(find(msk(xi2,yi2)==0))
        yil(I+1) = yil(I);  %i.e. when I get put within one grid box of land
        xil(I+1) = xil(I);  %I assume the iceberg don't move, until it doesn't happen anymore
    end
    %     else
    %         indic = dsearchn([LAT2,LON2],[yil(I+1),xil(I+1)]);
    %         [XI,YI]=find(indic==vec);
    %         if vel.mask(XI,YI)==1
    %             yil(I+1) = yil(I);  %i.e. when I get put within one grid box of land
    %             xil(I+1) = xil(I);  %I assume the iceberg don't move, until it doesn't happen anymore
    %         end
    %     end
end
% -----------------------------------------------------------------