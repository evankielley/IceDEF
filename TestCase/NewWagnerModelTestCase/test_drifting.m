% TEST CASE

% set dummy iterator ------------------------------------------------------
I = i;

% set initial conditions --------------------------------------------------
interpolate = false;
mob(I) = false;
mgrounded(I) = false;

% For regular grid only
YI = dsearchn(LAT,yil(I));
XI = dsearchn(LON,xil(I));
TI = dsearchn(T_ALL, I);       

ua = uaF(XI,YI,I);
va = vaF(XI,YI,I);
uw = uwF(XI,YI,I);
vw = vwF(XI,YI,I);
SST= sst(XI,YI,I);
    
% compute wind speed and "U tilde" at location (for given iceberg size)----
Ua = sqrt(ua^2+va^2);
UT = Ut(Ua,yil(I),S(l(I),w(I)));   %U tilde is \Lambda in the papers

% now compute analytic iceberg velocity solution---------------------------
if UT >= 0.1
    alpha = a_naive(UT);

else
    fprintf('%i Taylor approx used for alpha\n', I);
    alpha = a_taylor(UT);
end

if UT >= 0.6
    beta = b_naive(UT);
else
    fprintf('%i Taylor approx used for beta\n', I);
    beta = b_taylor(UT);
end
ui = uw - g*alpha*va + g*beta*ua;
vi = vw + g*alpha*ua + g*beta*va;

% iceberg translation (note that I'm converting from m to deg lat/lon)-----
dlon = ui*dtR;
dlat = vi*dtR;

yil(I+1) = yil(I) + dlat;
xil(I+1) = xil(I) + dlon/cos((yil(I+1)+yil(I))/2*pi/180);

% check you haven't gone out of bounds-------------------------------------
if xil(I+1)>maxLON || xil(I+1)<minLON || yil(I+1)>maxLAT || yil(I+1)<minLAT
    outofbound = 1;
    mob(I)=true;
    ob = ob+1;

else % now check you didn't send the iceberg on land-----------------------
    yi2(1) = find(LAT<=yil(I+1),1,'last');
    yi2(2) = find(LAT>yil(I+1),1,'first');
    xi2(1) = find(LON<=xil(I+1),1,'last');
    xi2(2) = find(LON>xil(I+1),1,'first');
    if any(find(msk(xi2,yi2)==0))
        fprintf('land ho');
        mgrounded(I) = true;
        yil(I+1) = yil(I);  %i.e. when I get put within one grid box of land
        xil(I+1) = xil(I);  %I assume the iceberg don't move, until it doesn't happen anymore
    end
end

% store into output vectors
mxi(I) = XI; myi(I) = YI;
mua(I) = Ua; mut(I) = UT;
uiv(I) = ui; viv(I) = vi;
uav(I) = ua; vav(I) = va;
uwv(I) = uw; vwv(I) = vw;
temp(I) = SST;
malpha(I) = alpha; mbeta(I) = beta;