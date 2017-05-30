def init_arrs1(trajnum, nt):
    XIL, YIL = np.empty((trajnum,numSteps))*np.nan, np.empty((trajnum,numSteps))*np.nan
    VOL, DVOL = np.empty((trajnum,numSteps))*np.nan, np.empty((trajnum,numSteps))*np.nan
    UI, VI = np.empty((trajnum,numSteps))*np.nan, np.empty((trajnum,numSteps))*np.nan
    UA, VA = np.empty((trajnum,numSteps))*np.nan, np.empty((trajnum,numSteps))*np.nan
    UW, VW = np.empty((trajnum,numSteps))*np.nan, np.empty((trajnum,numSteps))*np.nan
    TE = np.empty((trajnum,numSteps))*np.nan                            
    Memat, Mvmat = np.empty((trajnum,numSteps))*np.nan, np.empty((trajnum,numSteps))*np.nan
    Mbmat = np.empty((trajnum,numSteps))*np.nan 
    return XIL,YIL,VOL,DVOL,UI,VI,UA,VA,UW,VW,TE,Memat,Mvmat,Mbmat

def init_arrs2(lt):
    xil,yil = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    v,dv = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    uiv,viv = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    uav,vav = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    uwv,vwv = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    temp,Mev = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    Mvv,Mbv = np.empty(lt)*np.nan, np.empty(lt)*np.nan)
    return xil,yil,v,dv,uiv,viv,uav,vav,uwv,vwv,temp,Meb,Mev,Mbv
