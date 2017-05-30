import scipy.io as sio

def load_mat_vars(outloc):

    mrandoX = sio.loadmat(outloc + 'mrandox.mat')
    mrandoX = mrandoX['mrandox']; mrandoX = mrandoX[0,:]                                      
                                                                     
    mrandoY = sio.loadmat(outloc + 'mrandoy.mat')
    mrandoY = mrandoY['mrandoy']; mrandoY = mrandoY[0,:]                                      
                                                                     
    mts = sio.loadmat(outloc + 'mts.mat')
    mts = mts['mts']; mts = mts[0,:]                                              
                                                                     
    mt1= sio.loadmat(outloc + 'mt1.mat')
    mt1 = mt1['mt1']; mt1 = mt1[0,:]                                                      
                                                                     
    mtimestep = sio.loadmat(outloc + 'mtimestep.mat')
    mtimestep = mtimestep['mtimestep']; mtimestep = mtimestep[0,:]                                          
                                                                     
    mua = sio.loadmat(outloc + 'mua.mat')
    mua = mua['mua']; mua = mua[0,:]                                                      
                                                                     
    mva = sio.loadmat(outloc + 'mva.mat')
    mva = mva['mva']; mva = mva[0,:]                                                      
                                                                     
    muw = sio.loadmat(outloc + 'muw.mat')
    muw = muw['muw']; muw = muw[0,:]                                                      
                                                                     
    mvw = sio.loadmat(outloc + 'mvw.mat')
    mvw = mvw['mvw']; mvw = mvw[0,:]      

    mSST = sio.loadmat(outloc + 'mSST.mat')
    mSST = mSST['mSST']; mSST = mSST[0,:]                                                    
                                                                        
    mui = sio.loadmat(outloc + 'mui.mat')
    mui = mui['mui']; mui = mui[0,:]                                                      
                                                                        
    mvi = sio.loadmat(outloc + 'mvi.mat')
    mvi = mvi['mvi']; mvi = mvi[0,:]                                                      
                                                                        
    mXI = sio.loadmat(outloc + 'mXI.mat')
    mXI = mXI['mXI']; mXI = mXI[0,:]                                                      
                                                                        
    mYI = sio.loadmat(outloc + 'mYI.mat')
    mYI = mYI['mYI']; mYI = mYI[0,:] 

    mUT = sio.loadmat(outloc + 'mUT.mat')
    mUT = mUT['mUT']; mUT = mUT[0,:] 

    mUa = sio.loadmat(outloc + 'mUa.mat')
    mUa = mUa['mUa']; mUa = mUa[0,:] 

    mXIL = sio.loadmat(outloc + 'mXIL.mat')
    mXIL = mXIL['XIL']; mXIL = mXIL[0,:] 

    mYIL = sio.loadmat(outloc + 'mYIL.mat')
    mYIL = mYIL['YIL']; mYIL = mYIL[0,:] 

    ml = sio.loadmat(outloc + 'ml.mat')
    ml = ml['l']; ml = ml[0,:] 

    mw = sio.loadmat(outloc + 'mw.mat')
    mw = mw['w']; mw = mw[0,:] 

    mh = sio.loadmat(outloc + 'mh.mat')
    mh = mh['h']; mh = mh[0,:] 

    mv = sio.loadmat(outloc + 'mv.mat')
    mv = mv['v']; mv = mv[0,:] 

    return mrandoX,mrandoY,mts,mt1,mtimestep,mua,mva,muw,mvw,mSST,mui,mvi,mXI,mYI,mUT,mUa,mXIL,mYIL,ml,mw,mh,mv
