import os
from pathlib import Path
import numpy as np
import math

import prep_Coarse_read as readf
import prep_Coarse_upsc as upsc
import prep_Coarse_LGR as lgr
import prep_Coarse_write as writef



# =============================================================================
# User inputs
# =============================================================================
def createCoarseFile(path_curr, fact):
    # coarseLevels = [[5], 
    #                 [5], 
    #                 [5]] # scaling factors
    kz_ratio = 0.12
    bisection_tol = 1e-3
    
    # %% Define frequently used directory paths
    path_fine_model = path_curr
    path_pss_model = os.path.join(path_curr, 'PSS_Upsc') # path to the psuedo steady state upscaling model
    path_coarse_model = os.path.join(path_curr, 'Coarse_LGR')
    
    # %% Define variables for file names   
    PSS_DATA = os.path.join(path_pss_model, 'CO2_ECLIPSE.DATA')
    PSS_OUTPUT = os.path.join(path_pss_model,'CO2_ECLIPSE.RSM') # pss output
    PERM_FILE_TO_READ = os.path.join(path_pss_model,'PERMX.INC')
    PORO_FILE_TO_READ = os.path.join(path_pss_model,'PORO.INC')
    TWO_PHASE_FILE_SCHED = os.path.join(path_fine_model, 'base.sched')
    PSS_FILE_SCHED = os.path.join(path_pss_model, 'base.sched')
    COARSE_FILE_SCHED = os.path.join(path_coarse_model, 'base_C.sched')
    
    # %% Handling grid
    # Read from original fine_scale GPRS file: fine grid dimension and sizes
    grids, gridDims = readf.readDimsFromInput(PSS_DATA)
    totalGrids = grids[0] * grids[1] * grids[2]
    
    # Upscale grid and grid dimension
    # fact = [coarseLevels[coord][0] for coord in range(3)]
    cgrids, cgridDims = upsc.detUpsclGrids(grids, gridDims, fact)
    print(cgrids)
    print(cgridDims)
    
    # %% Handling reservoir properties 
    # Read results
    ssPresVec = readf.readBlockProps(PSS_OUTPUT, grids, item='BPR')
    ssDensity = readf.readBlockProps(PSS_OUTPUT, grids, item='BDENW')

    finePermVec = readf.readRockProps((PERM_FILE_TO_READ), totalGrids)
    finePoroVec = readf.readRockProps((PORO_FILE_TO_READ), totalGrids)
    fineTrans = upsc.fineScaleFlows(finePermVec, grids, gridDims)
    fineRates = upsc.fineScaleFlows_rate(ssPresVec, grids, gridDims, ssDensity, trans = 2, prevDat = fineTrans)
    
    #Upscale poro and transmissibility
    upPoroVec = upsc.upscaleProp(grids, cgrids, fact, finePoroVec) # porosity for output
    upPressVec = upsc.upscaleProp(grids, cgrids, fact, ssPresVec)
    upDensityVec = upsc.upscaleProp(grids, cgrids, fact, ssDensity)
    coarseRates = upsc.getCoarseFlux(fineRates, grids, cgrids, fact)
    
    gridGlobalT = upsc.Trans_upscale(upPressVec, gridDims, fact, cgrids, cgridDims, upDensityVec, trans = 3, prevDat = coarseRates)
    
    # calculate rate and BHP 
    TGeoC, geoAvgK = upsc.geometricTrans(fact, grids ,cgrids, cgridDims, finePermVec)
    
    finalT = upsc.correctTrans(gridGlobalT, TGeoC, wells = False) # corrected trans for output
    
    # %% LGR
    
    
    f= open("coarsePerm_log.txt","w+")
    f.close()
    
    # totalGrids_fine = totalGrids
    # totalGrids_coarse = cgrids[0]*cgrids[1]*cgrids[2]
    areas = {'x':gridDims[1]*gridDims[2],
              'y':gridDims[0]*gridDims[2],
              'z':gridDims[0]*gridDims[1]
            }
    
    distances = {'x':[gridDims[0]/2,cgridDims[0]/2],
                'y':[gridDims[1]/2,cgridDims[1]/2],
                'z':[gridDims[2]/2,cgridDims[2]/2]
                }

    
    m_trans = {'tranx':np.reshape(finalT[0], (cgrids[0],cgrids[1],cgrids[2]), 'F'), 
                'trany':np.reshape(finalT[1], (cgrids[0],cgrids[1],cgrids[2]), 'F'), 
                'tranz':np.reshape(finalT[2], (cgrids[0],cgrids[1],cgrids[2]), 'F')
                }
    
    m_finePerm = {'permx':np.reshape(finePermVec, (grids[0],grids[1],grids[2]), 'F'), 
                  'permy':np.reshape(finePermVec, (grids[0],grids[1],grids[2]), 'F'), 
                  'permz':np.reshape(finePermVec * kz_ratio, (grids[0],grids[1],grids[2]), 'F')}
    
    m_finePoro = {'poro':np.reshape(finePoroVec, (grids[0],grids[1],grids[2]), 'F')}
    
    # Get well info
    nameOfWells = readf.readWellsFromFine(PSS_FILE_SCHED)
    wellType, f_compdat, f_compdat_names, z_endpoints = readf.readWellTypeFromFine(PSS_FILE_SCHED, nameOfWells)
    f_welspecs = readf.readWellInfoFromFine(PSS_FILE_SCHED, 'WELSPECS')
    wellNumOfBlock = readf.genNumOfBlockPerWell(nameOfWells, f_compdat_names, wellType, z_endpoints)
    wellCoord = readf.readWellCoordFromFine(nameOfWells, f_compdat, wellType, wellNumOfBlock)
    wellCoordUpsc,wellCoordUpscAll = upsc.upscaleWellLoc(nameOfWells, wellCoord, fact)
    #print('wellCoordUpsc',wellCoordUpsc)
    #print('wellCoordUpscAll',wellCoordUpscAll)
    #print('f_compdat',f_compdat)
    #print('fact',fact)
    c_compdatl = lgr.generateWell(wellCoord, wellCoordUpsc,wellCoordUpscAll, nameOfWells, f_compdat, fact)
    c_welspecl = lgr.generateCoarseWelspecs(f_welspecs, wellType, c_compdatl, nameOfWells)
    print('c_welspec1', c_welspecl)
    lgr_entries, m_coarsePerm, m_coarsePerm_compute = lgr.generateLGR(wellCoordUpsc, nameOfWells, cgrids, fact, m_finePerm, m_finePoro, m_trans, areas, distances, bisection_tol)
    
    
    f = open("coarsePerm_log.txt", 'a')
    count_bisection = (m_coarsePerm_compute[0]==1).sum()+(m_coarsePerm_compute[1]==1).sum()+(m_coarsePerm_compute[2]==1).sum()
    count_geoavg = (m_coarsePerm_compute[0]==2).sum()+(m_coarsePerm_compute[1]==2).sum()+(m_coarsePerm_compute[2]==2).sum()
    f.write(str(count_bisection)+'\n')
    f.write(str(count_geoavg)+'\n')
    f.close()

    # %% Output
    # =============================================================================
    # We need to output the followings:
    # 1) Upscaled porosity
    # 2) Upscaled transmissibility X, Y, Z
    # 3) base.sched
    # 4) LGR
    # =============================================================================
    writef.mkTransFiles(finalT, path_coarse_model, nameBase = 'TRAN') 
    writef.mkTransFiles(m_coarsePerm, path_coarse_model, nameBase = 'PERM') 
    writef.mkSimpleInputFile(path_coarse_model, writef.PORO_NAME, 
                            writef.PORO_NAME, upPoroVec)
    
    f_lgr = os.path.join(path_coarse_model, "LGR.ECLIN")
    f= open(f_lgr,"w+")
    f.close()
    f_lgr_name = "LGR_names.txt"
    f= open(f_lgr_name,"w+")
    f.close()
    for well in nameOfWells:
        w_lgr = lgr_entries[well]['block_info']
        for block_info in w_lgr:
            f = open(f_lgr, 'a')
            f.write('CARFIN'+'\n')
            f.write(block_info['LGR_entry']+'\n')
            f.write('INCLUDE'+'\n')
            f.write(block_info['block_name']+'x.ECLIN /'+'\n')
            f.write('INCLUDE'+'\n')
            f.write(block_info['block_name']+'y.ECLIN /'+'\n')
            f.write('INCLUDE'+'\n')
            f.write(block_info['block_name']+'z.ECLIN /'+'\n')
            f.write('INCLUDE'+'\n')
            f.write(block_info['block_name']+'poro.ECLIN /'+'\n')
            f.write('ENDFIN'+'\n')
            f1 = open(f_lgr_name, 'a')
            f1.write(block_info['block_name']+'\n')
    f.close()        
    f1.close()
    
    f_lbgsat = os.path.join(path_coarse_model, "LBGSAT_C.ECLIN")
    f= open(f_lbgsat,"w+")
    f.close()
    f= open(f_lbgsat,'a')
    f.write('LBGSAT'+'\n')
    for well in nameOfWells:
        w_lgr = lgr_entries[well]['block_info']
        for block_info in w_lgr:
            for k in range(fact[2]):
                for j in range(fact[1]):
                    for i in range(fact[0]):
                        f.write(block_info['block_name']+' '+str(i+1)+' '+str(j+1)+' '+str(k+1)+' /'+'\n')
    f.write('/')
    f.close()
    
    f_lbpr = os.path.join(path_coarse_model, "LBPR_C.ECLIN")
    f= open(f_lbpr,"w+")
    f.close()
    f= open(f_lbpr,'a')
    f.write('LBPR'+'\n')
    for well in nameOfWells:
        w_lgr = lgr_entries[well]['block_info']
        for block_info in w_lgr:
            for k in range(fact[2]):
                for j in range(fact[1]):
                    for i in range(fact[0]):
                        f.write(block_info['block_name']+' '+str(i+1)+' '+str(j+1)+' '+str(k+1)+' /'+'\n')
    f.write('/')
    f.close()
            
    
    for well in nameOfWells:
        w_lgr = lgr_entries[well]['block_info']
        for block_info in w_lgr:
            f_name_x = os.path.join(path_coarse_model, block_info['block_name']+'x.ECLIN')
            f= open(f_name_x,"w+")
            f.close()
            f= open(f_name_x,"a")
            f.write('PERMX'+'\n')
            for perm in block_info['permx']:
                f.write(str(perm)+'\n')
            f.write('/'+'\n')
            f.close()
            
            f_name_y = os.path.join(path_coarse_model, block_info['block_name']+'y.ECLIN')
            f= open(f_name_y,"w+")
            f.close()
            f= open(f_name_y,"a")
            f.write('PERMY'+'\n')
            for perm in block_info['permy']:
                f.write(str(perm)+'\n')
            f.write('/'+'\n')
            f.close()
            
            f_name_z = os.path.join(path_coarse_model, block_info['block_name']+'z.ECLIN')
            f= open(f_name_z,"w+")
            f.close()
            f= open(f_name_z,"a")
            f.write('PERMZ'+'\n')
            for perm in block_info['permz']:
                f.write(str(perm)+'\n')
            f.write('/'+'\n')
            f.close()
            
            f_name_poro = os.path.join(path_coarse_model, block_info['block_name']+'poro.ECLIN')
            f= open(f_name_poro,"w+")
            f.close()
            f= open(f_name_poro,"a")
            f.write('PORO'+'\n')
            for poro in block_info['poro']:
                f.write(str(poro)+'\n')
            f.write('/'+'\n')
            f.close()
                

    f= open(COARSE_FILE_SCHED,"w+")
    f.close()
    
    udq_entry = ['UDQ','DEFINE RUDRM RWCD+RGCDI+RGCDM /','DEFINE FUDRMTOT FWCD+FGCDI+FGCDM /']
    f= open(COARSE_FILE_SCHED,"a")
    f.write('\n'.join(map(str, udq_entry))+'\n')
    f.write('/'+'\n\n')
    
    f.write('WELSPECL'+'\n')
    for welspecl in c_welspecl:
        f.write(welspecl+'\n')
    f.write('/'+'\n\n')
    
    
    compord_entry = ['COMPORD','I1 INPUT /']
    f.write('\n'.join(map(str, compord_entry))+'\n')
    f.write('/'+'\n\n')
    
    f.write('COMPDATL'+'\n')
    for well in nameOfWells:
        compdatls = c_compdatl[well]
        print(compdatls)
        for compdatl in compdatls:
            f.write(compdatl+'\n')
    f.write('/'+'\n\n')
    
    openFile = open(TWO_PHASE_FILE_SCHED)
    identifier = 'WELLSTRE'
    print_flag = False
    for line in openFile:
        temp = line.rstrip()
        if temp == identifier:
            print_flag = True
        if print_flag:
            f.write(temp+'\n')
            
            
    openFile.close()
    f.close()
    
    
    f= open('Coarse_well_blocks.txt',"w+")
    f.close()
    f= open('Coarse_well_blocks.txt')
    for well in wellCoordUpsc:
        for cord in well:
            f= open('Coarse_well_blocks.txt',"a")
            f.write(' '.join(map(str, cord))+'\n')
    f.close()
    
    

        
    
    
    
        
    
    
    
