# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 03:15:38 2023

@author: amyfzou
"""
import os
import numpy as np
import itertools
from prep_Coarse_upsc import geoMean
from sympy import *

###########################################################################
def getCoarseFineProp(b_coord, fact, m_finePerm, prop):
    finePerm = m_finePerm[prop]
    b_finePerm = finePerm[int(b_coord[0]*fact[0]-fact[0]):int(b_coord[0]*fact[0]), 
                          int(b_coord[1]*fact[1]-fact[1]):int(b_coord[1]*fact[1]), 
                          int(b_coord[2]*fact[2]-fact[2]):int(b_coord[2]*fact[2])]
    b_finePerm_reshape = np.reshape(b_finePerm, fact[0]*fact[1]*fact[2], 'F')
    
    return b_finePerm_reshape
###########################################################################
def check_neighbor_in_range(block, cgrids):
    check = False
    # check whether the neighbor block is still in domain
    if block[0] > 0 and block[0] <= cgrids[0] and block[1] > 0 and block[1] <= cgrids[1] and block[2] > 0 and block[2] <= cgrids[2]:
        check = True
        
    return check
###########################################################################
def check_neighbor_in_well(block, wellCoordUpsc):
    check = True
    # check whether the neighbor block is a well block
    for i in range(len(wellCoordUpsc)):
        w_blocks = wellCoordUpsc[i]
        if block in w_blocks:
            check = False

    return check
###########################################################################
def getTrans(b_coord, b_n_coord, b_n_tag, m_trans):
    if 'x' in b_n_tag:
        if '-' in b_n_tag:
            trans = m_trans['tranx'][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1]
        if '+' in b_n_tag:
            trans = m_trans['tranx'][b_coord[0]-1, b_coord[1]-1, b_coord[2]-1]
    if 'y' in b_n_tag:
        if '-' in b_n_tag:
            trans = m_trans['trany'][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1]
        if '+' in b_n_tag:
            trans = m_trans['trany'][b_coord[0]-1, b_coord[1]-1, b_coord[2]-1]
    if 'z' in b_n_tag:
        if '-' in b_n_tag:
            trans = m_trans['tranz'][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1]
        if '+' in b_n_tag:
            trans = m_trans['tranz'][b_coord[0]-1, b_coord[1]-1, b_coord[2]-1]
                    
    return trans
###########################################################################
def getInterfaceFinePerm(block_info, b_n_tag, fact):
    if 'x' in b_n_tag:
        m_fp = np.reshape(block_info['permx'], fact, 'F')
        if '-' in b_n_tag:
            finePerm = m_fp[0,:,:]
        if '+' in b_n_tag:
            finePerm = m_fp[fact[0]-1,:,:]
        finePerm_reshape = np.reshape(finePerm, fact[1]*fact[2], 'F')
    if 'y' in b_n_tag:
        m_fp = np.reshape(block_info['permy'], fact, 'F')
        if '-' in b_n_tag:
            finePerm = m_fp[:,0,:]
        if '+' in b_n_tag:
            finePerm = m_fp[:,fact[1]-1,:]
        finePerm_reshape = np.reshape(finePerm, fact[0]*fact[2], 'F')
    if 'z' in b_n_tag:
        m_fp = np.reshape(block_info['permz'], fact, 'F')
        if '-' in b_n_tag:
            finePerm = m_fp[:,:,0]
        if '+' in b_n_tag:
            finePerm = m_fp[:,:,fact[2]-1]
        finePerm_reshape = np.reshape(finePerm, fact[0]*fact[1], 'F')
    
    return finePerm_reshape

###########################################################################
def biSection_eval(interface_finePerm, A, dist_coarse, dist_fine, trans, k_c):
    val = 0
    for perm in interface_finePerm:
        val += 0.001127/((1/(A*2))*(dist_coarse/k_c+dist_fine/perm))
    val = val - trans
    
    return val

###########################################################################
def biSection(tag, interface_finePerm, trans, geoAvg, areas, distances, bisection_tol, N):
    factor = 0.2 # !!!!!!!! hyperparameter
    A = areas[tag[0]]
    dist_coarse = distances[tag[0]][1]
    dist_fine = distances[tag[0]][0]
    
    computation = 1
    
    
    # determine starting bracket
    tryNum = 1
    left = geoAvg
    right = geoAvg
    while True:
        left = left * (1 - factor)
        right = right * (1 + factor)

        left_eval = biSection_eval(interface_finePerm, A, dist_coarse, dist_fine, trans, left)
        right_eval = biSection_eval(interface_finePerm, A, dist_coarse, dist_fine, trans, right)
        
        if np.abs(left_eval) < bisection_tol:
            coarsePerm = left
            break
        
        if np.abs(right_eval) < bisection_tol:
            coarsePerm = right
            break
        
        if left_eval * right_eval < 0:
            break
        else:
            tryNum += 1
        
        if tryNum > 10:
            coarsePerm = geoAvg
            computation = 2
            break
    
    if computation == 1:
        mid = geoAvg
        n_iter = 0
        while True:
            left_eval = biSection_eval(interface_finePerm, A, dist_coarse, dist_fine, trans, left)
            right_eval = biSection_eval(interface_finePerm, A, dist_coarse, dist_fine, trans, right)
            mid_eval = biSection_eval(interface_finePerm, A, dist_coarse, dist_fine, trans, mid)
            
            if np.abs(mid_eval) < bisection_tol:
                error = np.abs(mid_eval)
                coarsePerm = mid
                break
            else:
                if n_iter > N:
                    coarsePerm = geoAvg
                    computation = 2
                    break
                
                # determine new bracket
                if mid_eval * left_eval <= 0:
                    right = mid
                elif mid_eval * right_eval <= 0:
                    left = mid
                mid = (left + right) / 2
                coarsePerm = mid
                
                n_iter += 1
      
    return coarsePerm, computation 
###########################################################################  
def generateLGR(wellCoordUpsc, nameOfWells, cgrids, fact, m_finePerm, m_finePoro, m_trans, areas, distances, bisection_tol):
    m_coarsePerm = [np.zeros((cgrids)) for dims in range(3)]
    m_coarsePerm_comp = [np.zeros((cgrids)) for dims in range(3)]
    lgr = dict()
    
    for i in range(len(nameOfWells)):
        well_coord = wellCoordUpsc[i]
        well_name = nameOfWells[i]
        # well_type = wellType[i]
        # well_ref_depth = float(f_welspecs[i].split()[4]) # this must be the 1st perf depth in fine scale. Took care of somewhere else
        lgr[well_name] = dict()
        lgr[well_name]['block_info'] = []
        for block in range(len(well_coord)):
            block_temp = dict()
            
            block_name = well_name + '_' + str(block+1)
            block_temp['block_name'] = block_name
            block_coord = well_coord[block]
            # if block == 0:
                # lgr[well_name]['well_first_region'] = block_name
                    
            block_coord_repeat = list(map(str, list(itertools.chain.from_iterable(itertools.repeat(a, 2) for a in block_coord))))
            block_temp['LGR_entry'] = (block_name + ' '+ ' '.join(block_coord_repeat) 
                                       + ' ' +str(fact[0])+ ' ' +str(fact[0]) + ' ' +str(fact[0])+ ' ' +str(fact[0]*fact[1]*fact[2])
                                       +' GLOBAL /')
            
            block_temp['permx'] = getCoarseFineProp(block_coord, fact, m_finePerm, prop = 'permx')
            block_temp['permy'] = getCoarseFineProp(block_coord, fact, m_finePerm, prop = 'permy')
            block_temp['permz'] = getCoarseFineProp(block_coord, fact, m_finePerm, prop = 'permz')
            block_temp['poro'] = getCoarseFineProp(block_coord, fact, m_finePoro, prop = 'poro')
            
            lgr[well_name]['block_info'].append(block_temp)
            
            # Calculate neighboring block permx, permy, permz based on trans
            block_neighbor =[[block_coord[0]+1, block_coord[1], block_coord[2]], 
                             [block_coord[0]-1, block_coord[1], block_coord[2]], 
                             [block_coord[0], block_coord[1]+1, block_coord[2]], 
                             [block_coord[0], block_coord[1]-1, block_coord[2]],
                             [block_coord[0], block_coord[1], block_coord[2]+1],
                             [block_coord[0], block_coord[1], block_coord[2]-1]
                            ]
            block_neighbor_tag = ['x+', 'x-', 'y+', 'y-', 'z+', 'z-']
            
            for n in range(len(block_neighbor)):
                b_n_coord = block_neighbor[n]
                b_n_tag = block_neighbor_tag[n]
                check_range = check_neighbor_in_range(b_n_coord, cgrids)
                check_not_wblock = check_neighbor_in_well(b_n_coord, wellCoordUpsc)
                if check_range == True and check_not_wblock == True:
                    b_n_trans = getTrans(block_coord, b_n_coord, b_n_tag, m_trans)
                    
                    if 'x' in b_n_tag:
                        b_n_finePerm = getCoarseFineProp(b_n_coord, fact, m_finePerm, prop = 'permx')
                    if 'y' in b_n_tag:
                        b_n_finePerm = getCoarseFineProp(b_n_coord, fact, m_finePerm, prop = 'permy')
                    if 'z' in b_n_tag:
                        b_n_finePerm = getCoarseFineProp(b_n_coord, fact, m_finePerm, prop = 'permz')
                    b_n_finePerm_geoAvg = geoMean(np.reshape(b_n_finePerm, fact, 'F'))
                    
                    interface_finePerm = getInterfaceFinePerm(block_temp, b_n_tag, fact)
                    
                    [coarsePerm,computation] = biSection(b_n_tag, interface_finePerm, b_n_trans, b_n_finePerm_geoAvg, areas, distances, bisection_tol, N=100)
                    
                    if 'x' in b_n_tag:
                        m_coarsePerm[0][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1] = coarsePerm
                        m_coarsePerm_comp[0][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1] = computation
                    if 'y' in b_n_tag:
                        m_coarsePerm[1][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1] = coarsePerm
                        m_coarsePerm_comp[1][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1] = computation
                    if 'z' in b_n_tag:
                        m_coarsePerm[2][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1] = coarsePerm
                        m_coarsePerm_comp[2][b_n_coord[0]-1, b_n_coord[1]-1, b_n_coord[2]-1] = computation
    
    m_coarsePerm_final = []
    for dim in range(3):
        m_coarsePerm_final.append(np.reshape(m_coarsePerm[dim],cgrids[0]*cgrids[1]*cgrids[2],'F'))
    
                    
    return lgr, m_coarsePerm_final, m_coarsePerm_comp
###########################################################################
def generateWell(wellCoord, wellCoordUpsc,wellCoordUpscAll, nameOfWells, f_compdat, fact):
        well_entry = dict()
        for i in range(len(nameOfWells)):
            well_coord = wellCoordUpsc[i]
            well_coord_all = wellCoordUpscAll[i]
            well_coord_fine = wellCoord[i]
            well_name = nameOfWells[i]
            well_index = []
            for compdat in f_compdat:
                if compdat.split()[0] == well_name:
                    well_index.append(compdat.split()[7])
            well_entry[well_name] = []
            for block in range(len(well_coord)):
                block_name = well_name + '_' + str(block+1)
                print(block_name)
                block_coord = well_coord[block]
                for n in range(len(well_coord_all)):
                    #print('well_coord_all',well_coord_all)
                    if block_coord == well_coord_all[n]:
                        temp_i = np.mod(well_coord_fine[n][0],fact[0])
                        if temp_i == 0:
                            temp_i = fact[0]
                        temp_j = np.mod(well_coord_fine[n][1],fact[1])
                        if temp_j == 0:
                            temp_j = fact[1]
                        temp_k = np.mod(well_coord_fine[n][2],fact[2])
                        if temp_k == 0:
                            temp_k = fact[2]
                            #print('temp_i',temp_i)
                            #print('temp_j',temp_j)
                            #print('temp_k',temp_k)
                        well_entry[well_name].append(' '.join([well_name, block_name, str(temp_i), str(temp_j), str(temp_k), str(temp_k),'1* 1*', str(well_index[n]),'0.3 3* /']))
                
        return well_entry
    
###########################################################################
def generateCoarseWelspecs(f_welspecs, wellType, c_compdatl, nameOfWells):
    c_welspecl = []
    for i in range(len(nameOfWells)):
        w_name = nameOfWells[i]
        w_type = wellType[i]
        w_compdatl = c_compdatl[w_name]
        for f_entry in f_welspecs:
            if w_name in f_entry:
                #print('f_entry', f_entry)
                entry_region = w_compdatl[0].split()[1]
                i_j_loc = ' '.join(w_compdatl[0].split()[2:4])
                #print(i_j_loc)
                #print('entry_region',entry_region)
                w_ref_depth = f_entry.split()[4]
                if w_type == 'vertical_1' or w_type == 'vertical_2' or w_type == 'deviated':
                    c_welspecl.append(w_name +' FIELD '+ entry_region + ' ' + i_j_loc + ' ' + w_ref_depth + ' GAS /')
                if w_type == 'horizontal':
                    c_welspecl.append(w_name +' FIELD '+ entry_region + ' ' + i_j_loc + ' 1* GAS /')
    return c_welspecl
    