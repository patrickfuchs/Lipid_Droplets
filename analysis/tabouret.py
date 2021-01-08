#!/usr/bin/python
#-*- coding: utf-8 -*-

#Amelie Bacle 
#25/01/2016

import numpy as np
import argparse 
import mdtraj as md
import sys
import operator
import scipy.spatial 
import random

MASK_FORK = [1, 0, 2]
MASK_CHAIR = [1, 0, 2]
MASK_TRIDENT = [3, 0, 0]
MASK_T = [0, 2, 1]
MASK_HAND = [0, 3, 0]
MASK_FENWICK = [1, 2, 0]

CONF_NAME = ["Trident", "Chair", "Fork", "T", "Hand", "Stacker", "Other"]


#Change this parameter to be more or less permissive on the conformation
thre = 0.3

def get_traj(traj, atom_select):
    top = traj.topology
    index_atom = top.select(atom_select)
    traj_atom = traj.atom_slice(index_atom)
    return traj_atom

def get_traj_residue(t, list_res):
    sel = "resname TO and ("
    for num_res in xrange(0, len(list_res)):
        if num_res == len(list_res)-1:
            sel += "resid %i "%list_res[num_res]
        else:
            sel += "resid %i or "%list_res[num_res]
    sel += ")"
    return get_traj(t, sel)
    
def get_resid(top, resname):
    table, bonds = top.to_dataframe()
    mask_res = table['resName'] == resname
    id_res = table['resSeq'][mask_res]
    #id_res = nAtom * nRes in length, need to set it to get unique values
    return list(set(id_res))
    
    
def get_vector(coord_start, coord_end):
    dx = coord_end[0] - coord_start[0]
    dy = coord_end[1] - coord_start[1]
    dz = coord_end[2] - coord_start[2]
    vec = [dx, dy, dz]
    norm = np.linalg.norm(vec)
    return vec /norm
    
def get_TO_vector(TO_coord, frame, nb_TO, model): 
    if model != "cg":
        #The First array (first chain), all frames, a TO in particular, the vector  
        sn1 = get_vector(TO_coord[frame, nb_TO*6 + 0, :], TO_coord[frame, nb_TO*6 + 1, :])
#        print TO_coord[frame, nb_TO*6 + 0, :], TO_coord[frame, nb_TO*6 + 1, :]
        sn2 = get_vector(TO_coord[frame, nb_TO*6 + 2, :], TO_coord[frame, nb_TO*6 + 3, :])
#        print TO_coord[frame, nb_TO*6 + 2, :], TO_coord[frame, nb_TO*6 + 3, :]
        sn3 = get_vector(TO_coord[frame, nb_TO*6 + 4, :], TO_coord[frame, nb_TO*6 + 5, :])
#        print TO_coord[frame, nb_TO*6 + 4, :], TO_coord[frame, nb_TO*6 + 5, :]
    else:
        #The First array (first chain), all frames, a TO in particular, the vector  
        sn1 = get_vector(TO_coord[frame, nb_TO*6 + 0, :], TO_coord[frame, nb_TO*6 + 3, :])
#        print TO_coord[frame, nb_TO*6 + 0, :], TO_coord[frame, nb_TO*6 + 3, :]
        sn2 = get_vector(TO_coord[frame, nb_TO*6 + 1, :], TO_coord[frame, nb_TO*6 + 4, :])
#        print TO_coord[frame, nb_TO*6 + 1, :], TO_coord[frame, nb_TO*6 + 4, :]
        sn3 = get_vector(TO_coord[frame, nb_TO*6 + 2, :], TO_coord[frame, nb_TO*6 + 5, :])
#        print TO_coord[frame, nb_TO*6 + 2, :], TO_coord[frame, nb_TO*6 + 5, :]
    return np.array([sn1, sn2, sn3])
    
def print_draw_VMD(vector):
    print "draw line {0 0 0} {%2.2f %2.2f %2.2f}" %(
            vector[0]*10, vector[1]*10, vector[2]*10)
        
    
def is_trident(mask):
    return mask == MASK_TRIDENT
    
def dist_to_trident(vec):
    return scipy.spatial.distance.euclidean(vec, [1, 1, 1])
    
def is_fork(mask, sn1sn3):
    return mask == MASK_FORK and sn1sn3 > 0 
    
def dist_to_fork(vec):
    return scipy.spatial.distance.euclidean(vec,  [-1, -1, 1])
    
def is_chair(mask, sn1sn3):
    return mask == MASK_CHAIR and sn1sn3 < 0
    
def dist_to_chair(vec):
    chair1 = [-1, 1, -1]
    chair2 = [1, -1, -1]
    d_chair1 = scipy.spatial.distance.euclidean(vec, chair1)
    d_chair2 = scipy.spatial.distance.euclidean(vec, chair2)
    return min([d_chair1, d_chair2])
    
def is_t(mask):
    return mask == MASK_T
    
def dist_to_t(vec):
    t1 = [-1, 0, 0]
    t2 = [0, -1, 0]
    t3 = [0, 0, -1]
    d_t1 = scipy.spatial.distance.euclidean(vec, t1)
    d_t2 = scipy.spatial.distance.euclidean(vec, t2)
    d_t3 = scipy.spatial.distance.euclidean(vec, t3)
    return min([d_t1, d_t2, d_t3])
    
def is_hand(mask):
    return mask == MASK_HAND
    
def dist_to_hand(vec):
    return scipy.spatial.distance.euclidean(vec, [0, 0, 0])
    
def is_fenwick(mask):
    return mask == MASK_FENWICK
    
def dist_to_fenwick(vec):
    fen1 = [1, 0, 0]
    fen2 = [0, 1, 0]
    fen3 = [0, 0, 1]
    d_fen1 = scipy.spatial.distance.euclidean(vec, fen1)
    d_fen2 = scipy.spatial.distance.euclidean(vec, fen2)
    d_fen3 = scipy.spatial.distance.euclidean(vec, fen3)
    return min([d_fen1, d_fen2, d_fen3])
    
def which_closer_conf(vec):
    d_trident = dist_to_trident(vec)
    d_chair = dist_to_chair(vec)
    d_fork = dist_to_fork(vec)
    d_t = dist_to_t(vec)
    d_hand = dist_to_hand(vec)
    d_fen = dist_to_fenwick(vec)
    dists_to_conf = [d_trident, d_chair, d_fork, d_t, d_hand, d_fen]
    min_index, min_value = min(enumerate(dists_to_conf), key=operator.itemgetter(1))
    return (min_index, min_value)
    
def is_POPC(top, popc_name):
    table, bonds = top.to_dataframe()
    residues_name = table['resName']
    for res in residues_name:
        if res == popc_name:
            return True
    return False
    
def get_z_coord(t, atom_name, res_name):
    carbon = get_traj(t, "resname %s and name %s" %(res_name, atom_name))
    z_carbon = carbon.xyz[:,:,2]
    return z_carbon
    
def get_inter(z_carbon):
    """From a trjactory returns the min and max z coordinate of the atom from a specific
    residue. 
    Returns : Two arrays of float, shape = nFrame. The max or min of the atom for each frame
    """
    max_inter = np.amax(z_carbon, axis=1)
    min_inter = np.amin(z_carbon, axis=1)
    return (max_inter, min_inter)

def last_atom_select(atom, sns):
    return atom == sns[-1][-1]
    
def print_conf(trident, chair, fork, t, hand, fenwick, nTO_total):
    pr = "Trident   %2.1f%% \n" %(trident / float(nTO_total) * 100)
    pr += "Chair   %2.1f%% \n" %(chair / float(nTO_total) * 100)
    pr += "Fork   %2.1f%% \n" %(fork / float(nTO_total) * 100)
    pr += "T   %2.1f%% \n" %(t / float(nTO_total) * 100)
    pr += "Right Hand   %2.1f%% \n" %(hand / float(nTO_total) * 100)
    pr += "Fenwick   %2.1f%% \n" %(fenwick / float(nTO_total) * 100)
    print pr
    
def print_inter(trident_inter_P, chair_inter_P, fork_inter_P, ti_inter_P, hand_inter_P,fenwick_inter_p):
    pr =  "%2.1f%% of the trident TO are at the interface\n" %trident_inter_P
    pr += "%2.1f%% of the chair TO are at the interface\n" %chair_inter_P
    pr += "%2.1f%% of the fork TO are at the interface\n" %fork_inter_P
    pr += "%2.1f%% of the T TO are at the interface\n" %ti_inter_P
    pr += "%2.1f%% of the right hand TO are at the interface\n" %hand_inter_P
    pr += "%2.1f%% of the fenwick TO are at the interface\n" %fenwick_inter_P
    print pr
    
def get_mask_class(inner1, inner2, inner3, thre):
    class1 = [is_class1(inner1, thre), is_class1(inner2, thre), is_class1(inner3, thre)]
    class2 = [is_class2(inner1, thre), is_class2(inner2, thre), is_class2(inner3, thre)]
    class3 = [is_class3(inner1, thre), is_class3(inner2, thre), is_class3(inner3, thre)]
    mask = [class1.count(True), class2.count(True), class3.count(True)]
    return mask

def is_class1(inner, thre): 
    return 1-thre < inner < 1
    
def is_class2(inner, thre):
    return -thre/2 < inner < thre/2
    
def is_class3(inner, thre):
    return -1 < inner < -1+thre

def get_frame_divided(t, nDiv):
    nFrames = t.n_frames
    length_div = nFrames / nDiv
    shuffle_frames = range(0, nFrames)
    random.shuffle(shuffle_frames)
    frame_div = []
    for div in range(0, nDiv):
        start = (div * length_div) 
        if div == nDiv - 1 :
            frame_div.append(shuffle_frames[start:])
        else :
            end = (div + 1) * length_div 
            frame_div.append(shuffle_frames[start:end])
    return frame_div

def which_block(frame, frame_blocks):
    for nBlock in xrange(0, len(frame_blocks)):
        if frame in frame_blocks[nBlock]:
            return nBlock
    
def write_csv(filename, conf_name, m_count, st_count, m_inter, st_inter):
    filin = open(filename, "w")
    nb_conf = len(conf_name)
    for i, conf in enumerate(conf_name):
        filin.write("%s;%2.1f;%2.1f\n" %(conf, m_count[i], st_count[i]))
    filin.write("\n")
    filin.write("Interface\n")
    for i, conf in enumerate(conf_name):
        filin.write("%s;%2.1f;%2.1f\n" %(conf, m_inter[i], st_inter[i]))
    filin.close()
#########################################################################################

parser = argparse.ArgumentParser(description = 'Arguments for triolein conformation clustering')
parser.add_argument('-f', action = 'store', dest = 'xtc', help = 'Trajectory file')
parser.add_argument('-s', action = 'store', dest = 'gro', help = 'Structure File')
parser.add_argument('-inter', action = 'store', dest = 'inter', help = 'Interface')
parser.add_argument('-d', action = 'store', dest = 'nBlocks', 
        help = 'Number of blocks for block averaging')
parser.add_argument('-m', action = 'store', dest = 'model', help = 'Interface')
parser.add_argument('-o', action = 'store', dest = 'output', help = 'Interface')
args = parser.parse_args()


if args.inter == "True" :
    inter = True
else :
    inter = False
    
nBlocks = int(args.nBlocks)

if args.model == "aa" : 
    t = md.load_xtc(args.xtc, top = args.gro)
    sns = [["C5", "C12"], ["C26", "C33"], ["C47", "C54"]]
    z_central_TO = get_z_coord(t, "C22", "TO")
elif args.model == "cg" : 
    t = md.load(args.xtc, top = args.gro)
    sns = [["ES1", "D2A"], ["ES2", "D2B"], ["ES3", "D2C"]]
    #sns = [["ES1", "C4A"], ["ES2", "C4B"], ["ES3", "C4C"]]
    z_central_TO = get_z_coord(t, "GLY", "TO")
elif args.model == "sdk" : 
    t = md.load(args.xtc, top = args.gro)
    sns = [["ETO1", "C13"], ["ETO2", "C23"], ["ETO3", "C33"]]
    #sns = [["ESTA", "C13"], ["ESTB", "C23"], ["ESTC", "C33"]]
    z_central_TO = get_z_coord(t, "GL", "TO")

#Only works for POPC and AA model
if args.model == "aa" :
    if is_POPC(t.top, "POPC"): 
        popc_name = "POPC"
        LD = True
    elif is_POPC(t.top, "POP"):
        popc_name = "POP"
        LD = True
    else :
        LD = False
elif args.model == "sdk":
    LD = is_POPC(t.top, "POPB")
elif args.model == "cg":
    LD = is_POPC(t.top, "POPC")

print "LD ? %s" %LD

if inter == True : 
#    print "Interface"
    if LD : 
        if args.model == "aa":
            z_carbon_inter = get_z_coord(t, "C50", popc_name)
        elif args.model == "cg" :
            z_carbon_inter = get_z_coord(t, "C4A", "POPC")
        elif args.model == "sdk" :
            print "C25 of POPBX"
            z_carbon_inter = get_z_coord(t, "C25", "POPB")
    else : 
        z_carbon_inter = z_central_TO
    (max_inter, min_inter) = get_inter(z_carbon_inter)


atom_select = "resname TO and ("   
for sn in sns:
    for atom in sn :
        if last_atom_select(atom, sns):
            atom_select += "name %s" %atom
        else:
            atom_select += "name %s or " %atom
        
atom_select += ")"

TO_coord = get_traj(t, atom_select).xyz


trident = 0
trident_inser = 0 
chair = 0
fork = 0
ti = 0
hand = 0 
fenwick = 0
if inter == True : 
    trident_inter = 0
    chair_inter = 0
    fork_inter = 0
    ti_inter = 0
    hand_inter = 0 
    fenwick_inter = 0 


count_confs2 = np.zeros((nBlocks,7))
confs_inter2 = np.zeros((nBlocks,7))
resid_conf = []

dist_to_conf2 = [[], [], [], [], [], []]

zs = []
nFrames = t.n_frames
blocks_frame = get_frame_divided(t, nBlocks)
idTO = get_resid(t.top, "TO")
nTO = len(idTO)
cpt = 0
#Shape 3(sns), nFrames, 3(x,y,z)
for frame in range(0, nFrames):
    block = which_block(frame, blocks_frame)
    trident_frame = []
    chair_frame = []
    fork_frame = []
    ti_frame = []
    hand_frame = []
    fenwick_frame = []
    for i, TO in enumerate(idTO):
        v_sns_TO = get_TO_vector(TO_coord, frame, i, args.model)
        sn1 = v_sns_TO[0]
        sn2 = v_sns_TO[1]
        sn3 = v_sns_TO[2]        
        sn1sn2 = np.inner(sn1, sn2)
        sn2sn3 = np.inner(sn2, sn3)
        sn1sn3 = np.inner(sn1, sn3)
        inner_sns = [sn1sn2, sn2sn3, sn1sn3]
        (conf, dist) = which_closer_conf(inner_sns)
        if dist < 0.7 :
            count_confs2[block][conf] += 1
            dist_to_conf2[conf].append(dist)
        else :
            count_confs2[block][-1] += 1
        if inter == True : 
            z = z_central_TO[frame, i]
            dist_z = min([max_inter[frame] - z, z - min_inter[frame]])
            #print dist_z
            if dist_z < 0.3 :
                cpt += 1
                if TO not in resid_conf :
                    resid_conf.append(TO)
                if dist < 0.7 :
                    confs_inter2[block][conf] += 1
                else:
                    confs_inter2[block][-1] += 1
                    #print "%s %i @ frame %i with dist of %2.1f from inter" \
                           # %(CONF_NAME[conf], i, frame, dist_z)
        #if dist < 0.1 :
            #print inner_sns, i
            #print v_sns_TO
            #print "%s %i @ frame %i with dist of %2.1f from inter" \
            #    %(CONF_NAME[conf], i+1001, frame, dist)


print resid_conf

nTO_total = float(nTO * nFrames) / nBlocks
count_P = count_confs2 / nTO_total * 100 
means_count_P = np.mean(count_P, axis = 0)
st_count_P = np.std(count_P, axis = 0)
#Reshape to broadcast the following division
tot_inter = np.sum(confs_inter2, axis = 1).reshape((nBlocks,1))
inter_P = confs_inter2 / tot_inter * 100 
means_inter_P = np.mean(inter_P, axis = 0)
st_inter_P = np.std(inter_P, axis = 0)
if args.model != "aa":
    write_csv(args.output, CONF_NAME, means_count_P, st_count_P, means_inter_P, st_inter_P)
# for i in xrange(0, len(dist_to_conf2)):
#     print "%s : %2.1f%% +/-  %2.1f%% total with %2.1f%% +/- %2.1f%% interface" \
#         %(CONF_NAME[i], means_count_P[i], st_count_P[i], means_inter_P[i], st_inter_P[i])  
#else: 
print "All trioleins"
print nTO_total 
print count_confs2
print "Triolein at the interface"
print cpt
print tot_inter
print confs_inter2       
     

# for nConf in xrange(0, len(dist_to_conf1)):
#     f = open("%s_dist.txt" %CONF_NAME[nConf], 'w')
#     for i in xrange(0, len(dist_to_conf1[nConf])):
#         dist_str = str(dist_to_conf1[nConf][i])
#         f.write("%s\n" %dist_str)
#     f.close()
    
