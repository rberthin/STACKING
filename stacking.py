import argparse
import numpy as np
from math import sqrt, acos

NMOLS = 10
NSTEP = 4200
TOL_ANGLE = 20.0     # tolerance angle en degré
TOL_DIST = 5.5    # tolerance distance en angstroms
RAD2DEG = 57.2958

parser = argparse.ArgumentParser(description = '  ')
parser.add_argument('-f', '--pos', default = 'trajectories.xyz',\
                        help = 'name of the file containing the position')
args = parser.parse_args()

print("Need 3 atoms A, B and C to build 2 non colinear vectors AB and AC") 
at1 = input("Name of atom A?\n")  # atom A
at2 = input("Name of atom B?\n")  # atom B
at3 = input("Name of atom C?\n")  # atom C
box = float(input("Box length in same unit as position\n"))
temp = input("To compute the center of mass, which atom to take ?\n")
AT_NAME = temp.split()
 
f1_out = open("number_stacking.dat", "w+")
f1_out.write("# Step    Number of pair\n")
f2_out = open("stacking_stats.dat", "w+")
f2_out.write("# mol1   mol2   distance   angle\n")

mass = []
# On lit le runtime pour avoir les masses de chaque atome
with open('runtime.inpt','r') as run:
    for line in run:
        if (line.lstrip()).startswith("name"):
            name = line.split()[1]
        if (line.lstrip()).startswith("mass"):
            if name in AT_NAME:
                mass.append(float(line.split()[1]))
sum_mass = sum(mass)

with open(args.pos, 'r') as f_in:
    for step in range(NSTEP):
#        f2_out.write("# STEP   "+str(step)+"\n")
        print('STEP', step)
        count_stacking = 0
        X = []
        Y = []
        Z = []
        A = []
        B = []
        C = []
        # On lit la trajectoire pour avoir les positions x y et z de tous les atomes
        line = f_in.readline()  # 2 ligne au début xyz
        ntot = int(line)
        f_in.readline()
        tempx = []
        tempy = []
        tempz = []
        for i in range(ntot):
            line = f_in.readline()
            name = line.split()[0]
            if name in AT_NAME:
                tempx.append(float(line.split()[1]))
                tempy.append(float(line.split()[2]))
                tempz.append(float(line.split()[3]))
                if len(tempx) == NMOLS:    
                    X.append(tempx)
                    Y.append(tempy)
                    Z.append(tempz)
                    tempx = []
                    tempy = []
                    tempz = []

            if (name == at1):
                A.append([ float(line.split()[1]), float(line.split()[2]), \
                           float(line.split()[3]) ]) 
            elif (name == at2):
                B.append([ float(line.split()[1]), float(line.split()[2]), \
                           float(line.split()[3]) ])
            elif (name == at3):
                C.append([ float(line.split()[1]), float(line.split()[2]), \
                           float(line.split()[3]) ])

        # On calcule le centre de masse de chaque molecule
        CMX = []
        CMY = []
        CMZ = []
        for j in range(10):
            cmx = 0
            cmy = 0
            cmz = 0
            for i in range(len(mass)):
                cmx = cmx + (X[i][j]*mass[i])
                cmy = cmy + (Y[i][j]*mass[i])
                cmz = cmz + (Z[i][j]*mass[i])
            CMX.append(cmx/sum_mass)
            CMY.append(cmy/sum_mass)
            CMZ.append(cmz/sum_mass)
        
        coeff_plan = []
    
        for n in range(len(A)):  
            # calcul des vecteurs AB et AC
            vec_AB = [ B[n][0] - A[n][0], B[n][1] - A[n][1], B[n][2] - A[n][2] ]
            vec_AC = [ C[n][0] - A[n][0], C[n][1] - A[n][1], C[n][2] - A[n][2] ]
        
            # produit vectoriel pour trouver le vecteur normal au plan
            temp = np.cross(vec_AB, vec_AC) 
        
            # d tq ax + by + cz + d = 0
            d = -( (temp[0]*A[n][0]) + (temp[1]*A[n][1]) + (temp[2]*A[n][2]) ) 
        
            # pour chaque plan les coefficients a b c et d
            coeff_plan.append([temp[0], temp[1], temp[2], d]) 
        
        DISTANCE_CM = []
    
        # calcul distance cm - cm
        for mol1 in range(len(CMX)):
            for mol2 in range(mol1+1, len(CMX)):
                xx = CMX[mol2]-CMX[mol1]
                yy = CMY[mol2]-CMY[mol1]
                zz = CMZ[mol2]-CMZ[mol1]
                if xx > box/2:
                    xx = xx - box
                elif xx < -box/2:
                    xx = xx + box
                elif yy > box/2:
                    yy = yy - box
                elif yy < -box/2:
                    yy = yy + box
                elif zz > box/2:
                    zz = zz - box
                elif zz < -box/2:
                    zz = zz + box
        
                distance = sqrt( xx**2 + yy**2 + zz**2 )
                if distance <= TOL_DIST:
                    DISTANCE_CM.append([mol1, mol2, distance])
                    # l'angle entre 2 plans est égal à l'angle formé par leurs vecteurs normaux
                    # calcul norme des vecteurs normaux
                    norme1 = sqrt( coeff_plan[mol1][0]**2 + coeff_plan[mol1][1]**2 + \
                             coeff_plan[mol1][2]**2 )
                    norme2 = sqrt( coeff_plan[mol2][0]**2 + coeff_plan[mol2][1]**2 + \
                             coeff_plan[mol2][2]**2 )
                    # calcul du produit scalaire entre les vecteurs
                    produit_scalaire = coeff_plan[mol1][0]*coeff_plan[mol2][0] + \
                                       coeff_plan[mol1][1]*coeff_plan[mol2][1] + \
                                       coeff_plan[mol1][2]*coeff_plan[mol2][2]

                    cos_theta = produit_scalaire / (norme1 * norme2)
                    theta = (acos(cos_theta))*RAD2DEG
                    if theta > 90:
                        theta = 180 - theta
                    if theta <= TOL_ANGLE:
                        count_stacking += 1
                        f2_out.write("{0}   {1}      {2}   {3:.4f}   {4:.4f}\n".format(
                                     step, mol1, mol2, distance, theta))
        f1_out.write("{0}               {1}\n".format(step, count_stacking))                    
f1_out.close()
f2_out.close()
