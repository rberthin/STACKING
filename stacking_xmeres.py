import numpy as np

NB_STEP = 4200 
NB_MOL = 10
stp, at1, at2, dist, angle = np.loadtxt('stacking_stats.dat', comments='#', unpack = True)

step = stp.tolist()
atom2 = []
count_step = []
count_xmere = [0, 0, 0]
dim = 0
tri = 0
for i in range(NB_STEP):
    print('STEP', i)
    nb = step.count(i)
    if nb != 0:
        index = step.index(i)
    #    print("index =", index)
        for j in range(index, index+nb):
            atom2.append(at2[j])
    #        print(i, atom2)
        for k in range(NB_MOL) :
            count_step.append(atom2.count(k))
    #        print(count_step)
        atom2 = []
        for n in range(len(count_step)):
            if count_step[n] == 1:
                count_xmere[1] += 1
                dim += 1
            elif count_step[n] == 2:
                count_xmere[2] += 1
                tri += 1
    count_xmere[0] += (NB_MOL -2*dim -3*tri)
    dim = 0
    tri = 0
    #    print(count_xmere[0], count_xmere[1], count_xmere[2]) 
    count_step = []
print("Results : \n")
print("Number of monomers :", count_xmere[0])
print("Average over the simulation :", count_xmere[0]/NB_STEP)
print("Number of dimers :", count_xmere[1])
print("Average over the simulation :", count_xmere[1]/NB_STEP)
print("Number of trimers :", count_xmere[2])
print("Average over the simulation :", count_xmere[2]/NB_STEP)
#print("mono =", count_xmere[0])
