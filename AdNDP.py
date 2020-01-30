import numpy as np

PNAt = 128
PNVl = 9
PNTot = 15
PNMax = 1024
PBSz = 1024

'''
! NAt - number of atoms
! NVal - number of valence electronic pairs
! NTot - number of core and valence electronic pairs
! NMax - the maximum number of NBO vectors to be found;
         defines the size of NBOVec , NBOOcc , NMax
! BSz - total number of basis functions;
        defines the size of the density matrix
'''

CMOFile = ''
AdNDPFile = ''
NAt = 0
NVal = 0
NTot = 0
BSz = 0
AtBs = np.zeros(PNAt, dtype=np.int32)
AtBsRng = np.zeros((PNAt, 2), dtype=np.int32)
DMNAO = np.zeros((PBSz, PBSz), dtype=np.float64)
Thr = np.zeros(PNAt, dtype=np.float64)

# (1,)
BOND = np.dtype([('nc', np.int32, (1,)),
                 ('ctr', np.int32, (PNAt,)),
                 ('occ', np.float64, (1,)),
                 ('vec', np.float64, (PBSz,)),
                 ('vocc', np.float64, (PBSz,))])


def calbond(b: BOND):

    global AtBsRng

    for i in range(b['nc'][0]):

        m = b['ctr'][i]
        n1 = AtBsRng[m, 0]
        n2 = AtBsRng[m, 1]

        b['vocc'][i] = b['vec'][n1:n2].dot(b['vec'][n1:n2]) * b['occ'][0]


def main():

    MnSrch = 0

    NBOOcc = np.zeros(PNMax, dtype=np.float64)
    NBOVec = np.zeros((PBSz, PNMax), dtype=np.float64)
    NBOVecAO = np.zeros((PBSz, PNMax), dtype=np.float64)

    DResid = 0.0

    NBOCtr = np.zeros((PNAt, PNMax), dtype=np.int32)
    NBOAmnt = 0

    NAOAO = np.zeros((PBSz, PBSz), dtype=np.float64)
    bond = np.zeros(PNMax, dtype=BOND)

    Input(NAOAO)


def AdNBO(MnSrch: int, NBOAmnt: int, DResid: int, bond):

    PrelOcc = np.zeros(PNMax, dtype=np.float64)
    PrelVec = np.zeros((PBSz, PNMax), dtype=np.float64)
    PrelCtr = np.zeros((PNAt, PNMax), dtype=np.int32)

    smode = 0
    Cnt = 0
    PP = 0
    IndS = 0
    IndF = 0
    NCtr = 0
    AtBlQnt = 0
    AtBl = np.zeros((100000, PNAt), dtype=np.int32)
    CBl = np.zeros(PNAt, dtype=np.int32)

    prebond = np.zeros(PNMax, dtype=BOND)
    b = np.zeros(1, dtype=BOND)

    threshold = 0.0
    vmax = 0.0
    DUMMY = np.zeros((PBSz, PBSz), dtype=np.float64)
    EiVal = np.zeros(PBSz, dtype=np.float64)
    EiVec = np.zeros((PBSz, PBSz), dtype=np.float64)

    with open("out_debug", 'w') as debug:

        print("Welcome to AdNDP-manual program v3.0, revised by Ran-wei")

        NBOAmnt = 0

        while True:
            while True:
                NCtr = int(input(
                    "Input the number of centers of the bond to be found (nc=1,2,3,...):"))
                if NCtr < 1 or NCtr > NAt:
                    print(
                        "Invalid nc number nc={:4d}! Reinit please!".format(NCtr))
                else:
                    break
            threshold = -1.0
            while True:
                threshold = float(input(
                    "Input the threshold value for the {:d}-center bonds(0.0~0.2):".format(NCtr)))
                if threshold < 0 or threshold > 2:
                    print(
                        "Invaild threshold value: {:f}! Reinit please!".format(threshold))
                else:
                    break
            while True:
                smode = int(
                    input("Input the searching mode: 1--direct searching; 2--given centers:"))
                if smode != 1 and smode != 2:
                    print("Invalid input {:d}".format(smode))
                else:
                    break

            if smode == 1:
                print(
                    "Searching for all possible {:d}-center bonds ......".format(NCtr))
                # Subsets(AtBl,...)

            if smode == 2:
                AtBlQnt = 1
                print("Input the given {:3d} centers:".format(NCtr))
                for i in range(NCtr):
                    AtBl[0][i] = int(input())

            print("NCtr = {:3d}, AtBlQnt = {:6d}".format(NCtr, AtBlQnt))

            for k in range(AtBlQnt):
                CBl[0:NCtr] = AtBl[k][0:NCtr]
                # BlockDMNAO()
                # EigenSystem()
                # EIgenSrt()
                for i range(BSz):
                    if EiVal[i] >= threshold:
                        PrelOcc[PP] = EiVal[i]
                        PrelVec[0:BSz][PP] = EiVec[0:BSz][PP]
                        PrelCtr[0:NCtr][PP] = CBl[0:NCtr]
                        PP += 1
                        debug.write("AtBl ")
                        for j in range(NCtr):
                            debug.write("{:3d}".format(CBl[j]))
                        debug.write("\nEiVal \n\n")
