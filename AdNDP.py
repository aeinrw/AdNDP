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

    # Input(NAOAO)


def AdNBO(MnSrch: int, NBOAmnt: int, DResid: int, bond):

    PrelOcc = np.zeros(PNMax, dtype=np.float64)
    PrelVec = np.zeros((PBSz, PNMax), dtype=np.float64)
    PrelCtr = np.zeros((PNAt, PNMax), dtype=np.int32)

    smode = 0
    #Cnt = 0
    PP = 0
    #IndS = 0
    #IndF = 0
    NCtr = 0
    AtBlQnt = 0
    AtBl = np.zeros((100000, PNAt), dtype=np.int32)
    CBl = np.zeros(PNAt, dtype=np.int32)

    prebond = np.zeros(PNMax, dtype=BOND)

    threshold = 0.0
    #vmax = 0.0
    #DUMMY = np.zeros((PBSz, PBSz), dtype=np.float64)
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
                for i in range(BSz):
                    if EiVal[i] >= threshold:
                        PrelOcc[PP] = EiVal[i]
                        PrelVec[0:BSz][PP] = EiVec[0:BSz][PP]
                        PrelCtr[0:NCtr][PP] = CBl[0:NCtr]
                        PP += 1
                        debug.write("AtBl ")
                        for j in range(NCtr):
                            debug.write("{:3d}".format(CBl[j]))
                        debug.write("\nEiVal {:f}\n\n".format(EiVal[i]))

            if PP <= 0:
                print("No bonds found with occ >= {:f}".format(threshold))
            else:
                # SortPerl()
                print(
                    "********\n\t{:5d} {:3d}-center bonds found with occ = {:f}:".format(PP, NCtr, threshold))
                for i in range(PP):
                    b = prebond[i]
                    print("[{:4d}]**** occ=({:.4f})".fotmat(i, b['occ']))
                    for j in range(b['nc'][0]):
                        print("{:3d}({:.3f})  ".format(
                            b['ctr'][j] + 1, b['vocc'][j]))

            print("Please input the bond indexes to be selected(end with -1):")

            while True:
                i = int(input())
                if i < 0:
                    break
                b = prebond[i]
                # DepleteDMNAO(b)
                bond[NBOAmnt] = b
                NBOAmnt += 1

            # TraceDMNAO()
            print(
                "\n*** {:5d} bonds found, Density residure = {:f}***".fotmat(NBOAmnt, DResid))
            print("0 --- end the AdNDP program\n")
            print("1 --- redo AdNDP search\n")
            print("*** Input you selection(0 or 1): ")

            while True:
                i = int(input())
                if i != 0 and i != 1:
                    print(
                        "\n*** invalid input, reinput your selection please (o or 1):")
                    continue
                else:
                    break

            if i == 1:
                continue
            else:
                break

    return MnSrch, NBOAmnt, DResid


def SortPrel(PrelOcc, PrelVec, PrelCtr, bond, PP, NCtr):

    Cnt1 = 0
    ThrOcc = 0.0
    PrelAccpt = np.zeros(PNMax, dtype=np.int32)

    print(" SortPrel OK")

    for j in range(PP):
        for i in range(PP):
            if PrelOcc[i] > ThrOcc:
                ThrOcc = PrelOcc[i]
                Cnt1 = i

        bond[j]['nc'] = NCtr
        i = Cnt1
        bond[j]['occ'] = PrelOcc[i]
        PrelOcc[i] = -1.0
        bond[j]['vec'][0:BSz] = PrelVec[0:BSz][i]
        bond[j]['ctr'][0:NCtr] = PrelCtr[0:NCtr][i]
        calbond(bond[j])
