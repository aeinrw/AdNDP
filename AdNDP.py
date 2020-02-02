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
    AdNBO(MnSrch, bond, NBOAmnt, DResid)


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
    DUMMY = np.zeros((PBSz, PBSz), dtype=np.float64)
    EiVal = np.zeros(PBSz, dtype=np.float64)
    EiVec = np.eye((PBSz, PBSz), dtype=np.float64)

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
                Subsets(AtBl,AtBlQnt,NCtr,NAt)

            if smode == 2:
                AtBlQnt = 1
                print("Input the given {:3d} centers:".format(NCtr))
                for i in range(NCtr):
                    AtBl[0][i] = int(input())

            print("NCtr = {:3d}, AtBlQnt = {:6d}".format(NCtr, AtBlQnt))

            for k in range(AtBlQnt):
                CBl[0:NCtr] = AtBl[k][0:NCtr]
                BlockDMNAO(CBl,NCtr,DUMMY)
                EigenSystem(DUMMY,BSz,EiVal,EiVec)
                EigenSrt(EiVal,EiVec,BSz)
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
                SortPerl(PrelOcc,PrelVec,PrelCtr,PP,NCtr,prebond)
                print(
                    "********\n\t{:5d} {:3d}-center bonds found with occ = {:f}:".format(PP, NCtr, threshold))
                for i in range(PP):
                    b = prebond[i]
                    print("[{:4d}]**** occ=({:.4f})".format(i, b['occ'][0]))
                    for j in range(b['nc'][0]):
                        print("{:3d}({:.3f})  ".format(
                            b['ctr'][j] + 1, b['vocc'][j]))

            print("Please input the bond indexes to be selected(end with -1):")

            while True:
                i = int(input())
                if i < 0:
                    break
                b = prebond[i]
                DepleteDMNAO(b)
                bond[NBOAmnt] = b
                NBOAmnt += 1

            DResid = DMNAO[BSz][BSz].trace()
            print(
                "\n*** {:5d} bonds found, Density residure = {:f}***".format(NBOAmnt, DResid))
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
    #PrelAccpt = np.zeros(PNMax, dtype=np.int32)

    print(" SortPrel OK")

    for j in range(PP):
        for i in range(PP):
            if PrelOcc[i] > ThrOcc:
                ThrOcc = PrelOcc[i]
                Cnt1 = i

        bond[j]['nc'][0] = NCtr
        i = Cnt1
        bond[j]['occ'][0] = PrelOcc[i]
        PrelOcc[i] = -1.0
        bond[j]['vec'][0:BSz] = PrelVec[0:BSz][i]
        bond[j]['ctr'][0:NCtr] = PrelCtr[0:NCtr][i]
        calbond(bond[j])


def DepleteDMNAO(b):

    global DMNAO

    for i in range(BSz):
        for j in range(BSz):
            DMNAO[i, j] -= b['occ'][0] * b['vec'][j] * b['vec'][i]


def BlockDMNAO(CBl, NCtr, DUMMY):
    flag = np.zeros(PBSz, dtype=np.bool)
    for l in range(NCtr):
        n1 = AtBsRng[CBl[l]][0]
        n2 = AtBsRng[CBl[l]][1]
        for i in range(n1, n2):
            flag[i] = True

    for i in range(BSz):
        for j in range(BSz):
            if (flag[i] and flag[j]):
                DUMMY[i][j] = DMNAO[i][j]
            else:
                DUMMY[i][j] = 0.0


def EigenSystem(a, n, d, v):

    b = np.zeros(n, dtype=np.float64)
    z = np.zeros(n, dtype=np.float64)

    for i in range(n):
        b[i] = a[i, i]
        d[i] = b[i]

    sm = 0.0
    for i in range(1, 101):
        sm = np.sum(np.triu(np.fabs(a[:n, :n]), 1))
        assert sm != 0.0, "sm = 0.0"

        #print("EigenSystem: loop= {:4d},sum= {:lf}".format(i, sm))

        if i < 4:
            tresh = 0.2 * sm / (n ** 2)
        else:
            tresh = 0.0

        for ip in range(n - 1):
            for iq in range((ip + 1), n):
                g = 100.0 * np.fabs(a[ip][iq])
                if i > 4 and (np.fabs(d[ip]) + g == np.fabs(d[ip])) and (np.fabs(d[iq]) + g == fabs(d[iq])):
                    a[ip][iq] = 0.0
                elif np.fabs(a[ip][iq] > tresh):
                    h = d[iq] - d[ip]
                    if np.fabs(h) + g == np.fabs(h):
                        t = a[ip][iq] / h
                    else:
                        theta = 0.5 * h / a[ip][iq]
                        t = 1.0 / (np.fabs(theta) +
                                   np.sqrt(1.0 + theta ** 2))
                        if theta < 0.0:
                            t = -t

                    c = 1.0 / np.sqrt(1 + t ** 2)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a[ip][iq]
                    z[ip] = z[ip] - h
                    z[iq] = z[iq] + h
                    d[ip] = d[ip] - h
                    d[iq] = d[iq] + h
                    a[ip][iq] = 0.0

                    for j in range(ip):
                        g = a[j][ip]
                        h = a[j][iq]
                        a[j][ip] = g - s * (h + g * tau)
                        a[j][iq] = h + s * (g - h * tau)

                    for j in range(ip+1, iqE):
                        g = a[ip][j]
                        h = a[j][iq]
                        a[ip][j] = g - s * (h + g * tau)
                        a[j][iq] = h + s * (g - h * tau)

                    for j in range(iq+1, n):
                        g = a[ip][j]
                        h = a[iq][j]
                        a[ip][j] = g - s * (h + g * tau)
                        a[iq][j] = h + s * (g - h * tau)

                    for j in range(n):
                        g = v[j][ip]
                        h = v[j][iq]
                        v[j][ip] = g - s * (h + g * tau)
                        v[j][iq] = h + s * (g - h * tau)

        for ip in range(n):
            b[ip] = b[ip] + z[ip]
            d[ip] = b[ip]
            z[ip] = 0.0

    print("\n****************\nToo many iterations in jacobi\n*************\n")


def EigenSrt(d, v, n):
    p = 0.0

    for i in range(n - 1):
        k = i
        p = d[i]
        for j in range(i + 1, n):
            if d[j] >= p:
                k = j
                p = d[j]

        if k != i:
            d[k] = d[i]
            d[i] = p
            for j in range(n):
                v[j, k], v[j, i] = v[j, i], v[j, k]


def Subsets(AtBl, NLn, NClmn, NAt):

    if NClmn == NAt:
        for j in range(NAt):
            AtBl[0, j] = j
        AtBl[0, NAt] = -1
        NLn = 1
        return NLn

    if NClmn == 1:
        for i in range(NAt):
            AtBl[i, 0] = i
            AtBl[i, 1] = -1
        NLn = NAt
        return NLn

    Done = 0
    NLn = 1

    for j in range(NClmn):
        AtBl[NLn - 1, j] = j

    while (True):
        for k in range(NClmn):
            if AtBl[NLn - 1, k] == (NAt - NClmn + k):
                Done = k + 1
                break

        if Done == 1:
            return NLn
        if Done == 0:
            NLn += 1
            for k in range(NClmn - 1):
                AtBl[NLn - 1, k] = AtBl[NLn - 2, k]
            if AtBl[NLn - 2, NClmn - 1] < NAt:
                AtBl[NLn - 1, NClmn - 1] = AtBl[NLn - 2, NClmn - 1] + 1
            else:
                NLn += 1
                for k in range(Done - 2):
                    AtBl[NLn - 1, k] = AtBl[NLn - 2, k]
                AtBl[NLn - 1, Done - 2] = AtBl[NLn - 2, Done - 2] + 1
                for k in range(Done - 1, NClmn):
                    AtBl[NLn - 1, k] = AtBl[NLn - 1, k - 1] + 1
                Done = 0


def Input(self):

    # 输入AdNDP.ini文件
    import configparser
    cf = configparser.ConfigParser()
    cf.read(self.__INIFile)
    sections = cf.sections()[0]
    nbofile = cf.get(sections, "nbofile")
    self.NAt = cf.getint(sections, "NAt")
    self.NVal = cf.getint(sections, "NVal")
    self.NTot = cf.getint(sections, "NTot")
    self.BSz = cf.getint(sections, "BSz")
    self.__AtBs[0:self.NAt] = np.array(
        cf.get(sections, "AtBs").split(','), dtype=np.int32)
    self.__Thr[0:self.NAt] = np.array(
        cf.get(sections, "Thr").split(','), dtype=np.float64)
    self.__CMOFile = cf.get(sections, "CMOFile")
    self.__ADNDPFile = cf.get(sections, "AdNDPFile")



if __name__ == '__main__':
    main()