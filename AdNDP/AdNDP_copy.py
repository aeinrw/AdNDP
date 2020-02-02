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

INIFile = './AdNDP.ini'

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
    #NBOVecAO = np.zeros((PBSz, PNMax), dtype=np.float64)

    DResid = 0.0

    NBOCtr = np.zeros((PNAt, PNMax), dtype=np.int32)
    NBOAmnt = 0

    NAOAO = np.zeros((PBSz, PBSz), dtype=np.float64)
    #bond = np.zeros(PNMax, dtype=BOND)

    Input(NAOAO)
    MnSrch,NBOAmnt,DResid = AdNBO(MnSrch,NBOOcc,NBOVec,NBOCtr,NBOAmnt,DResid)


def AdNBO(MnSrch,NBOOcc,NBOVec,NBOCtr,NBOAmnt,DResid):

    PrelOcc = np.zeros(PNMax, dtype=np.float64)
    PrelVec = np.zeros((PBSz, PNMax), dtype=np.float64)
    PrelCtr = np.zeros((PNAt, PNMax), dtype=np.int32)

    Cnt = 0
    PP = 0
    IndS = 0
    IndF = 0
    NCtr = 0
    AtBlQnt = 0
    AtBl = np.zeros((100000, PNAt), dtype=np.int32)
    CBl = np.zeros(PNAt, dtype=np.int32)

    #prebond = np.zeros(PNMax, dtype=BOND)

    threshold = 0.0
    vmax = 0.0
    DUMMY = np.zeros((PBSz, PBSz), dtype=np.float64)
    EiVal = np.zeros(PBSz, dtype=np.float64)
    EiVec = np.eye(PBSz, dtype=np.float64)

    with open("out_debug", 'w') as fp:

        print("Welcome to AdNDP-manual program v3.0, revised by Ran-wei")

        NBOAmnt = 0

        print("NAt= ", NAt)
        print(Thr[0]," and ",Thr[1]," and ",Thr[2])

        for NCtr in range(1,NAt+1):
            threshold = Thr[NCtr-1]
            if threshold == 0:
                print("skippint NCtr= ", NCtr)
                continue

            if (MnSrch == 0):
                AtBlQnt = Subsets(AtBl, AtBlQnt, NCtr, NAt)
            print("NCtr={0:3d},AtBlQnt={1:6d}".format(NCtr, AtBlQnt))
            PP = 0

            Cnt = 1
            while (1):
                for k in range(AtBlQnt):
                    CBl[:NCtr] = AtBl[k][:NCtr]
                    BlockDMNAO(CBl, NCtr, DUMMY)
                    EigenSystem(DUMMY, BSz, EiVal, EiVec)
                    EigenSrt(EiVal, EiVec, BSz)

                    vmax = np.max(EiVal)
                    imax = np.argmax(EiVal)

                    i = imax

                    if np.fabs(vmax - 2.0) <= threshold:
                        PrelOcc[PP] = EiVal[i]
                        PrelVec[:BSz][PP] = EiVal[:BSz][i]
                        PrelCtr[:NCtr][PP] = CBl[:NCtr]
                        PP += 1
                        fp.write("AtBl ")
                        for j in range(NCtr):
                            fp.write("{:3d}".format(CBl[j]))
                        fp.write("\nEival \n".format(EiVal[i]))

                Cnt += 1
                if PP > 0:
                    print("PP=", PP)
                    print("NCtr=",NCtr)
                    print(AtBl[:PP][:NCtr].shape)
                    print(PrelCtr[:NCtr][:PP].shape)
                    AtBl[:PP][:NCtr] = PrelCtr[:NCtr][:PP].T
                    AtBlQnt = PP
                    print("***AtBlQnt: {:5d}".format(AtBlQnt))

                    indF = SortPrel(PrelOcc, PrelVec, PrelCtr, PP, NBOOcc, NBOVec,NBOVec,NCtr,IndF)
                    DepleteDMNAO(NBOOcc, NBOVec, IndS, IndF)
                    DResid = DMNAO[:BSz][:BSz].trace()
                    IndS = IndF
                    NBOAmnt = IndF
                    print("PP={:4d},Indf={:4d},IndF={:4d},NBOAmnt={:4d}".format(
                        PP, IndS, IndF, NBOAmnt))
                    PP = 0
                else:
                    break
                if Cnt > 20:
                    break

    return MnSrch, NBOAmnt, DResid


def SortPrel(PrelOcc, PrelVec, PrelCtr, PP, NBOOcc,NBOVec,NBOCtr,NCtr,IndF):

    Cnt1 = 0
    ThrOcc = 0.0
    PrelAccpt = np.zeros(PNMax, dtype=np.int32)

    print(" SortPrel OK")
    for i in range(PP):
        if (int)(PrelOcc[i] * 100000) > (int)(ThrOcc * 100000):
            ThrOcc = (int)(PrelOcc[i] * 100000) / 100000.0
            Cnt1 = 0
            PrelAccpt[Cnt1] = i
            Cnt1 += 1
            print("A Ctr= ", PrelCtr[:NCtr][i] + 1)
            print(" ***ThrOcc= {:.5f}".format(ThrOcc))
        elif (int)(PrelOcc[i] * 100000) == (int)(ThrOcc * 100000):
            PrelAccpt[Cnt1] = i
            Cnt1 += 1
            print("A Ctr= ", PrelCtr[:NCtr][i] + 1)
            print(" ***ThrOcc= {:.5f}".format(ThrOcc))

    for j in range(Cnt1):
        i = PrelAccpt[j]
        NBOOcc[IndF] = PrelOcc[i]
        NBOVec[:BSz][IndF] = PrelVec[:BSz][i]
        NBOCtr[:NCtr][IndF] = PrelCtr[:NCtr][i]
        NBOCtr[NCtr][IndF] = -1
        IndF += 1

    return IndF


def DepleteDMNAO(NBOOcc,NBOVec,IndS,IndF):

    global DMNAO

    print("DeleteDMNAO OK")
    for l in range(IndS,IndF):
        for i in range(BSz):
            for j in range(BSz):
                DMNAO[i, j] -= NBOOcc[l]*NBOVec[j,l]*NBOVec[i,l]


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
                DUMMY[i, j] = DMNAO[i, j]
            else:
                DUMMY[i, j] = 0.0


def EigenSystem(a, n, d, v):

    b = np.zeros(n, dtype=np.float64)
    z = np.zeros(n, dtype=np.float64)

    for i in range(n):
        b[i] = a[i, i]
        d[i] = b[i]

    sm = 0.0
    for i in range(1, 101):
        sm = np.sum(np.triu(np.fabs(a[:n, :n]), 1))
        if sm == 0.0:
            return 0

        #print("EigenSystem: loop= {:4d},sum= {:f}".format(i, sm))

        if i < 4:
            tresh = 0.2 * sm / (n ** 2)
        else:
            tresh = 0.0

        for ip in range(n - 1):
            for iq in range((ip + 1), n):
                g = 100.0 * np.fabs(a[ip][iq])
                if i > 4 and (np.fabs(d[ip]) + g == np.fabs(d[ip])) and (np.fabs(d[iq]) + g == np.fabs(d[iq])):
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

                    for j in range(ip+1, iq):
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

    print(NLn,NClmn,NAt)

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
                #print("!!",NLn-1,Done-2)
                AtBl[NLn - 1, Done - 2] = AtBl[NLn - 2, Done - 2] + 1
                for k in range(Done - 1, NClmn):
                    AtBl[NLn - 1, k] = AtBl[NLn - 1, k - 1] + 1
                Done = 0

    return NLn


def Input(NAOAO):

    global NAt,NTot,BSz,AtBs,Thr,CMOFile,AdNDPFile

    import configparser
    cf = configparser.ConfigParser()
    cf.read(INIFile)
    sections = cf.sections()[0]
    nbofile = cf.get(sections, "nbofile")
    NAt = cf.getint(sections, "NAt")
    NVal = cf.getint(sections, "NVal")
    NTot = cf.getint(sections, "NTot")
    BSz = cf.getint(sections, "BSz")
    AtBs[0:NAt] = np.array(
        cf.get(sections, "AtBs").split(','), dtype=np.int32)
    Thr[0:NAt] = np.array(
        cf.get(sections, "Thr").split(','), dtype=np.float64)
    CMOFile = cf.get(sections, "CMOFile")
    ADNDPFile = cf.get(sections, "AdNDPFile")

    j = 0
    for i in range(NAt):
        AtBsRng[i][0] = j
        j += AtBs[i]
        AtBsRng[i][1] = j

    assert BSz == j, "Size of basis set error!! {:5d} != {:5d}".format(
        BSz, j)

    # ??nbo.log??
    # try:
    #     fp = open(nbofile, 'r')
    # except IOError as err:
    #     print('File Error:' + str(err))

    dtln = ""
    with open("./"+nbofile, 'r') as fp:
        while dtln[:20] != " NAO density matrix:":
            dtln = fp.readline()
        dtln = fp.readline()
        dtln = fp.readline()
        dtln = fp.readline()

        Resid = BSz
        k = 0
        Cnt = 0
        while Resid > 0:

            if Resid < 8:
                Cnt = Resid
            else:
                Cnt = 8

            for i in range(BSz):

                line = fp.readline()[17:-1].replace('  ', ' ').split(' ')
                if line[0] == '':
                    line.pop(0)
                DMNAO[i][8 * k:8 * k +
                         Cnt] = np.array(line, dtype=np.float64)
                #print(DMNAO[i][8 * k:8 * k + Cnt])

            dtln = fp.readline()
            dtln = fp.readline()
            dtln = fp.readline()

            Resid -= Cnt
            k += 1

        fp.seek(0)

        while dtln[:22] != " NAOs in the AO basis:":
            dtln = fp.readline()
        dtln = fp.readline()
        dtln = fp.readline()
        dtln = fp.readline()

        Resid = BSz
        k = 0
        Cnt = 0
        while Resid > 0:

            if Resid < 8:
                Cnt = Resid
            else:
                Cnt = 8

            for i in range(BSz):

                line = fp.readline()[17:-1].replace('  ', ' ').split(' ')
                if line[0] == '':
                    line.pop(0)
                NAOAO[i][8 * k:8 * k +
                         Cnt] = np.array(line, dtype=np.float64)
                #print(NAOAO[i][8 * k:8 * k + Cnt])

            dtln = fp.readline()
            dtln = fp.readline()
            dtln = fp.readline()

            Resid -= Cnt
            k += 1

        print("Read DMNAO NAOAO matrixs in nbofile OK!")


if __name__ == '__main__':
    main()
