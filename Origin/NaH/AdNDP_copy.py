import numpy as np

PNMax = 100

de=0

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

NAt = 0
NVal = 0
NTot = 0
BSz = 0
AtBs = None
AtBsRng = None
DMNAO = None
Thr = None
CMOFile = ''
AdNDPFile = ''

# (1,)
BOND = np.dtype([('nc', np.int32, (1,)),
                 ('ctr', np.int32, (NAt,)),
                 ('occ', np.float64, (1,)),
                 ('vec', np.float64, (BSz,)),
                 ('vocc', np.float64, (BSz,))])


def calbond(b: BOND):

    global AtBsRng

    for i in range(b['nc'][0]):

        m = b['ctr'][i]
        n1 = AtBsRng[m, 0]
        n2 = AtBsRng[m, 1]

        b['vocc'][i] = b['vec'][n1:n2].dot(b['vec'][n1:n2]) * b['occ'][0]


def main():
    NAOAO = np.zeros((21, 21), dtype=np.float64)
    Input(NAOAO)

    MnSrch = 0

    NBOOcc = np.zeros(PNMax, dtype=np.float64)
    NBOVec = np.zeros((BSz, PNMax), dtype=np.float64)
    # NBOVecAO = np.zeros((BSz, PNMax), dtype=np.float64)

    DResid = 0.0

    NBOCtr = np.zeros((NAt, PNMax), dtype=np.int32)
    NBOAmnt = 0

    #bond = np.zeros(PNMax, dtype=BOND)

    MnSrch,NBOAmnt,DResid = AdNBO(MnSrch,NBOOcc,NBOVec,NBOCtr,NBOAmnt,DResid)

    with open("result.txt","w") as fp:
        fp.write("--------------DMNAO----------------\n")
        for i in range(21):
            for j in range(21):
                fp.write("{:f} ".format(DMNAO[i][j]))
            fp.write("\n")
        fp.write("---------------NBOOcc--------------\n")
        for i in range(21):
            fp.write("{:f} ".format(NBOOcc[i]))
        fp.write('\n')
        fp.write("---------------NBOVec---------------\n")
        for i in range(21):
            for j in range(21):
                fp.write("{:f} ".format(NBOVec[i][j]))
            fp.write("\n")
        fp.write("---------------NBOCtr---------------\n")
        for i in range(NAt):
            for j in range(21):
                fp.write("{:d} ".format(NBOCtr[i][j]))
            fp.write("\n")
        fp.write("---------------Other---------------\n")
        fp.write("MnSrch={:d}\n".format(MnSrch))
        fp.write("NBOAmnt={:d}\n".format(NBOAmnt))
        fp.write("DResid={:f}\n".format(DResid))




def AdNBO(MnSrch,NBOOcc,NBOVec,NBOCtr,NBOAmnt,DResid):

    global DMNAO
    PrelOcc = np.zeros(PNMax, dtype=np.float64)
    PrelVec = np.zeros((BSz, PNMax), dtype=np.float64)
    PrelCtr = np.zeros((NAt, PNMax), dtype=np.int32)

    Cnt = 0
    PP = 0
    IndS = 0
    IndF = 0
    NCtr = 0
    AtBlQnt = 0
    AtBl = np.zeros((200, NAt), dtype=np.int32)
    CBl = np.zeros(NAt, dtype=np.int32)

    #prebond = np.zeros(PNMax, dtype=BOND)

    threshold = 0.0
    vmax = 0.0
    DUMMY = np.zeros((BSz, BSz), dtype=np.float64)
    EiVal = np.zeros(BSz, dtype=np.float64)
    EiVec = np.eye(BSz, dtype=np.float64)


    with open("out_debug", 'w') as fp:

        print("Welcome to AdNDP-manual program v3.0, revised by Ran-wei")

        NBOAmnt = 0

        for NCtr in range(1,NAt+1):
            threshold = Thr[NCtr-1]
            if threshold == 0.0:
                print("skippint NCtr= ", NCtr)
                continue

            if MnSrch == 0:
                AtBlQnt = Subsets(AtBl, AtBlQnt, NCtr, NAt)

            print("NCtr={:d},AtBlQnt={:d}".format(NCtr, AtBlQnt))

            PP = 0
            Cnt = 1
            while True:
                for k in range(AtBlQnt):
                    CBl[:NCtr] = AtBl[k][:NCtr]

                    BlockDMNAO(CBl, NCtr, DUMMY)



                    value,vector = np.linalg.eig(DUMMY)
                    EiVal,EiVec = np.real(value),np.real(vector)
                    imax = np.argmax(EiVal)
                    valmax = EiVal[imax]
                    vecmax = EiVec[:,imax]



                    if np.fabs(valmax - 2.0) <= threshold:
                        PrelOcc[PP] = valmax
                        PrelVec[:,PP] = vecmax
                        PrelCtr[:NCtr,PP] = CBl[:NCtr]
                        PP += 1
                        fp.write("AtBl ")
                        for j in range(NCtr):
                            fp.write("{:3d}".format(CBl[j]))
                        fp.write("\nEival {:f}\n".format(EiVal[imax]))

                Cnt += 1
                if PP > 0:
                    AtBl[:PP,:NCtr] = PrelCtr[:NCtr,:PP].T
                    AtBlQnt = PP
                    print("***AtBlQnt: {:5d}".format(AtBlQnt))
                    IndF = SortPrel(PrelOcc, PrelVec, PrelCtr, PP, NBOOcc, NBOVec,NBOVec,NCtr,IndF)
                    DepleteDMNAO(NBOOcc, NBOVec, IndS, IndF)

                    DMNAO -= vmax*np.outer(vecmax,vecmax)

                    DResid = DMNAO.trace()
                    print("DResid=",DResid)

                    IndS = IndF
                    NBOAmnt = IndF
                    print("PP={:d},Indf={:d},IndF={:d},NBOAmnt={:d},Cnt={:d}".format(
                        PP, IndS, IndF, NBOAmnt,Cnt))
                    PP = 0
                else:
                    break

                if Cnt > 20:
                    break

    return MnSrch, NBOAmnt, DResid


def SortPrel(PrelOcc, PrelVec, PrelCtr, PP, NBOOcc,NBOVec,NBOCtr,NCtr,IndF):


    Cnt1 = 1
    ThrOcc = 0.0
    PrelAccpt = np.zeros(PNMax, dtype=np.int32)



    print("SortPrel OK")
    PrelAccpt=np.argsort(-PrelOcc)
    PrelOcc_sorted = PrelOcc[PrelAccpt]

    while int(PrelOcc_sorted[Cnt1]*100000)==int(PrelOcc_sorted[Cnt1-1]*100000):
        Cnt1+=1

    # print("Ctr=",PrelCtr[:,PrelAccpt[0]]+1)
    # print("ThrOcc=",PrelOcc[PrelAccpt[0]])

    for j in range(Cnt1):
        i = PrelAccpt[j]
        NBOOcc[IndF] = PrelOcc[i]
        NBOVec[:,IndF] = PrelVec[:,i]
        NBOCtr[:NCtr,IndF] = PrelCtr[:NCtr,i]
        NBOCtr[NCtr,IndF] = -1
        IndF += 1

    return IndF

def DepleteDMNAO(NBOOcc,NBOVec,IndS,IndF):
    global DMNAO

    print("DeleteDMNAO OK")
    for l in range(IndS,IndF):
        a = NBOVec[:,l]
        DMNAO -= NBOOcc[l]*np.outer(a,a)

def BlockDMNAO(CBl, NCtr, DUMMY):
    flag = np.zeros(BSz, dtype=np.bool)
    for l in range(NCtr):
        n1 = AtBsRng[CBl[l],0]
        n2 = AtBsRng[CBl[l],1]
        print("-------->n1=",n1," n2=",n2)
        for i in range(n1, n2):
            flag[i] = True

    for i in range(BSz):
        for j in range(BSz):
            if (flag[i] and flag[j]):
                DUMMY[i,j] = DMNAO[i,j]
            else:
                DUMMY[i,j] = 0.0


def Subsets(AtBl, NLn, NClmn, NAt):

    if NClmn == NAt:
        for j in range(NAt):
            AtBl[0, j] = j
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

    global NAt,NVal,NTot,BSz,AtBs,AtBsRng,Thr,DMNAO,CMOFile,ADNDPFile

    import configparser
    cf = configparser.ConfigParser()
    cf.read(INIFile)
    sections = cf.sections()[0]

    NAt = cf.getint(sections, "NAt")
    NVal = cf.getint(sections, "NVal")
    NTot = cf.getint(sections, "NTot")
    BSz = cf.getint(sections, "BSz")

    nbofile = cf.get(sections, "nbofile")
    CMOFile = cf.get(sections, "CMOFile")
    ADNDPFile = cf.get(sections, "AdNDPFile")

    #AtBs = np.zeros(NAt, dtype=np.int32)
    AtBsRng = np.zeros((NAt, 2), dtype=np.int32)
    DMNAO = np.zeros((BSz, BSz), dtype=np.float64)
    #Thr = np.zeros(NAt, dtype=np.float64)

    AtBs = np.array(
        cf.get(sections, "AtBs").split(','), dtype=np.int32)
    Thr = np.array(
        cf.get(sections, "Thr").split(','), dtype=np.float64)

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
                DMNAO[i][8 * k:8 * k + Cnt] = np.array(line, dtype=np.float64)
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
                NAOAO[i][8 * k:8 * k + Cnt] = np.array(line, dtype=np.float64)
                #print(NAOAO[i][8 * k:8 * k + Cnt])

            dtln = fp.readline()
            dtln = fp.readline()
            dtln = fp.readline()

            Resid -= Cnt
            k += 1

        print("Read DMNAO NAOAO matrixs in nbofile OK!")

if __name__ == '__main__':
    prt=open("print.txt","w")
    main()
    prt.close()
