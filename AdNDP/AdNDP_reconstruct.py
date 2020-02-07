import numpy as np
from scipy.special import comb

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

    NBO,NBOAmnt,DResid = AdNBO()

    NBOOcc,NBOVec,NBOCtr = NBO


    #bond = np.zeros(PNMax, dtype=BOND)

    # MnSrch,NBOAmnt,DResid = AdNBO(MnSrch,NBOOcc,NBOVec,NBOCtr,NBOAmnt,DResid)

    with open("result.txt","w") as fp:
        fp.write("--------------DMNAO----------------\n")
        for i in range(BSz):
            for j in range(BSz):
                fp.write("{:f} ".format(DMNAO[i][j]))
            fp.write("\n")
        fp.write("---------------NBOOcc--------------\n")
        for i in range(BSz):
            fp.write("{:f} ".format(NBOOcc[i]))
        fp.write('\n')
        fp.write("---------------NBOVec---------------\n")
        for i in range(BSz):
            for j in range(BSz):
                fp.write("{:f} ".format(NBOVec[i][j]))
            fp.write("\n")
        fp.write("---------------NBOCtr---------------\n")
        for i in range(NAt):
            for j in range(BSz):
                fp.write("{:d} ".format(NBOCtr[i][j]))
            fp.write("\n")
        fp.write("---------------Other---------------\n")
        fp.write("NBOAmnt={:d}\n".format(NBOAmnt))
        fp.write("DResid={:f}\n".format(DResid))




def AdNBO():

    global DMNAO

    PrelOcc = np.zeros(BSz, dtype=np.float64)
    PrelVec = np.zeros((BSz, BSz), dtype=np.float64)
    PrelCtr = np.ones((NAt, BSz), dtype=np.int32)
    PrelCtr *= -1

    NBOOcc = np.zeros(BSz, dtype=np.float64)
    NBOVec = np.zeros((BSz, BSz), dtype=np.float64)
    NBOCtr = np.ones((NAt, BSz), dtype=np.int32)
    NBOCtr *= -1

    IndS=0
    IndF=0


    with open("out_debug.txt",'w') as debug:

        print("Welcome to AdNDP-manual program v3.0, revised by Ran-wei")
        NBOAmnt = 0

        for NCtr in range(1,NAt+1):
            threshold = Thr[NCtr-1]
            if threshold == 0.0:
                print("skippint NCtr=",NCtr)
                continue


            AtBl,AtBlQnt = Subsets(NCtr)



            print("NCtr={:d},AtBlQnt={:d}".format(NCtr, AtBlQnt))

            PP=0
            Counter = 1
            while True:
                for k in range(AtBlQnt):
                    CBl = AtBl[k,:] #当前的中心原子组合
                    DUMMY = BlockDMNAO(DMNAO,CBl)

                    value,vector = np.linalg.eig(DUMMY)
                    value,vector = np.real(value),np.real(vector)
                    imax = np.argmax(value)
                    EiVal = value[imax]
                    EiVec = vector[:,imax]

                    if np.fabs(EiVal - 2.0)<=threshold:
                        PrelOcc[PP] = EiVal
                        PrelVec[:,PP] = EiVec
                        PrelCtr[:NCtr,PP]=CBl[:NCtr]
                        PP+=1
                        debug.write("AtBl ")
                        for j in range(NCtr):
                            debug.write("{:d} ".format(CBl[j]))
                        debug.write("\n EiVal {:f}\n".format(EiVal))

                Counter+=1
                if PP>0:
                    AtBl[:PP,:NCtr] = PrelCtr[:NCtr,:PP].T
                    AtBlQnt = PP

                    #SortPrel
                    #-----------------------------
                    NBOOcc[IndF]=PrelOcc[0]
                    NBOVec[:,IndF]=PrelVec[:,0]
                    NBOCtr[:,IndF]=PrelCtr[:,0]
                    IndF+=1
                    #-----------------------------

                    for i in range(IndS,IndF):
                        a = NBOVec[:,i]
                        DMNAO -= NBOOcc[i]*np.outer(a,a)

                    debug.write("trace={:f}\n".format(DMNAO.trace()))

                    IndS=IndF
                    NBOAmnt = IndF
                    print("PP={:d},Indf={:d},IndF={:d},NBOAmnt={:d},Cnt={:d}".format(
                        PP, IndS, IndF, NBOAmnt,Counter))
                    PP=0
                else:
                    break

                if Counter>10:
                    break

    return (NBOOcc,NBOVec,NBOCtr),NBOAmnt,DMNAO.trace()






def Subsets(NCtr):

    global NAt
    AtBlQnt = int(comb(NAt,NCtr))
    AtBl = np.zeros((AtBlQnt,NCtr),dtype=np.int32)

    if NCtr == NAt:
        AtBl[0,:]=np.arange(NAt)
        return AtBl,AtBlQnt

    if NCtr == 1:
        AtBl[:,0]=np.arange(NAt)
        return AtBl,AtBlQnt

    Done = 0
    NLn = 1

    AtBl[NLn-1,:]=np.arange(NCtr)

    while (True):
        for k in range(NCtr):
            if AtBl[NLn - 1, k] == (NAt - NCtr + k):
                Done = k + 1
                break
        if Done == 1:
            return AtBl,AtBlQnt
        if Done == 0:
            NLn += 1
            AtBl[NLn - 1, :NCtr-1] = AtBl[NLn - 2, :NCtr-1]
            if AtBl[NLn - 2, NCtr - 1] < NAt:
                AtBl[NLn - 1, NCtr - 1] = AtBl[NLn - 2, NCtr - 1] + 1
        else:
            NLn += 1
            AtBl[NLn - 1, :Done-2] = AtBl[NLn - 2, :Done-2]
            AtBl[NLn - 1, Done - 2] = AtBl[NLn - 2, Done - 2] + 1
            for k in range(Done - 1, NCtr):
                AtBl[NLn - 1, k] = AtBl[NLn - 1, k - 1] + 1
            Done = 0

    return AtBl,AtBlQnt

def BlockDMNAO(DMNAO,CBl):

    DUMMY = np.zeros(DMNAO.shape,dtype=np.float64)

    flag = np.zeros(BSz,dtype=np.bool)

    for l in range(len(CBl)):
        n1 = AtBsRng[CBl[l],0]
        n2 = AtBsRng[CBl[l],1]
        flag[n1:n2] = True

    for i in range(BSz):
        for j in range(BSz):
            if flag[i] and flag[j]:
                DUMMY[i,j] = DMNAO[i,j]

    return DUMMY





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

    main()