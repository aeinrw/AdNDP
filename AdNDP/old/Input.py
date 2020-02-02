import numpy as np


class AdNDP:

    # 宏变量
    __PNAt = 128
    __PNVl = 9
    __PNTot = 15
    __PNMax = 4096
    __PBSz = 1024

    # 全局变量
    __INIFile = ""
    __CMOFile = ""
    __ADNDPFile = ""
    NAt = 0
    NVal = 0
    NTot = 0
    BSz = 0
    __AtBs = np.zeros(shape=__PNAt, dtype=np.int32)
    __AtBsRng = np.zeros(shape=(__PNAt, 2), dtype=np.int32)
    __DMNAO = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)
    __Thr = np.zeros(shape=__PNAt, dtype=np.float64)

    # main函数变量
    MnSrch = 0
    NBOOcc = np.zeros(shape=__PNMax, dtype=np.float64)
    NBOVec = np.zeros(shape=(__PBSz, __PNMax), dtype=np.float64)
    DResid = 0.0
    NBOCtr = np.zeros(shape=(__PNAt, __PNMax), dtype=np.int32)
    NBOAmnt = 0
    __NAOAO = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)

    def __init__(self, iniFile):

        self.__INIFile = iniFile

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

        j = 0
        for i in range(self.NAt):
            self.__AtBsRng[i][0] = j
            j += self.__AtBs[i]
            self.__AtBsRng[i][1] = j

        assert self.BSz == j, "Size of basis set error!! {:5d} != {:5d}".format(
            self.BSz, j)

        # 输入nbo.log文件
        # try:
        #     fp = open(nbofile, 'r')
        # except IOError as err:
        #     print('File Error:' + str(err))

        dtln = ""
        with open("../"+nbofile, 'r') as fp:
            while dtln[:20] != " NAO density matrix:":
                dtln = fp.readline()
            dtln = fp.readline()
            dtln = fp.readline()
            dtln = fp.readline()

            Resid = self.BSz
            k = 0
            Cnt = 0
            while Resid > 0:

                if Resid < 8:
                    Cnt = Resid
                else:
                    Cnt = 8

                for i in range(self.BSz):

                    line = fp.readline()[17:-1].replace('  ', ' ').split(' ')
                    if line[0] == '':
                        line.pop(0)
                    self.__DMNAO[i][8 * k:8 * k +
                                    Cnt] = np.array(line, dtype=np.float64)
                    #print(self.__DMNAO[i][8 * k:8 * k + Cnt])

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

            Resid = self.BSz
            k = 0
            Cnt = 0
            while Resid > 0:

                if Resid < 8:
                    Cnt = Resid
                else:
                    Cnt = 8

                for i in range(self.BSz):

                    line = fp.readline()[17:-1].replace('  ', ' ').split(' ')
                    if line[0] == '':
                        line.pop(0)
                    self.__NAOAO[i][8 * k:8 * k +
                                    Cnt] = np.array(line, dtype=np.float64)
                    #print(self.__NAOAO[i][8 * k:8 * k + Cnt])

                dtln = fp.readline()
                dtln = fp.readline()
                dtln = fp.readline()

                Resid -= Cnt
                k += 1

            print("Read DMNAO NAOAO matrixs in nbofile OK!")

    def AdNBO(self):

        PrelOcc = np.zeros(shape=self.__PNMax, dtype=np.float64)
        PrelVec = np.zeros(shape=(self.__PBSz, self.__PNMax), dtype=np.float64)
        PrelCtr = np.zeros(shape=(self.__PNAt, self.__PNMax), dtype=np.int32)

        Cnt = 0
        PP = 0
        IndS = 0
        IndF = 0
        NCtr = 0
        AtBlQnt = 0
        AtBl = np.zeros(shape=(100000, self.__PNAt), dtype=np.int32)
        CBl = np.zeros(shape=self.__PNAt, dtype=np.int32)
        threshold = 0.0
        vmax = 0.0
        DUMMY = np.zeros(shape=(self.__PBSz, self.__PBSz), dtype=np.float64)
        EiVal = np.zeros(shape=self.__PBSz, dtype=np.float64)
        EiVec = np.zeros(shape=(self.__PBSz, self.__PBSz), dtype=np.float64)

        with open("out_debug", 'w') as fp:

            for NCtr in range(self.NAt):

                threshold = self.__Thr[NCtr - 1]
                if threshold == 0:
                    print("skippint NCtr=", NCtr)
                    continue

                if (self.MnSrch == 0):
                    # Subsets()
                    pass
                print("NCtr={0:3d},AtBlQnt={1:6d}".format(NCtr, AtBlQnt))
                PP = 0

                Cnt = 1
                while (1):
                    for k in range(AtBlQnt):
                        CBl[:NCtr] = AtBl[k][:NCtr]
                        # BlockDMNAO()
                        # EigenSystem()
                        # EigenSrt()

                        vmax = np.max(EiVal)
                        imax = np.argmax(EiVal)

                        i = imax

                        if np.fabs(vmax - 2.0) <= threshold:
                            PrelOcc[PP] = EiVal[i]
                            PrelVec[:self.BSz][PP] = EiVal[:self.BSz][i]
                            PrelCtr[:NCtr][PP] = CBl[:NCtr]
                            PP += 1
                            fp.write("AtBl")
                            fp.write(CBl)
                            fp.write("\nEival \n", EiVal[i])

                    if PP > 0:
                        AtBl[:PP][:NCtr] = PrelCtr[:NCtr][:PP].T
                        AtBlQnt = PP
                        print("***AtBlQnt:", AtBlQnt)

                        # SortPrel()
                        # DepleteDMNAO()
                        DResid = self.__DMNAO[:self.BSz][:self.BSz].trace()
                        IndS = IndF
                        self.NBOAmnt = IndF
                        print("PP={:4d},Indf={:4d},IndF={:4d},NBOAmnt={:4d}".format(
                            PP, IndS, IndF, self.NBOAmnt))
                        PP = 0
                    else:
                        break
                    if Cnt > 20:
                        break

    # IndF应用型!!
    def __SortPrel(self, PP, PrelOcc, PrelCtr, NCtr, IndF, NBOOcc, NBOVec, PrelVec, NBOCtr):

        ThrOcc = 0.0
        Cnt1 = 0
        PrelAccpt = np.zeros(shape=self.__PNMax, dtype=np.int32)

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
            NBOVec[:self.BSz][IndF] = PrelVec[:self.BSz][i]
            NBOCtr[:NCtr][IndF] = PrelCtr[:NCtr][i]
            NBOCtr[NCtr][IndF] = -1
            IndF += 1

    def DepleteDMNAO(self, NBOOcc, NBOVec, IndS, IndF):
        print("DeleteDMNAO OK")
        for l in range(IndS, IndF, 1):
            for i in range(self.BSz):
                for i in range(self.BSz):
                    self.__DMNAO[i][j] -= NBOOcc[l] * \
                        NBOVec[i][l] * NBOVec[j][l]

    def BlockDMNAO(self, CBl, NCtr, DUMMY):
        flag = np.zeros(shape=self.__PBSz, dtype=np.bool)
        for l in range(NCtr):
            n1 = self.__AtBsRng[CBl[l]][0]
            n2 = self.__AtBsRng[CBl[l]][1]
            for i in range(n1, n2, 1):
                flag = True

        for i in range(self.BSz):
            for j in range(self.BSz):
                if (flag[i] and flag[j]):
                    DUMMY[i][j] = self.__DMNAO[i][j]
                else:
                    DUMMY[i][j] = 0.0

    def EigenSystem(self, DUMMY, EiVal, EiVec):
        # n=BSz
        b = np.zeros(shape=self.__PBSz, dtype=np.float64)
        z = np.zeros(shape=self.__PBSz, dtype=np.float64)

        b[:self.BSz] = np.array([DUMMY[i][i] for i in range(self.BSz)])
        EiVal[:self.BSz] = b[:self.BSz]

        sm = 0.0
        for i in range(1, 101):
            sm = np.sum(np.triu(np.fabs(DUMMY[:self.BSz, :self.BSz]), 1))
            assert sm != 0.0, "sm = 0.0"

            print("EigenSystem: loop= {:4d},sum= {:lf}".format(i, sm))

            if i < 4:
                tresh = 0.2 * sm / (self.BSz ** 2)
            else:
                tresh = 0.0

            for ip in range(self.BSz - 1):
                for iq in range((ip + 1), self.BSz):
                    g = 100.0 * np.fabs(DUMMY[ip][iq])
                    if i > 4 and (np.fabs(EiVal[ip]) + g == np.fabs(EiVal[ip])) and (np.fabs(EiVal[iq]) + g == fabs(EiVal[iq])):
                        DUMMY[ip][iq] = 0.0
                    elif np.fabs(DUMMY[ip][iq] > tresh):
                        h = EiVal[iq] - EiVal[ip]
                        if np.fabs(h) + g == np.fabs(h):
                            t = DUMMY[ip][iq] / h
                        else:
                            theta = 0.5 * h / DUMMY[ip][iq]
                            t = 1.0 / (np.fabs(theta) +
                                       np.sqrt(1.0 + theta ** 2))
                            if theta < 0.0:
                                t = -t

                        c = 1.0 / np.sqrt(1 + t ** 2)
                        s = t * c
                        tau = s / (1.0 + c)
                        h = t * DUMMY[ip][iq]
                        z[ip] = z[ip] - h
                        z[iq] = z[iq] + h
                        EiVal[ip] = EiVal[ip] - h
                        EiVal[iq] = EiVal[iq] + h
                        DUMMY[ip][iq] = 0.0

                        for j in range(ip):
                            g = DUMMY[j][ip]
                            h = DUMMY[j][iq]
                            DUMMY[j][ip] = g - s * (h + g * tau)
                            DUMMY[j][iq] = h + s * (g - h * tau)

                        for j in range(ip+1, iq):
                            g = DUMMY[ip][j]
                            h = DUMMY[j][iq]
                            DUMMY[ip][j] = g - s * (h + g * tau)
                            DUMMY[j][iq] = h + s * (g - h * tau)

                        for j in range(iq+1, self.BSz):
                            g = DUMMY[ip][j]
                            h = DUMMY[iq][j]
                            DUMMY[ip][j] = g - s * (h + g * tau)
                            DUMMY[iq][j] = h + s * (g - h * tau)

                        for j in range(self.BSz):
                            g = EiVec[j][ip]
                            h = EiVec[j][iq]
                            EiVec[j][ip] = g - s * (h + g * tau)
                            EiVec[j][iq] = h + s * (g - h * tau)

            for ip in range(self.BSz):
                b[ip] = b[ip] + z[ip]
                d[ip] = b[ip]
                z[ip] = 0.0

        print("\n****************\nToo many iterations in jacobi\n*************\n")

    def EigenSrt(self):

        p = 0.0

        # GIT

    def Subsets():
        pass
        # get

    def Output():
        pass
        # git

    def BasisChange():
        pass
        # git

    def NBOPlotMolden():
        pass
        # git
