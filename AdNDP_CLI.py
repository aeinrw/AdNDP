import numpy as np
import configparser
from scipy.special import comb
import time


class AdNDP:

    def __init__(self, nboFile):
        self.__nboFile = nboFile
        self.__readFile()

    def __readFile(self):
        print("读取文件中...")
        fp = open(self.__nboFile, 'r')

        dtln = fp.readline()
        while dtln[:9] != ' Charge =':
            dtln = fp.readline()

        NAt = 0
        dtln = fp.readline()
        while (dtln[:2]).lstrip() != '':
            NAt += 1
            dtln = fp.readline()

        while dtln[:20] != " NATURAL POPULATIONS":
            dtln = fp.readline()

        fp.readline()
        fp.readline()
        fp.readline()

        AtBs = []
        for i in range(NAt):
            dtln = fp.readline()
            atom = dtln[9:11]
            num = 0
            while dtln != '\n':
                num += 1
                dtln = fp.readline()
            AtBs.append(num)

        fp.close()

        self.NAt = NAt
        self.AtBs = np.array(AtBs)
        self.NBf = np.sum(AtBs)
        # self.Thr = np.array([0.1]*self.NAt)
        self.Thr = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        self.__AtBsRng = np.zeros((self.NAt, 2), dtype=np.int32)
        j = 0
        for i in range(self.NAt):
            self.__AtBsRng[i, 0] = j
            j += self.AtBs[i]
            self.__AtBsRng[i, 1] = j
        self.DMNAO = self.__readMatrix(self.__nboFile, "DMNAO")

        print("原子总数:{:d}".format(self.NAt))
        print("原子轨道总数:{:d}".format(self.NBf))
        print("每个原子的轨道总数:"+str(self.AtBs))
        print("占据数阈值" + str(self.Thr))

    def __readMatrix(self, nboFile, matrix):
        if matrix == 'DMNAO':
            string = " NAO density matrix:"
        elif matrix == "AONAO":
            string = " NAOs in the AO basi"

        fp = open(nboFile, 'r')
        Matrix = np.zeros((self.NBf, self.NBf), dtype=np.float64)

        dtln = ""
        while dtln[:20] != string:
            dtln = fp.readline()
        dtln = fp.readline()
        dtln = fp.readline()
        dtln = fp.readline()
        resid = self.NBf
        k = 0
        Cnt = 0
        while resid > 0:
            if resid < 8:
                Cnt = resid
            else:
                Cnt = 8
            for i in range(self.NBf):
                line = fp.readline()[17:].replace("  ", " ").split(" ")
                if line[0] == '':
                    line.pop(0)
                Matrix[i, 8*k:8*k+Cnt] = np.array(line, dtype=np.float64)
            dtln = fp.readline()
            dtln = fp.readline()
            dtln = fp.readline()
            resid -= Cnt
            k += 1

        print("读取" + matrix + "矩阵成功")
        fp.close()
        return Matrix

    def partition(self):

        print("分析中...")
        prelOcc = np.zeros(self.NBf, dtype=np.float64)
        prelVec = np.zeros((self.NBf, self.NBf), dtype=np.float64)
        prelCtr = np.ones((self.NAt, self.NBf), dtype=np.int32)
        prelCtr *= -1

        self.nboOcc = np.zeros(self.NBf, dtype=np.float64)
        self.nboVec = np.zeros((self.NBf, self.NBf), dtype=np.float64)
        self.nboCtr = np.ones((self.NAt, self.NBf), dtype=np.int32)
        self.nboCtr *= -1

        indS, indF = 0, 0
        nboAmnt = 0

        for NCtr in range(1, self.NAt+1):
            threshold = self.Thr[NCtr-1]
            if threshold == 0.0:
                print("跳过{:d}c-2e".format(NCtr))
                continue

            AtBl, AtBlQnt = self.__subsets(NCtr)
            print("搜索{:d}c-2e中......".format(NCtr))

            counter = 0

            PP = 0
            Cnt = 1
            while True:
                for k in range(AtBlQnt):
                    CBl = AtBl[k, :]
                    DUMMY = self.__blockDMNAO(CBl)
                    eigenVal, eigenVec = self.__eigenSystem(DUMMY)

                    if np.fabs(eigenVal - 2.0) <= threshold:
                        prelOcc[PP] = eigenVal
                        prelVec[:, PP] = eigenVec
                        prelCtr[:NCtr, PP] = CBl[:NCtr]
                        PP += 1
                Cnt += 1
                if PP > 0:
                    AtBl[:PP, :] = prelCtr[:NCtr, :PP].T
                    AtBlQnt = PP
                    indF = self.__sortPrel(prelOcc, prelVec, prelCtr, PP, indF)
                    for i in range(indS, indF):
                        counter += 1
                        a = self.nboVec[:, i]
                        self.DMNAO -= self.nboOcc[i]*np.outer(a, a)
                        print("    找到第{:d}个{:d}c-2e轨道".format(counter, NCtr))
                    indS = indF
                    nboAmnt = indF
                    PP = 0
                else:
                    break
                if Cnt > 20:
                    break

        self.nboAmnt = nboAmnt
        print("分析结束")
        return nboAmnt

    def output(self):
        print("分析结果如下")
        print("轨道总数为{:d}".format(self.nboAmnt))
        for i in range(self.nboAmnt):
            nc = np.sum(self.nboCtr[:, i] != -1)
            string = "第{:d}个轨道({:d}中心,{:.3f}电子):".format(
                i+1, nc, np.sum(self.nboOcc[i]))
            #print("bond {:3d} ({:2d} center,{:.3f} e):".format(i+1,nc,np.sum(self.nboOcc[i])))
            for j in range(nc):
                m = self.nboCtr[j, i]
                n1 = self.__AtBsRng[m, 0]
                n2 = self.__AtBsRng[m, 1]
                vocc = self.nboOcc[i]*(np.sum((self.nboVec[:, i]**2)[n1:n2]))
                string += "{:3d}({:.3f}) ".format(
                    self.nboCtr[j, i]+1, vocc)
                #print("{:3d}({:.3f}) ".format(self.nboCtr[j,i]+1,vocc))
            print(string)

        print("电子剩余:{:5.6f}".format(self.DMNAO.trace()))

    def nboPlotMolden(self):
        AONAO = self.__readMatrix(self.__nboFile, "AONAO")
        nboVecAO = AONAO@self.nboVec

        fp = open("adndp.log", 'w')
        resid = self.nboAmnt
        k = 0
        while resid > 0:
            fp.write("{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n"
                     .format(self.nboOcc[5*k], self.nboOcc[5*k+1], self.nboOcc[5*k+2], self.nboOcc[5*k+3], self.nboOcc[5*k+4]))
            fp.write("-"*50+'\n')
            for i in range(self.NBf):
                fp.write("{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n"
                         .format(nboVecAO[i][5*k], nboVecAO[i][5*k+1], nboVecAO[i][5*k+2], nboVecAO[i][5*k+3], nboVecAO[i][5*k+4]))
            k += 1
            resid -= 5
            fp.write('\n')

        fp.close()

    def __subsets(self, NCtr):
        AtBlQnt = int(comb(self.NAt, NCtr))
        AtBl = np.zeros((AtBlQnt, NCtr), dtype=np.int32)

        if NCtr == self.NAt:
            AtBl[0, :] = np.arange(self.NAt)
            return AtBl, AtBlQnt

        if NCtr == 1:
            AtBl[:, 0] = np.arange(self.NAt)
            return AtBl, AtBlQnt

        Done = 0
        NLn = 1

        AtBl[NLn-1, :] = np.arange(NCtr)

        while (True):
            for k in range(NCtr):
                if AtBl[NLn - 1, k] == (self.NAt - NCtr + k):
                    Done = k + 1
                    break
            if Done == 1:
                return AtBl, AtBlQnt
            if Done == 0:
                NLn += 1
                AtBl[NLn - 1, :NCtr-1] = AtBl[NLn - 2, :NCtr-1]
                if AtBl[NLn - 2, NCtr - 1] < self.NAt:
                    AtBl[NLn - 1, NCtr - 1] = AtBl[NLn - 2, NCtr - 1] + 1
            else:
                NLn += 1
                AtBl[NLn - 1, :Done-2] = AtBl[NLn - 2, :Done-2]
                AtBl[NLn - 1, Done - 2] = AtBl[NLn - 2, Done - 2] + 1
                for k in range(Done - 1, NCtr):
                    AtBl[NLn - 1, k] = AtBl[NLn - 1, k - 1] + 1
                Done = 0

        return AtBl, AtBlQnt

    def __blockDMNAO(self, CBl):

        DUMMY = np.zeros(self.DMNAO.shape, dtype=np.float64)
        flag = []

        for l in range(len(CBl)):
            n1 = self.__AtBsRng[CBl[l], 0]
            n2 = self.__AtBsRng[CBl[l], 1]
            flag.append((n1, n2))

        for i in flag:
            for j in flag:
                DUMMY[i[0]:i[1], j[0]:j[1]] = self.DMNAO[i[0]:i[1], j[0]:j[1]]

        return DUMMY

    def __eigenSystem(self, DUMMY):

        values, vectors = np.linalg.eig(DUMMY)
        values, vectors = np.real(values), np.real(vectors)
        imax = np.argmax(values)

        return values[imax], vectors[:, imax]

    def __sortPrel(self, prelOcc, prelVec, prelCtr, PP, indF):

        index = np.argsort(prelOcc[:PP])[::-1]
        temp = np.floor((prelOcc[index])*100000)
        Cnt = 1
        if len(temp) != 1:
            while temp[Cnt] == temp[Cnt-1]:
                Cnt += 1
                if Cnt == PP:
                    break
        for j in range(Cnt):
            i = index[j]
            self.nboOcc[indF] = prelOcc[i]
            self.nboVec[:, indF] = prelVec[:, i]
            self.nboCtr[:, indF] = prelCtr[:, i]
            indF += 1

        return indF


if __name__ == "__main__":
    Li5 = AdNDP("./example/Li5+.log")
    Li5.partition()
    Li5.output()
    Li5.nboPlotMolden()
