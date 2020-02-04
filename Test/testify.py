import numpy as np


BSz = 21
DMNAO = np.zeros(shape=(BSz, BSz), dtype=np.float64)
NAOAO = np.zeros(shape=(BSz, BSz), dtype=np.float64)

dtln=""
with open("./nbo.log", 'r') as fp:
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