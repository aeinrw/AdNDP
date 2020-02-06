import numpy as np
from scipy.special import comb

NCtr = 3
NAt = 6

Nu = comb(NAt,NCtr)


AtBl = np.zeros((int(Nu),NAt),dtype=np.int32)


def Subsets(AtBl,NClmn, NAt):

    NLn = 0
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

print(Subsets(AtBl,NCtr,NAt))

for i in range(int(Nu):
    print(AtBl[i,:])

