import numpy as np

# # 宏变量
# __PNAt = 128
# __PNVl = 9
# __PNTot = 15
# __PNMax = 4096
# __PBSz = 1024

# # 全局变量
# CMOFile = ""
# ADNDPFile = ""
# NAt = 0
# NVal = 0
# NTot = 0
# BSz = 0
# AtBs = np.zeros(shape=__PNAt, dtype=np.int32)
# AtBsRng = np.zeros(shape=(__PNAt, 2), dtype=np.int32)
# DMNAO = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)
# Thr = np.zeros(shape=__PNAt, dtype=np.float64)


# # main里面的变量
# MnSrch = 0
# NBOOcc = np.zeros(shape=__PNMax, dtype=np.float64)
# NBOVes = np.zeros(shape=(__PBSz, __PNMax), dtype=np.float64)
# DResid = 0.0
# NBOCtr = np.zeros(shape=(__PNAt, __PNMax), dtype=np.int32)
# NBOAmnt = 0
# NAOAO = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)

# # AdNBO里面的变量
# PrelOcc = np.zeros(shape=__PNMax, dtype=np.float64)
# PrelVec = np.zeros(shape=(__PBSz, __PNMax), dtype=np.float64)
# PrelCtr = np.zeros(shape=(__PNAt, __PNMax), dtype=np.int32)

# Cnt = 0
# PP = 0
# IndS = 0
# IndF = 0
# NCtr = 0
# AtBlQnt = 0
# AtBl = np.zeros(shape=(100000, __PNAt), dtype=np.int32)
# CBl = np.zeros(shape=__PNAt, dtype=np.int32)
# threshold = 0.0
# vmax = 0.0
# DUMMY = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)
# EiVal = np.zeros(shape=__PBSz, dtype=np.float64)
# EiVec = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)


class AdNDP:

    # 宏变量
    __PNAt = 128
    __PNVl = 9
    __PNTot = 15
    __PNMax = 4096
    __PBSz = 1024

    # 全局变量
    __CMOFile = ""
    __ADNDPFile = ""
    __NAt = 0
    __NVal = 0
    __NTot = 0
    __BSz = 0
    __AtBs = np.zeros(shape=__PNAt, dtype=np.int32)
    __AtBsRng = np.zeros(shape=(__PNAt, 2), dtype=np.int32)
    __DMNAO = np.zeros(shape=(__PBSz, __PBSz), dtype=np.float64)
    __Thr = np.zeros(shape=__PNAt, dtype=np.float64)

    def __init__(self):

        # 输入AdNDP.ini文件
        import configparser
        cf = configparser.ConfigParser()
        cf.read("../AdNDP.ini")
        sections = cf.sections()[0]
        nbofile = cf.get(sections, "nbofile")
        self.__NAt = cf.getint(sections, "NAt")
        self.__NVal = cf.getint(sections, "NVal")
        self.__NTot = cf.getint(sections, "NTot")
        self.__BSz = cf.getint(sections, "BSz")
        self.__AtBs[0:self.__NAt] = np.array(
            cf.get(sections, "AtBs").split(','), dtype=np.int32)
        self.__Thr[0:self.__NAt] = np.array(
            cf.get(sections, "Thr").split(','), dtype=np.float64)
        self.__CMOFile = cf.get(sections, "CMOFile")
        self.__ADNDPFile = cf.get(sections, "AdNDPFile")


