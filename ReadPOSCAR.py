
'''
author: Yilin Zhang
---
本脚本针对给定的·POSCAR·文件来进行信息提取与处理
具备的主要功能：
1.结构信息读取
* POSCAR的名称[str] getSystemName()
* 缩放系数[float] getScalingFactor()
* 晶格矢量矩阵[numpy] getLatticeMatrix()
* 所含元素及其数量[dic] getAtomsInfo()
* 坐标类型[str] getCoordinateType()
* 获取坐标位置信息[numpy] getAtomsPositionMatrix()
* 获取各元素的坐标矩阵[dic] getElementsPositionMatrix()

2.几何/非直接的信息读取
* 计算给定两个矢量的夹角(矢量需为笛卡尔坐标系)[float] calAngleBetween2Vectors(self, vector0, vector1); e.g. vector0 = [x0,y0,z0]
* 计算原子总数[int] calAtomsSum()
* 计算晶胞体积[float] calVolume()
* 计算晶格矢量的模[numpy] calLatticeNorm()

3.POSCAR的再输出
* 输出对应的POSCAR以后缀_D和_C来区分坐标[file] outputPOSCAR()

&& main()函数用于其他脚本的快速调用，会返还其中设置的所有信息，但一般用不到

other:
如需要几何变换操作可以查看SymmetricalOperation.py脚本，包含扩胞，旋转，镜面对称等操作

---
版本更新：2021/05/03
* 重新更新代码，更具备鲁棒性
* 本次在getAtomsPositionMatrix()/getElementsPositionMatrix()中新增功能:
- 选择获取direct的坐标还是笛卡尔坐标(默认给出笛卡尔坐标系)
* 新增POSCAR文件输出outputPOSCAR()
选择是否输出对应的POSCAR以后缀_D和_C来区分坐标
完成更新目标[2021/05/06]
---


'''
import os
import numpy as np
from collections import OrderedDict


class get_poscar_info():

    def __init__(self, path='./POSCAR', coor='None'):
        '''
        *coor:
        -None[default]:保持原来POSCAR
        -C:以笛卡尔坐标提取信息
        -D:以分数坐标提取信息
        '''
        self.path = path
        self.parameters = OrderedDict()
        self.coor = coor.upper()[0]
        # print(self.coor)
        f = open(self.path, 'r')
        self.lines = f.readlines()
        f.close()

    def __addParameters(self, key, value):
        self.parameters[key] = value
        if False:
            if 'Matrix' in key:
                print(str(key) + ':\n' + str(self.parameters[key]))
            else:
                print(str(key) + ':' + str(self.parameters[key]))

    def getSystemName(self):
        SystemName = str(self.lines[0].strip())
        self.__addParameters(key='SystemName', value=SystemName)
        return SystemName

    def getScalingFactor(self):
        ScalingFactor = np.array(
            str(self.lines[1]).strip().split()).astype(np.float)[0]
        self.__addParameters(key='ScalingFactor', value=ScalingFactor)
        return ScalingFactor

    def getLatticeMatrix(self):
        a = np.array(str(self.lines[2]).strip().split()).astype(np.float)
        b = np.array(str(self.lines[3]).strip().split()).astype(np.float)
        c = np.array(str(self.lines[4]).strip().split()).astype(np.float)
        LatticeMatrix = np.array([a, b, c])
        self.__addParameters(key='LatticeMatrix', value=LatticeMatrix)
        return LatticeMatrix

    def getAtomsInfo(self):
        AtomsInfo = OrderedDict()
        AtomsKeys = self.lines[5].strip().split()
        AtomsNumber = self.lines[6].strip().split()
        for i in range(len(AtomsKeys)):
            AtomsInfo[AtomsKeys[i]] = int(AtomsNumber[i])
        self.__addParameters(key='AtomsInfo', value=AtomsInfo)
        return AtomsInfo

    def getCoordinateType(self):
        CoordinateType = str(self.lines[7].strip())
        self.__addParameters(key='CoordinateType', value=CoordinateType)
        return CoordinateType

    def getAtomsPositionMatrix(self,):
        AtomsSum = self.calAtomsSum()
        AtomsPositionMatrix = np.zeros((AtomsSum, 3))
        for i in range(AtomsSum):
            AtomsPositionMatrix[i] = np.array(
                str(self.lines[i + 8]).strip().split())[0:3].astype(np.float)

        # 坐标转换
        if self.getCoordinateType().upper()[0] == 'D':
            if self.coor == 'C':
                AtomsPositionMatrix = np.dot(AtomsPositionMatrix, self.getLatticeMatrix())
        if self.getCoordinateType().upper()[0] == 'C':
            if self.coor == 'D':
                AtomsPositionMatrix = np.dot(AtomsPositionMatrix, np.linalg.inv(self.getLatticeMatrix()))

        self.__addParameters(key='AtomsPositionMatrix', value=AtomsPositionMatrix)

        return AtomsPositionMatrix

    def getElementsPositionMatrix(self,):
        AtomsInfo = self.getAtomsInfo()
        AtomsPositionMatrix = self.getAtomsPositionMatrix()

        ElementsPositionMatrix = OrderedDict()
        count = 0
        for key, value in AtomsInfo.items():
            ElementsPositionMatrix[key] = np.zeros((value, 3))
            for i in range(value):
                ElementsPositionMatrix[key][i] = AtomsPositionMatrix[i + count]
            count += value
        self.__addParameters(key='ElementsPositionMatrix',
                             value=ElementsPositionMatrix)
        return ElementsPositionMatrix

    def outputPOSCAR(self, vesta=True):
        def float2line(line):
            string = '%20.10f' % float(
                line[0]) + '%21.10f' % float(line[1]) + '%21.10f' % float(line[2]) + '\n'
            line = str(string)
            return line

        path = self.path[:-len(self.path.strip().split('/')[-1])]

        filedir = path + 'POSCAR_' + self.coor

        self.main()
        with open(filedir, 'w') as f:
            f.write(self.parameters['SystemName'] + '\n')
            f.write(str(self.parameters['ScalingFactor']) + '\n')
            LatticeMatrix = self.parameters['LatticeMatrix']
            for i in range(len(LatticeMatrix)):
                f.write(float2line(LatticeMatrix[i]))

            AtomsInfo = self.parameters['AtomsInfo']
            string = ['', '']
            for key in AtomsInfo.keys():
                string[0] += key + '    '
                string[1] += str(AtomsInfo[key]) + '    '
            string[0] += '\n'
            string[1] += '\n'
            f.write(string[0])
            f.write(string[1])
            if self.coor == 'C':
                f.write('Cartesian\n')
            elif self.coor == 'D':
                f.write('Direct\n')
            else:
                f.write(self.getCoordinateType() + '\n')
            ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy()
            for key in AtomsInfo.keys():
                arr = ElementsPositionMatrix[key]
                for i in range(len(arr)):
                    f.write(float2line(arr[i]))

        if vesta:
            os.chdir(path)
            os.system('VESTA ' + filedir)
        return filedir

    def calAngleBetween2Vectors(self, vector0, vector1):
        '''
        获取两个矢量的夹角
        '''
        angle = np.arccos(np.dot(vector0, vector1) /
                          (np.linalg.norm(vector0) * np.linalg.norm(vector1)))
        return angle

    def calLatticeMatrix_Transformation_2D(self):
        '''
        test,test,test!
        这部分还未添加到功能说明中
        '''
        LatticeMatrix = self.getLatticeMatrix()
        a = LatticeMatrix[0]
        b = LatticeMatrix[1]
        angle = self.calAngleBetween2Vectors(a, b)
        a_norm = np.linalg.norm(a)
        b_norm = np.linalg.norm(b)
        a = np.array([np.cos(angle / 2) * a_norm, np.sin(angle / 2) * a_norm, 0])
        b = np.array([np.cos(angle / 2) * b_norm, np.sin(angle / 2) * b_norm * (-1), 0])
        return a, b

    def calAtomsSum(self):
        AtomsInfo = self.getAtomsInfo()
        AtomsSum = 0
        for value in AtomsInfo.values():
            AtomsSum += value
        self.__addParameters(key='AtomsSum', value=AtomsSum)
        return AtomsSum

    def calVolume(self):
        """
        Get unit cell volume
        """
        sf = self.getScalingFactor()
        a = np.array(str(self.lines[2]).strip().split()).astype(np.float) * sf
        b = np.array(str(self.lines[3]).strip().split()).astype(np.float) * sf
        c = np.array(str(self.lines[4]).strip().split()).astype(np.float) * sf
        Volume = np.dot(np.cross(a, b), c)
        self.__addParameters(key='Volume', value=Volume)
        return Volume

    def direct2Cartesian(self):
        return np.dot(self.getAtomsPositionMatrix()[:, :], self.getLatticeMatrix())

    def calLatticeNorm(self,):
        '''
        获取晶格参数的模
        '''
        LatticeMatrix = self.getLatticeMatrix()
        LatticeNorm = []
        for i in range(3):
            LatticeNorm.append(np.linalg.norm(LatticeMatrix[i]))
        LatticeNorm = np.array(LatticeNorm).astype(float)

        return LatticeNorm

    def main(self):
        self.getSystemName()  # 获取体系名称
        self.getScalingFactor()  # 获取缩放系数
        self.getLatticeMatrix()  # 获取晶格参数矩阵
        self.getAtomsInfo()  # 获取原子(数量)信息
        self.getCoordinateType()  # 获取坐标类型
        self.getAtomsPositionMatrix()  # 获取原子位置矩阵
        self.getElementsPositionMatrix()  # 获取每种元素的位置矩阵

        self.calAtomsSum()  # 计算原子总数
        self.calVolume()  # 计算晶体面积

        # os.system('clear')  # 注释掉此行可以方便脚本检查
        # for key, value in self.parameters.items():
        #     if 'Matrix' in key:
        #         print(str(key) + ':\n' + str(value))
        #     else:
        #         print(str(key) + ':' + str(value))

        return self.parameters


if __name__ == '__main__':
    '''
    功能测试:
    1.从github下载POSCAR文件
    2.进行测试
    '''

    '''
    常用功能的测试
    '''
    gpi = get_poscar_info(path='./POSCAR')
    # print(gpi.getSystemName())
    # print(gpi.getScalingFactor())
    # print(gpi.getLatticeMatrix())
    # print(gpi.getAtomsInfo())
    # print(gpi.getCoordinateType())
    # print(gpi.getAtomsPositionMatrix())
    # print(gpi.getElementsPositionMatrix())
    gpi.outputPOSCAR()


'''
 1.结构信息读取
 * POSCAR的名称[str] getSystemName()
 * 缩放系数[float] getScalingFactor()
 * 晶格矢量矩阵[numpy] getLatticeMatrix()
 * 所含元素及其数量[dic] getAtomsInfo()
 * 坐标类型[str] getCoordinateType()
 * 获取坐标位置信息[numpy] getAtomsPositionMatrix()
 * 获取各元素的坐标矩阵[dic] getElementsPositionMatrix()

 2.几何/非直接的信息读取
 * 计算给定两个矢量的夹角(矢量需为笛卡尔坐标系)[float] calAngleBetween2Vectors(self, vector0, vector1); e.g. vector0 = [x0,y0,z0]
 * 计算原子总数[int] calAtomsSum()
 * 计算晶胞体积[float] calVolume()
 * 计算晶格矢量的模[numpy] calLatticeNorm()

 3.POSCAR的再输出
 * 输出对应的POSCAR以后缀_D和_C来区分坐标[file] outputPOSCAR()
'''
