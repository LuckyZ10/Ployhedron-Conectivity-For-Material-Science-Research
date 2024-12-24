'''
author: Yilin Zhang
---
本脚本针对给定的·POSCAR·文件来进行结构变换操作
具备的主要功能：
1.坐标变换
* coordinateTransfer()函数
- coor = 'D2C' Direct到Cartesian转换
- coor = 'C2D' Cartesian到Direct转换


---
版本更新：2021/05/06
* 重新更新代码，更具备鲁棒性
* 加入坐标转换的函数 coordinateTransfer():
- 把·Direct到Cartesian·和·Cartesian到Direct·的转换合并到coordinateTransfer()函数中
* 重写扩胞程序 simpleDIM()
- 正常的阔胞功能
- 加入阔胞取中心胞的功能
完成本次更新目标 2021/05/06
'''


from STRUCTURE.ReadPOSCAR import get_poscar_info
from collections import OrderedDict
import numpy as np
import os


class SOP():
    def __init__(self, path):
        self.readpath = path
        gpi = get_poscar_info(path=path)
        self.path = path[0:-len(path.strip().split('/')[-1])]
        self.parameters = gpi.main()
        self.process = []
        os.system('clear')

    def coordinateTransfer(self, coor='None'):
        '''
        coor:
        'None'[default]
        'D2C'
        'C2D'
        '''

        coor = coor.upper()[0: 3]
        if coor == 'Non':
            return 0

        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy()
        LatticeMatrix = self.parameters['LatticeMatrix'].copy()

        if self.parameters['CoordinateType'].upper()[0] == 'C':
            if coor == 'C2D':
                for key in ElementsPositionMatrix.keys():
                    ElementsPositionMatrix[key] = np.dot(ElementsPositionMatrix[key], np.linalg.inv(LatticeMatrix))
                self.process.append('Cartesian to Direct:YES')
                self.parameters['CoordinateType'] = 'Direct'

        if self.parameters['CoordinateType'].upper()[0] == 'D':
            if coor == 'D2C':
                for key in ElementsPositionMatrix.keys():
                    ElementsPositionMatrix[key] = np.dot(ElementsPositionMatrix[key], LatticeMatrix)
                self.process.append(' Direct to Cartesian:YES')
                self.parameters['CoordinateType'] = 'Cartesian'

        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

    def simpleDIM(self, dim=[3, 3, 3], get_center=False):
        '''
        简单的扩胞,[x,y,z]对应x,y,z方向上的扩胞倍数
        get_center = True 将返回中心胞的笛卡尔坐标系,并输出以_center为后缀的POSCAR文件[POSCAR_center]
        '''
        # 0.根据应用场景确定扩胞参数
        dim = np.array(dim)
        if get_center:
            dim = np.array([5, 5, 5])

        # 1.首先转为分数坐标系方便扩胞
        self.coordinateTransfer(coor='C2D')

        # 2.更新晶格矢量
        LatticeMatrix = self.parameters['LatticeMatrix']
        self.parameters['LatticeMatrix'] = [
            LatticeMatrix[0] * dim[0], LatticeMatrix[1] * dim[1], LatticeMatrix[2] * dim[2]]

        # 3.更新元素坐标矩阵
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy()
        ElementsPositionMatrix_copy = self.parameters['ElementsPositionMatrix'].copy()  # 用于后续的附加功能

        if self.parameters['CoordinateType'].upper()[0] == 'D':
            a = np.array([1, 0, 0])
            b = np.array([0, 1, 0])
            c = np.array([0, 0, 1])

        for key in ElementsPositionMatrix.keys():
            ElementsPositionMatrix_backup = ElementsPositionMatrix[key]
            for i in range(1, dim[0]):
                ElementsPositionMatrix[key] = np.vstack(
                    (ElementsPositionMatrix[key], ElementsPositionMatrix_backup + i * a))
            ElementsPositionMatrix_backup = ElementsPositionMatrix[key]
            for i in range(1, dim[1]):
                ElementsPositionMatrix[key] = np.vstack(
                    (ElementsPositionMatrix[key], ElementsPositionMatrix_backup + i * b))
            ElementsPositionMatrix_backup = ElementsPositionMatrix[key]
            for i in range(1, dim[2]):
                ElementsPositionMatrix[key] = np.vstack(
                    (ElementsPositionMatrix[key], ElementsPositionMatrix_backup + i * c))

            ElementsPositionMatrix[key] /= np.array(dim)
            self.parameters['AtomsInfo'][key] = len(ElementsPositionMatrix[key])

        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('simpleDIM:\n' + str(np.array(dim)))

        self.cutLattice([[0, 1], [0, 1], [0, 1]], flag=False)
        self.outputPOSCAR()
        # ----------- 额外的附加功能 ----------- #
        if get_center:
            ElementsPositionMatrix_copy = self.parameters['ElementsPositionMatrix'].copy()

            center_atoms = OrderedDict()
            self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix_copy
            self.cutLattice([[0.4, 0.599999], [0.4, 0.599999], [0.4, 0.599999]], flag=False)
            center_atoms = self.parameters['ElementsPositionMatrix']

            envir_atoms = OrderedDict()
            self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix_copy
            self.cutLattice([[0.2, 0.799999], [0.2, 0.799999], [0.2, 0.799999]], flag=False)
            envir_atoms = self.parameters['ElementsPositionMatrix']

            system_atoms = OrderedDict()
            self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix_copy
            self.cutLattice([[0, 0.999999], [0, 0.999999], [0, 0.999999]], flag=False)
            system_atoms = self.parameters['ElementsPositionMatrix']
            # print(envir_atoms)
            # for key in ElementsPositionMatrix_copy.keys():
            #     center_atom[key] = ElementsPositionMatrix_copy[key] + 2 * a + 2 * b + 2 * c
            #     center_atom[key] /= np.array(dim)
            #     center_atom[key] = np.dot(center_atom[key], self.parameters['LatticeMatrix'])

            for key in ElementsPositionMatrix_copy.keys():
                center_atoms[key] = np.dot(center_atoms[key], self.parameters['LatticeMatrix'])
                envir_atoms[key] = np.dot(envir_atoms[key], self.parameters['LatticeMatrix'])
                system_atoms[key] = np.dot(system_atoms[key], self.parameters['LatticeMatrix'])

            return center_atoms, envir_atoms, system_atoms

    def outputPOSCAR(self, vesta=False):
        def float2line(line):
            string = '%20.10f' % float(
                line[0]) + '%21.10f' % float(line[1]) + '%21.10f' % float(line[2]) + '\n'
            line = str(string)
            return line

        filedir = self.path

        f = open(filedir + 'SOP_step', 'w')

        f.write('\n!BEGIN\n')

        f.write('---------------------------------\n')
        for i in range(len(self.process)):
            f.write(str(i + 1) + '.' + self.process[i] + '\n')
            f.write('********************************\n')
        f.write('FINISHED!\n')
        f.write('---------------------------------\n\n')
        f.write('=================================\n')
        f.write('\nthe input file:' + self.readpath + '\n')
        pos = open(self.readpath)
        lines = pos.readlines()
        pos.close()

        for line in lines:
            f.write(line)
        f.close()

        with open(filedir + 'POSCAR_sop', 'w') as f:
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

            f.write(self.parameters['CoordinateType'] + '\n')
            ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy()
            for key in AtomsInfo.keys():
                arr = ElementsPositionMatrix[key]
                for i in range(len(arr)):
                    f.write(float2line(arr[i]))

        if vesta:
            os.chdir(self.path)
            os.system('VESTA ' + filedir + 'POSCAR_sop')

    def transformMatrixR(self, fixedAxis, angle):
        '''
        fixedAxis=[[0,0,0],[0,0,1]]:对应的旋转轴的起点和旋转轴的终点
        angle:绕轴旋转的角度
        注意：尽量用abc相等的包进行旋转，因为用的是分数坐标系(本问题要解决一下)
        '''
        angle = math.radians(angle)
        fixedAxis = np.array(fixedAxis)

        a = fixedAxis[0][0]
        b = fixedAxis[0][1]
        c = fixedAxis[0][2]

        vector = fixedAxis[1] - fixedAxis[0]

        u = vector[0] / np.linalg.norm(vector)
        v = vector[1] / np.linalg.norm(vector)
        w = vector[2] / np.linalg.norm(vector)

        cosA = np.cos(angle)
        sinA = np.sin(angle)

        K1 = [u * u + (v * v + w * w) * cosA, u * v * (1 - cosA) - w * sinA, u * w * (1 - cosA) +
              v * sinA, (a * (v * v + w * w) - u * (b * v + c * w)) * (1 - cosA) + (b * w - c * v) * sinA]
        K2 = [u * v * (1 - cosA) + w * sinA, v * v + (u * u + w * w) * cosA, v * w * (1 - cosA) -
              u * sinA, (b * (u * u + w * w) - v * (a * u + c * w)) * (1 - cosA) + (c * u - a * w) * sinA]
        K3 = [u * w * (1 - cosA) - v * sinA, v * w * (1 - cosA) + u * sinA, w * w + (u * u + v * v)
              * cosA, (c * (u * u + v * v) - w * (a * u + b * v)) * (1 - cosA) + (a * v - b * u) * sinA]
        K4 = [0, 0, 0, 1]
        return np.array([K1, K2, K3, K4])

    def move2point(self, refElement, point):
        '''
        将某个元素的第几个原子refElement移动到某个点point
        refElement=['Sn',1]
        '''
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        refElementPositionMatrix = ElementsPositionMatrix[refElement[0]
                                                          ][refElement[1] - 1]
        delta_r = np.array(point) - refElementPositionMatrix
        for key in ElementsPositionMatrix.keys():
            ElementsPositionMatrix[key] += delta_r

        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('move2point:\n' + '-refElement:'
                            + str(refElement) + '\n-point:' + str(point))

    def modifyingAtomicCoordinates(self, refElement, cor):
        '''
        将某个元素的第几个原子refElement坐标改为cor
        refElement=['Sn',1]
        '''
        LatticeMatrix = self.parameters['LatticeMatrix']
        if self.parameters['CoordinateType'] == 'Direct':
            a = b = c = 1.0
        else:
            a = np.linalg.norm(LatticeMatrix[0])
            b = np.linalg.norm(LatticeMatrix[1])
            c = np.linalg.norm(LatticeMatrix[2])
        basicvector = np.array([a, b, c])

        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        ElementsPositionMatrix[refElement[0]
                               ][refElement[1] - 1] = np.array(cor) * basicvector

        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('modifyingAtomicCoordinates:\n'
                            + '-refElement:' + str(refElement) + '\n-cor:' + str(cor))

    def delOneAtom(self, refElement):
        '''
        将某个元素的第几个原子refElement坐标改为cor
        refElement=['Sn',1]
        '''
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        ElementsPositionMatrix[refElement[0]] = np.delete(
            ElementsPositionMatrix[refElement[0]], refElement[1] - 1, axis=0)

        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        AtomsInfo = self.parameters['AtomsInfo'].copy()
        AtomsInfo[refElement[0]] -= 1
        self.parameters['AtomsInfo'] = AtomsInfo

        self.process.append('delAtom:\n'
                            + '-refElement:' + str(refElement))

    def rotationRoundFixedAxis(self, fixedAxis, angle):

        matrixR = self.transformMatrixR(fixedAxis, angle)

        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        for key in ElementsPositionMatrix.keys():
            matrix = ElementsPositionMatrix[key]
            matrix = np.insert(matrix, 3, [1], axis=1)
            for i in range(len(matrix)):
                matrix[i] = np.dot(matrix[i], matrixR.T)
            ElementsPositionMatrix[key] = matrix[:, 0:-1]
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('rotationRoundFixedAxis:\n' + '-fixedAxis:\n'
                            + str(np.array(fixedAxis)) + '\n-angle:' + str(angle))

        self.cutLattice([[0, 1], [0, 1], [0, 1]], flag=False)

        # LatticeMatrix=self.parameters['LatticeMatrix']
        # basicvector=np.array([LatticeMatrix[0][0],LatticeMatrix[1][1],LatticeMatrix[2][2],1])
        # vector=np.array(fixedAxis[1])-np.array(fixedAxis[0])
        # matrixR = transformMatrixR([[0,0,0],vector], angle)
        # basicvector = np.dot(basicvector, matrixR.T)
        # basicvector =basicvector[0:-1]
        # a=[basicvector[0],0,0]
        # b=[0,basicvector[1],0]
        # c=[0,0,basicvector[2]]
        # self.parameters['LatticeMatrix'] = np.array([a,b,c])

    def rotationRoundFixedPlane(self, fixedPlane, parallelAxis, angle, postionRange):
        '''
        fixedPlane:[[1,1,0],[0,1,0],[0,1,1]],平面由三点构成
        parallelAxis:平行轴[0，0，1]
        angle:绕平行轴转动的角度


        尚未写入的功能：选取某个范围绕着给平面的平行轴作个转动
        postionRange：[[0.5,1],[0.5,1],[0.5,1]]转动范围
        '''
        def getNormalVector(fixedPlane):
            fixedPlane = np.array(fixedPlane)
            vector1 = fixedPlane[1] - fixedPlane[0]
            vector2 = fixedPlane[2] - fixedPlane[0]
            normalVector = np.cross(vector1, vector2)
            normalVector = normalVector / np.linalg.norm(normalVector)
            return normalVector, fixedPlane[0]

        LatticeMatrix = self.parameters['LatticeMatrix']
        basicvector = np.array(
            [LatticeMatrix[0][0], LatticeMatrix[1][1], LatticeMatrix[2][2]])

        normalVector, inPoint = getNormalVector(
            np.array(fixedPlane) * basicvector)

        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        fixedAxis = np.zeros((2, 3))
        postionRange = np.array(postionRange)
        for key in ElementsPositionMatrix.keys():
            matrix = ElementsPositionMatrix[key]
            for i in range(len(matrix)):
                point = matrix[i]
                if (point[0] >= postionRange[0][0] and point[0] <= postionRange[0][1]):
                    if (point[1] >= postionRange[1][0] and point[1] <= postionRange[1][1]):
                        if (point[2] >= postionRange[2][0] and point[2] <= postionRange[2][1]):
                            point *= basicvector
                            dropfoot = point + normalVector * \
                                (normalVector * (inPoint - point))
                            fixedAxis = [dropfoot, dropfoot
                                         + np.array(parallelAxis)]
                            matrixR = self.transformMatrixR(fixedAxis, angle)
                            matrix[i] = np.dot(
                                np.insert(point, 3, [1], axis=0), matrixR.T)[:-1]
                            matrix[i] /= basicvector
            ElementsPositionMatrix[key] = matrix
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('rotationRoundFixedPlane:\n' + '-fixedPlane:\n' + str(np.array(
            fixedPlane)) + '\n-parallelAxis:' + str(np.array(parallelAxis)) + '\n-angle:' + str(angle) + '\n-postionRange:\n' + str(np.array(postionRange)))

        self.cutLattice([[0, 1], [0, 1], [0, 1]], flag=False)

    def cutLattice(self, cutRange, flag=True):
        '''
        cutRange:需要切割的范围
        '''
        self.coordinateTransfer(coor='C2D')
        LatticeMatrix = self.parameters['LatticeMatrix']
        if self.parameters['CoordinateType'].upper()[0] == 'D':
            a = b = c = 1.0
        else:
            a = np.linalg.norm(LatticeMatrix[0])
            b = np.linalg.norm(LatticeMatrix[1])
            c = np.linalg.norm(LatticeMatrix[2])
        basicvector = np.array([a, b, c])

        cutRange = np.array(cutRange)
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        for key in ElementsPositionMatrix.keys():
            matrix = ElementsPositionMatrix[key]
            for i in range(3):
                matrix = matrix[matrix[:, i]
                                >= cutRange[i][0] * basicvector[i], :]
                matrix = matrix[matrix[:, i]
                                < cutRange[i][1] * basicvector[i], :]
            ElementsPositionMatrix[key] = matrix
            self.parameters['AtomsInfo'][key] = len(
                ElementsPositionMatrix[key])
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix
        if self.parameters['CoordinateType'] == 'Cartesian':
            delta = np.array(cutRange[:, 1] - cutRange[:, 0])
            self.parameters['LatticeMatrix'] = [LatticeMatrix[0] * delta[0],
                                                LatticeMatrix[1] * delta[1], LatticeMatrix[2] * delta[2]]
        if flag:
            self.process.append('cutLattice:\n' + str(np.array(cutRange)))

    def localAtomicCoordinatesComponentReset(self, zoom, reset):
        '''
        zoom:需要修改的范围
        reset:需要修正的坐标分量['a', 0.]
        '''
        if reset[0] == 'a':
            axis = 0
        if reset[0] == 'b':
            axis = 1
        if reset[0] == 'c':
            axis = 2

        component = float(reset[1])

        LatticeMatrix = self.parameters['LatticeMatrix']

        zoom = np.array(zoom)
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )

        for key in ElementsPositionMatrix.keys():
            matrix = ElementsPositionMatrix[key]
            for i in range(len(matrix)):
                point = matrix[i]
                if (point[0] >= zoom[0][0] and point[0] <= zoom[0][1]):
                    if (point[1] >= zoom[1][0] and point[1] <= zoom[1][1]):
                        if (point[2] >= zoom[2][0] and point[2] <= zoom[2][1]):
                            matrix[i][axis] = component
            ElementsPositionMatrix[key] = matrix
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('localAtomicCoordinatesComponentReset:\n'
                            + str(np.array(zoom)) + '\nreset:' + str(np.array(reset)))

    def screen(self, screenPlane, justScreen=True):
        '''
        镜面对称：
        screenPlane:[]
        '''

        def getNormalVector(fixedPlane):
            fixedPlane = np.array(fixedPlane)
            vector1 = fixedPlane[1] - fixedPlane[0]
            vector2 = fixedPlane[2] - fixedPlane[0]
            normalVector = np.cross(vector1, vector2)
            normalVector = normalVector / np.linalg.norm(normalVector)
            return -normalVector, fixedPlane[0]

        normalVector, inPoint = getNormalVector(screenPlane)
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        for key in ElementsPositionMatrix.keys():
            matrix = ElementsPositionMatrix[key]
            screenpoint = np.zeros((len(matrix), 3))
            for i in range(len(matrix)):
                point = matrix[i]
                dropfoot = point + normalVector * \
                    (normalVector * (inPoint - point))
                screenpoint[i] = 2 * dropfoot - point
            ElementsPositionMatrix[key] = screenpoint
            if not justScreen:
                ElementsPositionMatrix[key] = np.unique(
                    np.vstack((matrix, screenpoint)), axis=0)

            self.parameters['AtomsInfo'][key] = len(
                ElementsPositionMatrix[key])
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('screen:\n' + str(np.array(screenPlane)))

    def autoMove2Zero(self):
        '''
        自动移动到原点
        '''
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )

        for key in ElementsPositionMatrix.keys():
            matrix = ElementsPositionMatrix[key]
            minA = min(matrix[:, 0])
            minB = min(matrix[:, 1])
            minC = min(matrix[:, 2])
            delta = np.array([minA, minB, minC])
            for key1 in ElementsPositionMatrix.keys():
                ElementsPositionMatrix[key1] -= delta
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('autoMove2Zero:YES')

        self.cutLattice([[0, 1], [0, 1], [0, 1]], flag=False)

    def autoMove2ZeroByElement(self, element):
        '''
        通过某个元素自动移动到原点
        '''
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )

        matrix = ElementsPositionMatrix[element]
        minA = min(matrix[:, 0])
        minB = min(matrix[:, 1])
        minC = min(matrix[:, 2])
        delta = np.array([minA, minB, minC])
        for key in ElementsPositionMatrix.keys():
            ElementsPositionMatrix[key] -= delta
        self.parameters['ElementsPositionMatrix'] = ElementsPositionMatrix

        self.process.append('autoMove2ZeroByElement:' + element)

        self.cutLattice([[0, 1], [0, 1], [0, 1]], flag=False)

    def autoCutLatticeByElement(self, element, cor):
        LatticeMatrix = self.parameters['LatticeMatrix']
        a = np.linalg.norm(LatticeMatrix[0])
        b = np.linalg.norm(LatticeMatrix[1])
        c = np.linalg.norm(LatticeMatrix[2])
        basicvector = np.array([a, b, c])

        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        matrix = ElementsPositionMatrix[element]

        delta = np.ones(3)
        scale = np.ones(3)

        corlist = {'a': 0, 'b': 1, 'c': 2}
        i = int(corlist[cor])
        delta[i] = max(matrix[:, i]) - min(matrix[:, i])
        if delta[i] != 0:
            scale[i] = delta[i] / basicvector[i]
        LatticeMatrix[i] *= scale[i]
        self.parameters['LatticeMatrix'] = LatticeMatrix

        self.process.append('autoMove2ZeroByElement:' + element + '_' + cor)

    def getAngle(self, refElement1, refElement2):
        '''
        选定某个元素的两个原子，来给出他们相对的夹角，前者作为参考系
        refElement1=['Sn',12]
        '''
        ElementsPositionMatrix = self.parameters['ElementsPositionMatrix'].copy(
        )
        LatticeMatrix = self.parameters['LatticeMatrix']
        basicvector = np.array(
            [LatticeMatrix[0][0], LatticeMatrix[1][1], LatticeMatrix[2][2]])
        refElement1 = np.array(
            ElementsPositionMatrix[refElement1[0]][refElement1[1] - 1])
        refElement2 = np.array(
            ElementsPositionMatrix[refElement2[0]][refElement2[1] - 1])
        vector = (refElement2 - refElement1) * basicvector
        normalVector = vector / np.linalg.norm(vector)

        gamma = np.degrees(np.arccos(normalVector[2]))
        lineXY = np.sqrt(normalVector[0]**2 + normalVector[1]**2)
        alpha = np.degrees(np.arccos(normalVector[0] / lineXY))
        beta = np.degrees(np.arccos(normalVector[1] / lineXY))

        Angle = np.array([alpha, beta, gamma])

        self.process.append('getAngle:\n' + '-refElement:\n' + str(refElement1) + '\n' + str(refElement2) +
                            '\n*Return:' + '\nalpha:' + str(Angle[0]) + '\nbeta :' + str(Angle[1]) + '\ngamma:' + str(Angle[2]))

        print(Angle)
        return Angle


if __name__ == '__main__':
    '''
    测试simpleDIM的选取中间胞的功能
    '''
    sop = SOP(path='./POSCAR')
    center_atoms = sop.simpleDIM(dim=[3, 3, 3], get_center=True)
    # for key in center_atoms.keys():
    #     print(key)
    #     print(center_atoms[key])
