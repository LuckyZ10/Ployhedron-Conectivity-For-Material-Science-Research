'''
author: Yilin Zhang
---
本脚本针对给定的·POSCAR·文件来进行多面体结构信息的提取,具备的主要功能：
1.多面体间的连接情况分析
* 获取中心胞和整体胞的原子位置__init__()进行后续操作
* 自动获取指定多面体mainPoly的顶点坐标并返回相应的多面体信息[dict] getMainPolyhedrons(main_element, NN_element, NN_element_number='auto') [main_element = 'Fe',NN_element = 'O']
* 自动获取mainPoly周边指定的多面体环境envirPloy[dict] getEnvirPolyhedrons(envirPloy) [envirPloy=[['Ni', 'O', 'auto'], ['Fe', 'O', 'auto']]]
* 统计单个多面体间的连接情况[numpy] getMatrixVEF_singleAndSingle(singlePloyA, singlePloyB) singlePloy = [[Fe_x,Fe_y,Fe_z,O1_x,O1_y,O1_z],...,[Fe_x,Fe_y,Fe_z,On_x,On_y,On_z]]
* 统计某种多面体与对应中间元素的多面体群之间的连接情况[dict] getMatrixVEF_singleAndElement(singlePloy, elementPloy)
* 统计主多面体与环境多面体之间的连接情况[dict] getMatrixVEF_mainAndEnvir(mainPloy, envirPloy)
* 混合整理主多面体与所有环境多面体的真实连接情况[dict] getMergeMatrixVEF(MatrixVEF, merge='ME')

2.后处理功能
* 计算多面体的体积 calPolyVolumes(poly)
* 计算原胞的体积 calCellVolume()
* 计算主多面体共面的情况[dict] calCommonFace(MatrixVEF)
* 计算主多面体共边的情况[dict] calCommonEdge(MatrixVEF)
* 计算主多面体共顶的复杂情况[vertex_dict, half_vertex_dict] calCommonVertex(MatrixVEF, envirPloy)

3.数据分析
* 统计连接方式[dict] statisticsConnection(face,edge,vertex,halfVertex) 默认输出统计矩阵mainPloy+'_SC'命名
- 共包括9种连接方式:|共面|共线|纯共顶(两个多面体)|三共顶|四共顶|半共顶(两个多面体共边)|二体连接+一个单体连接|连接一个共体|连接的真空层数|

---
版本更新：2022/10/26
---


'''

import numpy as np
import sys
from collections import OrderedDict
from STRUCTURE.SymmetricalOperation import SOP
from STRUCTURE.ReadPOSCAR import get_poscar_info
import matplotlib.pyplot as plt


class get_hedron_info():
    def __init__(self, path='./POSCAR'):
        '''初始化结构文件所在路径和结构文件名称'''
        self.POSCAR = path
        self.path = path[0:-len(path.strip().split('/')[-1])]
        sop = SOP(path=path)
        self.center_atoms, self.envir_atoms, self.system_atoms = sop.simpleDIM(get_center=True)

    def getMainPolyhedrons(self, main_element, NN_element, NN_element_number='auto', envir=False):
        # 1.统计main_element与顶点的NN_element原素
        atomA = self.center_atoms[main_element]
        atomB = self.envir_atoms[NN_element]
        if envir:
            atomA = self.envir_atoms[main_element]
            atomB = self.system_atoms[NN_element]
        atomA_number = len(atomA)
        atomB_number = len(atomB)

        # 2.创建循环矩阵
        cycle_1_C_PositionMatrix = np.repeat(atomA, len(atomB), axis=0)
        cycle_2_C_PositionMatrix = np.tile(atomB, (len(atomA), 1))

        # 3.计算所有键长
        bond_length = (cycle_1_C_PositionMatrix - cycle_2_C_PositionMatrix)**2
        bond_length = np.sqrt(np.sum(bond_length, axis=1))
        # 4.筛选出最近邻的NN原子坐标
        PositionMatrix = np.concatenate((cycle_1_C_PositionMatrix.reshape(atomA_number, atomB_number, 3),
                                         cycle_2_C_PositionMatrix.reshape(atomA_number, atomB_number, 3)), axis=2)
        bond_length = bond_length.reshape(atomA_number, atomB_number, 1)
        PositionMatrix = np.concatenate((PositionMatrix, bond_length), axis=2)

        # 5.排序操作
        for i in range(PositionMatrix.shape[0]):
            PositionMatrix_copy = PositionMatrix[i]
            PositionMatrix[i] = PositionMatrix_copy[PositionMatrix_copy[:, -1].argsort()]
        '''
        PositionMatrix 矩阵的每行的0，1，2元素对应main_element的坐标; 3，4，5对应NN_element的坐标，最后一个对应键长
        '''

        # 6.输出多面体
        dict = OrderedDict()
        # 尝试进行自动筛选NN_element_number
        NN_element_number = str(NN_element_number)
        if NN_element_number.upper()[0] == 'A':
            bond_length = PositionMatrix[:, :, -1]
            BL_changeRatio = (bond_length[:, 1:] - bond_length[:, 0:-1]) / bond_length[:, 0:-1]
            NN_element_number = np.array(BL_changeRatio.argmax(axis=1) + 1)
            # print(NN_element_number)

            for key in set(NN_element_number):
                dict[main_element + '_' + NN_element + str(key)] = PositionMatrix[np.where(NN_element_number == key), :key, :][0]
            '''测试'''
            # for i in range(bond_length.shape[0]):
            #     print(bond_length[i, 0:15])
            # print(bond_length[19, 0:15])
            # print(BL_changeRatio[19, :])
            # print(BL_changeRatio.argmax(axis=1))
            for i in range(len(BL_changeRatio)):
                plt.plot(np.arange(len(BL_changeRatio[i, :25])) + 1, BL_changeRatio[i, :25])
            plt.savefig(main_element + '_' + NN_element + '.png')
            plt.close()

            '''测试'''

            return dict

        # 手动设置NN_element_number的返回值设置
        dict = {main_element + '_' + NN_element + str(NN_element_number): PositionMatrix[:, :int(NN_element_number), :]}
        return dict

    def getEnvirPolyhedrons(self, envirPloy):
        dict = OrderedDict()
        for item in envirPloy:
            dict[item[0]] = self.getMainPolyhedrons(main_element=item[0], NN_element=item[1], NN_element_number=item[2], envir=True)
        return dict

    def getMatrixVEF_singleAndSingle(self, singlePloyA, singlePloyB):
        # 统计singleA与singleB的连接情况
        MatrixVEF = np.zeros((singlePloyA.shape[0], singlePloyB.shape[0], singlePloyA.shape[1]))

        for i in range(singlePloyA.shape[0]):
            main_atomA = singlePloyA[i, 0, 0:3]
            NN_atomA = singlePloyA[i, :, 3:6]
            for j in range(singlePloyB.shape[0]):
                main_atomB = singlePloyB[j, 0, 0:3]  # 判断是不是同一个多面体
                NN_atomB = singlePloyB[j, :, 3:6]
                if np.linalg.norm(main_atomA - main_atomB) < 0.001:
                    continue
                for k in range(len(NN_atomA)):
                    if any((NN_atomB[:] == NN_atomA[k]).all(1)):
                        # if beside_atomsj.any(beside_atomsi[k]):
                        MatrixVEF[i][j][k] += 1

        # 统除掉非接触的多面体结构
        MatrixVEF_connect = []
        for i in range(singlePloyA.shape[0]):
            number_commom_vertexs_each_polyhedron = MatrixVEF[i]
            flag = np.where(np.sum(number_commom_vertexs_each_polyhedron, axis=1) > 0)
            number_commom_vertexs_each_polyhedron = np.concatenate((number_commom_vertexs_each_polyhedron[flag], np.array(flag).T), axis=1)
            MatrixVEF_connect.append(number_commom_vertexs_each_polyhedron.astype(int))
        return MatrixVEF_connect

        return MatrixVEF.astype(int)

    def getMatrixVEF_singleAndElement(self, singlePloy, elementPloy):
        dict = OrderedDict()
        singlePloyA = singlePloy
        for key in elementPloy.keys():
            singlePloyB = elementPloy[key]
            dict[key] = self.getMatrixVEF_singleAndSingle(singlePloyA, singlePloyB)

        return dict

    def getMatrixVEF_mainAndEnvir(self, mainPloy, envirPloy):
        dict = OrderedDict()
        for key in mainPloy.keys():  # 提取mainPloy中的一个
            singlePloy = mainPloy[key]
            for element in envirPloy.keys():  # 提取envirPloy中的一种元素的elementPloy
                elementPloy = envirPloy[element]
                dict[key + '-' + element] = self.getMatrixVEF_singleAndElement(singlePloy, elementPloy)
        return dict

    def getMergeMatrixVEF(self, MatrixVEF, merge='ME'):
        if merge == 'ME':
            # print(MatrixVEF)
            # 重新整理关键词
            dict = OrderedDict()
            for key0 in MatrixVEF.keys():
                for key1 in MatrixVEF[key0].keys():
                    MatrixVEF_SS = np.array(MatrixVEF[key0][key1], dtype=object)
                    # 对MatrixVEF_SS后面贴多面体标签tag
                    MatrixVEF_SS_list = []
                    for i in range(len(MatrixVEF_SS)):
                        key_add = np.repeat([key1], MatrixVEF_SS[i].shape[0]).T
                        MatrixVEF_SS_i = np.vstack((MatrixVEF_SS[i].T, key_add)).T
                        MatrixVEF_SS_list.append(MatrixVEF_SS_i)

                    dict[key0 + '-' + key1] = np.array(MatrixVEF_SS_list, dtype=object)

            # 获取下mainPloy的名称
            mainPloy_name = []
            [mainPloy_name.append(key.strip().split('-')[0])for key in dict.keys()]
            mainPloy_name = set(mainPloy_name)

            # 收集与mainPloy相连的envirPoly
            mergeMatrixVEF = OrderedDict()
            for key0 in mainPloy_name:
                list = []
                for key in dict.keys():
                    key_mark = key.strip().split('-')[0]

                    if key0 == key_mark:
                        if dict[key].size > 0:
                            list.append(dict[key])
                SameMatrix = list[0]  # huode zishen duomianti lianjieqingkuang

                # TEST1

                SameMatrix_list = []
                for i in range(len(list[0])):
                    SameMatrix_list.append(SameMatrix[i])

                for i in range(1, len(list)):
                    DiffMatrix = np.array(list[i])

                    # if DiffMatrix.ndim < 3:
                    # if Unkownbug:
                    for j in range(len(list[0])):
                        SameMatrix_list[j] = np.vstack((SameMatrix_list[j], DiffMatrix[j]))
                    SameMatrix = np.array(SameMatrix_list, dtype=object)
                # TEST1

                # test2
                # for i in range(1, len(list)):
                #     DiffMatrix = np.array(list[i])
                #     # print(DiffMatrix.ndim)
                #     # print(SameMatrix)
                #     if DiffMatrix.ndim < 3:
                #         # if Unkownbug:
                #         for j in range(len(list[0])):
                #             SameMatrix[j] = np.vstack((SameMatrix[j], DiffMatrix[j]))
                #     else:
                #         SameMatrix = np.concatenate((SameMatrix, DiffMatrix), axis=1)
                # test2
                # print(SameMatrix)
                mergeMatrixVEF[key0] = SameMatrix
                # print(mergeMatrixVEF)
            return mergeMatrixVEF

    def calCommonFace(self, MatrixVEF):
        face_dict = OrderedDict()
        for key in MatrixVEF.keys():
            face_list = []
            for i in range(len(MatrixVEF[key])):
                MatrixVEF_i = MatrixVEF[key][i]
                # 拆分几何连接的finger 和对应的多面体ploys_tag
                finger_i = MatrixVEF_i[:, :-2].astype(int)
                # polys_tag_i = MatrixVEF_i[:, -2:]
                # 对公共拥有的点求和
                common_points_number_i = np.sum(finger_i, axis=1)
                # 筛选出共面
                common_points_number_i_copy = np.zeros((common_points_number_i.shape[0]))
                common_points_number_i_copy[np.where(common_points_number_i[:] > 2)] = 1
                face_list.append(common_points_number_i_copy)
            face_dict[key] = np.array(face_list, dtype=object)

        return face_dict

    def calCommonEdge(self, MatrixVEF):
        edge_dict = OrderedDict()
        for key in MatrixVEF.keys():
            edge_list = []
            for i in range(len(MatrixVEF[key])):
                MatrixVEF_i = MatrixVEF[key][i]
                # 拆分几何连接的finger 和对应的多面体ploys_tag
                finger_i = MatrixVEF_i[:, :-2].astype(int)
                # polys_tag_i = MatrixVEF_i[:, -2:]
                # 对公共拥有的点求和
                common_points_number_i = np.sum(finger_i, axis=1)
                # 筛选出共面
                common_points_number_i_copy = np.zeros((common_points_number_i.shape[0]))
                common_points_number_i_copy[np.where(common_points_number_i[:] == 2)] = 1
                edge_list.append(common_points_number_i_copy)
            edge_dict[key] = np.array(edge_list, dtype=object)
        return edge_dict

    def calCommonVertex(self, MatrixVEF, envirPloy):
        # 1.通过共边共面进行判定可能存在点特征的多面体原子
        def selectCommonPoint(MatrixVEF_i, poly_face_edge_i):
            # 拆分几何连接的finger 和对应的多面体ploys_tag
            finger_i = MatrixVEF_i[:, :-2].astype(int)

            flag = finger_i[poly_face_edge_i > 0]  # 如果顶点都存在则说明不存在点情况了

            # 创建可能的点特征统计矩阵
            possible_common_vertex = np.ones((finger_i.shape[1])).astype(int)

            if flag.size != 0:
                possible_common_vertex[np.sum(flag, axis=0) > 0] = 0
            return possible_common_vertex

        # 0.主要部分
        vacuum_Layer_dict = OrderedDict()
        vertex_dict = OrderedDict()
        half_vertex_dict = OrderedDict()
        # 获取共面和共边的情况
        face = ghi.calCommonFace(MatrixVEF)
        edge = ghi.calCommonEdge(MatrixVEF)
        for key in MatrixVEF.keys():
            # 共边共面的情况
            poly_face_edge = face[key] + edge[key]

            # vacuum_Layer
            vacuum_Layer = np.zeros((MatrixVEF[key].shape[0], 1))

            # 0索引为纯共顶(两个多面体)的，1索引为三共顶的点，2为四共顶的点
            vertex_matrix = np.zeros((MatrixVEF[key].shape[0], 3))

            # 0索引为半共顶(两个多面体共边)数，1索引为二体连接+一个单体连接，2索引为共体连接的
            half_vertex_matrix = np.zeros((MatrixVEF[key].shape[0], 3))

            # 1.通过共边共面进行判定可能存在点特征的多面体原子
            possible_common_vertex_matrix = []
            for i in range(len(MatrixVEF[key])):
                MatrixVEF_i = MatrixVEF[key][i]
                poly_face_edge_i = poly_face_edge[i]
                possible_common_vertex_matrix.append(selectCommonPoint(MatrixVEF_i, poly_face_edge_i))
            possible_common_vertex_matrix = np.array(possible_common_vertex_matrix)
            # print(possible_common_vertex_matrix)

            # 开始进行不同点特征的区分
            for i in range(possible_common_vertex_matrix.shape[0]):
                poly = MatrixVEF[key][i][:, 0:-2].astype(int)
                polys_tag = MatrixVEF[key][i][:, -2:]

                for j in range(possible_common_vertex_matrix.shape[1]):
                    # 设定共点的那行
                    point = np.zeros(possible_common_vertex_matrix.shape[1]).astype(int)
                    if possible_common_vertex_matrix[i][j] == 0:
                        continue
                    else:
                        point[j] = 1

                    # 统计共点的数据
                    point_flag = np.abs(poly - point)
                    # print(point_flag)
                    beside_polys_number = np.sum(np.sum(point_flag, axis=1) < 1)

                    beside_polys = polys_tag[np.sum(point_flag, axis=1) < 1]
                    # 合并周围的多面体的氧原子
                    # print(envirPloy['Al'])
                    if len(beside_polys) == 0:
                        vacuum_Layer[i][0] += 1
                        continue

                    merge_beside_polys = envirPloy[beside_polys[0][-1].split('_')[0]][beside_polys[0][-1]][int(beside_polys[0][0])][:, 3:6]
                    for k in range(1, len(beside_polys)):
                        beside_polys_k = envirPloy[beside_polys[k][-1].split('_')[0]][beside_polys[k][-1]][int(beside_polys[k][0])][:, 3:6]
                        merge_beside_polys = np.concatenate((merge_beside_polys, beside_polys_k), axis=0)

                    # 统计纯共点数和半共点数
                    extraPoints = len(merge_beside_polys) - len(np.unique(merge_beside_polys, axis=0))  # beside多面体间的多余点数(删除重复的只保留一个)
                    # print(extraPoints)
                    # print(beside_polys_number)
                    # 纯共顶统计
                    if (extraPoints + 1) == beside_polys_number:
                        # if beside_polys_number == 1:
                        #     vertex_matrix[0] += 1
                        # if beside_polys_number == 2:
                        #     vertex_matrix[1] += 1
                        # if beside_polys_number == 3:
                        #     vertex_matrix[2] += 1

                        vertex_matrix[i][beside_polys_number - 1] += 1
                        # print(vertex_matrix[i][beside_polys_number - 1])

                    # 统计半共点数
                    if (extraPoints + 1) > beside_polys_number:
                        beside_common_points = (extraPoints + 1) - beside_polys_number  # beside多面体间的共点数
                        if beside_polys_number == 2:
                            half_vertex_matrix[i][0] += 1
                        if beside_polys_number == 3:
                            if beside_common_points == 1:
                                half_vertex_matrix[i][1] += 1
                            if beside_common_points > 2:
                                half_vertex_matrix[i][2] += 1
            vacuum_Layer_dict[key] = vacuum_Layer
            vertex_dict[key] = vertex_matrix
            half_vertex_dict[key] = half_vertex_matrix
        return vacuum_Layer_dict, vertex_dict, half_vertex_dict

    def calPolyVolumes(self, poly):
        from scipy.spatial import ConvexHull
        poly_volumes_dict = OrderedDict()
        for key in poly.keys():
            hedron = poly[key]
            # print(hedron)
            ploy_volumes = np.zeros((hedron.shape[0], 4))

            for i in range(hedron.shape[0]):
                ploy_volumes[i, 0: 3] = hedron[i, 0, 0: 3]
                # print(ploy_volumes)
                if len(hedron[i, :, 3: 6]) < 3:
                    ploy_volumes[i, 3] = 0
                    print('In a straight line ')
                elif len(hedron[i, :, 3: 6]) == 3:
                    try:
                        ploy_volumes[i, 3] = ConvexHull(np.vstack((hedron[i, :, 3: 6], ploy_volumes[i, 0: 3]))).volume
                    except:
                        ploy_volumes[i, 3] = 0
                        print('In plane of triangle')
                else:
                    try:
                        ploy_volumes[i, 3] = ConvexHull(hedron[i, :, 3: 6]).volume
                    except:
                        ploy_volumes[i, 3] = 0
                        print('In plane of triangle')
                # print('ok')
            poly_volumes_dict[key] = ploy_volumes
        return poly_volumes_dict

    def calCellVolume(self,):
        gpi = get_poscar_info(path=self.POSCAR)
        return gpi.calVolume()

    def statisticsConnection(self, face, edge, vertex, halfVertex, vacuum_Layer, output=True):
        SC_dict = OrderedDict()
        data_files = []
        for key in face.keys():
            # 构建统计矩阵
            SC_matrix = np.zeros((face[key].shape[0], 9))
            # print(SC_matrix)
            # 统计共面数
            for i in range(len(face[key])):
                SC_matrix[i][0] = np.sum(face[key][i])
            # 统计共边数
            for i in range(len(edge[key])):
                SC_matrix[i][1] = np.sum(edge[key][i])

            # 统计vacuum_layer
            SC_matrix[:, 8] = vacuum_Layer[key].T
            # 统计纯共点数
            SC_matrix[:, 2:5] = vertex[key]
            # 统计半共点数
            SC_matrix[:, 5:8] = halfVertex[key]
            SC_dict[key] = SC_matrix

            if output:
                np.savetxt(key + '_SC', SC_matrix, fmt="%d")
                data_files.append(key + '_SC')
            # 统计平均数

            SC_matrix[:, 0:2] = SC_matrix[:, 0:2] / 2  # 共面，共边
            SC_matrix[:, 3] = SC_matrix[:, 3] / 3  # 三体共点平均
            SC_matrix[:, 4] = SC_matrix[:, 4] / 4  # 四体共点平均

            SC_matrix_mean = np.mean(SC_matrix, axis=0)
            print('包括9种连接方式:|共面|共线|纯共顶(两个多面体)|三共顶|四共顶|半共顶(两个多面体共边)|二体连接+一个单体连接|连接一个共体|连接的真空层数')
            print(key + ':' + 'per polyhedron')
            print(SC_matrix_mean)

        return SC_dict, data_files


if __name__ == "__main__":
    '''
    测试1：Fe8Ni4O16 没问题
    测试2：Cu8V8O28  没问题
    测试3：alpha.vasp 没问题
    POSCAR:
    '''

    '''
    genju keti xuyao zuo de DIY
    '''

    atoms_info = get_poscar_info(path='./POSCAR').getAtomsInfo()
    envirploy_list = []
    NN_atom = []
    main_atom = []
    for key in atoms_info.keys():
        if key == 'O' or key == 'F':
            NN_atom.append(key)
        else:
            main_atom.append(key)
    for keyN in NN_atom:
        for keyM in main_atom:
            envirploy_list.append([keyM, keyN, 'auto'])
    '''
    jin xing tong ji
    '''

    # unkownbug = int(sys.argv[1])
    # unkownbug = 0
    # if unkownbug == 0:
    #     unkownbug = False
    # if unkownbug == 1:
    #     unkownbug = True

    ghi = get_hedron_info(path='./POSCAR')
    data_file_summary = []

    mainPloyVolume_list = []

    for atom_pair in envirploy_list:
        print(atom_pair)
        atomM = atom_pair[0]
        atomN = atom_pair[1]
        mainPloy = ghi.getMainPolyhedrons(main_element=atomM, NN_element=atomN)  # 获取指定的多面体
        mainPloyVolumes = ghi.calPolyVolumes(poly=mainPloy)  # 计算多面体体积
        mainPloyVolume_list.append(mainPloyVolumes)

        envirPloy = ghi.getEnvirPolyhedrons(envirPloy=envirploy_list)
        #
        MatrixVEF_ME = ghi.getMatrixVEF_mainAndEnvir(mainPloy=mainPloy, envirPloy=envirPloy)
        mergeMatrixVEF = ghi.getMergeMatrixVEF(MatrixVEF=MatrixVEF_ME, merge='ME')  # Unkownbug
        # print(mergeMatrixVEF)
        face = ghi.calCommonFace(MatrixVEF=mergeMatrixVEF)
        edge = ghi.calCommonEdge(MatrixVEF=mergeMatrixVEF)
        # ghi.calCommonVertex(MatrixVEF=mergeMatrixVEF, envirPloy=envirPloy)
        vacuum_Layer, vertex, half_vertex = ghi.calCommonVertex(MatrixVEF=mergeMatrixVEF, envirPloy=envirPloy)
        SC_dict, data_files = ghi.statisticsConnection(face=face, edge=edge, vertex=vertex, halfVertex=half_vertex, vacuum_Layer=vacuum_Layer)
        data_file_summary.append(data_files)
    '''tezheng chuli'''

    '''获得体积参数'''
    cellVolume = ghi.calCellVolume()
    # mainPloyVolumeSum = np.sum(np.array(mainPloyVolume_list).astype(float))
    polyV = 0
    for poly in mainPloyVolume_list:
        for key in poly.keys():
            polyV += np.sum(poly[key], axis=0)[-1]

    '''获得基于体积的特征描述符'''
    from pymatgen.core.periodic_table import Element as EB
    cations_mass = 0
    anions_mass = 0
    atom_numer = 0
    for key, value in atoms_info.items():
        atom_numer += value
        if key == 'O' or key == 'F':
            anions_mass += EB(key).atomic_mass * int(value)
        else:
            cations_mass += EB(key).atomic_mass * int(value)

    atom_numer_perV = atom_numer / cellVolume
    cations_mass_perV = cations_mass / cellVolume
    anions_mass_perV = anions_mass / cellVolume
    total_masss_perV = (cations_mass + anions_mass) / cellVolume
    CA_ratio = cations_mass / anions_mass
    vaccum_ratio = 1 - polyV / cellVolume

    with open('SpaceFeature', 'w') as f:
        f.write('atom_numer_perV:  %f\n' % atom_numer_perV)
        f.write('cations_mass_perV:  %f\n' % cations_mass_perV)
        f.write('anions_mass_perV:  %f\n' % anions_mass_perV)
        f.write('total_masss_perV:  %f\n' % total_masss_perV)
        f.write('CA_ratio:  %f\n' % CA_ratio)
        f.write('vaccum_ratio:  %f\n' % vaccum_ratio)

    '''
    xiamian wei jingdian de jiaocheng
    '''
    # ghi = get_hedron_info(path='./POSCAR')
    # mainPloy = ghi.getMainPolyhedrons(main_element='Ag', NN_element='O')  # 获取指定的多面体
    # # mainPloyVolumes = ghi.calPolyVolumes(poly=mainPloy)  # 计算多面体体积
    # envirPloy = ghi.getEnvirPolyhedrons(envirPloy=[['Ag', 'O', 'auto']])
    #
    # MatrixVEF_ME = ghi.getMatrixVEF_mainAndEnvir(mainPloy=mainPloy, envirPloy=envirPloy)
    # mergeMatrixVEF = ghi.getMergeMatrixVEF(MatrixVEF=MatrixVEF_ME, merge='ME')
    # # print(mergeMatrixVEF)
    # face = ghi.calCommonFace(MatrixVEF=mergeMatrixVEF)
    # edge = ghi.calCommonEdge(MatrixVEF=mergeMatrixVEF)
    # # ghi.calCommonVertex(MatrixVEF=mergeMatrixVEF, envirPloy=envirPloy)
    # vertex, half_vertex = ghi.calCommonVertex(MatrixVEF=mergeMatrixVEF, envirPloy=envirPloy)
    # SC_dict = ghi.statisticsConnection(face=face, edge=edge, vertex=vertex, halfVertex=half_vertex)

    '''
    一些diy的后处理
    '''
    # # 计算对真空率
    # poly = ghi.calPolyVolumes(mainPloy)
    # polyV = 0
    # for key in poly.keys():
    #     polyV += np.sum(poly[key], axis=0)[-1]
    #
    # cellV = ghi.calCellVolume()
    # spcaeV = cellV - polyV
    #
    # # 统计连接数
    # connect = []
    # for key in SC_dict.keys():
    #     connect_mean = np.sum(SC_dict[key], axis=0)
    #     connect_mean[0:3] = connect_mean[0:3] / 2
    #     connect_mean[3] = connect_mean[3] / 3  # 三体共点平均
    #     connect_mean[4] = connect_mean[4] / 4  # 四体共点平均
    #     connect.append(connect_mean)
    # connect = np.array(connect).flatten()
    #
    # # 计算对外压缩率
    # print(connect / spcaeV)

    '''
    30172: 0.021 0.066 15.99
    30136: 0.025 0.063 16.20
    10148: 0.023 0.056 23.49
    10136: 0.025 0.049 25.60
      pri: 0.019 0.077 14.94
    '''
