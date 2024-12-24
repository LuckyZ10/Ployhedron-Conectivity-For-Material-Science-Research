# Ployhedron-Conectivity-For-Material-Science-Research
本脚本针对给定的·POSCAR·文件来进行多面体结构信息的提取,具备的主要功能：
1. 多面体间的连接情况分析
   * 获取中心胞和整体胞的原子位置init()进行后续操作
   * 自动获取指定多面体mainPoly的顶点坐标并返回相应的多面体信息[dict] 
     getMainPolyhedrons(main_element, NN_element, NN_element_number='auto') 
     [main_element = 'Fe',NN_element = 'O']
   * 自动获取mainPoly周边指定的多面体环境envirPloy[dict] getEnvirPolyhedrons(envirPloy) 
      [envirPloy=[['Ni', 'O', 'auto'], ['Fe', 'O', 'auto']]]
   * 统计单个多面体间的连接情况[numpy] getMatrixVEF_singleAndSingle(singlePloyA, singlePloyB) 
      singlePloy = [[Fe_x,Fe_y,Fe_z,O1_x,O1_y,O1_z],...,[Fe_x,Fe_y,Fe_z,On_x,On_y,On_z]]
   * 统计某种多面体与对应中间元素的多面体群之间的连接情况[dict] 
      getMatrixVEF_singleAndElement(singlePloy, elementPloy)
   * 统计主多面体与环境多面体之间的连接情况[dict] getMatrixVEF_mainAndEnvir(mainPloy, envirPloy)
   * 混合整理主多面体与所有环境多面体的真实连接情况[dict] getMergeMatrixVEF(MatrixVEF, merge='ME')
2. 后处理功能
   * 计算多面体的体积 calPolyVolumes(poly)
   * 计算原胞的体积 calCellVolume()
   * 计算主多面体共面的情况[dict] calCommonFace(MatrixVEF)
   * 计算主多面体共边的情况[dict] calCommonEdge(MatrixVEF)
   * 计算主多面体共顶的复杂情况[vertex_dict, half_vertex_dict] calCommonVertex(MatrixVEF, envirPloy)
3. 数据分析
  * 统计连接方式[dict] statisticsConnection(face,edge,vertex,halfVertex) 默认输出统计矩阵mainPloy+'_SC'命名
  * 共包括9种连接方式:|共面|共线|纯共顶(两个多面体)|三共顶|四共顶|半共顶(两个多面体共边)|二体连接+一个单体连接|连接一个共体|连接的真空层数|

