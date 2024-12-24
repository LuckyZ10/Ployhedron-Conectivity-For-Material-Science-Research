# Ployhedron-Conectivity-For-Material-Science-Research
This script extracts polyhedron structure information for a given POSCAR file. It has the following main functions:
1. Analysis of the connection between polyhedrons
   * Get the atomic positions of the central cell and the whole cell init () for subsequent operations
   * Automatically obtain the vertex coordinates of the specified polyhedron mainPoly and return the corresponding polyhedron information [dict]
     getMainPolyhedrons(main_element, NN_element, NN_element_number='auto') 
     [main_element = 'Fe',NN_element = 'O']
   * Automatically obtain the specified polyhedral environment around mainPoly. envirPloy[dict] getEnvirPolyhedrons(envirPloy) 
      [envirPloy=[['Ni', 'O', 'auto'], ['Fe', 'O', 'auto']]]
   * Statistics of the connection between individual polyhedrons. [numpy] getMatrixVEF_singleAndSingle(singlePloyA, singlePloyB) 
      singlePloy = [[Fe_x,Fe_y,Fe_z,O1_x,O1_y,O1_z],...,[Fe_x,Fe_y,Fe_z,On_x,On_y,On_z]]
   * Statistics on the connection between a certain polyhedron and the polyhedral group corresponding to the intermediate element. [dict] 
      getMatrixVEF_singleAndElement(singlePloy, elementPloy)
   * Statistics on the connection between the main polyhedron and the environmental polyhedron. [dict] getMatrixVEF_mainAndEnvir(mainPloy, envirPloy)
   * Mixedly sort out the real connection situation of the main polyhedron and all environmental polyhedrons. [dict] getMergeMatrixVEF(MatrixVEF, merge='ME')
2. Post-processing function
   * Calculate the volume of a polyhedron. calPolyVolumes(poly)
   * Calculate the volume of the primitive cell. calCellVolume()
   * Calculate the coplanar situation of the main polyhedron. [dict] calCommonFace(MatrixVEF)
   * Calculate the situation of common edges of the main polyhedron. [dict] calCommonEdge(MatrixVEF)
   * Calculate the complex situation of copoints of the main polyhedron. [vertex_dict, half_vertex_dict] calCommonVertex(MatrixVEF, envirPloy)
3. Data analytics
  * Count connection methods. [dict] statisticsConnection(face,edge,vertex,halfVertex) Default output statistical matrix mainPloy + '_SC' naming
  * A total of 9 connection methods are included: | Coplanar | Collinear | Pure cotopical (two polyhedra) | Tricotopical | Tetracotopical | Semicomponent (two polyhedra co-edge) | Disomeric connection + one monomer connection | Connect a community | Number of vacuum layers connected |

**Please cite the following paper if using the program:**
[Evaluating thermal expansion in fluorides and oxides: Machine learning predictions with connectivity descriptors](https://iopscience.iop.org/article/10.1088/1674-1056/accdca/meta)

 
