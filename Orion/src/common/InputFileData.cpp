#include "InputFileData.h"

InputFileData::InputFileData(std::string& filename,std::ofstream& logFile) :
  pointDispNormalSet(0),pointForceNormalSet(0),lineDispNormalSet(0),
  lineForceNormalSet(0),pointRotNormalSet(0),pointElectricNormalSet(0),
  pointMomentNormalSet(0),lineRotNormalSet(0),lineElectricNormalSet(0),
  lineMomentNormalSet(0) {

  using namespace std;

  // Set default values.
  setDefaultValues(logFile);

  // Read and store the input file.
  readInputFile(filename,logFile);

}

InputFileData::InputFileData(std::ofstream& logFile) :
  pointDispNormalSet(0),pointForceNormalSet(0),lineDispNormalSet(0),
  lineForceNormalSet(0),pointRotNormalSet(0),pointElectricNormalSet(0),
  pointMomentNormalSet(0),lineRotNormalSet(0),lineElectricNormalSet(0),
  lineMomentNormalSet(0) {

  using namespace std;

  // Set default values.
  setDefaultValues(logFile);

  std::string filename = "input.dat";

  // Read and store the input file.
  readInputFile(filename,logFile);

}

/**********************************************************************/
/**********************************************************************/
// Read and store the input file.
void InputFileData::readInputFile(std::string& filename,std::ofstream& logFile) {

  using namespace std;

//  string filename = "input.dat";

  // Open input file.
  ifstream inputFile(filename.c_str());
  if(!inputFile) {
    logFile<<"Can't open input file " << filename << "!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(inputFile.eof()) {
    logFile<<"Input file " << filename << " contains no data!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /**********************************************************************/
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  string name,token;
  double value,value1,value2,value3;

  logFile<<"###################################################"<<endl;
  logFile<<"############# input file data #####################"<<endl;

  /**********************************************************************/
  // Read material data.
  while(inputFile >> name && name != "END") {

    if(name == "FEM_GEOMETRY" || name == "MESHLESS_GEOMETRY"
       || name == "MODEL_INPUT" || name == "CALCULATION_CONTROL_INPUT"
       || name == "EQUATION_SOLVER_CONTROL_INPUT"
       || name == "POSTPROCESSING_INPUT" ||name == "graphs") {

      if(materials.size() < 1) {
	logFile<<"No material data has been stored while reading "
	       <<"input data!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      logFile<<name<<endl;

      break;
    }
    else if(name == "" || name == "MATERIAL_INPUT") {

      logFile<<name<<endl;

      continue;
    }
    else if(name == "Material") {

      inputFile >> value;
      materials.resize(materials.size()+1);
      logFile<<name<<" "<< value<<endl;

      // store macro ID
      materials[materials.size()-1]["macroID"] =
	materials.size();

    }
    else if(name == "microID") {

      inputFile >> value;

      // micro ID > 1
      if(value > 1) {
	materials.resize(materials.size()+1);

	materials[materials.size()-1]["macroID"] =
	  materials[materials.size()-2]["macroID"];

      }

      materials[materials.size()-1][name] = value;
      logFile<<name<<" "<< value<<endl;

    }
    else if(name == "Description") {

      logFile<<name<<" ";	//writes Description
      inputFile >> name;
      logFile<<name<<endl;	//writes e.g. "Linear-Hyper-Elastic(1a)"
    }

    else {
      inputFile >> value;
      materials[materials.size()-1][name] = value;

      logFile<<name<<" "<<value<<endl;		// writes e.g. "YoungsModulus 1000"
    }

  }

  // check whether each material has a default microID

  for(int i=0;i<materials.size();i++) {

    // micro ID
    map<string,double>::iterator p = materials[i].find("microID");

    if(p == materials[i].end())

      materials[i]["microID"] = 1;

    // micro volume ratio
    p = materials[i].find("microVolumeRatio");

    if(p == materials[i].end())

      materials[i]["microVolumeRatio"] = 1.0;

    // micro constitutive law ID
    p = materials[i].find("microConstitutiveLawID");

    if(p == materials[i].end())

      materials[i]["microConstitutiveLawID"] =
	materials[i]["constitutiveLawID"];

  }

  /**********************************************************************/
  // Read FEM-geometry.
  while(inputFile >> name && name != "END") {

    if(name == "MESHLESS_GEOMETRY"
       || name == "MODEL_INPUT"
       || name == "CALCULATION_CONTROL_INPUT"
       || name == "EQUATION_SOLVER_CONTROL_INPUT"
       || name == "POSTPROCESSING_INPUT"
       || name == "graphs"
       || name == "Graphs"
       || name == "Resultant-Internal-Traction"
       || name == "Resultant-Internal-Surface-Torque"
       || name == "BOUNDARY_CONDITIONS_INPUT")
      break;

    else if(name == "FEM_GEOMETRY" || name == "") {

      logFile<<name<<endl;

      continue;
    }
    else {
      inputFile >> value;

      if(name == "backGroundMeshElemType")

	backGroundMeshInfo["elemType"] = value;

      else if(name == "backGroundMeshElemOrder")

	backGroundMeshInfo["elemOrder"] = value;

      else if(name == "volumeGaussQuadratureOrder" ||
	      name == "surfaceGaussQuadratureOrder" ||
	      name == "lineGaussQuadratureOrder")

	backGroundMeshInfo[name] = value;

      else

	problemData[name] = value;


      logFile<<name<<" "<<value<<endl;

    }
  }

  // --------------------------------------------------------------------
  // set nodes and Gauss points per element
  map<string,double> params;
  intVector data(4);

  data[0] = (int)backGroundMeshInfo["elemType"];
  data[1] = (int)backGroundMeshInfo["elemOrder"];
  getFEMMeshData(data,params);

  int nodesPerElem = (int)params["nodesPerVolumeElement"];
  int nodesPerSurfaceElem = (int)params["nodesPerSurfaceElement"];
  int nodesPerLineElem = (int)params["nodesPerLineElement"];

  backGroundMeshInfo["nodesPerVolumeElement"] =
    params["nodesPerVolumeElement"];
  backGroundMeshInfo["nodesPerSurfaceElement"] =
    params["nodesPerSurfaceElement"];
  backGroundMeshInfo["nodesPerLineElement"] =
    params["nodesPerLineElement"];


  data[0] = (int)backGroundMeshInfo["elemType"];
  data[1] = (int)backGroundMeshInfo["volumeGaussQuadratureOrder"];
  data[2] = (int)backGroundMeshInfo["surfaceGaussQuadratureOrder"];
  data[3] = (int)backGroundMeshInfo["lineGaussQuadratureOrder"];

  getGaussQuadratureData(data,params);

  backGroundMeshInfo["gaussPointsPerVolumeElement"] =
    params["gaussPointsPerVolumeElement"];
  backGroundMeshInfo["gaussPointsPerSurfaceElement"] =
    params["gaussPointsPerSurfaceElement"];
  backGroundMeshInfo["gaussPointsPerLineElement"] =
    params["gaussPointsPerLineElement"];

  /**********************************************************************/
  // Read meshless geometry, common calculation and
  // post-processing.
  while(inputFile >> name && name != "END") {

    if(name == "graphs" ||
       name == "Graphs" ||
       name == "Resultant-Internal-Traction" ||
       name == "Resultant-Internal-Surface-Torque" ||
       name == "BOUNDARY_CONDITIONS_INPUT")
      break;

    else if(name == "MESHLESS_GEOMETRY"
            || name == "MODEL_INPUT"
            || name == "CALCULATION_CONTROL_INPUT"
            || name == "EQUATION_SOLVER_CONTROL_INPUT"
            || name == "POSTPROCESSING_INPUT"
            || name == "") {

      logFile<<name<<endl;

      continue;
    }
    else {

      inputFile >> value;
      problemData[name] = value;

      logFile<<name<<" "<<value<<endl;

    }

  }

  /**********************************************************************/
  // Read in the graph input for plotting at selected nodes.
  int coord,particle,type;

  int count = 0;
  int graphNum = 0;


  while(name != "END") {

    if(name == "Resultant-Internal-Traction" ||
       name == "Resultant-Internal-Surface-Torque" ||
       name == "BOUNDARY_CONDITIONS_INPUT") {

      if(graphs.size() < 1) {
	graphs.resize(1);
	graphs[0].push_back(1);//node
	graphs[0].push_back(1);//type
	graphs[0].push_back(1);//dof
	graphNum = 1;
      }

      break;
    }

    else if(name == "graphs" || name == "Graphs" || name == "") {
      inputFile >> name;
      continue;

    }
    else if(name == "node") {
      inputFile >> particle; // particle ID
      graphs.resize(graphs.size()+1);
      graphs[count].push_back(particle);

      inputFile >> name;
      continue;
    }
    else if(name == "type") {
      inputFile >> type; // type of graph
      inputFile >> name; // first dof

      if(type == 2)
	problemData["twistGraphPlot"] = 1.0;

      if(type == 3)
	problemData["strainGraphPlot"] = 1.0;

      else if(type == 4) {
	problemData["strainGraphPlot"] = 1.0;
	problemData["stressGraphPlot"] = 1.0;
      }
      else if(type == 7) {
	problemData["strainGraphPlot"] = 1.0;
	problemData["stressGraphPlot"] = 1.0;
      }

      if(type == 12) {
        problemData["activeStressGraphPlot"] = 1.0;
        graphNum++;
        graphs[count].push_back(type);
      }

      if(type == 13) {
        problemData["sarcomereLengthGraphPlot"] = 1.0;
        graphNum++;
        graphs[count].push_back(type);
      }

      if(type == 14) {
        problemData["contractileLengthGraphPlot"] = 1.0;
        graphNum++;
        graphs[count].push_back(type);
      }

      // load-deformation or stress strain
      if(name == "dof") {
	graphs[count].push_back(type);

	while(name == "dof") {
	  inputFile >> coord;
	  inputFile >> name;
	  graphs[count].push_back(coord);
	  graphNum++;

	}

      }
      else if((type == 1 || type == 3 ||
	       type == 4 || type == 5 ||
	       type == 6) &&
	      name != "dof") {
	logFile<<"Incorrect 'graphs' input data!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
      // else if// energy ....

      count++;
      continue;
    }

    else if(name == "Resultant-Internal-Traction" ||
	    name == "Resultant-Internal-Surface-Torque" ||
	    name == "BOUNDARY_CONDITIONS_INPUT")
      break;

    else {
      logFile<<"Incorrect 'graphs' input data!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }


  if((bool)problemData["plotSystemEnergy"] &&
     problemData["dynamicCalculation"] > 0)

    problemData["numberOfGraphs"] = graphNum + 1;

  else {

    problemData["numberOfGraphs"] = graphNum;
    problemData["plotSystemEnergy"] = 0;
  }

  // Write output in the log file.
  logFile<<"Graphs"<<endl;
  logFile<<"number of graphs: "<<problemData["numberOfGraphs"] <<endl;
  for(int i=0;i<graphs.size();i++) {
    int j = 0;
    while(j < graphs[i].size()) {
      if(j == 0) {
	logFile<<"node: "<<graphs[i][j]<<" ";
	j++;
	continue;
      }
      else if(j == 1) {
	logFile<<"graph type: "<<graphs[i][j]<<" ";
	j++;
	continue;
      }
      else if(j > 1) {
	logFile<<"dof "<<graphs[i][j]<<" ";
	j++;
	continue;
      }
    }
    logFile<<endl;
  }

  // --------------------------------------------------------------------
  // resultant force and torque on the specified surfaces
  int currentElem,currentSurface,currentVolume,currentCond,elem,volume;
  bool resultantInternalTraction;
  bool resultantInternalSurfaceTorque;

  intVector node((nodesPerElem > nodesPerSurfaceElem) ? nodesPerElem
		 : nodesPerSurfaceElem);

  while(name != "END") {

    if(name == "LOAD_CONDITIONS_INPUT" ||
       name == "BOUNDARY_CONDITIONS_INPUT")
      break;

    else if(name == "") {
      inputFile >> name;
      continue;
    }

    else {


      if(name == "Resultant-Internal-Traction") {
        resultantInternalTraction = true;

	while(inputFile >> name) {

	  if(name == "Surface") {
	    inputFile >> value;
	    currentSurface = resultantForceOnSurfaces.size();
	    resultantForceOnSurfaces.resize(resultantForceOnSurfaces.size()+1);
	    continue;
	  }

	  else if(name == "nodes") {

	    inputFile >> count;

	    resultantForceOnSurfaces[currentSurface].resize(count+1);

	    inputFile >> name >> resultantForceOnSurfaces[currentSurface][0]; // dof

	    for(int j=1;j<resultantForceOnSurfaces[currentSurface].size();j++)

	      inputFile >> resultantForceOnSurfaces[currentSurface][j];

	  }

	  else if(name == "")
	    continue;

	  else if(name == "Resultant-Internal-Surface-Torque" ||
		  name == "BOUNDARY_CONDITIONS_INPUT" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;

	  else {
	    logFile<<"Incorrect written 'Resultant-Internal-Traction' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }

      if(name == "Resultant-Internal-Surface-Torque") {
        resultantInternalSurfaceTorque = true;

	while(inputFile >> name) {

	  if(name == "Surface") {
	    inputFile >> value;
	    currentSurface = resultantTorqueOnSurfaces.size();
	    resultantTorqueOnSurfaces.resize(resultantTorqueOnSurfaces.size()+1);
	    continue;
	  }


	  else if(name == "nodes") {

	    inputFile >> count;

	    resultantTorqueOnSurfaces[currentSurface].resize(count+3);

	    inputFile >> name;

	    // torque axis
	    inputFile >> resultantTorqueOnSurfaces[currentSurface][0]
		      >> resultantTorqueOnSurfaces[currentSurface][1]
		      >> resultantTorqueOnSurfaces[currentSurface][2];

	    resultantTorqueOnSurfaces[currentSurface][0] =
	      problemData["rotatingBoundaryCondsRotVecX"];
	    resultantTorqueOnSurfaces[currentSurface][1] =
	      problemData["rotatingBoundaryCondsRotVecY"];
	    resultantTorqueOnSurfaces[currentSurface][2] =
	      problemData["rotatingBoundaryCondsRotVecZ"];

	    for(int j=3;j<resultantTorqueOnSurfaces[currentSurface].size();j++)

	      inputFile >> resultantTorqueOnSurfaces[currentSurface][j];

	  }

	  else if(name == "")
	    continue;

	  else if(name == "BOUNDARY_CONDITIONS_INPUT" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;

	  else {
	    logFile<<"Incorrect written 'Resultant-Internal-Traction' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }

    }

  }

  // ---------------------------------------------------------------------
  // pressure volume graph plotting
  int oldSize = graphs.size();

  if((bool)problemData["pressureVolumePlotting"]) {

    graphs.resize(graphs.size()+1);


    graphs[oldSize].push_back(1); // dummy node
    graphs[oldSize].push_back(7); // type
    graphs[oldSize].push_back(0); // DOF

    problemData["numberOfGraphs"] += 1;
  }

  oldSize = graphs.size();

  if((bool)problemData["pressureTimePlotting"]) {

    graphs.resize(graphs.size()+1);


    graphs[oldSize].push_back(1); // dummy node
    graphs[oldSize].push_back(15); // type
    graphs[oldSize].push_back(0); // DOF

    problemData["numberOfGraphs"] += 1;
  }


  oldSize = graphs.size();

  if((bool)problemData["volumeTimePlotting"]) {

    graphs.resize(graphs.size()+1);


    graphs[oldSize].push_back(1); // dummy node
    graphs[oldSize].push_back(16); // type
    graphs[oldSize].push_back(0); // DOF

    problemData["numberOfGraphs"] += 1;
  }

  // Write output in the log file.
  logFile<<"Load-Deformed-Volume Graphs"<<endl;
  logFile<<"number of graphs: "<<problemData["numberOfGraphs"]
	 <<endl;
  for(int i=oldSize;i<graphs.size();i++) {
    int j = 0;
    while(j < graphs[i].size()) {
      if(j == 0) {
	logFile<<"node: "<<graphs[i][j]<<" ";
	j++;
	continue;
      }
      else if(j == 1) {
	logFile<<"graph type: "<<graphs[i][j]<<" ";
	j++;
	continue;
      }
      else if(j > 1) {
	logFile<<"dof "<<graphs[i][j]<<" ";
	j++;
	continue;
      }
    }
    logFile<<endl;
  }


  // Write output in the log file.
  if (resultantInternalTraction == true){
    logFile<<"Resultant-Internal-Traction"<<endl;
    for(int s=0;s<resultantForceOnSurfaces.size();s++) {
      logFile<<"DOF: "<<resultantForceOnSurfaces[s][0]<<endl;
      logFile<<"nodes:"<<endl;
      for(int j=1;j<resultantForceOnSurfaces[s].size();j++)
        logFile<<resultantForceOnSurfaces[s][j]<<endl;
    }
  }
  if (resultantInternalSurfaceTorque == true){
    logFile<<"Resultant-Internal-Surface-Torque"<<endl;
    for(int s=0;s<resultantTorqueOnSurfaces.size();s++) {
      logFile<<"torque-axis: "<<resultantTorqueOnSurfaces[s][0]<<" "
	     <<resultantTorqueOnSurfaces[s][1]<<" "
	     <<resultantTorqueOnSurfaces[s][2]<<endl;
      logFile<<"nodes:"<<endl;
      for(int j=3;j<resultantTorqueOnSurfaces[s].size();j++)
        logFile<<resultantTorqueOnSurfaces[s][j]<<endl;
    }
  }


  /**********************************************************************/
  // Read the boundary conditions.
  double dof;
  dbVector normal(3);

  while(name != "END") {

    if(name == "LOAD_CONDITIONS_INPUT")
      break;

    else if(name == "BOUNDARY_CONDITIONS_INPUT" || name == "") {
      inputFile >> name;
      continue;
    }

    else {

      // ================================================================
      // deformation boundary conditions

      if(name == "Displacement-Point-Constraints") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointDispBoundConds.size();
	    pointDispBoundConds.resize(pointDispBoundConds.size()+1);
	    pointDispBoundConds[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointDispNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointDispBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "disp") {
	    inputFile >> dof >> name >> value;

	    pointDispBoundConds[currentCond].push_back(dof);
	    pointDispBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Displacement-Point-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Displacement-Line-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineDispBoundConds.size();
	    lineDispBoundConds.resize(lineDispBoundConds.size()+1);
	    lineDispBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineDispBoundConds[currentCond].push_back(node[j]);

	      break;

	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineDispBoundConds[currentCond].push_back(node[j]);

	      break;

	    default:
	      logFile<<"Variable 'nodesPerboundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "normal") {
	    lineDispNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineDispBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "disp") {
	    inputFile >> dof >> name >> value;

	    lineDispBoundConds[currentCond].push_back(dof);
	    lineDispBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Displacement-Line-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Displacement-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceDispBoundConds.size();
	    surfaceDispBoundConds.resize(surfaceDispBoundConds.size()+1);
	    surfaceDispBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceDispBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceDispBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceDispBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerBoundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "disp") {
	    inputFile >> dof >> name >> value;

	    surfaceDispBoundConds[currentCond].push_back(dof);
	    surfaceDispBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Displacement-Surface-Constraints' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

      // ----------------------------------------------------------------
      // rotation boundary conditions

      else if(name == "Rotation-Point-Constraints") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointRotBoundConds.size();
	    pointRotBoundConds.resize(pointRotBoundConds.size()+1);
	    pointRotBoundConds[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointRotNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointRotBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "rot") {
	    inputFile >> dof >> name >> value;

	    pointRotBoundConds[currentCond].push_back(dof);
	    pointRotBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Rotation-Point-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Rotation-Line-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineRotBoundConds.size();
	    lineRotBoundConds.resize(lineRotBoundConds.size()+1);
	    lineRotBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineRotBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineRotBoundConds[currentCond].push_back(node[j]);

	      break;

	    default:
	      logFile<<"Variable 'nodesPerboundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "normal") {
	    lineRotNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineRotBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "rot") {
	    inputFile >> dof >> name >> value;

	    lineRotBoundConds[currentCond].push_back(dof);
	    lineRotBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Rotation-Line-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Rotation-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceRotBoundConds.size();
	    surfaceRotBoundConds.resize(surfaceRotBoundConds.size()+1);
	    surfaceRotBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceRotBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceRotBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceRotBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerBoundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "rot") {
	    inputFile >> dof >> name >> value;

	    surfaceRotBoundConds[currentCond].push_back(dof);
	    surfaceRotBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Rotation-Surface-Constraints' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

      // ================================================================
      // electric boundary conditions

      else if(name == "Electric-Point-Constraints") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointElectricBoundConds.size();
	    pointElectricBoundConds.resize(pointElectricBoundConds.size()+1);
	    pointElectricBoundConds[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointElectricNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointElectricBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    pointElectricBoundConds[currentCond].push_back(dof);
	    pointElectricBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Electric-Point-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Electric-Line-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineElectricBoundConds.size();
	    lineElectricBoundConds.resize(lineElectricBoundConds.size()+1);
	    lineElectricBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineElectricBoundConds[currentCond].push_back(node[j]);

	      break;

	    default:
	      logFile<<"Variable 'nodesPerboundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "normal") {
	    lineElectricNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineElectricBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    lineElectricBoundConds[currentCond].push_back(dof);
	    lineElectricBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Electric-Line-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Electric-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceElectricBoundConds.size();
	    surfaceElectricBoundConds.resize(surfaceElectricBoundConds.size()+1);
	    surfaceElectricBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerBoundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    surfaceElectricBoundConds[currentCond].push_back(dof);
	    surfaceElectricBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Electric-Surface-Constraints' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Linear-Electric-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = linearSurfaceElectricBoundConds.size();
	    linearSurfaceElectricBoundConds.resize(linearSurfaceElectricBoundConds.size()+1);
	    linearSurfaceElectricBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		linearSurfaceElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		linearSurfaceElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		linearSurfaceElectricBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerBoundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    linearSurfaceElectricBoundConds[currentCond].push_back(dof);
	    linearSurfaceElectricBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "gradient") {
	    inputFile >> value1 >> value2 >> value3;

	    linearSurfaceElectricBoundConds[currentCond].push_back(value1);
	    linearSurfaceElectricBoundConds[currentCond].push_back(value2);
	    linearSurfaceElectricBoundConds[currentCond].push_back(value3);
	  }
	  else if(name == "x0") {
	    inputFile >> value1 >> value2 >> value3;

	    linearSurfaceElectricBoundConds[currentCond].push_back(value1);
	    linearSurfaceElectricBoundConds[currentCond].push_back(value2);
	    linearSurfaceElectricBoundConds[currentCond].push_back(value3);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Linear-Electric-Surface-Constraints' "
		   <<" input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

      // ================================================================
      // depolarisation boundary conditions

      else if(name == "Depolarisation-Time-Point-Constraints") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointDepolarisationBoundConds.size();
	    pointDepolarisationBoundConds.resize(pointDepolarisationBoundConds.size()+1);
	    pointDepolarisationBoundConds[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointDepolarisationNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointDepolarisationBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    pointDepolarisationBoundConds[currentCond].push_back(dof);
	    pointDepolarisationBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Electric-Point-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

      else if(name == "Depolarisation-Time-Line-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineDepolarisationBoundConds.size();
	    lineDepolarisationBoundConds.resize(lineDepolarisationBoundConds.size()+1);
	    lineDepolarisationBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineDepolarisationBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineDepolarisationBoundConds[currentCond].push_back(node[j]);

	      break;

	    default:
	      logFile<<"Variable 'nodesPerboundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "normal") {
	    lineDepolarisationNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineDepolarisationBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    lineDepolarisationBoundConds[currentCond].push_back(dof);
	    lineDepolarisationBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Depolarisation-Time-Line-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      else if(name == "Depolarisation-Time-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceDepolarisationBoundConds.size();
	    surfaceDepolarisationBoundConds.resize(surfaceDepolarisationBoundConds.size()+1);
	    surfaceDepolarisationBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceDepolarisationBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceDepolarisationBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceDepolarisationBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerBoundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    surfaceDepolarisationBoundConds[currentCond].push_back(dof);
	    surfaceDepolarisationBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Depolarisation-Time-Surface-Constraints' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}



      }
      // ================================================================
      // stress boundary conditions

      else if(name == "Stress-Point-Constraints") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointStressBoundConds.size();
	    pointStressBoundConds.resize(pointStressBoundConds.size()+1);
	    pointStressBoundConds[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointStressNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointStressBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "stress") {
	    inputFile >> dof >> name >> value;

	    pointStressBoundConds[currentCond].push_back(dof);
	    pointStressBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Stress-Point-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      else if(name == "Stress-Line-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineStressBoundConds.size();
	    lineStressBoundConds.resize(lineStressBoundConds.size()+1);
	    lineStressBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineStressBoundConds[currentCond].push_back(node[j]);

	      break;

	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineStressBoundConds[currentCond].push_back(node[j]);

	      break;

	    default:
	      logFile<<"Variable 'nodesPerLineElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "normal") {
	    lineStressNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineStressBoundConds[currentCond].push_back(normal[k]);
	  }
	  else if(name == "stress") {
	    inputFile >> dof >> name >> value;

	    lineStressBoundConds[currentCond].push_back(dof);
	    lineStressBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Stress-Line-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      else if(name == "Stress-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceStressBoundConds.size();
	    surfaceStressBoundConds.resize(surfaceStressBoundConds.size()+1);
	    surfaceStressBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceStressBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceStressBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceStressBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerboundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "stress") {
	    inputFile >> dof >> name >> value;

	    surfaceStressBoundConds[currentCond].push_back(dof);
	    surfaceStressBoundConds[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Electric-Surface-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Micro-Surface-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Stress-Surface-Constraints' "
		   <<"input data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }

      else if(name == "Micro-Surface-Constraints") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceMicroBoundConds.size();
	    surfaceMicroBoundConds.resize(surfaceMicroBoundConds.size()+1);
	    surfaceMicroBoundConds[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceMicroBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceMicroBoundConds[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceMicroBoundConds[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerBoundElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	    continue;
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    surfaceMicroBoundConds[currentCond].push_back(dof);
	    surfaceMicroBoundConds[currentCond].push_back(value);

	  }
	  else if(name == "")
	    continue;
	  else if(name == "Displacement-Point-Constraints" ||
		  name == "Displacement-Line-Constraints" ||
		  name == "Displacement-Surface-Constraints" ||
		  name == "Rotation-Point-Constraints" ||
		  name == "Rotation-Line-Constraints" ||
		  name == "Rotation-Surface-Constraints" ||
		  name == "Electric-Point-Constraints" ||
		  name == "Electric-Line-Constraints" ||
		  name == "Linear-Electric-Surface-Constraints" ||
		  name == "Micro-Point-Constraints" ||
		  name == "Micro-Line-Constraints" ||
		  name == "Stress-Point-Constraints" ||
		  name == "Stress-Line-Constraints" ||
		  name == "Stress-Surface-Constraints" ||
		  name == "Depolarisation-Time-Surface-Constraints" ||
		  name == "Depolarisation-Time-Point-Constraints" ||
		  name == "Depolarisation-Time-Line-Constraints" ||
		  name == "LOAD_CONDITIONS_INPUT")
	    break;
	  else {
	    logFile<<"Incorrect written 'Micro-Surface-Constraints' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

    }

  }

  // Write output in the log file.
  logFile<<"BOUNDARY_CONDITIONS"<<endl;
  logFile<<"Displacement-Point-Constraints"<<endl;
  for(int i=0;i<pointDispBoundConds.size();i++) {
    int j = 0;
    while(j < pointDispBoundConds[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointDispBoundConds[i][0]<<endl;
	if(pointDispNormalSet) {
	  logFile<<"normal: "<<pointDispBoundConds[i][1]<<" "
		 <<pointDispBoundConds[i][2]<<" "<<pointDispBoundConds[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      else {
	logFile<<"dof "<<pointDispBoundConds[i][j]
	       <<" amount "<<pointDispBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Displacement-Line-Constraints"<<endl;
  for(int i=0;i<lineDispBoundConds.size();i++) {
    int j = 0;
    while(j < lineDispBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineDispBoundConds[i][0]<<" nodes: "
		 <<lineDispBoundConds[i][1]<<" "<<lineDispBoundConds[i][2]<<endl;
	  if(lineDispNormalSet) {
	    logFile<<"normal: "<<lineDispBoundConds[i][3]<<" "
		   <<lineDispBoundConds[i][4]<<" "<<lineDispBoundConds[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineDispBoundConds[i][0]<<" nodes: "
		 <<lineDispBoundConds[i][1]<<" "<<lineDispBoundConds[i][2]<<" "
		 <<lineDispBoundConds[i][3]<<endl;
	  if(lineDispNormalSet) {
	    logFile<<"normal: "<<lineDispBoundConds[i][4]<<" "
		   <<lineDispBoundConds[i][5]<<" "<<lineDispBoundConds[i][6]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<lineDispBoundConds[i][j]
	       <<" amount "<<lineDispBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Displacement-Surface-Constraints"<<endl;
  for(int i=0;i<surfaceDispBoundConds.size();i++) {
    int j = 0;
    while(j < surfaceDispBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceDispBoundConds[i][0]<<" nodes: "
		 <<surfaceDispBoundConds[i][1]<<" "
		 <<surfaceDispBoundConds[i][2]<<" "
		 <<surfaceDispBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceDispBoundConds[i][0]<<" nodes: "
		 <<surfaceDispBoundConds[i][1]<<" "
		 <<surfaceDispBoundConds[i][2]<<" "
		 <<surfaceDispBoundConds[i][3]<<" "
		 <<surfaceDispBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceDispBoundConds[i][0]<<" nodes: "
		 <<surfaceDispBoundConds[i][1]<<" "
		 <<surfaceDispBoundConds[i][2]<<" "
		 <<surfaceDispBoundConds[i][3]<<" "
		 <<surfaceDispBoundConds[i][4]<<" "
		 <<surfaceDispBoundConds[i][5]<<" "
		 <<surfaceDispBoundConds[i][6]<<" "
		 <<surfaceDispBoundConds[i][7]<<" "
		 <<surfaceDispBoundConds[i][8]<<" "
		 <<surfaceDispBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<surfaceDispBoundConds[i][j]
	       <<" amount "<<surfaceDispBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Rotation-Point-Constraints"<<endl;
  for(int i=0;i<pointRotBoundConds.size();i++) {
    int j = 0;
    while(j < pointRotBoundConds[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointRotBoundConds[i][0]<<endl;
	if(pointRotNormalSet) {
	  logFile<<"normal: "<<pointRotBoundConds[i][1]<<" "
		 <<pointRotBoundConds[i][2]<<" "<<pointRotBoundConds[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      else {
	logFile<<"dof "<<pointRotBoundConds[i][j]
	       <<" amount "<<pointRotBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Rotation-Line-Constraints"<<endl;
  for(int i=0;i<lineRotBoundConds.size();i++) {
    int j = 0;
    while(j < lineRotBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineRotBoundConds[i][0]<<" nodes: "
		 <<lineRotBoundConds[i][1]<<" "<<lineRotBoundConds[i][2]<<endl;
	  if(lineRotNormalSet) {
	    logFile<<"normal: "<<lineRotBoundConds[i][3]<<" "
		   <<lineRotBoundConds[i][4]<<" "<<lineRotBoundConds[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineRotBoundConds[i][0]<<" nodes: "
		 <<lineRotBoundConds[i][1]<<" "<<lineRotBoundConds[i][2]
		 <<" "<<lineRotBoundConds[i][3]<<endl;
	  if(lineRotNormalSet) {
	    logFile<<"normal: "<<lineRotBoundConds[i][4]<<" "
		   <<lineRotBoundConds[i][5]<<" "<<lineRotBoundConds[i][5]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<lineRotBoundConds[i][j]
	       <<" amount "<<lineRotBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Rotation-Surface-Constraints"<<endl;
  for(int i=0;i<surfaceRotBoundConds.size();i++) {
    int j = 0;
    while(j < surfaceRotBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceRotBoundConds[i][0]<<" nodes: "
		 <<surfaceRotBoundConds[i][1]<<" "
		 <<surfaceRotBoundConds[i][2]<<" "
		 <<surfaceRotBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceRotBoundConds[i][0]<<" nodes: "
		 <<surfaceRotBoundConds[i][1]<<" "
		 <<surfaceRotBoundConds[i][2]<<" "
		 <<surfaceRotBoundConds[i][3]<<" "
		 <<surfaceRotBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceRotBoundConds[i][0]<<" nodes: "
		 <<surfaceRotBoundConds[i][1]<<" "
		 <<surfaceRotBoundConds[i][2]<<" "
		 <<surfaceRotBoundConds[i][3]<<" "
		 <<surfaceRotBoundConds[i][4]<<" "
		 <<surfaceRotBoundConds[i][5]<<" "
		 <<surfaceRotBoundConds[i][6]<<" "
		 <<surfaceRotBoundConds[i][7]<<" "
		 <<surfaceRotBoundConds[i][8]<<" "
		 <<surfaceRotBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<surfaceRotBoundConds[i][j]
	       <<" amount "<<surfaceRotBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  // --------------------------------------------------------------------
  logFile<<"Electric-Point-Constraints"<<endl;
  for(int i=0;i<pointElectricBoundConds.size();i++) {
    int j = 0;
    while(j < pointElectricBoundConds[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointElectricBoundConds[i][0]<<endl;
	if(pointElectricNormalSet) {
	  logFile<<"normal: "<<pointElectricBoundConds[i][1]<<" "
		 <<pointElectricBoundConds[i][2]<<" "<<pointElectricBoundConds[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      else {
	logFile<<"dof "<<pointElectricBoundConds[i][j]
	       <<" amount "<<pointElectricBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Electric-Line-Constraints"<<endl;
  for(int i=0;i<lineElectricBoundConds.size();i++) {
    int j = 0;
    while(j < lineElectricBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineElectricBoundConds[i][0]<<" nodes: "
		 <<lineElectricBoundConds[i][1]<<" "<<lineElectricBoundConds[i][2]<<endl;
	  if(lineElectricNormalSet) {
	    logFile<<"normal: "<<lineElectricBoundConds[i][3]<<" "
		   <<lineElectricBoundConds[i][4]<<" "<<lineElectricBoundConds[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineElectricBoundConds[i][0]<<" nodes: "
		 <<lineElectricBoundConds[i][1]<<" "<<lineElectricBoundConds[i][2]
		 <<" "<<lineElectricBoundConds[i][3]<<endl;
	  if(lineElectricNormalSet) {
	    logFile<<"normal: "<<lineElectricBoundConds[i][4]<<" "
		   <<lineElectricBoundConds[i][5]<<" "<<lineElectricBoundConds[i][5]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<lineElectricBoundConds[i][j]
	       <<" amount "<<lineElectricBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Electric-Surface-Constraints"<<endl;
  for(int i=0;i<surfaceElectricBoundConds.size();i++) {
    int j = 0;
    while(j < surfaceElectricBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceElectricBoundConds[i][0]<<" nodes: "
		 <<surfaceElectricBoundConds[i][1]<<" "
		 <<surfaceElectricBoundConds[i][2]<<" "
		 <<surfaceElectricBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceElectricBoundConds[i][0]<<" nodes: "
		 <<surfaceElectricBoundConds[i][1]<<" "
		 <<surfaceElectricBoundConds[i][2]<<" "
		 <<surfaceElectricBoundConds[i][3]<<" "
		 <<surfaceElectricBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceElectricBoundConds[i][0]<<" nodes: "
		 <<surfaceElectricBoundConds[i][1]<<" "
		 <<surfaceElectricBoundConds[i][2]<<" "
		 <<surfaceElectricBoundConds[i][3]<<" "
		 <<surfaceElectricBoundConds[i][4]<<" "
		 <<surfaceElectricBoundConds[i][5]<<" "
		 <<surfaceElectricBoundConds[i][6]<<" "
		 <<surfaceElectricBoundConds[i][7]<<" "
		 <<surfaceElectricBoundConds[i][8]<<" "
		 <<surfaceElectricBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<surfaceElectricBoundConds[i][j]
	       <<" amount "<<surfaceElectricBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Linear-Electric-Surface-Constraints"<<endl;
  for(int i=0;i<linearSurfaceElectricBoundConds.size();i++) {
    int j = 0;
    while(j < linearSurfaceElectricBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<linearSurfaceElectricBoundConds[i][0]<<" nodes: "
		 <<linearSurfaceElectricBoundConds[i][1]<<" "
		 <<linearSurfaceElectricBoundConds[i][2]<<" "
		 <<linearSurfaceElectricBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<linearSurfaceElectricBoundConds[i][0]<<" nodes: "
		 <<linearSurfaceElectricBoundConds[i][1]<<" "
		 <<linearSurfaceElectricBoundConds[i][2]<<" "
		 <<linearSurfaceElectricBoundConds[i][3]<<" "
		 <<linearSurfaceElectricBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<linearSurfaceElectricBoundConds[i][0]<<" nodes: "
		 <<linearSurfaceElectricBoundConds[i][1]<<" "
		 <<linearSurfaceElectricBoundConds[i][2]<<" "
		 <<linearSurfaceElectricBoundConds[i][3]<<" "
		 <<linearSurfaceElectricBoundConds[i][4]<<" "
		 <<linearSurfaceElectricBoundConds[i][5]<<" "
		 <<linearSurfaceElectricBoundConds[i][6]<<" "
		 <<linearSurfaceElectricBoundConds[i][7]<<" "
		 <<linearSurfaceElectricBoundConds[i][8]<<" "
		 <<linearSurfaceElectricBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof"<<linearSurfaceElectricBoundConds[i][j]<<" value "
	       <<linearSurfaceElectricBoundConds[i][j+1]
	       <<" grad: "<<linearSurfaceElectricBoundConds[i][j+2]<<" "
	       <<linearSurfaceElectricBoundConds[i][j+3]<<" "
	       <<linearSurfaceElectricBoundConds[i][j+4]
	       <<" x0:"<<linearSurfaceElectricBoundConds[i][j+5]<<" "
	       <<linearSurfaceElectricBoundConds[i][j+6]<<" "
	       <<linearSurfaceElectricBoundConds[i][j+7]<<endl;
	j+=8;
      }
    }
  }

  // --------------------------------------------------------------------
  logFile<<"Depolarisation-Time-Point-Constraints"<<endl;
  for(int i=0;i<pointDepolarisationBoundConds.size();i++) {
    int j = 0;
    while(j < pointDepolarisationBoundConds[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointDepolarisationBoundConds[i][0]<<endl;
	if(pointDepolarisationNormalSet) {
	  logFile<<"normal: "<<pointDepolarisationBoundConds[i][1]<<" "
		 <<pointDepolarisationBoundConds[i][2]<<" "<<pointDepolarisationBoundConds[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      else {
	logFile<<"dof "<<pointDepolarisationBoundConds[i][j]
	       <<" amount "<<pointDepolarisationBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Depolarisation-Time-Line-Constraints"<<endl;
  for(int i=0;i<lineDepolarisationBoundConds.size();i++) {
    int j = 0;
    while(j < lineDepolarisationBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineDepolarisationBoundConds[i][0]<<" nodes: "
		 <<lineDepolarisationBoundConds[i][1]<<" "<<lineDepolarisationBoundConds[i][2]<<endl;
	  if(lineDepolarisationNormalSet) {
	    logFile<<"normal: "<<lineDepolarisationBoundConds[i][3]<<" "
		   <<lineDepolarisationBoundConds[i][4]<<" "<<lineDepolarisationBoundConds[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineDepolarisationBoundConds[i][0]<<" nodes: "
		 <<lineDepolarisationBoundConds[i][1]<<" "<<lineDepolarisationBoundConds[i][2]
		 <<" "<<lineDepolarisationBoundConds[i][3]<<endl;
	  if(lineDepolarisationNormalSet) {
	    logFile<<"normal: "<<lineDepolarisationBoundConds[i][4]<<" "
		   <<lineDepolarisationBoundConds[i][5]<<" "<<lineDepolarisationBoundConds[i][5]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<lineDepolarisationBoundConds[i][j]
	       <<" amount "<<lineDepolarisationBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Depolarisation-Time-Surface-Constraints"<<endl;
  for(int i=0;i<surfaceDepolarisationBoundConds.size();i++) {
    int j = 0;
    while(j < surfaceDepolarisationBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceDepolarisationBoundConds[i][0]<<" nodes: "
		 <<surfaceDepolarisationBoundConds[i][1]<<" "
		 <<surfaceDepolarisationBoundConds[i][2]<<" "
		 <<surfaceDepolarisationBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceDepolarisationBoundConds[i][0]<<" nodes: "
		 <<surfaceDepolarisationBoundConds[i][1]<<" "
		 <<surfaceDepolarisationBoundConds[i][2]<<" "
		 <<surfaceDepolarisationBoundConds[i][3]<<" "
		 <<surfaceDepolarisationBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceDepolarisationBoundConds[i][0]<<" nodes: "
		 <<surfaceDepolarisationBoundConds[i][1]<<" "
		 <<surfaceDepolarisationBoundConds[i][2]<<" "
		 <<surfaceDepolarisationBoundConds[i][3]<<" "
		 <<surfaceDepolarisationBoundConds[i][4]<<" "
		 <<surfaceDepolarisationBoundConds[i][5]<<" "
		 <<surfaceDepolarisationBoundConds[i][6]<<" "
		 <<surfaceDepolarisationBoundConds[i][7]<<" "
		 <<surfaceDepolarisationBoundConds[i][8]<<" "
		 <<surfaceDepolarisationBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<surfaceDepolarisationBoundConds[i][j]
	       <<" amount "<<surfaceDepolarisationBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  // --------------------------------------------------------------------
  logFile<<"Stress-Point-Constraints"<<endl;
  for(int i=0;i<pointStressBoundConds.size();i++) {
    int j = 0;
    while(j < pointStressBoundConds[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointStressBoundConds[i][0]<<endl;
	if(pointStressNormalSet) {
	  logFile<<"normal: "<<pointStressBoundConds[i][1]<<" "
		 <<pointStressBoundConds[i][2]<<" "<<pointStressBoundConds[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      else {
	logFile<<"dof "<<pointStressBoundConds[i][j]
	       <<" amount "<<pointStressBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Stress-Line-Constraints"<<endl;
  for(int i=0;i<lineStressBoundConds.size();i++) {
    int j = 0;
    while(j < lineStressBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineStressBoundConds[i][0]<<" nodes: "
		 <<lineStressBoundConds[i][1]<<" "<<lineStressBoundConds[i][2]<<endl;
	  if(lineStressNormalSet) {
	    logFile<<"normal: "<<lineStressBoundConds[i][3]<<" "
		   <<lineStressBoundConds[i][4]<<" "<<lineStressBoundConds[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineStressBoundConds[i][0]<<" nodes: "
		 <<lineStressBoundConds[i][1]<<" "<<lineStressBoundConds[i][2]
		 <<" "<<lineStressBoundConds[i][3]<<endl;
	  if(lineStressNormalSet) {
	    logFile<<"normal: "<<lineStressBoundConds[i][4]<<" "
		   <<lineStressBoundConds[i][5]<<" "<<lineStressBoundConds[i][5]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<lineStressBoundConds[i][j]
	       <<" amount "<<lineStressBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Stress-Surface-Constraints"<<endl;
  for(int i=0;i<surfaceStressBoundConds.size();i++) {
    int j = 0;
    while(j < surfaceStressBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceStressBoundConds[i][0]<<" nodes: "
		 <<surfaceStressBoundConds[i][1]<<" "
		 <<surfaceStressBoundConds[i][2]<<" "
		 <<surfaceStressBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceStressBoundConds[i][0]<<" nodes: "
		 <<surfaceStressBoundConds[i][1]<<" "
		 <<surfaceStressBoundConds[i][2]<<" "
		 <<surfaceStressBoundConds[i][3]<<" "
		 <<surfaceStressBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceStressBoundConds[i][0]<<" nodes: "
		 <<surfaceStressBoundConds[i][1]<<" "
		 <<surfaceStressBoundConds[i][2]<<" "
		 <<surfaceStressBoundConds[i][3]<<" "
		 <<surfaceStressBoundConds[i][4]<<" "
		 <<surfaceStressBoundConds[i][5]<<" "
		 <<surfaceStressBoundConds[i][6]<<" "
		 <<surfaceStressBoundConds[i][7]<<" "
		 <<surfaceStressBoundConds[i][8]<<" "
		 <<surfaceStressBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<surfaceStressBoundConds[i][j]
	       <<" amount "<<surfaceStressBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }
  logFile<<"Micro-Surface-Constraints"<<endl;
  for(int i=0;i<surfaceMicroBoundConds.size();i++) {
    int j = 0;
    while(j < surfaceMicroBoundConds[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceMicroBoundConds[i][0]<<" nodes: "
		 <<surfaceMicroBoundConds[i][1]<<" "
		 <<surfaceMicroBoundConds[i][2]<<" "
		 <<surfaceMicroBoundConds[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceMicroBoundConds[i][0]<<" nodes: "
		 <<surfaceMicroBoundConds[i][1]<<" "
		 <<surfaceMicroBoundConds[i][2]<<" "
		 <<surfaceMicroBoundConds[i][3]<<" "
		 <<surfaceMicroBoundConds[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceMicroBoundConds[i][0]<<" nodes: "
		 <<surfaceMicroBoundConds[i][1]<<" "
		 <<surfaceMicroBoundConds[i][2]<<" "
		 <<surfaceMicroBoundConds[i][3]<<" "
		 <<surfaceMicroBoundConds[i][4]<<" "
		 <<surfaceMicroBoundConds[i][5]<<" "
		 <<surfaceMicroBoundConds[i][6]<<" "
		 <<surfaceMicroBoundConds[i][7]<<" "
		 <<surfaceMicroBoundConds[i][8]<<" "
		 <<surfaceMicroBoundConds[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      else {
	logFile<<"dof "<<surfaceMicroBoundConds[i][j]
	       <<" amount "<<surfaceMicroBoundConds[i][j+1]<<endl;
	j+=2;
      }
    }
  }

  /********************************************************************/
  // Read the loading conditions.
  while(name != "END") {

    if(name == "BOUNDARY_ELEMENTS" ||
       name == "END")
      break;

    else if(name == "LOAD_CONDITIONS_INPUT" || name == "") {
      inputFile >> name;
      continue;
    }

    else {

      if(name == "Point-Force-Loads") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointForceLoads.size();
	    pointForceLoads.resize(pointForceLoads.size()+1);
	    pointForceLoads[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointForceNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointForceLoads[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    pointForceLoads[currentCond].push_back(dof);
	    pointForceLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written Point-Force-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      // ----------------------------------------------------------------
      else if(name == "Line-Force-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineForceLoads.size();
	    lineForceLoads.resize(lineForceLoads.size()+1);
	    lineForceLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineForceLoads[currentCond].push_back(node[j]);

	      break;
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineForceLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerLineElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }
	  }
	  else if(name == "normal") {
	    lineForceNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineForceLoads[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    lineForceLoads[currentCond].push_back(dof);
	    lineForceLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written Line-Force-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      // ----------------------------------------------------------------
      else if(name == "Traction-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = tractionLoads.size();
	    tractionLoads.resize(tractionLoads.size()+1);
	    tractionLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		tractionLoads[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		tractionLoads[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		tractionLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerSurfaceElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    tractionLoads[currentCond].push_back(dof);
	    tractionLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Traction-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      // ----------------------------------------------------------------
      else if(name == "Surface-Pressure-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfacePressureLoads.size();
	    surfacePressureLoads.resize(surfacePressureLoads.size()+1);
	    surfacePressureLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfacePressureLoads[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfacePressureLoads[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfacePressureLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerSurfaceElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }
	  }
	  else if(name == "amount") {
	    inputFile >> value;
	    surfacePressureLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Surface-Pressure-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      // ----------------------------------------------------------------
      else if(name == "Body-Force-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = bodyForceLoads.size();
	    bodyForceLoads.resize(bodyForceLoads.size()+1);
	    bodyForceLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerElem) {
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		bodyForceLoads[currentCond].push_back(node[j]);

	      break;
	    case 8:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7];

	      for(int j=0;j<8;j++)
		bodyForceLoads[currentCond].push_back(node[j]);

	      break;
	    case 27:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8] >> node[9] >> node[10] >> node[11]
			>> node[12] >> node[13] >> node[14] >> node[15]
			>> node[16] >> node[17] >> node[18] >> node[19]
			>> node[20] >> node[21] >> node[22] >> node[23]
			>> node[24] >> node[25] >> node[26];

	      for(int j=0;j<27;j++)
		bodyForceLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    bodyForceLoads[currentCond].push_back(dof);
	    bodyForceLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Body-Force-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      // ----------------------------------------------------------------
      else if(name == "Point-Moment-Loads") {

	while(inputFile >> name) {

	  if(name == "node") {
	    inputFile >> particle;
	    currentCond = pointMomentLoads.size();
	    pointMomentLoads.resize(pointMomentLoads.size()+1);
	    pointMomentLoads[currentCond].push_back(particle);
	    continue;
	  }
	  else if(name == "normal") {
	    pointMomentNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      pointMomentLoads[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;

	    pointMomentLoads[currentCond].push_back(dof);
	    pointMomentLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written Point-Moment-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }
      // ----------------------------------------------------------------
      else if(name == "Line-Moment-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = lineMomentLoads.size();
	    lineMomentLoads.resize(lineMomentLoads.size()+1);
	    lineMomentLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerLineElem) {
	    case 2:
	      inputFile >> node[0] >> node[1];

	      for(int j=0;j<2;j++)
		lineMomentLoads[currentCond].push_back(node[j]);

	      break;
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		lineMomentLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerLineElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }
	  }
	  else if(name == "normal") {
	    lineMomentNormalSet = true;
	    inputFile >> normal[0] >> normal[1] >> normal[2];

	    for(int k=0;k<normal.size();k++)
	      lineMomentLoads[currentCond].push_back(normal[k]);
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    lineMomentLoads[currentCond].push_back(dof);
	    lineMomentLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written Line-Moment-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      // ----------------------------------------------------------------
      else if(name == "Surface-Moment-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceMomentLoads.size();
	    surfaceMomentLoads.resize(surfaceMomentLoads.size()+1);
	    surfaceMomentLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceMomentLoads[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceMomentLoads[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceMomentLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerSurfaceElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    surfaceMomentLoads[currentCond].push_back(dof);
	    surfaceMomentLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Surface-Moment-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }
      // ----------------------------------------------------------------
      else if(name == "Body-Moment-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = bodyMomentLoads.size();
	    bodyMomentLoads.resize(bodyMomentLoads.size()+1);
	    bodyMomentLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerElem) {
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		bodyMomentLoads[currentCond].push_back(node[j]);

	      break;
	    case 8:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7];

	      for(int j=0;j<8;j++)
		bodyMomentLoads[currentCond].push_back(node[j]);

	      break;
	    case 27:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8] >> node[9] >> node[10] >> node[11]
			>> node[12] >> node[13] >> node[14] >> node[15]
			>> node[16] >> node[17] >> node[18] >> node[19]
			>> node[20] >> node[21] >> node[22] >> node[23]
			>> node[24] >> node[25] >> node[26];

	      for(int j=0;j<27;j++)
		bodyMomentLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    bodyMomentLoads[currentCond].push_back(dof);
	    bodyMomentLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Body-Moment-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

      // ----------------------------------------------------------------
      else if(name == "Electric-Surface-Charge-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = surfaceElectricChargeLoads.size();
	    surfaceElectricChargeLoads.resize(surfaceElectricChargeLoads.size()+1);
	    surfaceElectricChargeLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerSurfaceElem) {
	    case 3:
	      inputFile >> node[0] >> node[1] >> node[2];

	      for(int j=0;j<3;j++)
		surfaceElectricChargeLoads[currentCond].push_back(node[j]);

	      break;
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		surfaceElectricChargeLoads[currentCond].push_back(node[j]);

	      break;
	    case 9:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8];

	      for(int j=0;j<9;j++)
		surfaceElectricChargeLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerSurfaceElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }
	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    surfaceElectricChargeLoads[currentCond].push_back(dof);
	    surfaceElectricChargeLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Surface-Charge-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Electric-Surface-Charge-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }

	}

      }

      // ----------------------------------------------------------------
      else if(name == "Electric-Body-Charge-Loads") {

	while(inputFile >> name) {

	  if(name == "Element") {
	    inputFile >> elem;
	    currentCond = bodyElectricChargeLoads.size();
	    bodyElectricChargeLoads.resize(bodyElectricChargeLoads.size()+1);
	    bodyElectricChargeLoads[currentCond].push_back(elem);
	    continue;
	  }
	  else if(name == "nodes") {

	    switch(nodesPerElem) {
	    case 4:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3];

	      for(int j=0;j<4;j++)
		bodyElectricChargeLoads[currentCond].push_back(node[j]);

	      break;
	    case 8:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7];

	      for(int j=0;j<8;j++)
		bodyElectricChargeLoads[currentCond].push_back(node[j]);

	      break;
	    case 27:
	      inputFile >> node[0] >> node[1] >> node[2] >> node[3]
			>> node[4] >> node[5] >> node[6] >> node[7]
			>> node[8] >> node[9] >> node[10] >> node[11]
			>> node[12] >> node[13] >> node[14] >> node[15]
			>> node[16] >> node[17] >> node[18] >> node[19]
			>> node[20] >> node[21] >> node[22] >> node[23]
			>> node[24] >> node[25] >> node[26];

	      for(int j=0;j<27;j++)
		bodyElectricChargeLoads[currentCond].push_back(node[j]);

	      break;
	    default:
	      logFile<<"Variable 'nodesPerElem' doesn't exist in class "
		     <<"InputFileData!"<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      break;
	    }

	  }
	  else if(name == "dof") {
	    inputFile >> dof >> name >> value;
	    bodyElectricChargeLoads[currentCond].push_back(dof);
	    bodyElectricChargeLoads[currentCond].push_back(value);
	  }
	  else if(name == "")
	    continue;
	  else if(name == "Point-Force-Loads" ||
		  name == "Line-Force-Loads" ||
		  name == "Body-Force-Loads" ||
		  name == "Traction-Loads"  ||
		  name == "Surface-Pressure-Loads" ||
		  name == "Point-Moment-Loads" ||
		  name == "Line-Moment-Loads" ||
		  name == "Surface-Moment-Loads" ||
		  name == "Body-Moment-Loads" ||
		  name == "Electric-Body-Charge-Loads" ||
		  name == "BOUNDARY_ELEMENTS" ||
		  name == "END")
	    break;
	  else {
	    logFile<<"Incorrect written 'Body-Force-Loads' input "
		   <<"data!"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }
	}

      }

      // ----------------------------------------------------------------
      else
	continue;

    }

  }

  // --------------------------------------------------------------------
  // Write output in the log file.
  logFile<<"LOAD_CONDITIONS"<<endl;
  logFile<<"Point-Force-Loads"<<endl;
  for(int i=0;i<pointForceLoads.size();i++) {
    int j = 0;
    while(j < pointForceLoads[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointForceLoads[i][0]<<endl;
	if(pointForceNormalSet) {
	  logFile<<"normal: "<<pointForceLoads[i][1]<<" "
		 <<pointForceLoads[i][2]<<" "<<pointForceLoads[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      logFile<<"dof "<<pointForceLoads[i][j]
	     <<" amount "<<pointForceLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Line-Force-Loads"<<endl;
  for(int i=0;i<lineForceLoads.size();i++) {
    int j = 0;
    while(j < lineForceLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineForceLoads[i][0]<<" nodes: "
		 <<lineForceLoads[i][1]<<" "<<lineForceLoads[i][2]<<endl;
	  if(lineForceNormalSet) {
	    logFile<<"normal: "<<lineForceLoads[i][3]<<" "
		   <<lineForceLoads[i][4]<<" "<<lineForceLoads[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineForceLoads[i][0]<<" nodes: "
		 <<lineForceLoads[i][1]<<" "<<lineForceLoads[i][2]<<" "
		 <<lineForceLoads[i][3]<<endl;
	  if(lineForceNormalSet) {
	    logFile<<"normal: "<<lineForceLoads[i][4]<<" "
		   <<lineForceLoads[i][5]<<" "<<lineForceLoads[i][6]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<lineForceLoads[i][j]
	     <<" amount "<<lineForceLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Traction-Loads"<<endl;
  for(int i=0;i<tractionLoads.size();i++) {
    int j = 0;
    while(j < tractionLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<tractionLoads[i][0]<<" nodes: "
		 <<tractionLoads[i][1]<<" "
		 <<tractionLoads[i][2]<<" "
		 <<tractionLoads[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<tractionLoads[i][0]<<" nodes: "
		 <<tractionLoads[i][1]<<" "
		 <<tractionLoads[i][2]<<" "
		 <<tractionLoads[i][3]<<" "
		 <<tractionLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<tractionLoads[i][0]<<" nodes: "
		 <<tractionLoads[i][1]<<" "
		 <<tractionLoads[i][2]<<" "
		 <<tractionLoads[i][3]<<" "
		 <<tractionLoads[i][4]<<" "
		 <<tractionLoads[i][5]<<" "
		 <<tractionLoads[i][6]<<" "
		 <<tractionLoads[i][7]<<" "
		 <<tractionLoads[i][8]<<" "
		 <<tractionLoads[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<tractionLoads[i][j]
	     <<" amount "<<tractionLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Surface-Pressure-Loads"<<endl;
  for(int i=0;i<surfacePressureLoads.size();i++) {
    int j = 0;
    while(j < surfacePressureLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfacePressureLoads[i][0]<<" nodes: "
		 <<surfacePressureLoads[i][1]<<" "
		 <<surfacePressureLoads[i][2]<<" "
		 <<surfacePressureLoads[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfacePressureLoads[i][0]<<" nodes: "
		 <<surfacePressureLoads[i][1]<<" "
		 <<surfacePressureLoads[i][2]<<" "
		 <<surfacePressureLoads[i][3]<<" "
		 <<surfacePressureLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfacePressureLoads[i][0]<<" nodes: "
		 <<surfacePressureLoads[i][1]<<" "
		 <<surfacePressureLoads[i][2]<<" "
		 <<surfacePressureLoads[i][3]<<" "
		 <<surfacePressureLoads[i][4]<<" "
		 <<surfacePressureLoads[i][5]<<" "
		 <<surfacePressureLoads[i][6]<<" "
		 <<surfacePressureLoads[i][7]<<" "
		 <<surfacePressureLoads[i][8]<<" "
		 <<surfacePressureLoads[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"amount "<<surfacePressureLoads[i][j]<<endl;
      j++;
    }
  }
  logFile<<"Body-Force-Loads"<<endl;

  for(int i=0;i<bodyForceLoads.size();i++) {
    int j = 0;
    while(j < bodyForceLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerElem) {
	case 4:
	  logFile<<"Element: "<<bodyForceLoads[i][0]<<" nodes: "
		 <<bodyForceLoads[i][1]<<" "
		 <<bodyForceLoads[i][2]<<" "
		 <<bodyForceLoads[i][3]<<" "
		 <<bodyForceLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 8:
	  logFile<<"Element: "<<bodyForceLoads[i][0]<<" nodes: "
		 <<bodyForceLoads[i][1]<<" "
		 <<bodyForceLoads[i][2]<<" "
		 <<bodyForceLoads[i][3]<<" "
		 <<bodyForceLoads[i][4]<<" "
		 <<bodyForceLoads[i][5]<<" "
		 <<bodyForceLoads[i][6]<<" "
		 <<bodyForceLoads[i][7]<<" "
		 <<bodyForceLoads[i][8]<<endl;
	  j+=9;
	  break;
	case 27:
	  logFile<<"Element: "<<bodyForceLoads[i][0]<<" nodes: "
		 <<bodyForceLoads[i][1]<<" "
		 <<bodyForceLoads[i][2]<<" "
		 <<bodyForceLoads[i][3]<<" "
		 <<bodyForceLoads[i][4]<<" "
		 <<bodyForceLoads[i][5]<<" "
		 <<bodyForceLoads[i][6]<<" "
		 <<bodyForceLoads[i][7]<<" "
		 <<bodyForceLoads[i][8]<<" "
		 <<bodyForceLoads[i][9]<<" "
		 <<bodyForceLoads[i][10]<<" "
		 <<bodyForceLoads[i][11]<<" "
		 <<bodyForceLoads[i][12]<<" "
		 <<bodyForceLoads[i][13]<<" "
		 <<bodyForceLoads[i][14]<<" "
		 <<bodyForceLoads[i][15]<<" "
		 <<bodyForceLoads[i][16]<<" "
		 <<bodyForceLoads[i][17]<<" "
		 <<bodyForceLoads[i][18]<<" "
		 <<bodyForceLoads[i][19]<<" "
		 <<bodyForceLoads[i][20]<<" "
		 <<bodyForceLoads[i][21]<<" "
		 <<bodyForceLoads[i][22]<<" "
		 <<bodyForceLoads[i][23]<<" "
		 <<bodyForceLoads[i][24]<<" "
		 <<bodyForceLoads[i][25]<<" "
		 <<bodyForceLoads[i][26]<<" "
		 <<bodyForceLoads[i][27]<<endl;
	  j+=28;
	  break;
	default:
	  logFile<<"Variable 'nodesPerElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<bodyForceLoads[i][j]
	     <<" amount "<<bodyForceLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Point-Moment-Loads"<<endl;
  for(int i=0;i<pointMomentLoads.size();i++) {
    int j = 0;
    while(j < pointMomentLoads[i].size()) {
      if(j == 0) {
	logFile<<"Node: "<<pointMomentLoads[i][0]<<endl;
	if(pointMomentNormalSet) {
	  logFile<<"normal: "<<pointMomentLoads[i][1]<<" "
		 <<pointMomentLoads[i][2]<<" "<<pointMomentLoads[i][3]<<" "
		 <<endl;
	  j+=4;
	}
	else
	  j++;
      }
      logFile<<"dof "<<pointMomentLoads[i][j]
	     <<" amount "<<pointMomentLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Line-Moment-Loads"<<endl;
  for(int i=0;i<lineMomentLoads.size();i++) {
    int j = 0;
    while(j < lineMomentLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerLineElem) {
	case 2:
	  logFile<<"Element: "<<lineMomentLoads[i][0]<<" nodes: "
		 <<lineMomentLoads[i][1]<<" "<<lineMomentLoads[i][2]<<endl;
	  if(lineMomentNormalSet) {
	    logFile<<"normal: "<<lineMomentLoads[i][3]<<" "
		   <<lineMomentLoads[i][4]<<" "<<lineMomentLoads[i][5]<<" "
		   <<endl;
	    j+=6;
	  }
	  else
	    j+=3;
	  break;
	case 3:
	  logFile<<"Element: "<<lineMomentLoads[i][0]<<" nodes: "
		 <<lineMomentLoads[i][1]<<" "<<lineMomentLoads[i][2]
		 <<" "<<lineMomentLoads[i][3]<<endl;
	  if(lineMomentNormalSet) {
	    logFile<<"normal: "<<lineMomentLoads[i][4]<<" "
		   <<lineMomentLoads[i][5]<<" "<<lineMomentLoads[i][6]<<" "
		   <<endl;
	    j+=7;
	  }
	  else
	    j+=4;
	  break;
	default:
	  logFile<<"Variable 'nodesPerLineElem' has got a incorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<lineMomentLoads[i][j]
	     <<" amount "<<lineMomentLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Surface-Moment-Loads"<<endl;
  for(int i=0;i<surfaceMomentLoads.size();i++) {
    int j = 0;
    while(j < surfaceMomentLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceMomentLoads[i][0]<<" nodes: "
		 <<surfaceMomentLoads[i][1]<<" "
		 <<surfaceMomentLoads[i][2]<<" "
		 <<surfaceMomentLoads[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceMomentLoads[i][0]<<" nodes: "
		 <<surfaceMomentLoads[i][1]<<" "
		 <<surfaceMomentLoads[i][2]<<" "
		 <<surfaceMomentLoads[i][3]<<" "
		 <<surfaceMomentLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceMomentLoads[i][0]<<" nodes: "
		 <<surfaceMomentLoads[i][1]<<" "
		 <<surfaceMomentLoads[i][2]<<" "
		 <<surfaceMomentLoads[i][3]<<" "
		 <<surfaceMomentLoads[i][4]<<" "
		 <<surfaceMomentLoads[i][5]<<" "
		 <<surfaceMomentLoads[i][6]<<" "
		 <<surfaceMomentLoads[i][7]<<" "
		 <<surfaceMomentLoads[i][8]<<" "
		 <<surfaceMomentLoads[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<surfaceMomentLoads[i][j]
	     <<" amount "<<surfaceMomentLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Body-Moment-Loads"<<endl;

  for(int i=0;i<bodyMomentLoads.size();i++) {
    int j = 0;
    while(j < bodyMomentLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerElem) {
	case 4:
	  logFile<<"Element: "<<bodyMomentLoads[i][0]<<" nodes: "
		 <<bodyMomentLoads[i][1]<<" "
		 <<bodyMomentLoads[i][2]<<" "
		 <<bodyMomentLoads[i][3]<<" "
		 <<bodyMomentLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 8:
	  logFile<<"Element: "<<bodyMomentLoads[i][0]<<" nodes: "
		 <<bodyMomentLoads[i][1]<<" "
		 <<bodyMomentLoads[i][2]<<" "
		 <<bodyMomentLoads[i][3]<<" "
		 <<bodyMomentLoads[i][4]<<" "
		 <<bodyMomentLoads[i][5]<<" "
		 <<bodyMomentLoads[i][6]<<" "
		 <<bodyMomentLoads[i][7]<<" "
		 <<bodyMomentLoads[i][8]<<endl;
	  j+=9;
	  break;
	case 27:
	  logFile<<"Element: "<<bodyMomentLoads[i][0]<<" nodes: "
		 <<bodyMomentLoads[i][1]<<" "
		 <<bodyMomentLoads[i][2]<<" "
		 <<bodyMomentLoads[i][3]<<" "
		 <<bodyMomentLoads[i][4]<<" "
		 <<bodyMomentLoads[i][5]<<" "
		 <<bodyMomentLoads[i][6]<<" "
		 <<bodyMomentLoads[i][7]<<" "
		 <<bodyMomentLoads[i][8]<<" "
		 <<bodyMomentLoads[i][9]<<" "
		 <<bodyMomentLoads[i][10]<<" "
		 <<bodyMomentLoads[i][11]<<" "
		 <<bodyMomentLoads[i][12]<<" "
		 <<bodyMomentLoads[i][13]<<" "
		 <<bodyMomentLoads[i][14]<<" "
		 <<bodyMomentLoads[i][15]<<" "
		 <<bodyMomentLoads[i][16]<<" "
		 <<bodyMomentLoads[i][17]<<" "
		 <<bodyMomentLoads[i][18]<<" "
		 <<bodyMomentLoads[i][19]<<" "
		 <<bodyMomentLoads[i][20]<<" "
		 <<bodyMomentLoads[i][21]<<" "
		 <<bodyMomentLoads[i][22]<<" "
		 <<bodyMomentLoads[i][23]<<" "
		 <<bodyMomentLoads[i][24]<<" "
		 <<bodyMomentLoads[i][25]<<" "
		 <<bodyMomentLoads[i][26]<<" "
		 <<bodyMomentLoads[i][27]<<endl;
	  j+=28;
	  break;
	default:
	  logFile<<"Variable 'nodesPerElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<bodyMomentLoads[i][j]
	     <<" amount "<<bodyMomentLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Electric-Surface-Charge-Loads"<<endl;
  for(int i=0;i<surfaceElectricChargeLoads.size();i++) {
    int j = 0;
    while(j < surfaceElectricChargeLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerSurfaceElem) {
	case 3:
	  logFile<<"Element: "<<surfaceElectricChargeLoads[i][0]<<" nodes: "
		 <<surfaceElectricChargeLoads[i][1]<<" "
		 <<surfaceElectricChargeLoads[i][2]<<" "
		 <<surfaceElectricChargeLoads[i][3]<<endl;
	  j+=4;
	  break;
	case 4:
	  logFile<<"Element: "<<surfaceElectricChargeLoads[i][0]<<" nodes: "
		 <<surfaceElectricChargeLoads[i][1]<<" "
		 <<surfaceElectricChargeLoads[i][2]<<" "
		 <<surfaceElectricChargeLoads[i][3]<<" "
		 <<surfaceElectricChargeLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 9:
	  logFile<<"Element: "<<surfaceElectricChargeLoads[i][0]<<" nodes: "
		 <<surfaceElectricChargeLoads[i][1]<<" "
		 <<surfaceElectricChargeLoads[i][2]<<" "
		 <<surfaceElectricChargeLoads[i][3]<<" "
		 <<surfaceElectricChargeLoads[i][4]<<" "
		 <<surfaceElectricChargeLoads[i][5]<<" "
		 <<surfaceElectricChargeLoads[i][6]<<" "
		 <<surfaceElectricChargeLoads[i][7]<<" "
		 <<surfaceElectricChargeLoads[i][8]<<" "
		 <<surfaceElectricChargeLoads[i][9]<<endl;
	  j+=10;
	  break;
	default:
	  logFile<<"Variable 'nodesPerSurfaceElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<surfaceElectricChargeLoads[i][j]
	     <<" amount "<<surfaceElectricChargeLoads[i][j+1]<<endl;
      j+=2;
    }
  }
  logFile<<"Electric-Body-Charge-Loads"<<endl;

  for(int i=0;i<bodyElectricChargeLoads.size();i++) {
    int j = 0;
    while(j < bodyElectricChargeLoads[i].size()) {
      if(j == 0) {
	switch(nodesPerElem) {
	case 4:
	  logFile<<"Element: "<<bodyElectricChargeLoads[i][0]<<" nodes: "
		 <<bodyElectricChargeLoads[i][1]<<" "
		 <<bodyElectricChargeLoads[i][2]<<" "
		 <<bodyElectricChargeLoads[i][3]<<" "
		 <<bodyElectricChargeLoads[i][4]<<endl;
	  j+=5;
	  break;
	case 8:
	  logFile<<"Element: "<<bodyElectricChargeLoads[i][0]<<" nodes: "
		 <<bodyElectricChargeLoads[i][1]<<" "
		 <<bodyElectricChargeLoads[i][2]<<" "
		 <<bodyElectricChargeLoads[i][3]<<" "
		 <<bodyElectricChargeLoads[i][4]<<" "
		 <<bodyElectricChargeLoads[i][5]<<" "
		 <<bodyElectricChargeLoads[i][6]<<" "
		 <<bodyElectricChargeLoads[i][7]<<" "
		 <<bodyElectricChargeLoads[i][8]<<endl;
	  j+=9;
	  break;
	case 27:
	  logFile<<"Element: "<<bodyElectricChargeLoads[i][0]<<" nodes: "
		 <<bodyElectricChargeLoads[i][1]<<" "
		 <<bodyElectricChargeLoads[i][2]<<" "
		 <<bodyElectricChargeLoads[i][3]<<" "
		 <<bodyElectricChargeLoads[i][4]<<" "
		 <<bodyElectricChargeLoads[i][5]<<" "
		 <<bodyElectricChargeLoads[i][6]<<" "
		 <<bodyElectricChargeLoads[i][7]<<" "
		 <<bodyElectricChargeLoads[i][8]<<" "
		 <<bodyElectricChargeLoads[i][9]<<" "
		 <<bodyElectricChargeLoads[i][10]<<" "
		 <<bodyElectricChargeLoads[i][11]<<" "
		 <<bodyElectricChargeLoads[i][12]<<" "
		 <<bodyElectricChargeLoads[i][13]<<" "
		 <<bodyElectricChargeLoads[i][14]<<" "
		 <<bodyElectricChargeLoads[i][15]<<" "
		 <<bodyElectricChargeLoads[i][16]<<" "
		 <<bodyElectricChargeLoads[i][17]<<" "
		 <<bodyElectricChargeLoads[i][18]<<" "
		 <<bodyElectricChargeLoads[i][19]<<" "
		 <<bodyElectricChargeLoads[i][20]<<" "
		 <<bodyElectricChargeLoads[i][21]<<" "
		 <<bodyElectricChargeLoads[i][22]<<" "
		 <<bodyElectricChargeLoads[i][23]<<" "
		 <<bodyElectricChargeLoads[i][24]<<" "
		 <<bodyElectricChargeLoads[i][25]<<" "
		 <<bodyElectricChargeLoads[i][26]<<" "
		 <<bodyElectricChargeLoads[i][27]<<endl;
	  j+=28;
	  break;
	default:
	  logFile<<"Variable 'nodesPerElem' has got a inccorrect value in "
		 <<"class InputFileData!"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	  break;
	}
      }
      logFile<<"dof "<<bodyElectricChargeLoads[i][j]
	     <<" amount "<<bodyElectricChargeLoads[i][j+1]<<endl;
      j+=2;
    }
  }

  /**********************************************************************/
  // all problem data (including default ones)

  logFile<<"*****************************************************"<<endl;
  logFile<<"******* input parameter (incl. default values) ******"<<endl;
  for(map<string,double>::const_iterator p = problemData.begin();
      p!=problemData.end();++p)
    logFile<<p->first<<" "<<p->second<<endl;

}

/***********************************************************************/
/***********************************************************************/
// Set default values.
void InputFileData::setDefaultValues(std::ofstream& logFile) {

  using namespace std;

  // Orion: General stuff
  problemData["influenceRangeFactor"] = 0.5;
  problemData["resultFileType"] = 2;
  problemData["PODEnergyLevel"] = 95;
  problemData["directInterpolation"] = 0;
  problemData["gridNodesType"] = 1;

  // Orion: MLS stuff
  problemData["interpolantionType"] = 2;
  problemData["MLSCalculationType"] = 2;
  problemData["MLSPolynomialDegree"] = 1;
  problemData["PODCalculationType"] = 1;
  problemData["shiftedBasisAlgorithm"] = 1;
  problemData["minMLSSupport"] = 1;

  problemData["dispMatrixRearrange"] = 1;

  // FEM stuff
  problemData["gaussPointsPerVolumeElement"] = 8;
  problemData["gaussPointsPerSurfaceElement"] = 4;
  problemData["gaussPointsPerLineElement"] = 2;
  problemData["increaseElemOrderOnBoundary"] = 0;
  problemData["lineToSurfaceElementConversion"] = 0;

  // preprocessing
  problemData["numberOfSortingPlanes"] = 0;
  problemData["dimensionlessBasisPolynom"] = 0;
  problemData["shapefunctionType"] = 2;
  problemData["windowfunctionNorming"] = 2;
  problemData["modifiedMLS"] = 0;
  problemData["supportComputationMode"] = 2;

  // maximum entropy stuff
  problemData["maxEntLocalPtcls"] = 8;
  problemData["MaxEntGamma"] = 1.8;
  problemData["target_zero"] = 1e-5;
  problemData["TolLag"] = 1e-6;

  // (1) location optimized (default) (2) computing load balanced
  problemData["parallelPtcleDistributingMode"] = 1;

  // particle influence radius stuff
  problemData["influenceRadiusMultiplier"] = 1;
  problemData["mininumInfluenceRadiusTolerance"] = 1.0e-03;
  problemData["plusMinusDirectionDependentRadius"] = 0;
  problemData["positiveNegativeRadiusRatio"] = 1.0e-03;

  problemData["radiusDistanceDepend"] = 1; // obsolete
  problemData["minDirectionalPtcleSupport"] = 1;
  problemData["minDirectPtcleSuppReduction"] = 1;
  problemData["influenceRadiusDeterminationAngle"] = 0;

  problemData["radiusDeterminationAlgorithm"] = 1;
  problemData["radiusDeterminationAlgorithm"] = 1;
  problemData["radiusElementDepend"] = 0;
  problemData["radiusWeightDepend"] = 0; // obsolete
  problemData["setPtcleRadsElemDependMotherElemsOnly"] = 0;

  problemData["x1InfluenceRadiusAddend"] = 0;
  problemData["x2InfluenceRadiusAddend"] = 0;
  problemData["x3InfluenceRadiusAddend"] = 0;

  problemData["x1InfluenceRadius"] = 0;
  problemData["x2InfluenceRadius"] = 0;
  problemData["x3InfluenceRadius"] = 0;

  problemData["balanceInfluenceRadii"] = 0;
  problemData["minPtcleSupportEnforcing"] = 0;
  problemData["ptcleSupportIncrease"] = 0;

  // --------------------------------------------------------------------
  // essential boundary condition enforcement
  problemData["boundaryEnforcementMethod"] = 2;
  problemData["penaltyParameter"] = 1.0e+15;

  /*********************************************************************/
  // model specific stuff

  problemData["problemType"] = 1;
  problemData["dynamicCalculation"] = 0;

  problemData["constitutiveLaw"] = 1;

  // updating of rotation field
  problemData["integrationPointRotationUpdating"] = 1;
  problemData["simoSpinorRotationUpdating"] = 1.0;
  problemData["spinorRotationFieldUpdating"] = 0.0;
  problemData["rotationUpdatingStepIncrement"] = 1.0e+06;

  problemData["dirichletBoundaryIntegration"] = 1;
  problemData["deformationDependentLoads"] = 0;

  // Nitsche method stuff
  problemData["nitscheParameter"] = problemData["penaltyParameter"];
  problemData["lineWidth"] = 1.0e-04;
  problemData["pointWidth"] = 1.0e-04;
  problemData["nitscheParameter"] = 0;
  problemData["nitscheScalingFactor"] = 1.0e+00;
  problemData["determineNitscheFactor"] = 0;
  problemData["boundaryEnforcementTolerance"] = 1.0e-08;
  problemData["boundaryEnforcementMaxTolerance"] = 1.0e-07;
  problemData["iteratingNitscheFactorDetermination"] = 0;
  problemData["nitschePenaltyOnly"] = 1;
  problemData["nitscheParameterIncreasement"] = 10;
  problemData["nitscheFactorUnsigned"] = 0;
  problemData["maxNitscheParameterDifference"] = 1.0e+06;

  // Generalized Continua stuff
  problemData["generalContinuumQuadratureOrder"] = 2.0;
  problemData["curvatureChangeTensor"] = 1;
  problemData["variableMicrospace"] = 0;
  problemData["microRotationDOF1"] = 0.0;
  problemData["microRotationDOF2"] = 0.0;
  problemData["microRotationDOF3"] = 0.0;

  problemData["micromorphicFullRank"] = 1.0;

  // electro-mechanical coupling
  problemData["withElectricalFreeSpaceTerm"] = 0.0;
  problemData["electroMechanicsDecoupled"] = 0.0;

  /*********************************************************************/
  // discrete equation system
  problemData["sparseKmatStorageScheme"] = 2;


  /*********************************************************************/
  // postprocessing
  problemData["plotMaxIncrement"] = 1000;
  problemData["plotMinIncrement"] = 1;
  problemData["constantPlotIncrement"] = 0;
  problemData["plotOnParticles"] = 1;
  problemData["plotOnGaussPoints"] = 0;
  problemData["plotDeformationFigure"] = 0;
  problemData["plotStrainDistribution"] = 1;
  problemData["plotPlasticStrainDistribution"] = 0;
  problemData["plotStressDistribution"] = 0;
  problemData["stressFieldScaling"] = 1.0;
  problemData["defGraphXScaling"] = 0.0;
  problemData["defGraphYScaling"] = 0.0;
  problemData["stressGraphPlot"] = 0.0;
  problemData["strainGraphPlot"] = 0.0;
  problemData["pressureVolumePlotting"] = 0;
  problemData["volumeTimePlotting"] = 0;
  problemData["pressureTimePlotting"] = 0;
  problemData["computeDeformedVolumes"] = 0;
  problemData["trueStressPlotting"] = 1;
  problemData["plotSurfaceNormalDistribution"] = 1;
  problemData["plotCurrentMesh"] = 0;

  // loading vs. twist graphs
  problemData["twistGraphPlot"] = 0.0;
  // active stress vs. time graphs
  problemData["activeStressGraphPlot"] = 0.0;

  problemData["plotFibreDirectionDistribution"] = 0.0;
  problemData["propagationSpeedPlot"] = 0.0;

  /*********************************************************************/
  // calculation control
  problemData["extendCalculationID"] = 0;
  problemData["extendTimeIncrement"] = 0;
  problemData["calculationInitialConditionsID"] = 0;

  problemData["loadingPeriod"] = 1.0e+16;
  problemData["timeIncrement"] = 1.0e-02;
  problemData["maxTimeIncrement"] = problemData["timeIncrement"]*50.0;
  problemData["minTimeIncrement"] = problemData["timeIncrement"]*0.001;
  problemData["maxCalculationTime"] = 1.0e+16;


  problemData["neumannLoadingFunction"] =  0;
  problemData["dirichletLoadingFunction"] = 0;

  problemData["calculationControlMode"] =  1;
  problemData["KmatUpdatingWhileIterating"] = 1;// works only for a very
                                                // small loading increment
  problemData["deallocateParallelKmat"] = 0;
  problemData["calcStiffnessTangentDeterminant"] = 0;

  problemData["numberOfProcsPerNode"] = 1;
  problemData["PCSubMatSizeLimit"] = 1.0e+06;
  problemData["subBlockKSPDefinition"] = 0;
  problemData["subBlockSolverPreconditioner"] = 5;
  problemData["gmresRestartIterationNum"] = 200;
  problemData["tangentSymmetric"] = 0;
  problemData["tangentZeroPivotShift"] = 0;
  problemData["directSolver"] = 3; // 1: PETSc, 2: MUMPS, 3: SuperLU

  //for cuda
  problemData["cudaAssembly"] = 0;
  problemData["numThreads"] = 1;

  // displacement controlled simulation
  problemData["displacementControlParticle"] = 1;
  problemData["displacementControlDegreeOfFreedom"] = 1;
  problemData["displacementControlMethod"] = 1;
  //problemData["displacementIncrementScaling"] = 1.0;


  problemData["constantNeumannCalculationIncrement"] = 0;
  problemData["constantDirichletCalculationIncrement"] = 0;

  problemData["neumannCalculationIncrement"] = 1.0;
  problemData["dirichletCalculationIncrement"] = 1.0;

  problemData["maxNeumannLoadingFactor"] =  1;
  problemData["maxDirichletLoadingFactor"] = 1;
  problemData["maxInverseLoadingFactor"] = 1;

  problemData["minNeumannLoadingFactor"] =
    -1.0*problemData["maxNeumannLoadingFactor"];

  problemData["stepMultiplier"] = 0.2;
  problemData["startPerturbation"] = 0;
  problemData["determineEQSystemMinMaxValues"] = 0;
  problemData["equationSystemMultiplier"] = 1;

  problemData["convergenceType"] = 0.0;

  // surface deformation in a certain direction is not respected by
  // computing surface normal
  problemData["fixedPressureDOF"] = -1;

  problemData["rotatingNeumannBoundaryConditions"] = 0;
  problemData["rotatingDirichletBoundaryConditions"] = 0;
  problemData["rotatingBoundaryCondsRotVecX"] = 1.0;
  problemData["rotatingBoundaryCondsRotVecY"] = 0;
  problemData["rotatingBoundaryCondsRotVecZ"] = 0;
  problemData["rotatingBoundaryCondsFixPointX"] = 0;
  problemData["rotatingBoundaryCondsFixPointY"] = 0;
  problemData["rotatingBoundaryCondsFixPointZ"] = 0;

  problemData["postProcessingOnly"] = 0;

  /*********************************************************************/
  // debug stuff
  problemData["checkShapefunctions"] = 0;

  /*********************************************************************/
  // cardiac mechanics stuff

  //controls loading in cardiac simulations
  //this overrides the setup loading
  //0: not cardiac loading
  problemData["cardiacLoadingType"] = 0;

  //0 - if time dependent active tension
  //1 - if steady active tension
  problemData["steadyActiveTension"] = 1;

  //adaptive time stepping type
  problemData["adaptiveTimeSteppingType"] = 1;

  // cardiac simulation parameters
  problemData["endSystolicPressure"] = 0.25e3;
  problemData["endSystolicTime"] = 1250e-3;
  problemData["endDiastolicPressure"] = 0.0;
  problemData["endDiastolicTime"] = 0.0;

  problemData["preDiastolicTime"] = 0;
  problemData["ejectionPressure"] = 5.0e3;

  problemData["maxVolumeIterations"] = 4;
  problemData["volumeTolerance"] = 1.0e-08;
  problemData["constantPressureVolumeGradient"] = 0.0;

  // ---------------------------------------------------------------------
  //for ellipsoidal geometries
  problemData["muEpi"] = 0;
  problemData["muEndo"] = 0;
  problemData["prolateC"] = 0;
  problemData["p1"] = 0;
  problemData["p2"] = 0;
  problemData["p3"] = 0;
  problemData["c0Max"] = 0;
  problemData["c0Min"] = 0;
  problemData["k0"] = 0;

  // ---------------------------------------------------------------------
  // active tension stuff

  problemData["cardiacMechanicsProblem"] = 0;

  //for active tension plotting
  problemData["plotActiveTensionDistribution"] = 0;
  problemData["activeTensionFieldScaling"] = 1;

  problemData["plotContractileLengthDistribution"] = 0;
  problemData["plotSarcomereLengthDistribution"] = 0;

  // Guccione active stress model
  problemData["TMax"] = 135.7*1000; // tension Pa
  problemData["Ca0Max"] = 4.35;
  problemData["Ca0"] = 4.35;
  problemData["m"] = 1.0489;
  problemData["b"] = -1.429;
  problemData["B"] = 4.75;
  problemData["l0"] = 1.58;
  problemData["lR"] = 1.85; //average value
  problemData["lRMin"] = 1.78;
  problemData["lRMax"] = 1.91;
  problemData["maxLoadingFactor"] = 1;

  // ---------------------------------------------------------------------
  // eikonal equation

  problemData["plotDepolarisationTimeDistribution"] = 0;
  problemData["depolarisationTimeFieldScaling"] = 1;
  problemData["propagationSpeedPlot"] = 0;

  //kerckhoff constants
  //   problemData["a6"] = 2.0;
  //   problemData["a7"] = 1.5;
  //   problemData["T0"] = 180*1000; //in Pa
  //   problemData["Ea"] = 20;
  //   problemData["v0"] = 7.5;
  //   problemData["l0"] = 1.9;
  //   problemData["tr"] = 0.075;
  //   problemData["td"] = 0.075;
  //   problemData["b"] = 0.15;
  //   problemData["ld"] = -0.4;

  //   problemData["depolTime"] = 0;

  // two model simulation stuff
  problemData["twoModelCalculation"] = 0;

  /**********************************************************************/
  // deformed Lagrangian computation
  problemData["inverseMaterialComputation"] = 0;
  problemData["inverseProblemType"] = 0;


  /**********************************************************************/
  // inverse material characterization

  problemData["inverseProblemType"] = 0;
  problemData["ExperimentalDataType"] =  0;
  problemData["lambda"]= 0;
  problemData["v"]= 0;
  problemData["epsilonOne"]= 0;
  problemData["epsilonTwo"]= 0;
  problemData["epsilonTwo"]= 0;
  problemData["JacobianParameterIncrement"]= 0;
  problemData["searchType"] =0;
  problemData["convergenceTester"] =0;
  problemData["sigma"] =0;
  problemData["alpha"] =0;
  problemData["delta"] =0;

#ifdef _commonDebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"************ Stored input default values ************"<<endl;
  for(map<string,double>::const_iterator p = problemData.begin();
      p!=problemData.end();++p)
    logFile<<p->first<<" "<<p->second<<endl;
#endif

}

/***********************************************************************/
/***********************************************************************/
void InputFileData::setValue(const char* name,double value) {

  using namespace std;

  problemData[name] = value;
}

/***********************************************************************/
/***********************************************************************/
bool InputFileData::setValue(std::string name,double value) {

  using namespace std;

  bool flag = false;

  map<string,double>::iterator p = problemData.find(name);

  if(p != problemData.end()) {

    p->second = value;
    flag = true;

  }

  return flag;
}

/************************************************************************/
/************************************************************************/
double InputFileData::getValue(const char* name) {

  using namespace std;

  map<string,double>::iterator p = problemData.find(name);

  // Key was found.
  if(p != problemData.end())
    return (double)p->second;

  // Key was not found -> looking for it in the default values.
  else {

    p = defaultProblemData.find(name);

    if(p != defaultProblemData.end()) {
      cerr<<"Input parameter '"<< name <<"' hasn't been found - "
	  <<" default value '"<<(double)p->second<<"' is used!"<< endl;
      return (double)p->second;
    }
    else {
      cerr<<"Input parameter '"<< name <<"' hasn't been found!"
	  <<" Check input file!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

}

/************************************************************************/
/************************************************************************/
double InputFileData::getMatValue(int matID,const char* name) {

  using namespace std;

  if(matID >= materials.size()) {
    cerr<<"Material set no. "<<matID+1<<" does not exist!"
	<<" Check input file!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  map<string,double>::iterator p = materials[matID].find(name);

  if(p != materials[matID].end())
    return (double)p->second;
  else {
    cerr<<"Material input parameter '"<< name <<"' haven't been found!"
	<<" Check input file!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/************************************************************************/
/************************************************************************/
bool InputFileData::setMatValue(int matID,std::string name,double value) {

  using namespace std;

  bool flag = false;

  map<string,double>::iterator p = materials[matID].find(name);

  if(p != materials[matID].end()) {

    p->second = value;
    flag = true;

  }

  return flag;
}

/***********************************************************************/
/***********************************************************************/
// assemble arrays needed to be compute material data via inverse
// modelling
void InputFileData::constructMaterialDataVectors(int& matID){

  using namespace std;

  std::string Max = "Max";
  std::string Min = "Min";
  for (int i=0; i<paramNames.size(); i++){
    maxName.push_back(((string)string(paramNames[i]) + Max));
    minName.push_back(((string)string(paramNames[i]) + Min));
    paramMaxNames.push_back(maxName[i].c_str());
    paramMinNames.push_back(minName[i].c_str());
    paramMaxNamesCopy.push_back(maxName[i].c_str());
    paramMinNamesCopy.push_back(minName[i].c_str());
  }

  if (materialDataSets.size() != matID+1)
    materialDataSets.resize(matID+1);

  // resize the materialDataMax vector if necessary

  if (materialDataMax.size() != matID+1)
    materialDataMax.resize(matID+1);

  // resize the materialDataMin vector if necessary

  if (materialDataMin.size() != matID+1)
    materialDataMin.resize(matID+1);

  /********************set materialDataSets vector********************/
  int n=0;
  while(n < paramNames.size()) {
    //remove parameter from optimization list if not to be optimized for
    if ((getValue(paramMinNames[n])) == 1 ||
	(getValue(paramMaxNames[n])) == 1){
      paramMaxNames.erase(paramMaxNames.begin()+n);
      paramMinNames.erase(paramMinNames.begin()+n);
      paramNames.erase(paramNames.begin()+n);
    }
    else{
      n++;
    }
  }

  //set materialDataSets vector
  for (int i=0;i<paramNames.size();i++){
    //add required parameters to optimization list
    materialDataSets[matID][paramNames[i]] =
      getMatValue(matID,paramNames[i]);
  }
  /********************set materialDataMax vector********************/
  int m =0;
  while (m<paramMaxNamesCopy.size() && m<maxName.size()){
    //if the parameter is not to be optimized for then remove its range from the list
    if ((getValue(paramMaxNamesCopy[m]))== 1){
      paramMaxNamesCopy.erase(paramMaxNamesCopy.begin()+m);
      maxName.erase(maxName.begin()+m);
    }
    else{
      m++;
    }
  }

  for (int i=0;i<paramMaxNamesCopy.size();i++){
    //for the case of unconstrained optimization
    if ((getValue(paramMaxNamesCopy[i]))== 0){
      setValue(maxName[i],10000000000);
      materialDataMax[matID][paramMaxNamesCopy[i]] =
	getValue(paramMaxNamesCopy[i]);
    }
    //else add the desrired parameter range to the parameter max range vector
    else {
      materialDataMax[matID][paramMaxNamesCopy[i]] =
	getValue(paramMaxNamesCopy[i]);
    }
  }
  /********************set materialDataMin vector********************/
  int k=0;
  while (k<paramMinNamesCopy.size() && k<minName.size()){
    //if the parameter is not to be optimized for then remove its range from the list
    if ((getValue(paramMinNamesCopy[k])) == 1){
      paramMinNamesCopy.erase(paramMinNamesCopy.begin()+k);
      minName.erase(minName.begin()+k);
    }
    else{
      k++;
    }
  }

  for (int i=0;i<paramMinNamesCopy.size();i++){
    //for the case of unconstrained optimization
    if ((getValue(paramMinNamesCopy[i])) == 0){
      setValue(minName[i],-10000000000);
      materialDataMin[matID][paramMinNamesCopy[i]] =
	getValue(paramMinNamesCopy[i]);
    }
    //else add the desrired parameter range to the parameter min range vector
    else {
      materialDataMin[matID][paramMinNamesCopy[i]] =
	getValue(paramMinNamesCopy[i]);
    }
  }


#ifdef _inverseDebugMode_
  if(rank == 0) {
    cout<<"####################################################"<<endl
      cout<<"**************** extracted material data ***********"<<endl;
    for(map<string,double>::iterator q = materialDataSets[matID].begin();
        q!=materialDataSets[matID].end(); ++q)
      cout<<q->first<<" "<<q->second<<endl;
    cout<<"####################################################"<<endl
      cout<<"*********** extracted max material data ***********"<<endl;
    for(map<string,double>::iterator q = materialDataMax[matID].begin();
        q!=materialDataMax[matID].end(); ++q)
      cout<<q->first<<" "<<q->second<<endl;
    cout<<"####################################################"<<endl;
    cout<<"*********** extracted min material data ***********"<<endl;
    for(map<string,double>::iterator q = materialDataMin[matID].begin();
        q!=materialDataMin[matID].end(); ++q)
      cout<<q->first<<" "<<q->second<<endl;
  }

  logFile<<"####################################################"<<endl;
  logFile<<"**************** extracted material data ***********"<<endl;
  for(map<string,double>::iterator q = materialDataSets[matID].begin();
      q!=materialDataSets[matID].end(); ++q)
    logFile<<q->first<<" "<<q->second<<endl;
  logFile<<"####################################################"<<endl;
  logFile<<"************ extracted max material data ***********"<<endl;
  for(map<string,double>::iterator q = materialDataMax[matID].begin();
      q!=materialDataMax[matID].end(); ++q)
    logFile<<q->first<<" "<<q->second<<endl;
  logFile<<"####################################################"<<endl;
  logFile<<"*********** extracted min material data ************"<<endl;
  for(map<string,double>::iterator q = materialDataMin[matID].begin();
      q!=materialDataMin[matID].end(); ++q)
    logFile<<q->first<<" "<<q->second<<endl;
#endif

}

/***********************************************************************/
void InputFileData::resizeMaterialDataVectors(){

  using namespace std;

  // resize the materialDataSets vector if necessary
  materialDataSets.resize(materialDataSets.size() + 1);

  // resize the materialDataMax vector if necessary
  materialDataMax.resize(materialDataMax.size() + 1);

  // resize the materialDataMin vector if necessary
  materialDataMin.resize(materialDataMin.size() + 1);

}

/***********************************************************************/
void InputFileData::constructActiveMaterialDataVectors(int& matID){

  using namespace std;

  // resize the materialDataSets vector if necessary
  materialDataSets.resize(materialDataSets.size() + 1);

  // resize the materialDataMax vector if necessary
  materialDataMax.resize(materialDataMax.size() + 1);

  // resize the materialDataMin vector if necessary
  materialDataMin.resize(materialDataMin.size() + 1);

  std::string Max = "Max";
  std::string Min = "Min";
  //construct the paramMaxNames and paramMinNames vector
  for (int i=0; i<paramNames.size(); i++){
    maxName.push_back(((string)string(paramNames[i]) + Max));
    minName.push_back(((string)string(paramNames[i]) + Min));
    paramMaxNames.push_back(maxName[i].c_str());
    paramMinNames.push_back(minName[i].c_str());
    paramMaxNamesCopy.push_back(maxName[i].c_str());
    paramMinNamesCopy.push_back(minName[i].c_str());
  }

  /********************set materialDataSets vector********************/
  int n=0;
  while(n < paramNames.size()) {
    //remove parameter from optimization list if not to be optimized for
    if ((getValue(paramMinNames[n])) == 1 ||
	(getValue(paramMaxNames[n])) == 1){
      paramMaxNames.erase(paramMaxNames.begin()+n);
      paramMinNames.erase(paramMinNames.begin()+n);
      paramNames.erase(paramNames.begin()+n);
    }
    else{
      n++;;
    }
  }

  //set materialDataSets vector
  for (int i=0;i<paramNames.size();i++){
    //add required parameters to optimization list
    materialDataSets[materialDataSets.size() - 1][paramNames[i]] =
      getMatValue(matID,paramNames[i]);
  }

  /********************set materialDataMax vector********************/
  int m =0;
  while (m<paramMaxNamesCopy.size() && m<maxName.size()){
    //if the parameter is not to be optimized for then remove its range from the list
    if ((getValue(paramMaxNamesCopy[m]))== 1){
      paramMaxNamesCopy.erase(paramMaxNamesCopy.begin()+m);
      maxName.erase(maxName.begin()+m);
    }
    else{
      m++;
    }
  }

  for (int i=0;i<paramMaxNamesCopy.size();i++){
    //for the case of unconstrained optimization
    if ((getValue(paramMaxNamesCopy[i]))== 0){
      setValue(maxName[i],10000000000);
      materialDataMax[materialDataMax.size() - 1][paramMaxNamesCopy[i]] =
	getValue(paramMaxNamesCopy[i]);
    }
    //else add the desrired parameter range to the parameter max range vector
    else {
      materialDataMax[materialDataMax.size() - 1][paramMaxNamesCopy[i]] =
	getValue(paramMaxNamesCopy[i]);
    }
  }

  /********************set materialDataMin vector********************/
  int k=0;
  while (k<paramMinNamesCopy.size() && k<minName.size()){
    //if the parameter is not to be optimized for then remove its range from the list
    if ((getValue(paramMinNamesCopy[k])) == 1){
      paramMinNamesCopy.erase(paramMinNamesCopy.begin()+k);
      minName.erase(minName.begin()+k);
    }
    else{
      k++;
    }
  }

  for (int i=0;i<paramMinNamesCopy.size();i++){
    //for the case of unconstrained optimization
    if ((getValue(paramMinNamesCopy[i])) == 0){
      setValue(minName[i],-10000000000);
      materialDataMin[materialDataMin.size() - 1][paramMinNamesCopy[i]] =
	getValue(paramMinNamesCopy[i]);
    }
    //else add the desrired parameter range to the parameter min range vector
    else {
      materialDataMin[materialDataMin.size() - 1][paramMinNamesCopy[i]] =
	getValue(paramMinNamesCopy[i]);
    }
  }

}
/***********************************************************************/
void InputFileData::constructSystemicMaterialDataVectors(){

  using namespace std;

  std::string Max = "Max";
  std::string Min = "Min";
  for (int i=0; i<systParamNames.size(); i++){
    systMaxName.push_back(((string)string(systParamNames[i]) + Max));
    systMinName.push_back(((string)string(systParamNames[i]) + Min));
    systParamMaxNames.push_back(systMaxName[i].c_str());
    systParamMinNames.push_back(systMinName[i].c_str());
    systParamMaxNamesCopy.push_back(systMaxName[i].c_str());
    systParamMinNamesCopy.push_back(systMinName[i].c_str());
  }

  int n=0;
  while(n < systParamNames.size()) {
    //remove parameter from optimization list if not to be optimized for
    if ((getValue(systParamMinNames[n])) == 1 ||
	(getValue(systParamMaxNames[n])) == 1){
      paramMaxNames.erase(systParamMaxNames.begin()+n);
      paramMinNames.erase(systParamMinNames.begin()+n);
      paramNames.erase(systParamNames.begin()+n);
    }
    else{
      n++;
    }
  }

  //set materialDataSets vector

  for (int i=0;i<systParamNames.size();i++){
    //add required parameters to optimization list
    materialDataSets[materialDataSets.size() - 1][systParamNames[i]] =
      getValue(systParamNames[i]);
  }


  int m =0;
  while (m<systParamMaxNamesCopy.size() && m<systMaxName.size()){
    //if the parameter is not to be optimized for then remove its range from the list
    if ((getValue(systParamMaxNamesCopy[m]))== 1){
      systParamMaxNamesCopy.erase(systParamMaxNamesCopy.begin()+m);
      systMaxName.erase(systMaxName.begin()+m);
    }
    else{
      m++;
    }
  }

  for (int i=0;i<systParamMaxNamesCopy.size();i++){
    //for the case of unconstrained optimization
    if ((getValue(systParamMaxNamesCopy[i]))== 0){
      setValue(systMaxName[i],10000000000);
      materialDataMax[materialDataMax.size() - 1][systParamMaxNamesCopy[i]] =
	getValue(systParamMaxNamesCopy[i]);
    }
    //else add the desrired parameter range to the parameter max range vector
    else {
      materialDataMax[materialDataMax.size() - 1][systParamMaxNamesCopy[i]] =
	getValue(systParamMaxNamesCopy[i]);
    }
  }

  int k=0;
  while (k<systParamMinNamesCopy.size() && k<systMinName.size()){
    //if the parameter is not to be optimized for then remove its range from the list
    if ((getValue(systParamMinNamesCopy[k])) == 1){
      systParamMinNamesCopy.erase(systParamMinNamesCopy.begin()+k);
      systMinName.erase(systMinName.begin()+k);
    }
    else{
      k++;
    }
  }

  for (int i=0;i<systParamMinNamesCopy.size();i++){
    //for the case of unconstrained optimization
    if ((getValue(systParamMinNamesCopy[i])) == 0){
      setValue(systMinName[i],-10000000000);
      materialDataMin[materialDataMin.size() - 1][systParamMinNamesCopy[i]] =
	getValue(systParamMinNamesCopy[i]);
    }
    //else add the desrired parameter range to the parameter min range vector
    else {
      materialDataMin[materialDataMin.size() - 1][systParamMinNamesCopy[i]] =
	getValue(systParamMinNamesCopy[i]);
    }
  }
}
/***********************************************************************/
/***********************************************************************/
void InputFileData::clearArrays(std::ofstream &logFile){

  using namespace std;

  backGroundMeshInfo.clear();

  materials.clear();
  problemData.clear();
  defaultProblemData.clear();

  resizeArray(graphs,0);
  linePlotGraphData.clear();

  resizeArray(microLengths,0);

  pointDispNormalSet = false;
  pointForceNormalSet = false;
  lineDispNormalSet = false;
  lineForceNormalSet = false;
  pointRotNormalSet = false;
  pointMomentNormalSet = false;
  lineRotNormalSet = false;
  lineMomentNormalSet = false;
  pointElectricNormalSet = false;
  lineElectricNormalSet = false;
  pointMicroNormalSet = false;
  lineMicroNormalSet = false;
  pointStressNormalSet = false;
  lineStressNormalSet = false;
  pointDepolarisationNormalSet = false;
  lineDepolarisationNormalSet = false;

  // Displacement boundary conditions
  resizeArray(pointDispBoundConds,0);
  resizeArray(lineDispBoundConds,0);
  resizeArray(surfaceDispBoundConds,0);

  // Rotation boundary conditions
  resizeArray(pointRotBoundConds,0);
  resizeArray(lineRotBoundConds,0);
  resizeArray(surfaceRotBoundConds,0);

  // Electric boundary conditions
  resizeArray(pointElectricBoundConds,0);
  resizeArray(lineElectricBoundConds,0);
  resizeArray(surfaceElectricBoundConds,0);
  resizeArray(linearSurfaceElectricBoundConds,0);

  // Depolarisation boundary conditions
  resizeArray(pointDepolarisationBoundConds,0);
  resizeArray(lineDepolarisationBoundConds,0);
  resizeArray(surfaceDepolarisationBoundConds,0);

  // Micro boundary conditions
  resizeArray(pointMicroBoundConds,0);
  resizeArray(lineMicroBoundConds,0);
  resizeArray(surfaceMicroBoundConds,0);

  // stress boundary conditions
  resizeArray(pointStressBoundConds,0);
  resizeArray(lineStressBoundConds,0);
  resizeArray(surfaceStressBoundConds,0);

  // load conditions
  resizeArray(pointForceLoads,0);
  resizeArray(lineForceLoads,0);
  resizeArray(tractionLoads,0);
  resizeArray(surfacePressureLoads,0);
  resizeArray(bodyForceLoads,0);

  resizeArray(pointMomentLoads,0);
  resizeArray(lineMomentLoads,0);
  resizeArray(surfaceMomentLoads,0);
  resizeArray(bodyMomentLoads,0);

  resizeArray(surfaceElectricChargeLoads,0);
  resizeArray(bodyElectricChargeLoads,0);

  resizeArray(surfaceElems,0);

  resizeArray(resultantForceOnSurfaces,0);
  resizeArray(resultantTorqueOnSurfaces,0);

}

