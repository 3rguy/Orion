#include "InputFileData.h"

InputFileData::InputFileData(std::ofstream& logFile) {

  using namespace std;

  // Set default simulation parameters.
  setDefaultSimulationParameters(logFile);

  std::string filename = "input.dat";

  // Read and store the input file.
  readInputFile(filename,logFile);

  // Set default material parameters.
  setDefaultMaterialParameters(logFile);

}

InputFileData::InputFileData(std::string& filename,std::ofstream& logFile) {

  using namespace std;

  // Set default simulation parameters.
  setDefaultSimulationParameters(logFile);

  // Read and store the input file.
  readInputFile(filename,logFile);

  // Set default material parameters.
  setDefaultMaterialParameters(logFile);

}

/**********************************************************************/
/**********************************************************************/
// Read and store the input file.
void InputFileData::readInputFile(std::string& filename,std::ofstream& logFile) {

  using namespace std;

  // Open input file.
  ifstream inputFile(filename.c_str());
  if( !inputFile) {
    logFile << "Can't open input file " << filename << "!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(inputFile.eof()) {
    logFile << "Input file " << filename << " contains no data!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /**********************************************************************/
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  string line,name,tag,token;
  double value,value1,value2,value3;

  logFile << "###################################################" << endl;
  logFile << "############# input file data #####################" << endl;

  /**********************************************************************/
  // Read material data.
  while(inputFile >> name && name != "END") {

    if(name == "FEM_GEOMETRY" || name == "MESHLESS_GEOMETRY"
      || name == "MODEL_INPUT" || name == "SIMULATION_CONTROL_INPUT" || name == "CALCULATION_CONTROL_INPUT"
      || name == "EQUATION_SOLVER_CONTROL_INPUT"
      || name == "POSTPROCESSING_INPUT" || name == "GRAPHS_INPUT"
      || name == "RESULTANT_REACTION_INPUT") {

      if(materials.size() < 1) {
        logFile << "No material data has been stored while reading "
            << "input data!" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      logFile << name << endl;

      break;
    }
    else if(name == "" || name == "MATERIAL_INPUT") {

      logFile << name << endl;

      continue;
    }
    else if(name == "Material") {

      inputFile >> value;
      materials.resize(materials.size() + 1);
      logFile << name << " " << value << endl;

      // store macro ID
      materials[materials.size() - 1]["macroID"] = materials.size();

    }
    else if(name == "microID") {

      inputFile >> value;

      // micro ID > 1
      if(value > 1) {
        materials.resize(materials.size() + 1);

        materials[materials.size() - 1]["macroID"] = materials[materials.size()
          - 2]["macroID"];

      }

      materials[materials.size() - 1][name] = value;
      logFile << name << " " << value << endl;

    }
    else if(name == "Description") {
      
      logFile << name << " ";	//writes Description
      inputFile >> name;
      logFile << name << endl;	//writes e.g. "Linear-Hyper-Elastic(1a)"
    }
    
    else {
      inputFile >> value;
      materials[materials.size() - 1][name] = value;
      
      logFile << name << " " << value << endl;// writes e.g. "YoungsModulus 1000"
    }
    
  }

  /**********************************************************************/
  // Read FEM-geometry.
  while(inputFile >> name && name != "END") {

    if(name == "MESHLESS_GEOMETRY" || name == "MODEL_INPUT"
      || name == "SIMULATION_CONTROL_INPUT"  || name == "CALCULATION_CONTROL_INPUT"
      || name == "EQUATION_SOLVER_CONTROL_INPUT"
      || name == "POSTPROCESSING_INPUT" || name == "GRAPHS_INPUT"
      || name == "RESULTANT_REACTION_INPUT"
      || name == "BOUNDARY_CONDITIONS_INPUT") break;

    else if(name == "FEM_GEOMETRY" || name == "") {

      logFile << name << endl;

      continue;
    }
    else {
      inputFile >> value;

      if(name == "backGroundMeshElemType")

      backGroundMeshInfo["elemType"] = value;

      else if(name == "backGroundMeshElemOrder")

      backGroundMeshInfo["elemOrder"] = value;

      else if(name == "volumeGaussQuadratureOrder"
        || name == "surfaceGaussQuadratureOrder"
        || name == "lineGaussQuadratureOrder")

      backGroundMeshInfo[name] = value;
      
      else

      problemData[name] = value;

      logFile << name << " " << value << endl;

    }
  }

  // --------------------------------------------------------------------
  // set nodes and Gauss points per element
  map<string,double> params;
  intVector data(4);

  data[0] = (int) backGroundMeshInfo["elemType"];
  data[1] = (int) backGroundMeshInfo["elemOrder"];
  getFEMMeshData(data,params);

  int nodesPerElem = (int) params["nodesPerVolumeElement"];
  int nodesPerSurfaceElem = (int) params["nodesPerSurfaceElement"];
  int nodesPerLineElem = (int) params["nodesPerLineElement"];

  backGroundMeshInfo["nodesPerVolumeElement"] = params["nodesPerVolumeElement"];
  backGroundMeshInfo["nodesPerSurfaceElement"] =
    params["nodesPerSurfaceElement"];
  backGroundMeshInfo["nodesPerLineElement"] = params["nodesPerLineElement"];

  data[0] = (int) backGroundMeshInfo["elemType"];
  data[1] = (int) backGroundMeshInfo["volumeGaussQuadratureOrder"];
  data[2] = (int) backGroundMeshInfo["surfaceGaussQuadratureOrder"];
  data[3] = (int) backGroundMeshInfo["lineGaussQuadratureOrder"];

  getGaussQuadratureData(data,params);

  backGroundMeshInfo["gaussPointsPerVolumeElement"] =
    params["gaussPointsPerVolumeElement"];
  backGroundMeshInfo["gaussPointsPerSurfaceElement"] =
    params["gaussPointsPerSurfaceElement"];
  backGroundMeshInfo["gaussPointsPerLineElement"] =
    params["gaussPointsPerLineElement"];

  /**********************************************************************/
  // Read meshless geometry, common simulation and 
  // post-processing.
  while(name != "GRAPHS_INPUT" && name != "RESULTANT_REACTION_INPUT"
    && name != "BOUNDARY_CONDITIONS_INPUT" && name != "END") {

    while(getline(inputFile,line) && line.empty());

    stringstream ss(line);
    ss >> name;

    if(name == "GRAPHS_INPUT" || name == "RESULTANT_REACTION_INPUT"
      || name == "BOUNDARY_CONDITIONS_INPUT" || name == "END") break;

    else if(name == "MESHLESS_GEOMETRY" || name == "MODEL_INPUT"
      || name == "SIMULATION_CONTROL_INPUT"  || name == "CALCULATION_CONTROL_INPUT"
      || name == "EQUATION_SOLVER_CONTROL_INPUT"
      || name == "POSTPROCESSING_INPUT") {

      logFile << "*****************************************************"
          << endl;
      logFile << name << endl;

      continue;
    }
    else {

      ss >> value;
      problemData[name] = value;

      logFile << name << " " << value << endl;

    }

  }

  /**********************************************************************/
  // Read in the graph input for plotting at selected nodes.
  int coord,particle,type;

  int counter = 0;
  Graph dummyGraph;

  while(name != "RESULTANT_REACTION_INPUT"
    && name != "BOUNDARY_CONDITIONS_INPUT" && name != "END") {

    while(getline(inputFile,line) && line.empty());
    stringstream ss(line);
    ss >> name;

    if(name == "RESULTANT_REACTION_INPUT" || name == "BOUNDARY_CONDITIONS_INPUT"
      || name == "END") break;

    else if(name == "GRAPHS_INPUT" && ss.rdstate() == ios::failbit) {

      //logFile << "*****************************************************" << endl;
      //logFile << name << endl;
      continue;
    }

    // ------------------------------------------------------------------
    // graphs input (DOFs must be at the end)
    //
    // load_deformation
    // loadConditionID 1
    // dirichletConditionID 1
    // node 12
    // dof 1
    // dof 3

    // allocate a new graph
    else if(dummyGraph.isType(name)) {

      counter = graphs.size();
      graphs.resize(graphs.size() + 1);

      graphs[counter].getID() = counter;
      graphs[counter].getType() = name;

    }
    // read in graph parameters/properties
    else {
      ss >> value;

      // check whether value could be read
      if(ss.rdstate() == ios::failbit) {
        logFile << "In InputFileData::readInputFile 'GRAPHS_INPUT' incorrect!"
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);

      }
      // another DOF at the same particle -> create a new graph
      else if(name == "DOF" && graphs[counter].getDOF() != 0) {
        counter = graphs.size();
        resizeArray(graphs,graphs.size() + 1);

        graphs[counter] = graphs[counter - 1]; // copy everything
        graphs[counter].setParam(name,value); // overwrite DOF

      }
      else if(name == "ID") continue; // ignored
      else {
        graphs[counter].setParam(name,value);

        // DOFs must be the last parameters read in.
        if(name != "DOF" && graphs[counter].getDOF() != 0) {
          logFile
              << "In InputFileData::readInputFile 'GRAPHS_INPUT' incorrect.\n"
              << "DOFs must be the last graph parameters read in." << endl;
          MPI_Abort(MPI_COMM_WORLD,1);
        }

      }

    }

  }

//  // no graphs defined
//  if(graphs.size() == 0) {
//    counter = graphs.size();
//    graphs.resize(1);
//    graphs[counter].getLoadingConditionID() = 1;
//    graphs[counter].getID() = counter;
//    graphs[counter].getNode() = 1;
//    graphs[counter].getType() = "load_deformation";
//    graphs[counter].getDOF() = 1;
//  }

  // ---------------------------------------------------------------------
  // enable some other graphs

  vector<string> otherGraphsTypes;

  if((bool) problemData["plotSystemEnergyGraph "]) pushBackVector(
      otherGraphsTypes,(string) "energy_time");

  // default cardiac mechanics graphs
  if((bool) problemData["cardiacMechanicsProblem"]) {

    if((bool) problemData["plotPressureVolumeGraph"]) pushBackVector(
        otherGraphsTypes,(string) "pressure_volume");

    if((bool) problemData["plotPressureTimeGraph"]) pushBackVector(
        otherGraphsTypes,(string) "load_time");

    if((bool) problemData["plotVolumeTimeGraph"]) pushBackVector(
        otherGraphsTypes,(string) "volume_time");

  }

  for(int i = 0;i < otherGraphsTypes.size();i++) {

    bool typeSet = false;

    for(int j = 0;j < graphs.size();j++) {

      if(graphs[j].getType() == otherGraphsTypes[i]) {
        typeSet = true;
        break;
      }

    }

    if( !typeSet && otherGraphsTypes[i] == "pressure_volume") {
      counter = graphs.size();
      graphs.resize(graphs.size() + 1);
      graphs[counter].getLoadingConditionID() = 0;
      graphs[counter].getDirichletControlID() = 0;
      graphs[counter].getID() = counter;
      graphs[counter].getType() = otherGraphsTypes[i];
    }
    else if( !typeSet && otherGraphsTypes[i] == "load_time") {
      counter = graphs.size();
      graphs.resize(graphs.size() + 1);
      graphs[counter].getLoadingConditionID() = 0;
      graphs[counter].getID() = counter;
      graphs[counter].getType() = otherGraphsTypes[i];
    }
    else if( !typeSet && otherGraphsTypes[i] == "volume_time") {
      counter = graphs.size();
      graphs.resize(graphs.size() + 1);
      graphs[counter].getDirichletControlID() = 0;
      graphs[counter].getID() = counter;
      graphs[counter].getType() = otherGraphsTypes[i];
    }

  }
  
  // Write output in the log file.
  logFile << "**************************************************" << endl;
  logFile << "GRAPHS_INPUT" << endl;
  for(int i = 0;i < graphs.size();i++) {
    int& ID = graphs[i].getID();
    string& type = graphs[i].getType();
    logFile << i << ".) " << type << " ID=" << ID << ": " << endl;
    logFile << "abscissa-" << graphs[i].getAbscissaLabel() << "scaling="
        << graphs[i].getAbscissaScaling() << endl;
    logFile << "ordinate-" << graphs[i].getOrdinateLabel() << "scaling="
        << graphs[i].getOrdinateScaling() << endl;
    if(graphs[i].getLoadingConditionID() != -1) logFile << "loadCondID="
        << graphs[i].getLoadingConditionID() << endl;
    if(graphs[i].getDirichletConditionID() != -1) logFile << "DirichletCondID="
        << graphs[i].getDirichletConditionID() << endl;
    if(graphs[i].getDirichletControlID() != -1) logFile << "DirichletControlID="
        << graphs[i].getDirichletControlID() << endl;
    logFile << "node=" << graphs[i].getNode() << endl;
    logFile << "DOF=" << graphs[i].getDOF() << endl;
  }

  /*********************************************************************/
  // resultant force and torque on the specified surfaces
  //int currentElem,currentSurface,currentVolume,elem,volume;
  bool resultantInternalTraction = false;
  bool resultantInternalSurfaceTorque = false;

  int currentElem,ID,dof,count;
  int usedDims = 3;

  Condition* condition;
  intVector* condElemNodes;

  Condition dummyCond;

  while(name != "LOAD_CONDITIONS_INPUT" && name != "BOUNDARY_CONDITIONS_INPUT"
    && name != "END") {

    if(dummyCond.isReactionType(name)) {

      while(getline(inputFile,line) && line.empty());
      stringstream ss2(line);
      ss2 >> tag >> ID; // ID
      ID -= 1;

      if(ID < resultantReactions.size()
        && resultantReactions[ID].getID() != 0) {
        logFile << "In InputFileData::readInputFile incorrectly written '"
            << name << "' input data, resultant reaction condition ID=" << ID
            << " is invalid!\n"
            << " Resultant reaction conditions need to have unique and consecutive IDs."
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }
      else if(ID >= resultantReactions.size()) resultantReactions.resize(
          ID + 1);

      condition = &resultantReactions[ID];
    }

    else if(name != "END" && name != "BOUNDARY_CONDITIONS_INPUT"
      && !dummyCond.isReactionType(name)) {

      while(getline(inputFile,line) && line.empty())
        ;
      stringstream ss(line);
      ss >> name;
      continue;
    }
    else if(name == "BOUNDARY_CONDITIONS_INPUT" || name == "END") break;
    else {
      logFile << "In InputFileData::readInputFile resultant reaction '" << name
          << "' is not supported." << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    // ================================================================
    // read-in the properties of the current condition
    int& condID = ( *condition).getID();
    string& condType = ( *condition).getType();

    vector<ConditionElement>& condElems = ( *condition).getElements();
    blVector& condDOFs = ( *condition).getConditionDOFs();

    conditionSet[name] = true;
    condID = ID;
    condType = name;

    currentElem = 0;

    // condition parameters

    while(getline(inputFile,line)) {

      if(line.empty()) continue;

      stringstream ss2(line);
      ss2 >> tag;

      if(( *condition).isReactionType(tag) || tag == "BOUNDARY_CONDITIONS_INPUT"
        || tag == "LOAD_CONDITIONS_INPUT" || tag == "END") {
        name = tag;
        break; // new reaction or end of reaction block
      }
      else if(tag == "Element") {
        currentElem = condElems.size();
        condElems.resize(currentElem + 1);
        int& elemID = condElems[currentElem].getID();
        ss2 >> elemID;
      }
      else if(tag == "nodes") {

        intVector& elemNodes = condElems[currentElem].getNodes();

        // loop over all nodes of current element
        while(ss2 >> particle && ss2.rdstate() != ios::failbit)
          elemNodes.push_back(particle);

      }

      // format e.g. "DOF <dof> value <value>
      else if(tag == "DOF") {

        ss2 >> dof;

        if(dof <= condDOFs.size()) {
          condDOFs[dof - 1] = true;
        }
        else {
          resizeArray(condDOFs,dof);
          condDOFs[dof - 1] = true;
        }

      }

    }

  }

  /**********************************************************************/
  // Read the Dirichlet and loading conditions.
  while(name != "END") {

    if(dummyCond.isDirichletType(name)) {

      while(getline(inputFile,line) && line.empty());
      stringstream ss2(line);
      ss2 >> tag >> ID; // ID
      ID -= 1;

      if(ID >= dirichletConditions.size()) dirichletConditions.resize(
          ID + 1);
      condition = &dirichletConditions[ID];

      if(dirichletConditions[ID].getType() != "none") {
        logFile << "In InputFileData::readInputFile condition '" << name
            << "' has not a unique ID." << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }
    }

    else if(dummyCond.isLoadingType(name)) {

      while(getline(inputFile,line) && line.empty());
      stringstream ss2(line);
      ss2 >> tag >> ID; // ID
      ID -= 1;

      if(ID >= loadingConditions.size()) loadingConditions.resize(ID + 1);
      condition = &loadingConditions[ID];

      if(loadingConditions[ID].getType() != "none") {
        logFile << "In InputFileData::readInputFile condition '" << name
            << "' has not a unique ID." << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

    }
    else if(dummyCond.isDirichletControlType(name)) {

      while(getline(inputFile,line) && line.empty());
      stringstream ss2(line);
      ss2 >> tag >> ID; // ID
      ID -= 1;

      if(ID >= dirichletControlConditions.size()) dirichletControlConditions.resize(ID + 1);
      condition = &dirichletControlConditions[ID];

      if(dirichletControlConditions[ID].getType() != "none") {
        logFile << "In InputFileData::readInputFile condition '" << name
            << "' has not a unique ID." << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

    }
    else if(name != "END" && !dummyCond.isDirichletType(name)
      && !dummyCond.isLoadingType(name)
      && !dummyCond.isDirichletControlType(name)) {

      while(getline(inputFile,line) && line.empty());
      stringstream ss(line);
      ss >> name;
      continue;
    }
    else if(name == "END") break;
    else {
      logFile << "In InputFileData::readInputFile condition '" << name
          << "' is not supported." << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // ================================================================
    // read-in the properties of the current condition
    int& condID = ( *condition).getID();
    string& condType = ( *condition).getType();

    vector<ConditionElement>& condElems = ( *condition).getElements();
    dbVector& condNormal = ( *condition).getSurfaceNormal();
    blVector& condDOFs = ( *condition).getConditionDOFs();
    dbVector& condValues = ( *condition).getConditionValues();
    allocateArray(condNormal,usedDims);

    conditionSet[name] = true;
    condID = ID;
    condType = name;

    currentElem = 0;

    // condition parameters

    while(getline(inputFile,line)) {

      if(line.empty()) continue;

      stringstream ss2(line);
      ss2 >> tag;

      if(( *condition).isDirichletType(tag) || ( *condition).isLoadingType(tag)
        || ( *condition).isDirichletControlType(tag)
        || tag == "LOAD_CONDITIONS_INPUT" || tag == "END") {
        name = tag;
        break; // new condition or end of Dirichlet conditions
      }
      else if(tag == "Element") {
        currentElem = condElems.size();
        condElems.resize(currentElem + 1);
        int& elemID = condElems[currentElem].getID();
        ss2 >> elemID;
      }
      else if(tag == "node") {
        ss2 >> particle;
        ( *condition).getNodes(logFile).push_back(particle);
        continue;
      }
      else if(tag == "nodes") {

        intVector& elemNodes = condElems[currentElem].getNodes();

        // loop over all nodes of current element
        while(ss2 >> particle && ss2.rdstate() != ios::failbit)
          elemNodes.push_back(particle);

      }
      else if(tag == "normal") {
        ss2 >> condNormal[0] >> condNormal[1] >> condNormal[2];
      }

      // format e.g. "DOF <dof> value <value>
      else if(tag == "DOF") {

        ss2 >> dof >> tag >> value;

        if(dof <= condDOFs.size()) {
          condDOFs[dof - 1] = true;
          condValues[dof - 1] = value;
        }
        else {
          resizeArray(condDOFs,dof);
          resizeArray(condValues,dof);
          condDOFs[dof - 1] = true;
          condValues[dof - 1] = value;
        }

      }
      // format e.g. value <value>
      else if(tag == "amount") {

        ss2 >> value;

        resizeArray(condDOFs,1);
        resizeArray(condValues,1);
        condDOFs[0] = true;
        condValues[0] = value;

      }
      // condition application function
      else if(tag == "function") {

        ss2 >> value;  // function ID
        ( *condition).getFunction() = value;
        dbVector& functionParams = ( *condition).getFunctionParams();

        // loop over all function parameters: a0, a1, etc.
        if(functionParams.size() == 0) {

          while(ss2 >> tag >> value && ss2.rdstate() != ios::failbit)

            pushBackVector(functionParams,value);

        }
        // already stored
        else while(ss2 >> tag >> value && ss2.rdstate() != ios::failbit);

      }
      // some other parameters
      else {
        ss2 >> value;
        ( *condition).setParam(tag,value);
      }

    }

  }

  // Write output in the log file.
  logFile << "*****************************************************" << endl;
  logFile << "BOUNDARY_CONDITIONS_INPUT" << endl;
  for(int i = 0;i < dirichletConditions.size();i++) {
    logFile << "---------------------------------------------------" << endl;
    Condition& condition = dirichletConditions[i];
    logFile << condition.getType() << "ID=" << condition.getID() << endl;
    logFile << "controlMode=" << condition.getControlMode() << endl;
    logFile << "function=" << condition.getFunction() << endl;
    for(int j = 0;j < condition.getFunctionParams().size();j++)
      logFile << "a" << j << "=" << condition.getFunctionParams()[j] << " ";
    logFile << endl;
    logFile << "factor=" << condition.getFactor() << endl;
    logFile << "previousFactor=" << condition.getPreviousFactor() << endl;
    logFile << "maxFactor=" << condition.getMaxFactor() << endl;
    logFile << "minFactor=" << condition.getMinFactor() << endl;
    dbVector& condNormal = condition.getSurfaceNormal();
    for(int k = 0;k < condNormal.size();k++)
      logFile << "normal[" << k << "]=" << condNormal[k] << endl;
    blVector& condDOFs = condition.getConditionDOFs();
    dbVector& condValues = condition.getConditionValues();
    for(int k = 0;k < condDOFs.size();k++) {
      if(condDOFs[k]) logFile << "condition[" << k << "]: DOF=" << k
          << " value=" << condValues[k] << endl;
    }
    vector<ConditionElement>& condElems = condition.getElements();
    for(int k = 0;k < condElems.size();k++) {
      logFile << "Element " << condElems[k].getID() << endl;
      logFile << "nodes: ";
      intVector& elemNodes = condElems[k].getNodes();
      for(int l = 0;l < elemNodes.size();l++)
        logFile << elemNodes[l] << " ";
      logFile << endl;
    }
  }
  logFile << "======================================================" <<endl;
  for(int i = 0;i < dirichletControlConditions.size();i++) {
    logFile << "---------------------------------------------------" << endl;
    Condition& condition = dirichletControlConditions[i];
    logFile << condition.getType() << "ID=" << condition.getID() << endl;
    logFile << "controlMode=" << condition.getControlMode() << endl;
    logFile << "function=" << condition.getFunction() << endl;
    for(int j = 0;j < condition.getFunctionParams().size();j++)
      logFile << "a" << j << "=" << condition.getFunctionParams()[j] << " ";
    logFile << endl;
    logFile << "factor=" << condition.getFactor() << endl;
    logFile << "previousFactor=" << condition.getPreviousFactor() << endl;
    logFile << "maxFactor=" << condition.getMaxFactor() << endl;
    logFile << "minFactor=" << condition.getMinFactor() << endl;
    dbVector& condNormal = condition.getSurfaceNormal();
    for(int k = 0;k < condNormal.size();k++)
      logFile << "normal[" << k << "]=" << condNormal[k] << endl;
    blVector& condDOFs = condition.getConditionDOFs();
    dbVector& condValues = condition.getConditionValues();
    for(int k = 0;k < condDOFs.size();k++) {
      if(condDOFs[k]) logFile << "condition[" << k << "]: DOF=" << k
          << " value=" << condValues[k] << endl;
    }
    vector<ConditionElement>& condElems = condition.getElements();
    for(int k = 0;k < condElems.size();k++) {
      logFile << "Element " << condElems[k].getID() << endl;
      logFile << "nodes: ";
      intVector& elemNodes = condElems[k].getNodes();
      for(int l = 0;l < elemNodes.size();l++)
        logFile << elemNodes[l] << " ";
      logFile << endl;
    }
  }
  logFile << "*****************************************************" << endl;
  logFile << "LOADING_CONDITIONS" << endl;
  for(int i = 0;i < loadingConditions.size();i++) {
    logFile << "---------------------------------------------------" << endl;
    Condition& condition = loadingConditions[i];
    logFile << condition.getType() << "ID=" << condition.getID() << endl;
    logFile << "function=" << condition.getFunction() << endl;
    logFile << "factor=" << condition.getFactor() << endl;
    logFile << "previousFactor=" << condition.getPreviousFactor() << endl;
    logFile << "maxFactor=" << condition.getMaxFactor() << endl;
    logFile << "minFactor=" << condition.getMinFactor() << endl;
    if((bool) problemData["cardiacMechanicsProblem"] && condition.getType() == "Surface-Pressure-Loading") {
      logFile << "endDiastolicPressure=" << condition.getEndDiastolicPressure() << endl;
      logFile << "endIVCPressure=" << condition.getEndIVCPressure() << endl;
    }
    dbVector& condNormal = condition.getSurfaceNormal();
    for(int k = 0;k < condNormal.size();k++)
      logFile << "normal[" << k << "]=" << condNormal[k] << endl;
    blVector& condDOFs = condition.getConditionDOFs();
    dbVector& condValues = condition.getConditionValues();
    for(int k = 0;k < condDOFs.size();k++) {
      if(condDOFs[k]) logFile << "condition[" << k << "]: DOF=" << k
          << " value=" << condValues[k] << endl;
    }
    vector<ConditionElement>& condElems = condition.getElements();
    for(int k = 0;k < condElems.size();k++) {
      logFile << "Element " << condElems[k].getID() << endl;
      logFile << "nodes: ";
      intVector& elemNodes = condElems[k].getNodes();
      for(int l = 0;l < elemNodes.size();l++)
        logFile << elemNodes[l] << " ";
      logFile << endl;
    }
  }

  // =====================================================================
  // check condition numbering
  blVector usedConditionIDs(dirichletConditions.size());

  for(int i = 0;i < dirichletConditions.size();i++) {
    Condition& condition = dirichletConditions[i];
    if(condition.getID() > dirichletConditions.size()) break;
    else if(condition.getID() > -1)
      usedConditionIDs[condition.getID()]=true;
  }

  for(int i=0;i<usedConditionIDs.size();i++) {

    if(!usedConditionIDs[i]) {
      logFile << "In InputFileData::readInputFile 'BOUNNDARY CONDITION' "
      << "input data  don't have unique and consecutive IDs."<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }
  resizeArray(usedConditionIDs,dirichletControlConditions.size());
  clearArray(usedConditionIDs);

  for(int i = 0;i < dirichletControlConditions.size();i++) {
    Condition& condition = dirichletControlConditions[i];
    if(condition.getID() > dirichletControlConditions.size()) break;
    else if(condition.getID() > -1)
      usedConditionIDs[condition.getID()]=true;
  }

  for(int i=0;i<usedConditionIDs.size();i++) {

    if(!usedConditionIDs[i]) {
      logFile << "In InputFileData::readInputFile 'BOUNNDARY CONDITION' "
      << "input data  don't have unique and consecutive IDs."<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }
  resizeArray(usedConditionIDs,loadingConditions.size());
  clearArray(usedConditionIDs);

  for(int i = 0;i < loadingConditions.size();i++) {
    Condition& condition = loadingConditions[i];
    if(condition.getID() > loadingConditions.size()) break;
    else if(condition.getID() > -1)
      usedConditionIDs[condition.getID()]=true;
  }

  for(int i=0;i<usedConditionIDs.size();i++) {

    if(!usedConditionIDs[i]) {
      logFile << "In InputFileData::readInputFile 'LOADING CONDITION' "
      << "input data don't have unique and consecutive IDs."<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

  /**********************************************************************/
  // all problem data (including default ones)
  logFile << "####################################################" << endl;
  logFile << "******* input parameter (incl. default values) ******" << endl;
  for(map<string,double>::const_iterator p = problemData.begin();
      p != problemData.end();++p)
    logFile << p->first << " " << p->second << endl;

}

/***********************************************************************/
/***********************************************************************/
// Set default simulation parameters.
void InputFileData::setDefaultSimulationParameters(std::ofstream& logFile) {

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
  problemData["minMLSSupport"] = 2;

  problemData["dispMatrixRearrange"] = 1;


  //calc control mode 6 stuff
  problemData["cardiacPhase"] = 1;

  // FEM stuff
  problemData["gaussPointsPerVolumeElement"] = 8;
  problemData["gaussPointsPerSurfaceElement"] = 4;
  problemData["gaussPointsPerLineElement"] = 2;
  problemData["increaseElemOrderOnBoundary"] = 0;
  problemData["lineToSurfaceElementConversion"] = 0;

  /**********************************************************************/
  // preprocessing
  problemData["numberOfSortingPlanes"] = 0;
  problemData["dimensionlessBasisPolynom"] = 0;
  problemData["shapefunctionType"] = 2;
  problemData["windowFunctionShape"] = 1;

  problemData["modifiedMLS"] = 0;
  problemData["supportComputationMode"] = 2;

  // --------------------------------------------------------------------
  // asymmetric shape functions
  problemData["windowfunctionNorming"] = 1;
  problemData["normedBasisPolynom"] = 1;

  // --------------------------------------------------------------------
  // maximum entropy stuff
  problemData["maxEntLocalPtcls"] = 8;
  problemData["MaxEntGamma"] = 1.8;
  problemData["target_zero"] = 1e-5;
  problemData["TolLag"] = 1e-6;
  
  // --------------------------------------------------------------------
  // (1) location optimized (default) (2) computing load balanced
  problemData["parallelPtcleDistributingMode"] = 1;

  // particle influence radius stuff
  problemData["influenceRadiusMultiplier"] = 1;
  problemData["mininumInfluenceRadiusTolerance"] = 1.0e-03;
  problemData["plusMinusDirectionDependentRadius"] = 0;
  problemData["positiveNegativeRadiusRatio"] = 1.0e-03;

  problemData["radiusDistanceDepend"] = 1; // obsolete
  problemData["minDirectionalPtcleSupport"] = 2;
  problemData["minDirectPtcleSuppReduction"] = 2;
  problemData["influenceRadiusDeterminationAngle"] = 0;

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
  problemData["dynamicSimulation"] = 0;

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
  problemData["plotMaxIncrement"] = 20;
  problemData["plotMinIncrement"] = 1;
  problemData["constantPlotIncrement"] = 0;
  problemData["plotOnParticles"] = 1;
  problemData["plotOnGaussPoints"] = 0;

  problemData["stressFieldScaling"] = 1.0;

  problemData["computeDeformedStructuralVolumes"] = 0;

  problemData["trueStressPlotting"] = 1;

  problemData["plotCurrentMesh"] = 0; // unknown reference config.

  problemData["plotDeformationDistribution"] = 0;
  problemData["plotStrainDistribution"] = 1;
  problemData["plotPlasticStrainDistribution"] = 0;
  problemData["plotStressDistribution"] = 0;
  problemData["plotCoupleStressDistribution"] = 0;
  problemData["plotSurfaceNormalDistribution"] = 1;
  problemData["plotLineTangentDistribution"] = 0;
  problemData["plotFibreDirectionDistribution"] = 0.0;

  problemData["twistGraphPlot"] = 0.0;
  problemData["plotSystemEnergyGraph"] = 0.0;

  /*********************************************************************/
  // simulation control
  problemData["restartFileID"] = 0;
  problemData["restartConditionIncrementScaling"] = 1.0;

  problemData["timeIncrement"] = 1.0e-02;
  problemData["maxTimeIncrement"] = problemData["timeIncrement"] * 50.0;
  problemData["minTimeIncrement"] = problemData["timeIncrement"] * 0.001;
  problemData["maxSimulationTime"] = 1.0e+16;

  problemData["KmatUpdatingWhileIterating"] = 1; // works only for a very small loading increment
  problemData["deallocateParallelKmat"] = 0;
  problemData["calcStiffnessTangentDeterminant"] = 0;
  problemData["symmetrizeKmat"] = 0;

  problemData["constantLoadingRate"] = 1;
  problemData["pressureControl"] = 0;

  problemData["startPerturbation"] = 0;
  problemData["determineEQSystemMinMaxValues"] = 0;
  problemData["equationSystemMultiplier"] = 1;

  problemData["convergenceType"] = 0.0;
  problemData["divergenceSimulationReset"] = 1.0;
  problemData["fastConvergenceStepNumTol"] = 4.0;
  problemData["slowConvergenceStepNumTol"] =
    problemData["fastConvergenceStepNumTol"] + 1;
  problemData["slowConvergenceStepNumTolWhenYielded"] =
    problemData["slowConvergenceStepNumTol"];

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
  // PETSc equation solver settings
  problemData["numberOfProcsPerNode"] = 1;
  problemData["PCSubMatSizeLimit"] = 1.0e+06;
  problemData["subBlockKSPDefinition"] = 0;
  problemData["subBlockSolverPreconditioner"] = 5;
  problemData["gmresRestartIterationNum"] = 200;
  problemData["tangentSymmetric"] = 0;
  problemData["tangentZeroPivotShift"] = 0;
  problemData["directSolver"] = 3; // 1: PETSc, 2: MUMPS, 3: SuperLU

  /*********************************************************************/
  //for cuda
  problemData["cudaAssembly"] = 0;
  problemData["numThreads"] = 1;

  /*********************************************************************/
  // debug stuff
  problemData["checkShapefunctions"] = 0;

  /*********************************************************************/
  // cardiac mechanics stuff
  // controls loading in cardiac simulations
  // this overrides the setup loading
  problemData["cardiacMechanicsProblem"] = 0;
  problemData["plotPressureVolumeGraph"] = 0;
  problemData["plotVolumeTimeGraph"] = 0;
  problemData["plotPressureTimeGraph"] = 0;

  // 0 - if time dependent active tension
  // 1 - if steady active tension
  problemData["steadyActiveTension"] = 1;
  
  // adaptive simulation increment method
  problemData["adaptiveConditionIncrementType"] = 1;

  // cardiac simulation parameters
  problemData["diastolePressureControlled"] = 1;
  problemData["endSystolicPressure"] = 0.25e3;
  problemData["ejectionPressure"] = 5.0e3;
  problemData["endSystolicTime"] = 1250e-3;

  problemData["maxVolumeIterations"] = 10;
  problemData["volumeTolerance"] = 1.0e-08;

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

  //for active tension plotting
  problemData["plotActiveTensionDistribution"] = 0;
  problemData["activeTensionFieldScaling"] = 1;

  problemData["plotContractileLengthDistribution"] = 0;
  problemData["plotSarcomereLengthDistribution"] = 0;

  // Guccione active stress model
//  problemData["TMax"] = 135.7 * 1000; // tension Pa
//  problemData["Ca0Max"] = 4.35;
//  problemData["Ca0"] = 4.35;
//  problemData["m"] = 1.0489;
//  problemData["b"] = -1.429;
//  problemData["B"] = 4.75;
//  problemData["l0"] = 1.58;
//  problemData["lR"] = 1.85; //average value
//  problemData["lRMin"] = 1.78;
//  problemData["lRMax"] = 1.91;
//  problemData["maxLoadingFactor"] = 1;

  // ---------------------------------------------------------------------
  // eikonal equation
  problemData["plotPropagationSpeedDistribution"] = 0.0;
  problemData["plotDepolarisationTimeDistribution"] = 0;
  problemData["depolarisationTimeFieldScaling"] = 1;

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

  // deformed Lagrangian computation
  problemData["inverseProblemType"] = 0;

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
// Set default material parameters.
void InputFileData::setDefaultMaterialParameters(std::ofstream& logFile) {

  using namespace std;

  // check whether each material has a default microID

  for(int i = 0;i < materials.size();i++) {

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

    materials[i]["microConstitutiveLawID"] = materials[i]["constitutiveLawID"];

  }

  // --------------------------------------------------------------------
  // Mises plasticity: set default material weakness coordinates 

  for(int i = 0;i < materials.size();i++) {

    // check whether the material is Cauchy-based viscoplasticity or
    // logarithmic stretch tensor-based von Mises plasticity
    map<string,double>::iterator p = materials[i].find("constitutiveLawID");

    map<string,double>::iterator q1 = materials[i].find(
        "materialWeaknessCoordX");
    map<string,double>::iterator q2 = materials[i].find(
        "materialWeaknessCoordY");
    map<string,double>::iterator q3 = materials[i].find(
        "materialWeaknessCoordZ");
    map<string,double>::iterator q4 = materials[i].find(
        "materialWeaknessSigmaYreduction");
    map<string,double>::iterator q5 = materials[i].find(
        "materialWeaknessRadius");

    if(p != materials[i].end() && (p->second == 7 || p->second == 16)
      && (q1 == materials[i].end() || q2 == materials[i].end()
        || q3 == materials[i].end() || q4 == materials[i].end()
        || q5 == materials[i].end())) {

      materials[i]["materialWeaknessCoordX"] = 0.0;
      materials[i]["materialWeaknessCoordY"] = 0.0;
      materials[i]["materialWeaknessCoordZ"] = 0.0;
      materials[i]["materialWeaknessSigmaYreduction"] = 1.0;
      materials[i]["materialWeaknessRadius"] = 0.0;
    }

    p = materials[i].find("microConstitutiveLawID");

    if(p != materials[i].end() && (p->second == 7 || p->second == 16)
      && (q1 == materials[i].end() || q2 == materials[i].end()
        || q3 == materials[i].end() || q4 == materials[i].end()
        || q5 == materials[i].end())) {

      materials[i]["materialWeaknessCoordX"] = 0.0;
      materials[i]["materialWeaknessCoordY"] = 0.0;
      materials[i]["materialWeaknessCoordZ"] = 0.0;
      materials[i]["materialWeaknessSigmaYreduction"] = 1.0;
      materials[i]["materialWeaknessRadius"] = 0.0;
    }

  }

  //#ifdef _commonDebugMode_
  logFile << "*****************************************************" << endl;
  logFile << "******* material input incl. default settings *******" << endl;
  for(int i = 0;i < materials.size();i++) {
    logFile << "MATERIAL " << i << ":" << endl;
    for(map<string,double>::const_iterator r = materials[i].begin();
        r != materials[i].end();++r)
      logFile << r->first << " " << r->second << endl;
  }
  //#endif

}

/***********************************************************************/
/***********************************************************************/
// delete not used conditions and adjust the array allocation
void InputFileData::adjustConditionAllocation(
    InputFileData* InputData,std::map<std::string,double>& modelData,
    std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int) modelData["usedDimensions"];
  int rotDOF = (int) modelData["rotationDegreesOfFreedom"];
  int defDOF = (int) modelData["deformationDegreesOfFreedom"];
  int microDOF = (int) modelData["microDegreesOfFreedom"];
  int stressDOF = (int) modelData["stressDegreesOfFreedom"];
  int electricDOF = (int) modelData["electricDegreesOfFreedom"];
  int depolarisationDOF = (int) modelData["depolarisationDegreesOfFreedom"];

  /**********************************************************************/
  // delete not used conditions
  intVector newOrder(dirichletConditions.size());

  int m = 0;
  int n = dirichletConditions.size() - 1;

  for(int i = 0;i < dirichletConditions.size();i++) {
    Condition& condition = dirichletConditions[i];

    if(rotDOF == 0
      && (condition.getType() == "Rotation-Point-Constraint"
        || condition.getType() == "Rotation-Line-Constraint"
        || condition.getType() == "Rotation-Surface-Constraint")) {
      newOrder[n] = i;
      n--;
    }

    else if(stressDOF == 0
      && (condition.getType() == "Stress-Point-Constraint"
        || condition.getType() == "Stress-Line-Constraint"
        || condition.getType() == "Stress-Surface-Constraint")) {
      newOrder[n] = i;
      n--;
    }

    else if(electricDOF == 0
      && (condition.getType() == "Electric-Point-Constraint"
        || condition.getType() == "Electric-Line-Constraint"
        || condition.getType() == "Electric-Surface-Constraint")) {
      newOrder[n] = i;
      n--;
    }

    else if(microDOF == 0
      && (condition.getType() == "Micro-Point-Constraint"
        || condition.getType() == "Micro-Line-Constraint"
        || condition.getType() == "Micro-Surface-Constraint")) {
      newOrder[n] = i;
      n--;
    }

    else if(depolarisationDOF == 0
      && (condition.getType() == "Depolarisation-Time-Point-Constraint"
        || condition.getType() == "Depolarisation-Time-Line-Constraint"
        || condition.getType() == "Depolarisation-Time-Surface-Constraint")) {
      newOrder[n] = i;
      n--;
    }

    else {

      newOrder[m] = i;
      m++;
    }

  }

  reorderVector(dirichletConditions,newOrder);
  resizeArray(dirichletConditions,m);

  // ------------------------------

  resizeArray(newOrder,loadingConditions.size());
  m = 0;
  n = loadingConditions.size() - 1;

  for(int i = 0;i < loadingConditions.size();i++) {
    Condition& condition = loadingConditions[i];

    if(rotDOF == 0
      && (condition.getType() == "Point-Moment-Loading"
        || condition.getType() == "Surface-Moment-Loading"
        || condition.getType() == "Body-Moment-Loading")) {
      newOrder[n] = i;
      n--;
    }
    else if(electricDOF == 0
      && (condition.getType() == "Electric-Surface-Charge-Loading"
        || condition.getType() == "Electric-Body-Charge-Loading")) {
      newOrder[n] = i;
      n--;
    }
    else {

      newOrder[m] = i;
      m++;
    }

  }

  reorderVector(loadingConditions,newOrder);
  resizeArray(loadingConditions,m);

  /*********************************************************************/
  // adjust the array allocation
  for(int i = 0;i < dirichletConditions.size();i++) {
    Condition& condition = dirichletConditions[i];

    if(condition.getType() == "Displacement-Point-Constraint") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Rotation-Point-Constraint") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Micro-Point-Constraint") {
      resizeArray(condition.getConditionDOFs(),microDOF);
      resizeArray(condition.getConditionValues(),microDOF);
    }
    else if(condition.getType() == "Stress-Point-Constraint") {
      resizeArray(condition.getConditionDOFs(),stressDOF);
      resizeArray(condition.getConditionValues(),stressDOF);
    }
    else if(condition.getType() == "Electric-Point-Constraint") {
      resizeArray(condition.getConditionDOFs(),electricDOF);
      resizeArray(condition.getConditionValues(),electricDOF);
    }
    else if(condition.getType() == "Depolarisation-Time-Point-Constraint") {
      resizeArray(condition.getConditionDOFs(),1);
      resizeArray(condition.getConditionValues(),1);
    }
    else if(condition.getType() == "Displacement-Line-Constraint") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Rotation-Line-Constraint") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Micro-Line-Constraints") {
      resizeArray(condition.getConditionDOFs(),microDOF);
      resizeArray(condition.getConditionValues(),microDOF);
    }
    else if(condition.getType() == "Electric-Line-Constraint") {
      resizeArray(condition.getConditionDOFs(),electricDOF);
      resizeArray(condition.getConditionValues(),electricDOF);
    }
    else if(condition.getType() == "Stress-Line-Constraints") {
      resizeArray(condition.getConditionDOFs(),stressDOF);
      resizeArray(condition.getConditionValues(),stressDOF);
    }
    else if(condition.getType() == "Depolarisation-Time-Line-Constraint") {
      resizeArray(condition.getConditionDOFs(),1);
      resizeArray(condition.getConditionValues(),1);
    }
    else if(condition.getType() == "Displacement-Surface-Constraint") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Rotation-Surface-Constraint") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Electric-Surface-Constraint") {
      resizeArray(condition.getConditionDOFs(),electricDOF);
      resizeArray(condition.getConditionValues(),electricDOF);
    }
    else if(condition.getType() == "Micro-Surface-Constraint") {
      resizeArray(condition.getConditionDOFs(),microDOF);
      resizeArray(condition.getConditionValues(),microDOF);
    }
    else if(condition.getType() == "Stress-Surface-Constraint") {
      resizeArray(condition.getConditionDOFs(),stressDOF);
      resizeArray(condition.getConditionValues(),stressDOF);
    }
    else if(condition.getType() == "Depolarisation-Time-Surface-Constraint") {
      resizeArray(condition.getConditionDOFs(),1);
      resizeArray(condition.getConditionValues(),1);
    }

  }
  for(int i = 0;i < loadingConditions.size();i++) {
    Condition& condition = loadingConditions[i];

    if(condition.getType() == "Point-Force-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Line-Force-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Traction-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Surface-Pressure-Loading") {
      resizeArray(condition.getConditionDOFs(),1);
      resizeArray(condition.getConditionValues(),1);
    }
    else if(condition.getType() == "Body-Force-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Point-Moment-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Line-Moment-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Surface-Moment-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Body-Moment-Loading") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Electric-Surface-Charge-Loading") {
      resizeArray(condition.getConditionDOFs(),electricDOF);
      resizeArray(condition.getConditionValues(),electricDOF);
    }
    else if(condition.getType() == "Electric-Body-Charge-Loading") {
      resizeArray(condition.getConditionDOFs(),electricDOF);
      resizeArray(condition.getConditionValues(),electricDOF);
    }
    else if(condition.getType() == "Elastic-Point-Force") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Elastic-Line-Force") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }
    else if(condition.getType() == "Elastic-Surface-Force") {
      resizeArray(condition.getConditionDOFs(),usedDims);
      resizeArray(condition.getConditionValues(),usedDims);
    }

  }

}

/************************************************************************/
/************************************************************************/
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
  if(p != problemData.end()) return (double) p->second;

  // Key was not found -> looking for it in the default values. 
  else {

    p = defaultProblemData.find(name);

    if(p != defaultProblemData.end()) {
      cerr << "Input parameter '" << name << "' hasn't been found - "
          << " default value '" << (double) p->second << "' is used!" << endl;
      return (double) p->second;
    }
    else {
      cerr << "Input parameter '" << name << "' hasn't been found!"
          << " Check input file!" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

}

/************************************************************************/
/************************************************************************/
double InputFileData::getMatValue(int matID,const char* name) {

  using namespace std;

  if(matID >= materials.size()) {
    cerr << "Material set no. " << matID + 1 << " does not exist!"
        << " Check input file!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  map<string,double>::iterator p = materials[matID].find(name);

  if(p != materials[matID].end()) return (double) p->second;
  else {
    cerr << "Material input parameter '" << name << "' haven't been found!"
        << " Check input file!" << endl;
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
void InputFileData::clearArrays(std::ofstream &logFile) {

  using namespace std;

  backGroundMeshInfo.clear();

  materials.clear();
  problemData.clear();
  defaultProblemData.clear();

  graphs.resize(0);
  linePlotGraphData.clear();

  resizeArray(microLengths,0);

  vector<Condition>(0,Condition()).swap(dirichletConditions);
  vector<Condition>(0,Condition()).swap(loadingConditions);

  vector<Condition>(0,Condition()).swap(resultantReactions);

}

