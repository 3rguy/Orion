// Stores all properties of a loading or Dirichlet condition applied to a
// set of elements or nodes.
//
// Note: The current magnitude of the condition is computed as
// p(t) = lambda(t)*p0 where p0 denotes the reference value and lambda
// the current condition factor. This factor is incremented in time.

#include "Graph.h"

Graph::Graph() :
    ID(0), type("load_deformation"), node(0), dof(0), plotCounter(0),
    abscissaScaling(0), ordinateScaling(0), loadingConditionID(-1),
    dirichletConditionID(-1), dirichletControlID(-1), resultantReactionID(-1) {

  using namespace std;

  pushBackVector(admissibleTypes,(string) "load_deformation");
  pushBackVector(admissibleTypes,(string) "load_twist");
  pushBackVector(admissibleTypes,(string) "stress_strain");
  pushBackVector(admissibleTypes,(string) "devStress_devStrain");
  pushBackVector(admissibleTypes,(string) "traction_displacement");
  pushBackVector(admissibleTypes,(string) "deformation_time");
  pushBackVector(admissibleTypes,(string) "load_time");
  pushBackVector(admissibleTypes,(string) "energy_time");
  pushBackVector(admissibleTypes,(string) "pressure_volume");
  pushBackVector(admissibleTypes,(string) "volume_time");
  pushBackVector(admissibleTypes,(string) "fibreStress_time");
  pushBackVector(admissibleTypes,(string) "activeTension_time");
  pushBackVector(admissibleTypes,(string) "sarcomereLength_time");
  pushBackVector(admissibleTypes,(string) "contractileElementLength_time");

}

Graph::~Graph() {

//  if((*outputFile).is_open())
//    (*outputFile).close();
  //delete[] outputFile;

}

/***********************************************************************/
/***********************************************************************/
// set the axis labels of the graph
void Graph::setLabels(std::map<std::string,double>& calcData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile) {


  using namespace std;

  bool extend = (bool) calcData["restartFileID"];

  // set label of axes
  char* cstr = new char[type.length() + 1];
  strcpy(cstr,type.c_str());

  char* p = strtok(cstr,"_");

  if(p != 0) ordinateLabel = p;
  else {
    logFile << "In Graph::initPtcleGraph axis initialization failed." << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  p = strtok(NULL," ");
  if(p != 0) abscissaLabel = p;
  else {
    logFile << "In Graph::initPtcleGraph axis initialization failed." << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  delete cstr;

#ifdef _postProcDebugMode_
  logFile<<"*****************************************************"<<endl;
  logFile<<"**************** Graph::setLabels *******************"<<endl;
  logFile<<type<<" "<<ID<<": "<<ordinateLabel<<" vs. "
      <<abscissaLabel<<endl;
#endif

}

/***********************************************************************/
/***********************************************************************/
// initialize the graph illustrating a result at particle (open graph-file etc.)
void Graph::initPtcleGraph(std::map<std::string,double>& calcData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile) {


  using namespace std;

  bool extend = (bool) calcData["restartFileID"];

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//  // set label of axes
//  char* cstr = new char[type.length() + 1];
//  strcpy(cstr,type.c_str());
//
//  char* p = strtok(cstr,"_");
//
//  if(p != 0)
//    ordinateLabel = p;
//  else {
//    logFile << "In Graph::initPtcleGraph axis initialization failed." << endl;
//    MPI_Abort(MPI_COMM_WORLD,1);
//  }
//
//  p = strtok(NULL," ");
//  if(p != 0)
//    abscissaLabel = p;
//  else {
//    logFile << "In Graph::initPtcleGraph axis initialization failed." << endl;
//    MPI_Abort(MPI_COMM_WORLD,1);
//  }
//
//  delete cstr;
//
//
//#ifdef _postProcDebugMode_
//  logFile<<"*****************************************************"<<endl;
//  logFile<<"*********** Graph::initPtcleGraph *******************"<<endl;
//  logFile<<type<<" "<<ID<<": "<<ordinateLabel<<" vs. "
//      <<abscissaLabel<<endl;
//#endif

  /*********************************************************************/
  // flag graph type to be active
  if(abscissaLabel == "twist" || abscissaLabel == "strain"
    || abscissaLabel == "stress" || abscissaLabel == "activeTension"
    || abscissaLabel == "sarcomereLength"
    || abscissaLabel == "contractileElementLength_time"
    || abscissaLabel == "energy") {

    string name = abscissaLabel + (string) "Plotting";
    calcData[name] = 1;

#ifdef _postProcDebugMode_
  logFile<<name<<endl;
#endif
  }

  if(ordinateLabel == "twist" || ordinateLabel == "strain"
    || ordinateLabel == "stress" || ordinateLabel == "activeTension"
    || ordinateLabel == "sarcomereLength"
    || ordinateLabel == "contractileElementLength_time"
    || ordinateLabel == "energy") {

    string name = ordinateLabel + (string) "Plotting";
    calcData[name] = 1;

#ifdef _postProcDebugMode_
  logFile<<name<<endl;
#endif
  }

  /*********************************************************************/
  // open graph output file

  if(rank == 0) {

    outputFile = new ofstream;
    string fileName;

    // graph for a specific node
    if(node != 0 && dof != 0)  {
      fileName = type + (string) "-node"
      + convertIntToString((int) node) + (string) "-DOF" + (char) (48 + dof)
      + (string) ".grf";

      // Write the graph headline if the current simulation is no restart.
      if( !extend) {
        ( *outputFile).open(fileName.c_str());
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);

        ( *outputFile) << "# Graf: \"" << type << " at particle " << node
            << " dof " << dof << "\"" << endl;
        ( *outputFile) << "#" << endl;
        ( *outputFile) << "# " << abscissaLabel << ": \" \" " << ordinateLabel
            << ": \" \"" << endl;
      }
      else {
        ( *outputFile).open(fileName.c_str(),ios_base::app);
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);
      }

    }
    // graph for condition applied to parts/entire structure -
    // loading vs. Dirichlet control condition
    else if(loadingConditionID != -1 && dirichletControlID != -1) {
      fileName = type + (string) "-loadID_"
            + convertIntToString((int) loadingConditionID)
            + (string) "-dirichletID_"
            + convertIntToString((int) dirichletControlID)
            + (string) ".grf";

      // Write the graph headline if the current simulation is no restart.
      if( !extend) {
        ( *outputFile).open(fileName.c_str());
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);

        ( *outputFile) << "# Graf: \"" << type << " loadID " << loadingConditionID
            << " dirichletID " << dirichletControlID << "\"" << endl;
        ( *outputFile) << "#" << endl;
        ( *outputFile) << "# " << abscissaLabel << ": \" \" " << ordinateLabel
            << ": \" \"" << endl;
      }
      else {
        ( *outputFile).open(fileName.c_str(),ios_base::app);
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);
      }

    }
    // graph for condition applied to parts/entire structure -
    // resultant loading vs. Dirichlet control condition
    else if(resultantReactionID != -1 && dirichletControlID != -1) {
      fileName = type + (string) "-resultantID_"
            + convertIntToString((int) resultantReactionID)
            + (string) "-dirichletID_"
            + convertIntToString((int) dirichletControlID)
            + (string) ".grf";
      // Write the graph headline if the current simulation is no restart.
      if( !extend) {
        ( *outputFile).open(fileName.c_str());
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);

        ( *outputFile) << "# Graf: \"" << type << " loadID " << loadingConditionID
            << " dirichletID " << dirichletControlID << "\"" << endl;
        ( *outputFile) << "#" << endl;
        ( *outputFile) << "# " << abscissaLabel << ": \" \" " << ordinateLabel
            << ": \" \"" << endl;
      }
      else {
        ( *outputFile).open(fileName.c_str(),ios_base::app);
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);
      }

    }
    // graph for condition applied to parts/entire structure -
    // loading condition
    else if(loadingConditionID != -1) {
      fileName = type + (string) "-loadID_"
            + convertIntToString((int) loadingConditionID)
            + (string) ".grf";

      // Write the graph headline if the current simulation is no restart.
      if( !extend) {
        ( *outputFile).open(fileName.c_str());
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);

        ( *outputFile) << "# Graf: \"" << type << " loadID " << loadingConditionID
             << "\"" << endl;
        ( *outputFile) << "#" << endl;
        ( *outputFile) << "# " << abscissaLabel << ": \" \" " << ordinateLabel
            << ": \" \"" << endl;
      }
      else {
        ( *outputFile).open(fileName.c_str(),ios_base::app);
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);
      }

    }
    // graph for condition applied to parts/entire structure -
    // Dirichlet control condition
    else if(dirichletControlID != -1) {
      fileName = type + (string) "-dirichletID_"
            + convertIntToString((int) dirichletControlID)
            + (string) ".grf";

      // Write the graph headline if the current simulation is no restart.
      if( !extend) {
        ( *outputFile).open(fileName.c_str());
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);

        ( *outputFile) << "# Graf: \"" << type << " dirichletID " << dirichletControlID << "\"" << endl;
        ( *outputFile) << "#" << endl;
        ( *outputFile) << "# " << abscissaLabel << ": \" \" " << ordinateLabel
            << ": \" \"" << endl;
      }
      else {
        ( *outputFile).open(fileName.c_str(),ios_base::app);
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);
      }

    }
    // some other graph
    else {
      fileName = type + (string) ".grf";

      // Write the graph headline if the current simulation is no restart.
      if( !extend) {
        ( *outputFile).open(fileName.c_str());
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);

        ( *outputFile) << "# Graf: \"" << type << "\"" << endl;
        ( *outputFile) << "#" << endl;
        ( *outputFile) << "# " << abscissaLabel << ": \" \" " << ordinateLabel
            << ": \" \"" << endl;
      }
      else {
        ( *outputFile).open(fileName.c_str(),ios_base::app);
        ( *outputFile).precision(12);
        ( *outputFile).setf(ios_base::scientific,ios_base::floatfield);
      }

    }

  }

}

/***********************************************************************/
/***********************************************************************/
// check whether specified graph type is admissible
bool Graph::isType(std::string& t) {

  using namespace std;

  bool flag=false;

  if (find(admissibleTypes.begin(),admissibleTypes.end(),t) != admissibleTypes.end())
    flag=true;

  return flag;
}

/***********************************************************************/
/***********************************************************************/
// add a new value pair to the graph
void Graph::addValuePair(dbVector& values,std::ofstream& logFile) {

  using namespace std;
  pushBackVector(valuePairs,values);

}

/***********************************************************************/
/***********************************************************************/
// write a chunk of value pairs to file
void Graph::writeValuePairs(std::ofstream& logFile) {

  using namespace std;

#ifdef _postProcDebugMode_
  logFile<<"write graph("<<ID<<") values: "<<endl;
#endif

  // loop over all values pairs which are not written to file yet
  for(int i = plotCounter;i < valuePairs.size();i++) {

    // x-coordinate
    if(abscissaScaling == 0)

    (*outputFile) << fabs(valuePairs[i][0]) << " ";

    else

    (*outputFile) << valuePairs[i][0] * abscissaScaling << " ";

    // y-coordinate
    if(ordinateScaling == 0)

    (*outputFile) << fabs(valuePairs[i][1]) << endl;

    else (*outputFile) << valuePairs[i][1] * ordinateScaling<< endl;

    plotCounter++;

#ifdef _postProcDebugMode_
    logFile<<valuePairs[i][0]<<", "<<valuePairs[i][1]<<endl;
#endif

  }

}

/***********************************************************************/
/***********************************************************************/
// check whether graph is plotting loading
bool Graph::isLoadingPlotting() {

  if(abscissaLabel == "load" || ordinateLabel == "load")
    return true;
  else
    return false;
}

/***********************************************************************/
/***********************************************************************/
// check whether graph is plotting Dirichlet control conditions
bool Graph::isDirichletControlConditionPlotting() {

  if(abscissaLabel == "volume" || ordinateLabel == "volume")
    return true;
  else
    return false;
}

/***********************************************************************/
/***********************************************************************/
// assign to the graph a loading or internal reaction
void Graph::setLoadingCondition(Condition& condition) {

  loadingCondition = &condition;
}

/***********************************************************************/
/***********************************************************************/
// assign to the graph a Dirichlet control condition
void Graph::setDirichletControlCondition(Condition& condition) {

  dirichletControlCondition = &condition;
}




