/// Stores all properties of a graph illustrating a result at particle,
/// along line of particles etc.
///
/// admissible types: load-deformation - note: input variable 'loadingConditionID' must be specified
/// admissible types: load_twist - note: input variable 'loadingConditionID' must be specified
/// admissible types: stress_strain
/// admissible types: devStress_devStrain
/// admissible types: traction_displacement
/// admissible types: deformation_time
/// admissible types: load_time - note: input variable 'loadingConditionID' must be specified
/// admissible types: energy_time
/// admissible types: pressure_volume - note: input variable 'loadingConditionID' and 'dirichletControlID' must be specified
/// admissible types: volume_time - note: input variable 'dirichletControlID' must be specified
/// admissible types: fibreStress_time
/// admissible types: activeTension_time
/// admissible types: sarcomereLength_time
/// admissible types: contractileElementLength_time


#ifndef Graph_h_
#define Graph_h_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "Condition.h"

class Graph {

  public:

    Graph();
    ~Graph();

    /// initialize the graph plotting the result at a particle (open graph-file etc.)
    void initPtcleGraph(std::map<std::string,double>& calcData,
                        std::map<std::string,double>& modelData,
                        std::ofstream& logFile);

    /// check whether specified graph type is admissible
    bool isType(std::string& t);

    /// add a new value pair to the graph
    void addValuePair(dbVector& valuePair,std::ofstream& logFile);

    /// write a chunk of value pairs to file
    void writeValuePairs(std::ofstream& logFile);

    /// check whether graph is plotting loading
    bool isLoadingPlotting();

    /// assign to the graph a loading or internal reaction
    void setLoadingCondition(Condition& condition);

    // assign to the graph a Dirichlet control condition
    void setDirichletControlCondition(Condition& condition);

    bool setParam(std::string name,std::string value) {
      using namespace std;
      bool info = false;
      if(name == "type") {
        type = value;
        info = true;
      }

      return info;
    };

    bool setParam(std::string name,double value) {
      using namespace std;
      bool info = false;
      if(name == "ID") {
        ID = value;
        info = true;
      }
      else if(name == "node") {
        node = value;
        info = true;
      }
      else if(name == "DOF") {
        dof = value;
        info = true;
      }
      else if(name == "loadingConditionID") {
        loadingConditionID = value-1;
        info = true;
      }
      else if(name == "dirichletConditionID") {
        dirichletConditionID = value-1;
        info = true;
      }
      else if(name == "dirichletControlID") {
        dirichletControlID = value-1;
        info = true;
      }
      else if(name == "resultantReactionID") {
        resultantReactionID = value-1;
        info = true;
      }
      else if(name == "abscissaScaling") {
        abscissaScaling = value;
        info = true;
      }
      else if(name == "ordinateScaling") {
        ordinateScaling = value;
        info = true;
      }

      return info;
    };



    int& getID() { return ID; };

    std::string& getType() { return type; };
    int& getNode() { return node; };
    int& getDOF() { return dof; };
    int& getLoadingConditionID() { return loadingConditionID; };
    int& getDirichletConditionID() { return dirichletConditionID; };
    int& getResultantReactionID() { return resultantReactionID; };
    int& getDirichletControlID() { return dirichletControlID; };

    std::string& getAbscissaLabel() { return abscissaLabel; };
    std::string& getOrdinateLabel() { return ordinateLabel; };
    double& getAbscissaScaling() { return abscissaScaling; };
    double& getOrdinateScaling() { return ordinateScaling; };

    Condition* getLoadingCondition() { return loadingCondition; };
    Condition* getDirichletControlCondition() { return dirichletControlCondition; };

  private:

    int ID,node,dof,loadingConditionID,dirichletConditionID,resultantReactionID,
        dirichletControlID;
    std::string type,abscissaLabel,ordinateLabel;
    std::vector<std::string> admissibleTypes;

    int plotCounter;
    double abscissaScaling,ordinateScaling;
    dbMatrix valuePairs;


    std::ofstream* outputFile;

    Condition* loadingCondition;
    Condition* dirichletControlCondition;

};

#endif
