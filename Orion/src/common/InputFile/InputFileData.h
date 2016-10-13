/// Read and store all input data.

#ifndef InputFileData_h_
#define InputFileData_h_

#include "commonFunctions.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "AnisotropyCondition.h"
#include "commonTypedefs.h"
#include "Condition.h"
#include "defs.h"
#include "mpi.h"
#include "Graph.h"
#include "ResultantReactionCondition.h"

class InputFileData {

 public:
  InputFileData(std::ofstream& logFile);
  InputFileData(std::string& filename,std::ofstream& logFile);

  ~InputFileData() {};
  
  void adjustConditionAllocation(InputFileData* InputData,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile);

  void adjustSimulationOutputFrequency(
        double adjustmentFactor,std::map<std::string,double>& calcData,
        std::map<std::string,double>& modelData,std::ofstream& logFile);

  /**********************************************************************/
  /// various simulation input data

  void setValue(const char* name,double value);
  bool setValue(std::string name,double value);
  double getValue(const char* name);
  std::map<std::string,double>& getProblemData() { return problemData; };

  /**********************************************************************/
  /// material data input
  bool setMatValue(int matID,std::string name,double value);
  double getMatValue(int matID,const char* name);
  std::vector<std::map<std::string,double> >& getMaterials() { return materials; };
  int getNumberOfMats() { return materials.size(); };

  /**********************************************************************/
  /// graphs input

  std::vector<Graph>& getGraphs() { return graphs; };
  std::vector<std::map<std::string,double> >& getLinePlotGraphData() 
    { return linePlotGraphData; };

  std::map<std::string,double>& getBackGroundMeshInfo() {
    return backGroundMeshInfo; };

  /**********************************************************************/
  /// Return applied loading and Dirichlet conditions

  std::vector<Condition>& getLoadingConditions() { return loadingConditions; };
  std::vector<Condition>& getDirichletConditions() { return dirichletConditions; };
  std::map<std::string,bool>& getConditionSet() { return conditionSet; };

  /// particle-displacement- or cavity-volume-control conditions
  std::vector<Condition>& getDirichletControlConditions() {
    return dirichletControlConditions;
  };

  /// surfaces where the resultant force and torque needs to be computed
  std::vector<ResultantReactionCondition>& getResultantReactions() { return resultantReactions; };

  /// surface anisotropy information
  std::map<std::string,std::vector<AnisotropyCondition> >& getAnisotropyConditions() {
    return anisotropyConditions; };


  /**********************************************************************/
  /**********************************************************************/

 private:

  /// Read and store the input file.
  void readInputFile(std::string& filename,std::ofstream& logFile);

  /// Set default simulation parameters.
  void setDefaultSimulationParameters(std::ofstream& logFile);

  /// Set default material parameters.
  void setDefaultMaterialParameters(std::ofstream& logFile);

  void clearArrays(std::ofstream& logFile);

  std::map<std::string,double> backGroundMeshInfo;
  std::vector<std::map<std::string,double> > materials;
  std::map<std::string,std::vector<AnisotropyCondition> > anisotropyConditions;

  /**********************************************************************/
  /// various simulation input data

  std::map<std::string,double> problemData;
  std::map<std::string,double> defaultProblemData;

  std::vector<Graph> graphs;
  std::vector<std::map<std::string,double> > linePlotGraphData;

  std::vector<std::vector<double> > microLengths;


  /**********************************************************************/
  /// Loading and Dirichlet conditions

  std::vector<Condition> loadingConditions,dirichletConditions,
        dirichletControlConditions;
  std::map<std::string,bool> conditionSet;

  /// surfaces where the resultant force and torque needs to be computed
  std::vector<ResultantReactionCondition> resultantReactions;

};

#endif 
