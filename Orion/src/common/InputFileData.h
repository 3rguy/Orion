/*
 * InputFileData.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef INPUTFILEDATA_H_
#define INPUTFILEDATA_H_

#include "commonFunctions.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "commonTypedefs.h"
#include "defs.h"
#include "mpi.h"

class InputFileData {

 public:
  InputFileData(std::string& filename,std::ofstream& logFile);
  InputFileData(std::ofstream& logFile);

  ~InputFileData() {};

  void setValue(const char* name,double value);
  bool setValue(std::string name,double value);
  double getValue(const char* name);
  std::map<std::string,double>& getProblemData() { return problemData; };

  bool setMatValue(int matID,std::string name,double value);
  double getMatValue(int matID,const char* name);
  std::vector<std::map<std::string,double> >& getMaterials() {return materials; };
  std::vector<std::map<std::string,double> >& getMaterialDataSets() {return materialDataSets; };
  std::vector<std::map<std::string,double> >& getMaterialDataMax() {return materialDataMax; };
  std::vector<std::map<std::string,double> >& getMaterialDataMin() {return materialDataMin; };
  std::vector<const char*>& getParamNames(){return paramNames;};
  std::vector<const char*>& getSystParamNames(){return systParamNames;};
  std::vector<const char*>& getParamMaxNames(){return paramMaxNames;};
  std::vector<const char*>& getParamMinNames(){return paramMinNames;};
  std::vector<const char*>& getSystParamMaxNames(){return systParamMaxNames;};
  std::vector<const char*>& getSystParamMinNames(){return systParamMinNames;};
  std::vector<const char*>& getParamMaxNamesCopy(){return paramMaxNamesCopy;};
  std::vector<const char*>& getParamMinNamesCopy(){return paramMinNamesCopy;};
  std::vector<const char*>& getSystParamMaxNamesCopy(){return systParamMaxNamesCopy;};
  std::vector<const char*>& getSystParamMinNamesCopy(){return systParamMinNamesCopy;};
  std::vector<std::string>& getMaxName(){return maxName;};
  std::vector<std::string>& getMinName(){return minName;};
  std::vector<std::string>& getSystMaxName(){return systMaxName;};
  std::vector<std::string>& getSystMinName(){return systMinName;};
  dbMatrix& getMaterialParameters() {return materialParameters; };
  std::vector<std::vector<std::string> >& getMaterialParameterNames() {return materialParameterNames;};
  int getNumberOfMats() { return materials.size(); };

  void constructMaterialDataVectors(int& matID);
  void constructActiveMaterialDataVectors(int& matID);
  void constructSystemicMaterialDataVectors();
  void resizeMaterialDataVectors();

  intMatrix& getGraphs() { return graphs; };
  std::vector<std::map<std::string,double> >& getLinePlotGraphData()
    { return linePlotGraphData; };

  std::map<std::string,double>& getBackGroundMeshInfo() {
    return backGroundMeshInfo; };


  // Retrun all applied boundary constraints
  dbMatrix& getPointDispBoundConds() { return pointDispBoundConds; };
  dbMatrix& getLineDispBoundConds() { return lineDispBoundConds; };
  dbMatrix& getSurfaceDispBoundConds() { return surfaceDispBoundConds; };

  dbMatrix& getPointRotBoundConds() { return pointRotBoundConds; };
  dbMatrix& getLineRotBoundConds() { return lineRotBoundConds; };
  dbMatrix& getSurfaceRotBoundConds() { return surfaceRotBoundConds; };

  dbMatrix& getPointElectricBoundConds() { return pointElectricBoundConds; };
  dbMatrix& getLineElectricBoundConds() { return lineElectricBoundConds; };
  dbMatrix& getSurfaceElectricBoundConds() { return surfaceElectricBoundConds; };
  dbMatrix& getLinearSurfaceElectricBoundConds() { return linearSurfaceElectricBoundConds; };

  dbMatrix& getPointDepolarisationBoundConds() { return pointDepolarisationBoundConds; };
  dbMatrix& getLineDepolarisationBoundConds() { return lineDepolarisationBoundConds; };
  dbMatrix& getSurfaceDepolarisationBoundConds() { return surfaceDepolarisationBoundConds; };

  dbMatrix& getPointMicroBoundConds() { return pointMicroBoundConds; };
  dbMatrix& getLineMicroBoundConds() { return lineMicroBoundConds; };
  dbMatrix& getSurfaceMicroBoundConds() { return surfaceMicroBoundConds; };

  dbMatrix& getPointStressBoundConds() { return pointStressBoundConds; };
  dbMatrix& getLineStressBoundConds() { return lineStressBoundConds; };
  dbMatrix& getSurfaceStressBoundConds() { return surfaceStressBoundConds; };

  // Return all applied loads.
  dbMatrix& getPointForceLoads() { return pointForceLoads; };
  dbMatrix& getLineForceLoads() { return lineForceLoads; };
  dbMatrix& getTractionLoads() { return tractionLoads; };
  dbMatrix& getSurfacePressureLoads() { return surfacePressureLoads; };
  dbMatrix& getBodyForceLoads() { return bodyForceLoads; };

  dbMatrix& getPointMomentLoads() { return pointMomentLoads; };
  dbMatrix& getLineMomentLoads() { return lineMomentLoads; };
  dbMatrix& getSurfaceMomentLoads() { return surfaceMomentLoads; };
  dbMatrix& getBodyMomentLoads() { return bodyMomentLoads; };

  dbMatrix& getSurfaceElectricChargeLoads()
    { return surfaceElectricChargeLoads; };
  dbMatrix& getBodyElectricChargeLoads()
    { return bodyElectricChargeLoads; };

  // surfaces where the resultant force and torque needs to be computed
  dbMatrix& getResultantForceOnSurfaces()
    { return resultantForceOnSurfaces; };
  dbMatrix& getResultantTorqueOnSurfaces()
    { return resultantTorqueOnSurfaces; };


 private:

  // Read and store the input file.
  void readInputFile(std::string& filename,std::ofstream& logFile);


  void clearArrays(std::ofstream& logFile);

  // Set default values.
  void setDefaultValues(std::ofstream& logFile);
  std::map<std::string,double> backGroundMeshInfo;
  std::vector<std::map<std::string,double> > materials;
  std::vector<std::map<std::string,double> > materialDataSets;
  std::vector<std::map<std::string,double> > materialDataMax;
  std::vector<std::map<std::string,double> > materialDataMin;
  std::vector<const char*> paramNames;
  std::vector<const char*> systParamNames;
  std::vector<const char*> paramMaxNames;
  std::vector<const char*> paramMinNames;
  std::vector<const char*> systParamMaxNames;
  std::vector<const char*> systParamMinNames;
  std::vector<const char*> paramMaxNamesCopy;
  std::vector<const char*> paramMinNamesCopy;
  std::vector<const char*> systParamMaxNamesCopy;
  std::vector<const char*> systParamMinNamesCopy;
  std::vector<std::string> maxName;
  std::vector<std::string> minName;
  std::vector<std::string> systMaxName;
  std::vector<std::string> systMinName;
  dbMatrix materialParameters;
  std::vector<std::vector<std::string> > materialParameterNames;

  std::map<std::string,double> problemData;
  std::map<std::string,double> defaultProblemData;

  intMatrix graphs;
  std::vector<std::map<std::string,double> > linePlotGraphData;

  std::vector<std::vector<double> > microLengths;

  bool pointDispNormalSet,pointForceNormalSet,lineDispNormalSet,
    lineForceNormalSet,pointRotNormalSet,pointMomentNormalSet,
    lineRotNormalSet,lineMomentNormalSet;
  bool pointElectricNormalSet,lineElectricNormalSet;
  bool pointMicroNormalSet,lineMicroNormalSet;
  bool pointStressNormalSet,lineStressNormalSet;
  bool pointDepolarisationNormalSet,lineDepolarisationNormalSet;

  // Displacement boundary conditions
  dbMatrix pointDispBoundConds;
  dbMatrix lineDispBoundConds;
  dbMatrix surfaceDispBoundConds;

  // Rotation boundary conditions
  dbMatrix pointRotBoundConds;
  dbMatrix lineRotBoundConds;
  dbMatrix surfaceRotBoundConds;

  // Electric boundary conditions
  dbMatrix pointElectricBoundConds;
  dbMatrix lineElectricBoundConds;
  dbMatrix surfaceElectricBoundConds;
  dbMatrix linearSurfaceElectricBoundConds;

  // Depolarisation boundary conditions
  dbMatrix pointDepolarisationBoundConds;
  dbMatrix lineDepolarisationBoundConds;
  dbMatrix surfaceDepolarisationBoundConds;

  // Micro boundary conditions
  dbMatrix pointMicroBoundConds;
  dbMatrix lineMicroBoundConds;
  dbMatrix surfaceMicroBoundConds;

  // stress boundary conditions
  dbMatrix pointStressBoundConds;
  dbMatrix lineStressBoundConds;
  dbMatrix surfaceStressBoundConds;

  // load conditions
  dbMatrix pointForceLoads;
  dbMatrix lineForceLoads;
  dbMatrix tractionLoads;
  dbMatrix surfacePressureLoads;
  dbMatrix bodyForceLoads;

  dbMatrix pointMomentLoads;
  dbMatrix lineMomentLoads;
  dbMatrix surfaceMomentLoads;
  dbMatrix bodyMomentLoads;

  dbMatrix surfaceElectricChargeLoads;
  dbMatrix bodyElectricChargeLoads;

  dbMatrix surfaceElems;

  // surfaces where the resultant force and torque needs to be computed
  dbMatrix resultantForceOnSurfaces;
  dbMatrix resultantTorqueOnSurfaces;

};

#endif /* INPUTFILEDATA_H_ */
