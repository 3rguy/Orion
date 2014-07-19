/*
 * FEMGeometryExt.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef FEMGEOMETRYEXT_H_
#define FEMGEOMETRYEXT_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "commonFunctions.h"
#include "defs.h"
#include "ElementTemplateX.h"
#include "FEMElementX.h"
#include "InputFileData.h"
#include "LineElementTemplatesX.h"
#include "ParticleX.h"
#include "SurfaceElementTemplatesX.h"
#include "VolumeElementTemplatesX.h"

#include "FEMElementExt.h"
#include "ParticleExt.h"

// Seska lib
#include "FEMGeometry.h"
#include "LineElementTemplates.h"

#include "mpi.h"
#include "petscvec.h"
#include "petscsys.h"

class FEMGeometryExt {

  public:

    FEMGeometryExt(InputFileData* InputData,
		std::map<std::string,double>& modelData,
		std::string& meshFileName,
		std::ofstream& logFile);

    FEMGeometryExt(InputFileData* InputData,
    		std::map<std::string,double>& modelData,
    		std::string& meshFileName,std::string& inputFileName,
    		std::ofstream& logFile);

    ~FEMGeometryExt();


    /*!************************************************************************/

	void readMeshFile(InputFileData* InputData,
			std::map<std::string, double>& modelData, std::string& meshFileName,
			std::ofstream& logFile);

	/*!************************************************************************/

	int getNodesNum() { return particles.size(); };
	std::vector<ParticleExt>& getNodesVec() { return particles; };

	std::vector<FEMElementExt>& getNodesElemsVec() { return nodesElements; };

	std::vector<FEMElementExt>& getSurfaceNodesElemsVec()
	      { return surfaceNodesElems; };

	/*!************************************************************************/

	void combiningSurfaceVolumeElems(InputFileData* InputData,
			std::map<std::string, double>& modelData, std::ofstream& logFile);

	intMatrix decompVolumeToSurfaceElems(intVector& volNodes, int& elemType,
			InputFileData* InputData, std::map<std::string, double>& modelData,
			std::ofstream& logFile);

	intMatrix decompTetraToTrianElems(intVector& volNodes,
			InputFileData* InputData, std::map<std::string, double>& modelData,
			std::ofstream& logFile);

	intMatrix decompHexaToQuadElems(intVector& volNodes,
			InputFileData* InputData, std::map<std::string, double>& modelData,
			std::ofstream& logFile);

	intVector surfaceIDGenerator(intMatrix surfaceNodesList,
			InputFileData* InputData, std::map<std::string, double>& modelData,
			std::ofstream& logFile);

	int storeSurfaceElem(intVector& surfaceNodes, InputFileData* InputData,
			std::map<std::string, double>& modelData, std::ofstream& logFile);

	/*!************************************************************************/
	void FEMGeoDataSetup(InputFileData* InputData,
			std::map<std::string, double>& modelData,
    		std::string& meshFileName, std::string& inputFileName,
			std::ofstream& logFile);

	FEMGeometry* getFEMGeoData(void){return FEMGeoData;};
	InputFileData* getFEMInputData(void){return FEMInputData;};

	/*!************************************************************************/
	int findSurfElemInSurfList(intVector& surfaceNodes,
			InputFileData* InputData, std::map<std::string, double>& modelData,
			std::ofstream& logFile);

	int findPointInGeometry(dbVector& pointCoord, InputFileData* InputData,
			std::map<std::string, double>& modelData, std::ofstream& logFile);

	bool findPointInVolumeElem(dbVector& pointCoord, intVector& surfaceElems,
			intVector& volumeNodes, InputFileData* InputData,
			std::map<std::string, double>& modelData, std::ofstream& logFile);

	int findPtcleIDInPtcleList(int& ptcleID, std::ofstream& logFile);

	int findSurfaceElemIDInList(int& surfaceID, std::ofstream& logFile);

	int findVolumeElemIDInVolList(int& volumeID, std::ofstream& logFile);

	intVector findAdjacentVolElemsOfVolumeElement(int& volID,
			std::ofstream& logFile);

	dbVector caclSurfaceNormal(dbMatrix& surfacePointsCoords,
			InputFileData* InputData, std::ofstream& logFile);

	dbMatrix setupLocalBasis(dbMatrix& surfacePointsCoords,
			InputFileData* InputData, std::map<std::string, double>& modelData,
			std::ofstream& logFile);

	dbVector normaliseVec(dbVector& vec, std::ofstream& logFile);

	void writeMeshFile(InputFileData* InputData, std::ofstream& logFile);

	void printVolumePtclsDetails(InputFileData* InputData, std::ofstream& logFile);

  private:

	FEMGeometry* FEMGeoData;
	InputFileData* FEMInputData;

//	std::vector<ParticleX> particles;
//	std::vector<FEMElementX> nodesElements;
//	std::vector<FEMElementX> surfaceNodesElems;

//	std::vector<ParticleExt> particleExt;
//	std::vector<FEMElementExt> nodesElementExt;
//	std::vector<FEMElementExt> surfaceNodesElemExt;

	std::vector<ParticleExt> particles;
	std::vector<FEMElementExt> nodesElements;
	std::vector<FEMElementExt> surfaceNodesElems;


	std::vector<std::map<int,int> > lineElemTemplateID;
	std::vector<std::map<int,int> > surfaceElemTemplateID;
	std::vector<std::map<int,int> > volumeElemTemplateID;
	std::vector<ElementTemplateX*> lineElemTemplates;
	std::vector<ElementTemplateX*> surfaceElemTemplates;
	std::vector<ElementTemplateX*> volumeElemTemplates;

    std::vector<GaussPointSetX*> volumeGaussPtTemplates;
    std::vector<GaussPointSetX*> surfaceGaussPtTemplates;
    std::vector<GaussPointSetX*> lineGaussPtTemplates;

};

#endif /* FEMGEOMETRYEXT_H_ */
