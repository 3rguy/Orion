/*
 * GridNodes.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef GRIDNODES_H_
#define GRIDNODES_H_


#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "defs.h"
#include "Data.h"
#include "Node.h"
#include "InputFileData.h"
#include "Interpolation.h"

#include "FEMGeometryExt.h"
#include "ParticleExt.h"
#include "MeshlessApproximation.h"

class GridNodes {
public:
//	GridNodes(std::vector<Data>& iDataList, InputFileData* InputData,
//			ofstream& logFile);

	GridNodes(Data& iData, InputFileData* InputData, ofstream& logFile);

	GridNodes(InputFileData* InputData, ofstream& logFile);

	~GridNodes() {}	;

	void setNumDimensions(int val){nDims = val;}
	int getNumDimensions(void){return nDims;}

	void addNode(Node myNode);
	Node& getNodeId(int val){return nodalList[findNode(val)];}

	Node& getNodeEntry(int val){return nodalList[val];}
	int findNode(int val);

	void setGridNodes(InputFileData* InputData, ofstream& logFile);

	void generateGridNodes(InputFileData* InputData, ofstream& logFile);

	void importGridNodes(InputFileData* InputData, ofstream& logFile);

	std::vector<Node>& getNodalList(){return nodalList;};

//	void calcInterpolants(Data& iData, InputFileData* InputData,
//			ofstream& logFile);

	void findSuppParticles(dbVector& coord, std::vector<ParticleExt>& ptclList,
			intVector& sPtcls, ofstream& logFile);

//	void interpolateDisp(Data& iData, InputFileData* InputData,
//			ofstream& logFile);

	dbMatrix calcDispOnGridPoint(Data& iData, InputFileData* InputData,
			ofstream& logFile);

	dbMatrix interpolateResultOnGridPoint(Data& iData,
			InputFileData* InputData,ofstream& logFile);

	void interpolantSetup(Data& iData,InputFileData* InputData,
			ofstream& logFile);

	void setSupportingParticles(Data& iData, InputFileData* InputData,
			ofstream& logFile);

	void setSupportingParticles_two(FEMGeometryExt* FEMDataExt,
			std::vector<Node>& nodesVec, InputFileData* InputData,
			ofstream& logFile);

	void setSupportingParticles_(FEMGeometryExt* FEMDataExt,
			std::vector<Node>& nodesVec, InputFileData* InputData,
			ofstream& logFile);

	void findSupportingParticles(std::map<std::string, oPType>& modelData,
			InputFileData* InputData, ofstream& logFile);

	void filterSupportingParticles(std::vector<Node>& nodesVec,
			vector<Particle>& ptclList, InputFileData* InputData,
			ofstream& logFile);

	void setInterpolantsOnNodes(FEMGeometryExt* FEMData, InputFileData* InputData,
			ofstream& logFile);

	dbMatrix setInterpolantsOnParticles(Data& iData, InputFileData* InputData,
			ofstream& logFile);

	void calcMatrixFieldMLSApproximants(dbVector& iPoint, dbMatrix& coords,
			dbVector& radiusVec, dbVector& interpolants,
			InputFileData* InputData, ofstream& logFile);

	void interpolateNodalResult(Data& iData, InputFileData* InputData,
			ofstream& logFile);

//	dbMatrix calcResultOnGridPoint(Data& iData, InputFileData* InputData,
//			ofstream& logFile);

//	dbMatrix initCalcResultOnParticles(Data& iData, InputFileData* InputData,
//				ofstream& logFile);

	void initCalcResultOnParticles(Data& iData, InputFileData* InputData,
					ofstream& logFile);

	dbMatrix calcResultOnParticles(Data& iData, dbMatrix& interpolantsList,
			InputFileData* InputData, ofstream& logFile);

	void setSupportingNodes(Data& iData, InputFileData* InputData,
			ofstream& logFile);

	void interpolateParticleResult(Data& iData, dbMatrix interpolantsList,
			InputFileData* InputData, ofstream& logFile);

	void interpolateParticleResult_(Data& iData, dbMatrix interpolantsList,
				InputFileData* InputData, ofstream& logFile);

	void setResultOnNodes(Data& iData, dbMatrix dispMatrix,
			InputFileData* InputData, ofstream& logFile);

	void resetNodes();
	void resetNodesStepDOFMat();

	dbMatrix assembleNodalResultMatrix(Data& iData, InputFileData* InputData,
			ofstream& logFile);

	dbMatrix assemblePtclResultMatrix(Data& iData, InputFileData* InputData,
				ofstream& logFile);

	void saveNodesToFile(ofstream& logFile);
	void saveSelectedNodesToFile(intVector& nodesVec,ofstream& logFile);

	FEMGeometryExt* getGridGeometry(ofstream& logFile);
	void setGridGeometry(FEMGeometryExt* p){gridGeometry=p;};
	void delVolumeElement(intVector& volElem, InputFileData* InputData,
			ofstream& logFile);

	void clearNodalList();

private:
	int nDims;
	std::vector<Node> nodalList;

	FEMGeometryExt* gridGeometry;
	MeshlessApproximation* MeshlessData;
};

#endif /* GRIDNODES_H_ */
