#ifndef CardiacFibreGenerator_h_
#define CardiacFibreGenerator_h_

#include "float.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include <functional>   // std::plus
#include <algorithm>	// std::transform

#include "petsc.h"

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "Particle.h"

#include "FEMGeometry.h"
#include "FEMElement.h"

using namespace std;

/*
 Paper: Jonathan Wong and Ellen Kuhl,"Generating fibre orientation maps in "
 "human heart models using Poisson Interpolation", Computer Methods in
 Biomechanics and Biomedical Engineering, 2012, 1-10.

 Objective: To compute the anisotropy directions of any surface defined using
 the epi-cardial and endo-cardial helix angle.

 Note : The "surfaceNormal" dbMatrix in the Particle class is defined by a
 set of vectors (in a row-wise manner) as follows:
 [0] n, averaged nodal normal
 [1] c, circumferential direction
 [2] n_cz, outward pointing normal
 [3] s, sheet normal
 [4] p, fibre direction projection
 [5] f, fibre direction
 [6] m, orthogonal fibre direction to f and s
 */

class CardiacFibreGenerator {

public:

	CardiacFibreGenerator(InputFileData* InputData, std::ofstream& logFile);
	~CardiacFibreGenerator() {};

	dbVector normaliseVec(dbVector& vec, std::ofstream& logFile);

	void readMeshFile(vector<Particle>& ptcls, vector<FEMElement>& surfaceElems,
			dbVector& surfaceAngles, InputFileData* InputData,
			std::ofstream& logFile);

	void generateFEMresFile(InputFileData* InputData, std::ofstream& logFile);

	void calcAllSurfaceNormals(vector<FEMElement>& sElems,
			InputFileData* InputData, std::ofstream& logFile);

	dbVector calcSurfaceNormal(dbMatrix& surfacePointsCoords,
			InputFileData* InputData, std::ofstream& logFile);

	void calcNodalNormal(vector<Particle>& ptcls, vector<FEMElement>& sElems,
			intVector& sNodes, InputFileData* InputData,
			std::ofstream& logFile);

	intVector compileSurfaceNodes(vector<Particle>& ptcls,
			vector<FEMElement>& elemVec, InputFileData* InputData,
			std::ofstream& logFile);

	void calcAveragedNodalNormals(vector<Particle>& ptcls,
			vector<FEMElement>& elemVec, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void calcNodalCircumDirections(vector<Particle>& ptcls, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void calcOutwardNormal(vector<Particle>& ptcls, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void calcSheetNormal(vector<Particle>& ptcls, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void calcFibreDirectProjection(vector<Particle>& ptcls, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void writingResultToFile(vector<Particle>& ptcls, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void calcFibreDirection(vector<Particle>& ptcls, intVector& sNodes,
			InputFileData* InputData, std::ofstream& logFile);

	void calcOrthogonalFibreDirection(vector<Particle>& ptcls,
			intVector& sNodes, InputFileData* InputData,
			std::ofstream& logFile);

	void saveResultsToFile(vector<Particle>& ptcls, InputFileData* InputData,
			std::ofstream& logFile);

	void saveResultsToFile_res_format(std::ofstream& femResFile,
			std::string& resultName, dbMatrix& resultMatrix, double& stepValue,
			InputFileData* InputData, ofstream& logFile);

	void saveResultsToFile_flexible(vector<Particle>& ptcls,
			InputFileData* InputData, std::ofstream& logFile);

	void saveResultsToFile_res_format_flexible(std::string& fileName,
			std::vector<std::string>& resultName,
			std::vector<dbMatrix>& resultMatrix, dbVector& stepValue,
			InputFileData* InputData, ofstream& logFile);

	void saveResultsToFile_res_format_flexible_resultTypes(
			std::string& fileName, std::vector<std::string>& vectorResultName,
			std::vector<dbMatrix>& vectorResultMatrix,
			std::vector<std::string>& scalarResultName,
			std::vector<dbVector>& scalarResultVector, dbVector& stepValue,
			InputFileData* InputData, ofstream& logFile);

	void eliminateSpecificSurface(intVector& surfElems, int surfID,
			vector<FEMElement>& elemVec, InputFileData* InputData,
			std::ofstream& logFile);

	void printVector(intVector& A, const char* msg, std::ofstream& logFile);
	void printVector(dbVector& A, const char* msg, std::ofstream& logFile);

	dbVector z; // longitudinal vector
	double pi;

};

#endif
