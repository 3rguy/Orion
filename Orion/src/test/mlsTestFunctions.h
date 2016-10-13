/*
 * mlsTestFunctions.h
 *
 *  Created on: Jul 25, 2014
 *      Author: rama
 */

#ifndef MLSTESTFUNCTIONS_H_
#define MLSTESTFUNCTIONS_H_

#define PI 3.14159265

#include <vector>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <string>

#include <time.h>
#include <chrono>
#include <ctime>

#include "commonTypedefs.h"
#include "commonFunctions.h"

#include "Data.h"
#include "defs.h"
#include "Interpolation.h"
#include "DataContainer.h"
#include "Particle.h"
#include "ErrorCalc.h"


// *****************************************************************************
// *****************************************************************************
void mlsTest(InputFileData* InputData,ofstream& logFile);

void mlsTest_one(InputFileData* InputData,ofstream& logFile);

void mlsTest_two(InputFileData* InputData,ofstream& logFile);

void create_domain(vector<Particle>& ptcls,DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile);

dbMatrix read_nodes(DataContainer* MLSData, InputFileData* InputData,
		ofstream& logFile);

dbMatrix create_nodes(DataContainer* MLSData,InputFileData* InputData,
		ofstream& logFile);

dbMatrix create_randomNodes(DataContainer* MLSData,InputFileData* InputData,
		ofstream& logFile);

void create_samplingPoints(vector<Particle>& ptcls,DataContainer* MLSData,
		InputFileData* InputData,ofstream& logFile);

void particlesSelection(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

void setupFunction(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

double calcFunction(dbVector& coords,DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile);

bool isInInfluenceZone(Particle& ptcl_A, Particle& ptcl_B, double& influenceRadius,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

//double calcDistance(dbVector& A, dbVector B);

void findSupports(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

void calcInterpolants(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

void calcApproxPolynomialFunction(vector<Particle>& ptcls,
		vector<Particle>& ptclList,DataContainer* MLSData,InputFileData* InputData,
		ofstream& logFile);

void errorCalculation(vector<Particle>& ptcls, vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

void plotShapeFunctions(vector<Particle>& ptcls, vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile);

void plot1DShapeFunctionsOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile);

void plot1DErrorOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile);


void plot2DShapeFunctionsOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile);

void plot2DErrorOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile);

#endif /* MLSTESTFUNCTIONS_H_ */
