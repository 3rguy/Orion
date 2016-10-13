/*
 * testFunctions.h
 *
 *  Created on: Jul 25, 2014
 *      Author: rama
 */

#ifndef TESTFUNCTIONS_H_
#define TESTFUNCTIONS_H_

#include <iostream>
#include <limits>
#include <stdlib.h>

#include <time.h>
#include <chrono>
#include <ctime>

#include "defs.h"
#include "Data.h"
#include "Interpolation.h"
#include "GridNodesTest.h"
#include "PODCalc.h"
#include "PODICalc.h"
#include "mlsTestFunctions.h"

#define PI 3.14159265

// *****************************************************************************
// *****************************************************************************
void testFunctions(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void findPointInPolygonTest(InputFileData* InputData, ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void PODCalcTest_readFile(dbMatrix& fullMatrix, InputFileData* InputData,
		ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void PODCalcTest(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void readResultFileTest(InputFileData* InputData,ofstream& logFile);
void readResultFileTest_M1(InputFileData* InputData,ofstream& logFile);
void readResultFileTest_M2(InputFileData* InputData,ofstream& logFile);
//void readResultFileTest_M3(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void readResultFileTest_writeResultToFile(string& resultName,dbVector& resultStepVector,
		dbMatrix& resultMatrix,int& numDofs,InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void supportingPtclsTest(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void calcVolume(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void calcVolumeTwo(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
void printSnapshotPOMs(InputFileData* InputData,ofstream& logFile);

// *****************************************************************************
// *****************************************************************************
//void pcl_test(InputFileData* InputData,ofstream& logFile);

dbMatrix createRandomMatrix(int row,int col,ofstream& logFile);

#endif /* TESTFUNCTIONS_H_ */
