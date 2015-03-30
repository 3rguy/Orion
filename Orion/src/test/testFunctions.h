/*
 * testFunctions.h
 *
 *  Created on: Jul 25, 2014
 *      Author: rama
 */

#ifndef TESTFUNCTIONS_H_
#define TESTFUNCTIONS_H_

#define PI 3.14159265

//#include "commonTypedefs.h"
//#include "commonFunctions.h"
//#include "fortranFunctions.h"

#include <limits>

//#include "Database.h"
#include "Data.h"
//#include "PODCalc.h"
//#include "FEMGeometryX.h"
//#include "ParticleX.h"
//#include "FEMElementX.h"

#include "defs.h"
#include "ROMCalc.h"
#include "Interpolation.h"
#include "testFunctions.h"
#include "mlsTestFunctions.h"

#include <time.h>
#include <chrono>
#include <ctime>

#include <iostream>
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_types.h>

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

#endif /* TESTFUNCTIONS_H_ */
