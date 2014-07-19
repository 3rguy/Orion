/*
 * debugFunctions.h
 *
 *  Created on: Jul 17, 2014
 *      Author: rama
 */

#ifndef DEBUGFUNCTIONS_H_
#define DEBUGFUNCTIONS_H_

#include "float.h"
#include <fstream>
#include <iostream>
#include "mpi.h"
#include <vector>

#include "AsymShapeFunc.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "EFGShapeFunc.h"
#include "fortranFunctions.h"
#include "MaxEntShapeFunc.h"
#include "OrthoShapeFunc.h"
#include "Particle.h"
#include "RKPMShapeFunc.h"
#include "ShepardShapeFunc.h"
#include "WindowFuncAsym.h"
#include "WindowFunctionSets.h"

// Test the Sansour algorithm to update the rotational degrees of
// freedom.
void testUpdate(InputFileData* InputData,
		std::ofstream& logFile);
void calcRotTens(dbVector& rotations,dbMatrix& rotationTens,
		 std::ofstream& logFile);

// Test routine for continuum mechanics' stuff.
void testContinuumMechanics(std::ofstream& logFile);

// Plot a cubic spline and its first and second order derivations.
void testSpline(std::ofstream& logFile);

// Plot the RKPM or EFG shape functions and its first and second order
// derivations(test particle 0,62,124!).
void testShapes(InputFileData* InputData,std::vector<Particle>& ptcls,
		std::map<std::string,double>& modelData,
		std::ofstream& logFile,PetscViewer& viewerSEQ);


// Test the curve fitting property of a MLS shape functions and its
// first and second order derivations.
void testMLS(InputFileData* InputData,std::vector<Particle>& ptcls,
	     std::map<std::string,double>& modelData,
	     std::ofstream& logFile,PetscViewer& viewerSEQ);

// Test the curve fitting property of a MaxEnt shape functions and its
// first within the domain and on its boundary.
void testMaxEnt(InputFileData* InputData,std::vector<Particle>& ptcls,
		std::map<std::string,double>& modelData,
		std::ofstream& logFile,PetscViewer& viewerSEQ);


// test the modified boundary collocation method
void testBoundCollocation(InputFileData* InputData,
			  std::vector<Particle>& ptcls,
			  intMatrix& boundaryPtcls,
			  intVector& newIdx,
			  intVector& newDOFID,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile,
			  PetscViewer& viewerMPI,
			  PetscViewer& viewerSEQ);

// test Carlo's 'dexpo' fortran routine
void testDexpo(dbMatrix& IN,dbMatrix& OUT,bool seskaStorage,std::ofstream& logFile);

#endif /* DEBUGFUNCTIONS_H_ */
