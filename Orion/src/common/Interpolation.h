/*
 * Interpolation.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "float.h"
#include <map>
#include "mpi.h"
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "InputFileData.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"

using namespace std;

// Calculate the Moving Least Square Approximants.
class Interpolation {
public:

	Interpolation(dbVector& iPoint,dbMatrix& coords,dbVector& radiusVec,
			dbVector& interpolants,InputFileData* InputData, ofstream& logFile);
	~Interpolation(){};

	/*!************************************************************************/
	//! Calculate the Moving Least Square interpolants
	dbVector MLSCalc(dbVector& iPoint,dbMatrix& coords,dbVector& radiusVec,
			InputFileData* InputData,ofstream& logFile);

	/*!************************************************************************/
	//! Function to setup the Pascal polynomial's basis
	dbVector PascalBasisCalc(dbVector& coord,int& orderPoly,ofstream& logFile);

	/*!************************************************************************/
	//! Calculate the Cubic Spline Weight function
	double cubicSplineWgtCalc(dbVector& coordOne,dbVector& coordTwo,
			dbVector& radii,ofstream& logFile);

	/*!************************************************************************/
	//! Find the neighbours of a particular point
	intVector findNeighbours(dbVector& iPoint,
				dbMatrix& coords,dbVector& radiusVec,ofstream& logFile);


	/*!************************************************************************/
	/*!************************************************************************/
	//! Debugging Functions
	void MLSUnitTest(ofstream& logFile);
	void MLSUnitTest_meshMultiDim(double& length,double& nodal_dist,
			int& ndim,dbMatrix& x,int& numCoords, intVector SBM_dim,
			dbVector SBM_coeff, ofstream& logFile);

};


#endif /* INTERPOLATION_H_ */
