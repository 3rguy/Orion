/*
 * PETScFunctions.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef PETSCFUNCTIONS_H_
#define PETSCFUNCTIONS_H_


#include "float.h"
#include <fstream>
#include <iostream>
#include <map>
#include "mpi.h"
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "defs.h"
#include "fortranFunctions.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"


// Create a dense sequential PETSc coefficient matrix.
void createDenseSequentialPETScMat(int rowNum,int colNum,
				   std::map<std::string,bool>& mOptions,
				   Mat& mat,std::ofstream& logFile);

// Create a dense parallel PETSc coefficient matrix.
void createDenseParallelPETScMat(int localRowNum,int localColNum,
				 int globalRowNum,int globalColNum,
				 std::map<std::string,bool>& mOptions,
				 Mat& mat,std::ofstream& logFile);

// Create a sparse parallel PETSc coefficient matrix.
void createSparseParallelPETScMat(int localRowNum,int localColNum,
				  int globalRowNum,int globalColNum,
				  int maxDiagRowNonzeros,
				  intVector& diagEntries,
				  int maxOffDiagRowNonzeros,
				  intVector& offDiagEntries,
				  std::map<std::string,bool>& mOptions,
				  Mat& mat,std::ofstream& logFile);

// Create a parallel PETSc vector.
void createParallelPETScVec(int localRowNum,int globalRowNum,Vec& vec);


// Create a parallel PETSc solver object.
void createParallelPETScSolver(Mat& mat,
			       std::map<std::string,double>& problemData,
			       std::vector<std::string> solvingMethod,
			       KSP& ksp,PC& pc,
			       std::map<std::string,double>& calcData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile);

void writeKSPConvergenceReason(KSPConvergedReason reason,
			       int solverIterations,
			       std::ofstream& logFile);

// check whether PETSc equation solve was successful
bool checkPETScConvergence(KSP& ksp,int& solverIterations,
			   std::ofstream& logFile);

// Assemble a vector containing all global entries of a parallel PETSc
// vector.
void getGlobalParallelPETScVec(int localSize,int globalSize,Vec& vec,
			       dbVector& globalVec);

void destroyPETScSolver(KSP& ksp);

void destroyPETScMat(Mat& mat);

void destroyPETScVec(Vec& vec);

// compute the conditioning number of a sparse parallel matrix
void checkMatrixConditioning(Mat& mat,
			     std::map<std::string,double>& problemData,
			     std::ofstream& logFile);

// compute the conditioning number of a sparse parallel matrix
void checkMatrixConditioning(KSP& ksp,
			     std::map<std::string,double>& problemData,
			     std::ofstream& logFile);


#endif /* PETSCFUNCTIONS_H_ */
