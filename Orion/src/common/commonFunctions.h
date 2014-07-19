/*
 * commonFunctions.h
 *
 *  Created on: 03 Jul 2014
 *      Author: ritesh
 */

#ifndef COMMONFUNCTIONS_H_
#define COMMONFUNCTIONS_H_

#include "float.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include "mpi.h"
#include <unistd.h>
#include <vector>
#include <unistd.h>

#include "commonTypedefs.h"
#include "defs.h"
#include "fortranFunctions.h"

#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"


#include <stdio.h>
#include <assert.h>

// ======================================================================
// Ritesh:
#include <string>

/***********************************************************************/
// Clear arrays
void clearArray(blVector& v);
//void clearArray(intVector& v);
//void clearArray(dbVector& v);
template<typename T>
void clearArray(std::vector<T>& v) {
	using namespace std;

	v.assign(v.size(), 0);
}


//void clearArray(intMatrix& T);
//void clearArray(dbMatrix& T);
template<typename T>
void clearArray(std::vector<std::vector<T> >& v) {
	using namespace std;

	for (int i = 0; i < v.size(); i++)
		v[i].assign(v[i].size(), 0);

}

//void clearArray(dbMatrix3& T);
//void clearArray(intMatrix3& T);
template<typename T>
void clearArray(std::vector<std::vector<std::vector<T> > >& mat) {

  using namespace std;

  for(int i=0;i<mat.size();i++)
    for(int j=0;j<mat[i].size();j++)
    	mat[i][j].assign(mat[i][j].size(),0);

}

void clearArray(dbMatrix4& T);
void clearArray(dbMatrix5& T);

// set all entries of a map<string,bool> to false
void clearMap(std::map<std::string,bool>& data,
	      std::ofstream& logFile);

void clearMap(std::map<std::string,bool>& data);


/***********************************************************************/
// Quicksort functions: left = first valid entry index (e.g. left = 0)
//                      right = last valid entry index (e.g.
//                                               right = vec.size() - 1

// Quicksort a int vector.
void sortIntVector(intVector& vec,int left,int right);

// Quicksort a double vector.
void sortDoubleVector(dbVector& vec,int left,int right);

// Quicksort a integer array together with its indices.
void sortValuesIdx(intVector& values,intVector& idx,
		   int left,int right);

// Quicksort a double array together with its indices.
void sortValuesIdx(dbVector& values,intVector& idx,
		   int left,int right);

// Rerrange the elements of a vector
// vOrder[new] = old
template<typename T>
void reorderVector(std::vector<T>& vA,std::vector<int> vOrder) {

  using namespace std;

  assert(vA.size() == vOrder.size());

  // for all elements to put in place
  for(int i=0;i<vA.size()-1;++i) {

    // while the element i is not yet in place
    while(i != vOrder[i]) {

       // swap it with the element at its final place
       int alt = vOrder[i];
       swap(vA[i],vA[alt] );
       swap(vOrder[i],vOrder[alt]);
     }
   }
};

// Sort the particles in 3 dimensions.
void sortParticles(dbMatrix& ptcls,intVector& idx,int sortCoord,
		   int numOfPlanes,std::ofstream& logFile);

/***********************************************************************/
// Search for a certain integer-vector entry and return the position.
int findIntVecPos(int& entry,int startIdx,int endIdx,
		  intVector& vec);

int findDoubleVecPos(double& entry,int startIdx,int endIdx,
		  dbVector& vec);

// Search in a matrix column for a certain integer-matrix entry
// and return the position.
int findIntMatPos(int& entry,int startIdx,int endIdx,int vecPos,
		  intMatrix& mat);

// Search for a certain row in an integer matrix and return the position.
int findIntMatPos(intVector& entries,int startIdx,int endIdx,
		  std::string& mode,intMatrix& mat);

// Search for a certain row in double matrix and return the position.
int findDoubleMatPos(dbVector& entries,int startIdx,int endIdx,
		     std::string mode,double tol,dbMatrix& mat,
		     std::ofstream& logFile);

/***********************************************************************/
// Compare two integer vectors.
bool compareIntVecs(std::string& mode,intVector& small,intVector& big,
		    std::ofstream& logFile);

// remove redundant entries of a integer vector
void removeRedundantEntries(std::vector<int>& vec,int startIdx,int& endIdx,
			    std::ofstream& logFile);

// remove redundant entries in a certain column of a matrix
void removeRedundantEntries(std::vector<std::vector<int> >& mat,int idx);

// find max value in a vector
int findMaxValue(std::vector<int>& vec);

// find min value in a vector
int findMinValue(std::vector<int>& vec);

/***********************************************************************/
// Calculate scalar product of two vectors.
void scalarProduct(dbVector& a,dbVector& b,double& scalar,
		   std::ofstream& logFile);

// Calculate scalar product of two matrices or tensors.
void scalarProduct(dbMatrix& A,dbMatrix& B,double& scalar,bool aTrans,
		   bool bTrans,std::ofstream& logFile);


// Calculate inner product of matrix or tensor T and vector u.
void innerTensorProduct(dbMatrix& T,dbVector& u,dbVector& result,
			bool tTrans,std::ofstream& logFile);


// Calculate inner product of matrix or tensor T and tensor S.
void innerTensorProduct(dbMatrix& T,dbMatrix& S,dbMatrix& Result,
			bool tTrans,bool sTrans,std::ofstream& logFile);


// Calculate cross product of two vectors u and v .
void crossProduct(dbVector& u,dbVector& v,dbVector& result);

// Calculate cross product of two vectors u and v .
void crossProduct(dbVector& u,dbVector& v,dbVector& result,double& det,
		  std::ofstream& logFile);

/***********************************************************************/
// Set the Kronecker symbol.
dbMatrix getKroneckerSymbol(int size);

// Calculate the permutations unequal zero if the size of used arrays
// is 'size'.
intMatrix getPermutations(int size);

// Determine the missing indices which build permutations for
// a given permutation index.
intMatrix3 getPermutationMissingIdx(int dim,int pos);

// Test routine for kontinuum mechanics' stuff.
void testContinuumMechanics(std::ofstream& logFile);

// Return Pi
double getPi();

/************************************************************************/
// Calculate the determinante of a single precision matrix.
int calcDetSingle(dbMatrix& Amat,double& det,std::ofstream& logFile);

// Calculate the determinante of a square double precision matrix.
int calcDetDouble(dbMatrix& Amat, double& det);
int calcDetDouble(dbMatrix& Amat, double& det, std::ofstream& logFile);

// Calculate the determinante of a small double precision matrix.
void calcDetDoubleSmall(dbMatrix& Amat, double& det, std::ofstream& logFile);
// Calculate the determinante of a double precision tensor.
void calcDetDoubleTens(dbMatrix& Amat, double& det, std::ofstream& logFile);

// Calculate the determinante of a dense square double precision matrix
// making use of the Cayley-Hamilton theorem
int calcDetDoubleDense(dbMatrix& Amat,double& det,std::ofstream& logFile);
int calcDetDoubleDense(dbMatrix& Amat,double& det);

// Calculate the determinante of a sparse square double precision matrix
// using LAPACK routines.
int calcDetDoubleSparse(dbMatrix& Amat,double& det,std::ofstream& logFile);
int calcDetDoubleSparse(dbMatrix& Amat,double& det);

// Calculate the single precision inverse of a matrix.
void calcInvSingle(dbMatrix& Amat,dbMatrix& inverse,
		   std::ofstream& logFile);

// Calculate the double precision inverse of a matrix.
void calcInvDouble(dbMatrix& Amat,dbMatrix& inverse,
		   std::ofstream& logFile);

// Invert 3x3 matrix without using fortran.
void calcInv3x3(dbMatrix& A, dbMatrix& inv, std::ofstream& logFile);

// Calculate the double precision inverse of a dense square matrix
// making use of the Cayley-Hamilton theorem
int calcInvDoubleDense(dbMatrix& Amat,dbMatrix& inverse,
		       std::ofstream& logFile);

// Calculate the double precision inverse of a sparse square matrix using
// LAPACK routines.
int calcInvDoubleSparse(dbMatrix& Amat,dbMatrix& inverse,
			std::ofstream& logFile);

// Check the Positive Definiteness of a rank 4 Tensor
void posDefCheck4(dbMatrix4& C4,dbVector& det,std::ofstream& logFile);

// Check the Positive Definiteness of a Matrix or rank 2 Tensor
void posDefCheck(dbMatrix& C4,dbVector& det,std::ofstream& logFile);

/***********************************************************************/
// Compute a vector norm.
double computeNorm(dbVector& aVector,int type,std::ofstream& logFile);

void calcL1Norm(dbVector& oldGlobalForceVec,double& residualNorm,
		std::ofstream& logFile);

void calcL2Norm(dbVector& oldGlobalForceVec,double& residualNorm,
		std::ofstream& logFile);

void calcLInftyNorm(dbVector& oldGlobalForcVec,double& residualNorm,
		    std::ofstream& logFile);

// Compute a matrix norm.
double computeNorm(dbMatrix& aMatrix,int type,std::ofstream& logFile);

/***********************************************************************/
void getFEMMeshData(intVector& data,std::map<std::string,double>& params);
void getGaussQuadratureData(intVector& data,
			    std::map<std::string,double>& params);

/***********************************************************************/
// Auxiliary functions to dealt with symmetric tensors.

// converts the indices of a symmetric matrix to the indices of a
// cooresponding vector (reduced)
intMatrix matrixToVector(int size);

// converts the indices of a general matrix to the indices of a
// cooresponding vector (full)
intMatrix matrixToVectorFull(int size);

// converts the indices of a vector to the pairs of indices of a
// cooresponding symmetric matrix
intMatrix vectorToMatrix(int size);
intMatrix vectorToMatrix(intVector& idx);

// Get the symmetric matrix which is stored as a vector.
dbMatrix convertVectorToMatrix(int size,dbVector& vec);

/***********************************************************************/
// Computation of roots of polynom 3rd degree ax^3 + bx^2 + cx + d
void polynom3rdRoots(dbVector& coefficients,dbVector& roots,
		     std::ofstream& logFile);

dbVector quadraticPolynom(double x);
dbVector dQuadraticPolynom(double x);
dbVector d2QuadraticPolynom(double x);

dbVector cubicPolynom(double x);
dbVector dCubicPolynom(double x);
dbVector d2CubicPolynom(double x);
dbVector d3CubicPolynom(double x);

dbVector quarticPolynom(double x);
dbVector dQuarticPolynom(double x);
dbVector d2QuarticPolynom(double x);
dbVector d3QuarticPolynom(double x);

/***********************************************************************/
// Compute all eigenvalues of a real symmetric matrix.
int calcEigenValues(dbMatrix& aMatrix,dbVector& eigenvalues,
		    std::ofstream& logFile);

// Compute all eigenvalues of a real general matrix.
int calcEigenValues(dbMatrix& aMatrix,dbMatrix& eigenvalues,
		    std::ofstream& logFile);

// Compute all eigenvalues of a generalized eigenvalue problem
// consisting of real general matrices.
int calcEigenValues(dbMatrix& aMatrix,dbMatrix& bMatrix,
		    dbMatrix& eigenvalues,std::ofstream& logFile);

// Compute the maximum eigenvalue of a real general matrix.
int computeMaxEigenValue(dbMatrix& aMatrix,double& maxEigenvalue,
			 std::ofstream& logFile);

// Compute all eigenvalues and eigenvectors of a real symmetric matrix.
int calcEigenValuesVectors(dbMatrix& aMatrix,dbVector& eigenvalues,
			   dbMatrix& eigenvectors,std::ofstream& logFile);

// Compute all eigenvalues and eigenvectors of a real general matrix.
int calcEigenValuesVectors(dbMatrix& aMatrix,dbMatrix& eigenvalues,
			   dbMatrix& eigenvectors,std::ofstream& logFile);

// Compute the singular value decomposition of a real general matrix
// A = USV^T.
int computeGeneralMatrixSVD(dbMatrix& A,dbVector& S,dbMatrix& U,
			    dbMatrix& V,std::ofstream& logFile);

// Compute the natural logarithm of a diagonalizable symmetric matrix.
int calcMatrixLn(dbMatrix& A,dbMatrix& lnA,std::ofstream& logFile);

// Compute the natural logarithm of a diagonalizable nonsymmetric matrix.
int calcNonsymMatrixLn(dbMatrix& A,dbMatrix& lnA,std::ofstream& logFile);

// Compute the natural logarithm of a diagonalizable non-symmetric
// dense matrix using a Taylor series with an initial value A0 = 1
// (i.e. ||A|| approximately 1)
int taylorApproxDenseMatrixLn(dbMatrix& A,dbMatrix& lnA,
			      std::ofstream& logFile);

// Compute the natural logarithm of a diagonalizable non-symmetric
// dense matrix (Cheng, Higham, Kenney and Laub 2001, p. 1118).
int highamApproxDenseMatrixLn(dbMatrix& A,dbMatrix& lnA,
			      std::ofstream& logFile);

// Compute the exponential of a diagonalizable non-symmetric
// dense matrix using a Taylor series
int taylorApproxDenseMatrixExp(dbMatrix& A,dbMatrix& expA,
			       std::ofstream& logFile);

// Compute the exponential of a diagonalizable non-symmetric matrix.
int calcNonsymMatrixExp(dbMatrix& A,dbMatrix& expA,std::ofstream& logFile);

// Compute the exponential of a matrix and its partial derivatives with
// respect to this matrix
//void calcExponential(dbMatrix& A,dbMatrix& expA,dbMatrix& dexpA,
//		     dbMatrix3& d2expA,int derivOrder,
//		     std::ofstream& logFile);


/***********************************************************************/
// Vector allocation routines

//void allocateArray(intVector& vec,int idx);
//void allocateArray(dbVector& vec,int idx);
//void allocateArray(blMatrix& mat,int idx1);
//void allocateArray(dbMatrix& mat,int idx1);
//void allocateArray(dbMatrix3& mat,int idx);
//void allocateArray(dbMatrix4& mat,int idx);
template<typename T>
void allocateArray(std::vector<T>& vec, int idx) {

	using namespace std;

	vec.resize(idx);

}


//void allocateArray(blMatrix& mat,int idx1,int idx2);
//void allocateArray(intMatrix& mat,int idx1,int idx2);
//void allocateArray(dbMatrix& mat,int idx1,int idx2);
//void allocateArray(ldbMatrix& mat,int idx1,int idx2);
template<typename T>
void allocateArray(std::vector<std::vector<T> >& mat,int idx1,int idx2) {

  using namespace std;

  mat.resize(idx1);

  for(int i=0;i<idx1;i++)
    mat[i].resize(idx2);

}

//void allocateArray(intMatrix3& mat,int idx1,int idx2,int idx3);
//void allocateArray(dbMatrix3& mat,int idx1,int idx2,int idx3);
template<typename T>
void allocateArray(std::vector<std::vector<std::vector<T> > >& mat,
		int idx1,int idx2,int idx3) {

  using namespace std;

  mat.resize(idx1);

  for(int i=0;i<idx1;i++) {
    mat[i].resize(idx2);

    for(int j=0;j<idx2;j++) {
      mat[i][j].resize(idx3);

    }

  }

}

void allocateArray(intMatrix3& outMat,intMatrix3& inMat);
void allocateArray(dbMatrix4& mat,int idx1,int idx2,int idx3,int idx4);
void allocateArray(dbMatrix5& mat,int idx1,int idx2,int idx3,int idx4,int idx5);


//void pushBackVector(intVector& vec,int n);
//void pushBackVector(dbVector& vec,double n);
template<typename T>
void pushBackVector(std::vector<T>& vec, T n){

	  using namespace std;

	  vec.resize(vec.size()+1);
	  vec[vec.size()-1] = n;
}

//void pushBackVector(intMatrix& mat,intVector& vec);
//void pushBackVector(dbMatrix& mat,dbVector& vec);
template<typename T>
void pushBackVector(std::vector<std::vector<T> >& mat,std::vector<T>& vec){
	using namespace std;

	mat.resize(mat.size()+1);
	mat[mat.size()-1] = vec;
}

void pushBackVector(intMatrix& oldMat,intMatrix& deltaMat);

/************************************************************************/
// Array resize routines which shrink and fit also its reserved capacity
//void resizeArray(intVector& vec,int size);
//void resizeArray(dbVector& vec,int size);
template<typename T>
void resizeArray(std::vector<T>& vec,int size) {
  vec.resize(size);
  std::vector<T>(vec).swap(vec);
}

//void resizeArray(intMatrix& mat,int size);
//void resizeArray(dbMatrix& mat,int size);
template<typename T>
void resizeArray(std::vector<std::vector<T> >& mat,int size) {
  mat.resize(size);
  std::vector<std::vector<T> >(mat).swap(mat);
}

void resizeArray(dbMatrix& mat,int size1,int size2);
void resizeArray(dbMatrix3& mat,int size);
void resizeArray(dbMatrix4& mat,int size);

/***********************************************************************/
// convert a integer number to a string
std::string convertIntToString(int number);

// Check if a double precision number equals to zero within numerical
// precision limits.
bool checkZero(double& value);

// Determine the twist of a point given by its reference configuration
// and its current configuration around a given axis going through
// the origin.
//
// Buggy for larger rotations
double calcPlacementTwist(dbVector& reference,dbVector& current,
			  dbVector& rotAxis,std::ofstream& logFile);

// Determine the twist of a point given by its origin 'reference' and its
// image 'current' around a given axis going through a certain point
// (fix point).
//
// buggy?!
double calcPlacementTwist(dbVector& reference,dbVector& current,
			  dbVector& axis,dbVector& fixPoint,
			  std::ofstream& logFile);

// Project a vector onto a plane defined by its normal vector.
void projectVecToPlane(dbVector& a,dbVector& normal,
		       dbVector& b,std::ofstream& logFile);

// Calculate an orthogonal basis to a given vector v
void calcOrthogonalBasis(dbVector& vec,dbMatrix& basis,
			 std::ofstream& logFile);

// Calculate rotation tensor.
void calcRotationMatrix(dbVector& rotVec,dbMatrix& rTens,
			std::ofstream& logFile);

// Determine the rotation matrix corresponding the rotation of the
// circumferential tangent vector around the cylinder's axis
// so that the tangent goes through the given position
//
// axis    ... cylinder axis (x1,x2,x3 Cartesian axis)
// X       ... position vector
// R       ... rotation matrix (output)
// dR      ... spatial derivatives of rotation matrix (output)
//
void setCylinderMetrics(dbVector& direction,dbVector& X,dbMatrix& R,
			dbMatrix3& dR,std::ofstream& logFile);

// compute the distance between two points
double calcDistance(dbVector& P,dbVector& Q);

// populate a symmetric matrix with upper triangle inputted
void populateSymmetricMatrix(dbMatrix& sym, int size);

/************************************************************************/

// Determine an equal distribution of an array over the processors.
void getParallelArrayDistribution(int globalNum,int& startIdx,
				  int& endIdx,std::ofstream& logFile);


/************************************************************************/

// returns the process' physical memory usage in Mbyte.
int process_RAM_usage();

// returns the process' memory usage in Mbyte.
int process_mem_usage();

// returns total physical memory of one node in Mbyte
int getTotalRAM();

// returns the free physical memory of one node in Mbyte
int getFreeRAM();

// check the RAM usage
void checkRAMusage(std::ofstream& logFile);

// check the current local RAM usage in Mbyte
void checkRAMusage(int numberOfProcsPerNode,std::ofstream& logFile);

int convertByteToMByte(long int value);
int convertkByteToMByte(long int value);

/************************************************************************/
// Geometric information

//convert from cartesian to prolate spheroidal coordinates.
void cartToProSpheroid(dbVector& cartCoords,
                       dbVector& proSpheroidCoords, std::map<std::string,double>& modelData,
                       std::ofstream& logFile);

//convert from prolate spheroidal to cartesian coordinates.
void proSpheroidToCart(dbVector& cartCoords,
                       dbVector& proSpheroidCoords, std::map<std::string,double>& modelData);

//get prolate spheroidal unit basis vectors at a point.
void getProSpheroidVecs(dbVector& proSpheroidCoords,
                     dbMatrix& unitVectors, std::map<std::string,double>& modelData);


//convert from cartesian to prolate spheroidal coordinates more accurately
void cart2ProSpheroid(dbVector& cartCoords,
                       dbVector& proSpheroidCoords, std::map<std::string,double>& modelData,
                       std::ofstream& logFile);

//combine n sorted vectors of differing lengths into a single vector
//recursive function
//rows represent each vector
void combineSortedVecs(intMatrix& inputMat);

// combines sorted vectors inputVec1 and inputVec2 and puts the sorted vector into inputVec1
void combine2SortedVecs(intVector& inputVec1, intVector& inputVec2);

// ========================================================================
// ========================================================================
// Ritesh common functions

// Display contents of a matrix
//void printMatrix(dbMatrix& matrix);
//void printMatrix(dbMatrix& matrix, std::ofstream& logFile);
//void printMatrix(dbMatrix& matrix, std::string& msg);
void printMatrix(dbMatrix& matrix, const char* msg, std::ofstream& logFile);
void printMatrix(intMatrix& matrix, const char* msg, std::ofstream& logFile);

// Display contents of a vector
//void printVector(dbVector& A);
//void printVector(dbVector& A, std::ofstream& logFile);
//void printVector(dbVector& A, std::string& msg);
void printVector(dbVector& A, const char* msg, std::ofstream& logFile);

// Display contents of a vector
//void printVector(intVector& A);
//void printVector(intVector& A, std::ofstream& logFile);
//void printVector(intVector& A, std::string& msg);
void printVector(intVector& A, const char* msg, std::ofstream& logFile);




#endif /* COMMONFUNCTIONS_H_ */
