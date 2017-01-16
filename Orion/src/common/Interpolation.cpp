#include "Interpolation.h"

using namespace std;

Interpolation::Interpolation(dbVector& iPoint, dbMatrix& coords,
		dbVector& radiusVec, dbVector& interpolants, InputFileData* InputData,
		ofstream& logFile) {

	int choice = InputData->getValue("MLSCalculationType");

	switch (choice) {
	case 1:
		MLSUnitTest(InputData,logFile);
		break;

	case 2:
		interpolants = MLSCalc(iPoint, coords, radiusVec, InputData, logFile);
		break;

	default:
		logFile << "ERROR: Interpolation type = " << choice << " does not exist"
				<< endl;
		cout << "ERROR: Interpolation type = " << choice << " does not exist"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}
;

///*!****************************************************************************/
////! Calculate the Moving Least Square interpolants
//dbVector Interpolation::MLSCalc(dbVector iPoint, dbMatrix coords,
//		dbVector radiusVec, InputFileData* InputData, ofstream& logFile) {
//
//#ifdef _MLSDebugMode_
//	logFile << "********* In Interpolation::MLSCalc *********" << endl;
//	printVector(iPoint, "iPoint", logFile);
//	printMatrix(coords, "Coords", logFile);
//#endif
//
//	removeRedundantDim(iPoint, coords, radiusVec, InputData, logFile);
//
//	int numCoords = coords.size();
//	int ndim = iPoint.size();
//
//	bool activate_SBM;
//	if (InputData->getValue("shiftedBasisAlgorithm") == 1)
//		// Activate the shifted basis correction to enhance MLS stability
//		activate_SBM = true;
//	else
//		activate_SBM = false;
//
//	//! ========================================================================
//	//! Search for the neighbours of the iPoint
//
////	intVector neighbours = findNeighbours(iPoint, coords, radiusVec,
////			logFile);
////	dbVector neighbours_db(neighbours.begin(), neighbours.end());
////	msg = "Neighbours"; printVector(neighbours_db, msg, logFile);
//
//	// Selected all coordinates as neighbours
//	intVector neighbours(coords.size());
//	for (int i = 0; i < coords.size(); i++) {
//		neighbours[i] = i;
//	}
//
//	//! ========================================================================
//	//! Setup the basis polynomial, calculate the weight function
//	//! and the B_matrix
//	dbMatrix basisPolyMat, B_matrix;
//	dbVector weightVec;
//	int degPoly = InputData->getValue("MLSPolynomialDegree");
//	for (int i = 0; i < neighbours.size(); i++) {
//
//		// Weight function
//		double weight = cubicSplineWgtCalc(iPoint, coords[neighbours[i]],
//				radiusVec, logFile);
//		weightVec.push_back(weight);
//
//#ifdef _MLSDebugMode_
//		printVector(radiusVec, "radiusVec", logFile);
//#endif
//
//		//Basis function
//		dbVector nCoords(ndim);
//		if (activate_SBM == true) {
//			for (int j = 0; j < ndim; j++) {
//				nCoords[j] = (coords[neighbours[i]][j] - iPoint[j])
//						/ radiusVec[j];
//			}
//
//#ifdef _MLSDebugMode_
//			printVector(nCoords, "Shifted Coords", logFile);
//#endif
//
//		} else {
//			for (int j = 0; j < ndim; j++) {
//				nCoords[j] = coords[neighbours[i]][j];
//			}
//		}
//
//		dbVector basisVec = PascalBasisCalc(nCoords, degPoly, logFile);
//		basisPolyMat.push_back(basisVec);
//
//
//
//		// B_matrix setup: w*P
//		B_matrix.push_back(dbVector(0));
//		for (int j = 0; j < basisVec.size(); j++) {
//			B_matrix[i].push_back(weight * basisVec[j]);
//		}
//	}
//
//	if(basisPolyMat[0].size()-1 > numCoords){
//		logFile << "In Interpolation::MLSCalc, bad interpolants calculations."
//				"Number of Basis polynomial unknowns is higher than the number of supports."
//				<< endl;
//	}
//
//	//! ========================================================================
//	//! Assemble the Moment Matrix : w*P*(P^T)
//	dbMatrix momentMat(basisPolyMat[0].size(),
//			dbVector(basisPolyMat[0].size(), 0));
//	for (int i = 0; i < neighbours.size(); i++) {
//		for (int j = 0; j < basisPolyMat[i].size(); j++) {
//			for (int k = 0; k < basisPolyMat[i].size(); k++) {
//				momentMat[j][k] += weightVec[i] * basisPolyMat[i][j]
//						* basisPolyMat[i][k];
//			}
//		}
//	}
//
//#ifdef _MLSDebugMode_
//	printVector(weightVec, "Weight Vector", logFile);
//	printMatrix(basisPolyMat, "Basis Polynomial Matrix", logFile);
//	printMatrix(B_matrix, "B_matrix", logFile);
//	printMatrix(momentMat, "Moment Matrix", logFile);
//#endif
//
//	//! ========================================================================
//	//! Calculate the inverse of the moment matrix
//	dbMatrix inverseMomentMat(momentMat.size(), dbVector(momentMat[0].size()));
//	calcInvDouble(momentMat, inverseMomentMat, logFile);
//
//	//! ========================================================================
//	//! Calculate the basis of iPoint
//	dbVector iPointBasisVec;
//	if (activate_SBM == true) {
//		dbVector newiPoint(ndim, 0);
//		iPointBasisVec = PascalBasisCalc(newiPoint, degPoly, logFile);
//	} else {
//		iPointBasisVec = PascalBasisCalc(iPoint, degPoly, logFile);
//	}
//
//	//! ========================================================================
//	//! Calculate the shape functions
//	dbVector interpolants(neighbours.size(), 0);
//	for (int i = 0; i < neighbours.size(); i++) {
//		for (int j = 0; j < iPointBasisVec.size(); j++) {
//			for (int k = 0; k < iPointBasisVec.size(); k++) {
//				interpolants[i] += iPointBasisVec[j] * inverseMomentMat[j][k]
//						* B_matrix[i][k];
//			}
//		}
//	}
//
//#ifdef _MLSDebugMode_
//	printMatrix(inverseMomentMat, "Inverse Moment Matrix", logFile);
//	printVector(iPointBasisVec, "iPointBasisVec", logFile);
//	printVector(interpolants, "Interpolants", logFile);
//
//	double intSum = 0;
//	for (int i = 0; i < interpolants.size(); i++)
//		intSum += interpolants[i];
//
//	logFile << "Sum of interpolants: " << intSum << endl;
//#endif
//
//	return interpolants;
//
//}
//;

/*!****************************************************************************/
//! Calculate the Moving Least Square interpolants
dbVector Interpolation::MLSCalc(dbVector iPoint, dbMatrix coords,
		dbVector radiusVec, InputFileData* InputData, ofstream& logFile) {

#ifdef _MLSDebugLiteMode_
	logFile << "********* In Interpolation::MLSCalc *********" << endl;
	printVector(iPoint, "iPoint", logFile);
	printMatrix(coords, "Coords", logFile);
	printVector(radiusVec, "radiusVec", logFile);
#endif

	if(InputData->getValue("shiftedBasisType") != 0)
		shiftedBasisCalc(iPoint, coords, radiusVec, InputData, logFile);

	removeRedundantDim(iPoint, coords, radiusVec, InputData, logFile);

	int numCoords = coords.size();
	int ndim = iPoint.size();

#ifdef _MLSDebugMode_
	if(InputData->getValue("shiftedBasisAlgorithm") == 1)
		logFile << "shiftedBasisAlgorithm activated" << endl;
	else
		logFile << "shiftedBasisAlgorithm de-activated" << endl;
#endif


	//! ========================================================================
	//! Search for the neighbours of the iPoint

	// Selected all coordinates as neighbours
	intVector neighbours(coords.size());
	for (int i = 0; i < coords.size(); i++) {
		neighbours[i] = i;
	}

	//! ========================================================================
	//! Setup the basis polynomial, calculate the weight function
	//! and the B_matrix
	dbMatrix basisPolyMat, B_matrix;
	dbVector weightVec = weightCalc(iPoint, neighbours, coords, radiusVec,
			 	 	 	 	 	 	 	 	 	 	 	 	InputData, logFile);

	int degPoly = InputData->getValue("MLSPolynomialDegree");
	for (int i = 0; i < neighbours.size(); i++) {

		dbVector nCoords = coords[i];

		// Basis function
		dbVector basisVec = PascalBasisCalc(nCoords, degPoly, logFile);
		basisPolyMat.push_back(basisVec);

		// B_matrix setup: w*P
		B_matrix.push_back(dbVector(0));
		for (int j = 0; j < basisVec.size(); j++) {
			B_matrix[i].push_back(weightVec[i] * basisVec[j]);
		}
	}

	// Basic check if moment matrix will be invertible
	if (basisPolyMat[0].size() - 1 > numCoords) {
		logFile << "In Interpolation::MLSCalc, bad interpolants calculations."
					"Number of Basis polynomial unknowns is higher than the"
					" number of supports." << endl;
	}

	//! ========================================================================
	//! Assemble the Moment Matrix : w*P*(P^T)
	dbMatrix momentMat(basisPolyMat[0].size(),
			dbVector(basisPolyMat[0].size(), 0));
	for (int i = 0; i < neighbours.size(); i++) {
		for (int j = 0; j < basisPolyMat[i].size(); j++) {
			for (int k = 0; k < basisPolyMat[i].size(); k++) {
				momentMat[j][k] += weightVec[i] * basisPolyMat[i][j]
						* basisPolyMat[i][k];
			}
		}
	}

#ifdef _MLSDebugMode_
	printVector(weightVec, "Weight Vector", logFile);
	printMatrix(basisPolyMat, "Basis Polynomial Matrix", logFile);
	printMatrix(B_matrix, "B_matrix", logFile);
	printMatrix(momentMat, "Moment Matrix", logFile);
#endif

	//! ========================================================================
	//! Calculate the inverse of the moment matrix
	dbMatrix inverseMomentMat(momentMat.size(), dbVector(momentMat[0].size()));
	calcInvDouble(momentMat, inverseMomentMat, logFile);

	//! ========================================================================
	//! Calculate the basis of iPoint
	dbVector iPointBasisVec = PascalBasisCalc(iPoint, degPoly, logFile);

	//! ========================================================================
	//! Calculate the shape functions
	dbVector interpolants(neighbours.size(), 0);
	for (int i = 0; i < neighbours.size(); i++) {
		for (int j = 0; j < iPointBasisVec.size(); j++) {
			for (int k = 0; k < iPointBasisVec.size(); k++) {
				interpolants[i] += iPointBasisVec[j] * inverseMomentMat[j][k]
						* B_matrix[i][k];
			}
		}
	}

#ifdef _MLSDebugMode_
	printMatrix(inverseMomentMat, "Inverse Moment Matrix", logFile);
	printVector(iPointBasisVec, "iPointBasisVec", logFile);
#endif

#ifdef _MLSDebugLiteMode_
	printVector(interpolants, "Interpolants", logFile);

	double intSum = 0;
	for (int i = 0; i < interpolants.size(); i++)
		intSum += interpolants[i];

	logFile << "Sum of interpolants: " << intSum << endl;
#endif

	return interpolants;

}
;

/*!****************************************************************************/
//!
void Interpolation::shiftedBasisCalc(dbVector& iPoint, dbMatrix& coords,
		dbVector& radiusVec, InputFileData* InputData, ofstream& logFile) {

	// Determine the num of dimensions
	int ndim = coords[0].size();

	int option = InputData->getValue("shiftedBasisType");

	switch (option) {

	case 1: {
		dbVector lowestCoords = coords[0];

		// find the lowest value in each direction
		for (int i = 1; i < coords.size(); i++) {
			for (int j = 0; j < ndim; j++) {
				if (lowestCoords[j] > coords[i][j]) {
					lowestCoords[j] = coords[i][j];
				}
			}
		}

		// normalise the coords with respect to the lowest value
		for (int i = 1; i < coords.size(); i++) {
			for (int j = 0; j < ndim; j++) {
				coords[i][j] = coords[i][j] / lowestCoords[j];
			}
		}

		// normalise ipoint
		for (int j = 0; j < ndim; j++) {
			iPoint[j] = iPoint[j] / lowestCoords[j];
		}

		break;
	}

	case 2: {

		double shiftedBasisCoef = InputData->getValue("shiftedBasisCoef");
//		logFile << "shiftedBasisCoef used: " << shiftedBasisCoef << endl;

		for (int i = 0; i < coords.size(); i++) {
			for (int j = 0; j < ndim; j++) {
				coords[i][j] = ((coords[i][j] - iPoint[j]) / radiusVec[j])
															+ shiftedBasisCoef;
//				logFile << "coords[i][j] = " << coords[i][j] << endl;
			}
		}

		dbVector newiPoint(ndim, 0);
		for (int j = 0; j < ndim; j++) {
			iPoint[j] = newiPoint[j] + shiftedBasisCoef;
			radiusVec[j] = 1;
		}

		break;
	}

	default:

		logFile << "shiftedBasisType is invalid" << endl;
		cout << "shiftedBasisType is invalid" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);

	}

#ifdef _MLSDebugMode_
	logFile << "********* In Interpolation::shiftedBasisCalc *********" << endl;
	printVector(iPoint, "iPoint", logFile);
	printMatrix(coords, "Coords", logFile);
	printVector(radiusVec, "radiusVec", logFile);
#endif

}

/*!****************************************************************************/
//! Function to setup the Pascal polynomial's basis
dbVector Interpolation::PascalBasisCalc(dbVector& coord, int& orderPoly,
		ofstream& logFile) {

	// Determine the num of dimensions
	int ndim = coord.size();

	// Account for the fact that the degree of a polynomial starts from "0"
	// Note: orderPoly++ is not recommended. If PascalBasisCalc is inside loop
	// orderPoly increase for each loop
	int degPoly = orderPoly + 1;

	intMatrix idMat(degPoly, intVector(ndim));
	dbVector previousVec;
	dbVector pascalBasisVec;

	for (int i = 0; i < degPoly; i++) {
		dbVector basisVec;

		if (i == 0) {
			basisVec.push_back(1.0);
		} else if (i == 1) {
			for (int j = 0; j < ndim; j++) {
				basisVec.push_back(coord[j]);
				idMat[i][j] = j;
			}
		} else if (i > 1) {
			for (int j = 0; j < ndim; j++) {
				double selectDim = coord[j];
				idMat[i][j] = basisVec.size();
				int startIndex = idMat[i - 1][j];

				for (int k = startIndex; k < previousVec.size(); k++) {
					basisVec.push_back((previousVec[k] * selectDim));
				}
			}
		}

		previousVec.clear();
		previousVec = basisVec;

		for (int j = 0; j < basisVec.size(); j++)
			pascalBasisVec.push_back(basisVec[j]);
	}

	return pascalBasisVec;
}
;

/*!****************************************************************************/
//! Calculate the Weight function
dbVector Interpolation::weightCalc(dbVector& iPoint, intVector& neighbours,
		dbMatrix& coords, dbVector& radiusVec, InputFileData* InputData,
		ofstream& logFile) {

	dbVector weightVec(neighbours.size(), 1);

	int choice = InputData->getValue("MLSWeightFunc");

	switch (choice) {
	// *************************************************************************
	// Cubic spline
	case 0:
		for (int i = 0; i < neighbours.size(); i++) {
			weightVec[i] = cubicSplineWgtCalc(iPoint, coords[neighbours[i]],
					radiusVec, logFile);
		}
		break;

	// *************************************************************************
	// Gaussian
	case 1:
		for (int i = 0; i < neighbours.size(); i++) {
			weightVec[i] = gaussianWgtCalc(iPoint, coords[neighbours[i]],
					radiusVec, logFile);
		}
		break;

	// *************************************************************************
	// Regularised
	case 2: {
		double weightSum = 0;
		for (int i = 0; i < neighbours.size(); i++) {
			weightVec[i] = regularizedWgtCalc(iPoint, coords[neighbours[i]],
					radiusVec, logFile);
			weightSum += weightVec[i];
		}

		for (int i = 0; i < neighbours.size(); i++) {
			weightVec[i] = weightVec[i] / weightSum;
		}

		break;
	}

	// *************************************************************************
	// Regularised modified
	case 3: {

		int ndim = iPoint.size();

		double weight = 0;
		for (int d = 0; d < ndim; d++) {

			dbVector weightDimVec(weightVec.size(), 0);
			double weightSum = 0;

			for (int i = 0; i < neighbours.size(); i++) {
				weightDimVec[i] = regularizedWgtCalc_mod(iPoint,
						coords[neighbours[i]], radiusVec, d, logFile);
				weightSum += weightDimVec[i];
			}

			for (int i = 0; i < neighbours.size(); i++) {
				weightVec[i] *= weightDimVec[i] / weightSum;
			}

		}

		break;
	}
	default:
		logFile << "ERROR: In Interpolation::weightCalc, MLSWeightFunc = "
				<< choice << "doesn't exist" << endl;
		cout << "ERROR: In Interpolation::weightCalc, MLSWeightFunc = "
				<< choice << "doesn't exist" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return weightVec;

}


/*!****************************************************************************/
//! Calculate the Cubic Spline Weight function
double Interpolation::cubicSplineWgtCalc(dbVector& coordOne, dbVector& coordTwo,
		dbVector& radii, ofstream& logFile) {

	int ndim = coordOne.size();
	double weight = 1;

	for (int i = 0; i < ndim; i++) {

		double r, w;

		r = abs(coordTwo[i] - coordOne[i]) / radii[i];
		if (r > 0.5) {
			w = (4.0 / 3.0) - (4.0 * r) + (4.0 * pow(r, 2))
					- ((4.0 / 3.0) * pow(r, 3));
		} else if (r <= 0.5) {
			w = (2.0 / 3.0) - (4.0 * pow(r, 2)) + (4.0 * pow(r, 3));
		}

		weight *= w;
	}

	return weight;

}
;

/*!****************************************************************************/
//! Calculate the Gaussian Weight function
double Interpolation::gaussianWgtCalc(dbVector& coordOne, dbVector& coordTwo,
		dbVector& radii, ofstream& logFile) {

	int ndim = coordOne.size();
	double weight = 1;
	double alpha = 0.4;

	for (int i = 0; i < ndim; i++) {

		double r, w;

		r = abs(coordTwo[i] - coordOne[i]) / radii[i];

//		w = (exp(-1*pow((r/alpha),2)) - exp(1/pow(alpha,2)))/(1-exp(-1/pow(alpha,2)));
		w = exp(-1*pow((r/alpha),2));


		weight *= w;
	}

	return weight;

}
;

/*!****************************************************************************/
//! Calculate the Cubic Spline Weight function
double Interpolation::regularizedWgtCalc(dbVector& coordOne, dbVector& coordTwo,
		dbVector& radii, ofstream& logFile) {

	int ndim = coordOne.size();
	double weight = 1;

	for (int i = 0; i < ndim; i++) {

		double r, w;

		r = abs(coordTwo[i] - coordOne[i]) / radii[i];

		w = pow(r,-2)-1;

		weight *= w;
	}

	return weight;

}
;

/*!****************************************************************************/
//! Calculate the Cubic Spline Weight function
double Interpolation::regularizedWgtCalc_mod(dbVector& coordOne, dbVector& coordTwo,
		dbVector& radii, int dim, ofstream& logFile) {

	double weight = 0;

		double r, w;

		r = abs(coordTwo[dim] - coordOne[dim]) / radii[dim];

		weight = r - 1;

	return weight;

}
;

/*!****************************************************************************/
//! Find the neighbours of a particular point
intVector Interpolation::findNeighbours(dbVector& iPoint, dbMatrix& coords,
		dbVector& radiusVec, ofstream& logFile) {

	int ndim = coords[0].size();
	intVector neighbourCheckVec(coords.size(), 0);

	for (int i = 0; i < ndim; i++) {
		for (int j = 0; j < coords.size(); j++) {
			if (abs(iPoint[i] - coords[j][i]) / radiusVec[i] <= 1) {
				neighbourCheckVec[j]++;
			}
		}
	}

	// If value in neighbourCheckVec = ndim, then node is inside influence zone
	intVector neighbourList;
	for (int k = 0; k < neighbourCheckVec.size(); k++) {
		if (neighbourCheckVec[k] == ndim)
			neighbourList.push_back(k);
	}

	return neighbourList;

}

/*!****************************************************************************/
//!
void Interpolation::removeRedundantDim(dbVector& iPoint, dbMatrix& coords,
			dbVector& radiusVec, InputFileData* InputData, ofstream& logFile){

	intVector radDimToDelete;
	for(int i=0; i < radiusVec.size(); i++){
		if(fabs(radiusVec[i]) < 1e-14){
			radDimToDelete.push_back(i);
		}
	}

	if(radiusVec.size() == radDimToDelete.size()){
		logFile << "ERROR: In Interpolation::removeRedundantDim, all radii are zero" << endl;
		cout << "ERROR: In Interpolation::removeRedundantDim, all radii are zero" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else if(radDimToDelete.size() > 0){
		for(int i = radDimToDelete.size()-1 ; i > -1 ; --i){
			int delEntry  = radDimToDelete[i];

			//delete radiusVec entry
			radiusVec.erase(radiusVec.begin()+delEntry);

			//delete coords entry
			for(int j = 0; j < coords.size(); j++){
				coords[j].erase(coords[j].begin()+delEntry);
			}

			//delete iPoint entry
			iPoint.erase(iPoint.begin()+delEntry);
		}

		printVector(iPoint,"New iPoint",logFile);
		printMatrix(coords,"New coords",logFile);
		printVector(radiusVec,"New radiusVec",logFile);
	}



}

/*!****************************************************************************/
/*!****************************************************************************/
//! Debugging
void Interpolation::MLSUnitTest(InputFileData* InputData,ofstream& logFile) {

	double length = 10, nodal_dist = 2, radiusFactor = 0.35;
	int numCoords, ndim = 3;
	dbMatrix coords;

	bool activate_SBM = true;
	intVector SBM_dim(1);
	SBM_dim[0] = ndim;
	dbVector SBM_coeff(1);
	SBM_coeff[0] = 1000000000;

	//! ========================================================================
	//! Generate nodal points over the domain
	MLSUnitTest_meshMultiDim(length, nodal_dist, ndim, coords, numCoords,
			SBM_dim, SBM_coeff, logFile);
	printMatrix(coords, "Nodal Coordinates", logFile);

	//! ========================================================================
	//! Define the iPoint(point of interest) and the radius in each dimension
	dbVector iPoint(ndim), radiusVec(ndim);
	for (int i = 0; i < ndim; i++) {

		// Choosing the middle point
		iPoint[i] = length * 0.25;

		// Set radius as a factor of the nodal distance
		radiusVec[i] = ((radiusFactor * length) / nodal_dist) * nodal_dist;

		for (int j = 0; j < SBM_dim.size(); j++)
			radiusVec[SBM_dim[j] - 1] *= SBM_coeff[j];
	}
	printVector(iPoint, "iPoint:", logFile);
	printVector(radiusVec, "radiusVec:", logFile);

	//! ========================================================================
	//! Search for the neighbours of the iPoint
	intVector neighbours = findNeighbours(iPoint, coords, radiusVec, logFile);
	dbVector neighbours_db(neighbours.begin(), neighbours.end());
	printVector(neighbours_db, "Neighbours", logFile);

	//! ========================================================================
	//! Setup the basis polynomial, calculate the weight function and the B_matrix
	dbMatrix basisPolyMat, B_matrix;
	dbVector weightVec;
	int degPoly = 1;
	for (int i = 0; i < neighbours.size(); i++) {

		// Weight function
		double weight = cubicSplineWgtCalc(iPoint, coords[neighbours[i]],
				radiusVec, logFile);
		weightVec.push_back(weight);

		//Basis function
		dbVector nCoords(ndim);
		if (activate_SBM == true) {
			for (int j = 0; j < ndim; j++) {
				nCoords[j] = (coords[neighbours[i]][j] - iPoint[j])
						/ radiusVec[j];
			}

		} else {
			for (int j = 0; j < ndim; j++) {
				nCoords[j] = coords[neighbours[i]][j];
			}
		}

		dbVector basisVec = PascalBasisCalc(nCoords, degPoly, logFile);
		basisPolyMat.push_back(basisVec);

		// B_matrix setup: w*P
		B_matrix.push_back(dbVector(0));
		for (int j = 0; j < basisVec.size(); j++) {
			B_matrix[i].push_back(weight * basisVec[j]);
		}

	}
	printVector(weightVec, "Weight Vector", logFile);
	printMatrix(basisPolyMat, "Basis Polynomial Matrix", logFile);
	printMatrix(B_matrix, "B_matrix", logFile);

	//! ========================================================================
	//! Assemble the Moment Matrix : w*P*(P^T)
	dbMatrix momentMat(basisPolyMat[0].size(),
			dbVector(basisPolyMat[0].size(), 0));
	for (int i = 0; i < neighbours.size(); i++) {
		for (int j = 0; j < basisPolyMat[i].size(); j++) {
			for (int k = 0; k < basisPolyMat[i].size(); k++) {
				momentMat[j][k] += weightVec[i] * basisPolyMat[i][j]
						* basisPolyMat[i][k];
			}
		}
	}

	printMatrix(momentMat, "Moment Matrix:", logFile);

	//! ========================================================================
	//! Calculate the inverse of the moment matrix
	dbMatrix inverseMomentMat(momentMat.size(), dbVector(momentMat[0].size()));

	calcInvDouble(momentMat, inverseMomentMat, logFile);

	printMatrix(inverseMomentMat, "Inverse Moment Matrix", logFile);

	//! ========================================================================
	//! Calculate the shape functions
	dbVector iPointBasisVec;
	if (activate_SBM == true) {
		dbVector newiPoint(ndim, 0);
		iPointBasisVec = PascalBasisCalc(newiPoint, degPoly, logFile);
	} else {
		iPointBasisVec = PascalBasisCalc(iPoint, degPoly, logFile);
	}

	dbVector interpolants(neighbours.size(), 0);

	for (int i = 0; i < neighbours.size(); i++) {
		for (int j = 0; j < iPointBasisVec.size(); j++) {
			for (int k = 0; k < iPointBasisVec.size(); k++) {
				interpolants[i] += iPointBasisVec[j] * inverseMomentMat[j][k]
						* B_matrix[i][k];
			}
		}
	}

	printVector(iPointBasisVec, "iPointBasisVec", logFile);
	printVector(interpolants, "Interpolants", logFile);

	double intSum = 0;
	for (int i = 0; i < interpolants.size(); i++)
		intSum += interpolants[i];

	logFile << "Sum of interpolants: " << intSum << endl;

}

/*!****************************************************************************/
//! Coordinate generator
//! Points are defined along the columns of coords and their dimensional
//! coordinates along the rows
void Interpolation::MLSUnitTest_meshMultiDim(double& length, double& nodal_dist,
		int& ndim, dbMatrix& coords, int& numCoords, intVector SBM_dim,
		dbVector SBM_coeff, ofstream& logFile) {

	// Define the coordinate list in a specific dimension
	dbVector coordList;
	for (int i = 0; (i) * nodal_dist <= length; i++)
		coordList.push_back(i * nodal_dist);

	printVector(coordList, "coordList", logFile);

	dbMatrix x;
	x.push_back(coordList);
	if (ndim > 1) {
		for (int i = 1; i < ndim; i++) {
			int dimLookup = i + 1;
			int position = findIntVecPos(dimLookup, 0, SBM_dim.size(), SBM_dim);
			logFile << "Position: " << position << endl;

			dbMatrix xBlock = x;
			int xBlockSize = xBlock[0].size();
			x.clear();
			x.resize(i + 1);

			for (int j = 0; j < coordList.size(); j++) {

				// Select a particular coordinate
				double coordSelect = coordList[j];

				// Find if dimension needs to be factored
				dbMatrix x_add = xBlock;
				if (position != -1) {
					logFile << "Dim factored" << endl;

					x_add.push_back(
							dbVector(xBlockSize,
									SBM_coeff[position] * coordSelect));

					//dbMatrix x_add = [SBM_coeff*coordSelect*ones(1,xBlockSize);xBlock];
				} else {
					// Repeat the previous block for each selected coordinate
					//x_add = [coordSelect*ones(1,xBlockSize);xBlock];
					x_add.push_back(dbVector(xBlockSize, coordSelect));
				}

				printMatrix(x_add, "x_add", logFile);

				//x = [x x_add];
				for (int k = 0; k < x_add.size(); k++) {
					for (int l = 0; l < x_add[k].size(); l++) {
						x[k].push_back(x_add[k][l]);
					}
				}

				printMatrix(x, "x", logFile);

			}
		}
	}

	// Transposing the coordinate matrix
	coords.resize(x[0].size());
	for (int i = 0; i < coords.size(); i++) {
		coords[i].resize(ndim);
		for (int j = 0; j < ndim; j++)
			coords[i][j] = x[j][i];
	}

	numCoords = coords.size();
}

