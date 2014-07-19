#include "Interpolation.h"

using namespace std;

Interpolation::Interpolation(dbVector& iPoint, dbMatrix& coords,
		dbVector& radiusVec, dbVector& interpolants, InputFileData* InputData,
		ofstream& logFile) {

	int choice = InputData->getValue("MLSCalculationType");

	switch (choice) {
	case 1:
		MLSUnitTest(logFile);
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

/*!****************************************************************************/
//! Calculate the Moving Least Square interpolants
dbVector Interpolation::MLSCalc(dbVector& iPoint, dbMatrix& coords,
		dbVector& radiusVec, InputFileData* InputData, ofstream& logFile) {

#ifdef _MLSDebugMode_
	logFile << "********* In Interpolation::MLSCalc *********" << endl;
	printVector(iPoint, "iPoint", logFile);
	printMatrix(coords, "Coords", logFile);
#endif

	int numCoords = coords.size();
	int ndim = iPoint.size();

	bool activate_SBM;
	if (InputData->getValue("shiftedBasisAlgorithm") == 1)
		// Activate the shifted basis correction to enhance MLS stability
		activate_SBM = true;
	else
		activate_SBM = false;

	//! ========================================================================
	//! Search for the neighbours of the iPoint

//	intVector neighbours = findNeighbours(iPoint, coords, radiusVec,
//			logFile);
//	dbVector neighbours_db(neighbours.begin(), neighbours.end());
//	msg = "Neighbours"; printVector(neighbours_db, msg, logFile);

	// Selected all coordinates as neighbours
	intVector neighbours(coords.size());
	for (int i = 0; i < coords.size(); i++) {
		neighbours[i] = i;
	}

	//! ========================================================================
	//! Setup the basis polynomial, calculate the weight function
	//! and the B_matrix
	dbMatrix basisPolyMat, B_matrix;
	dbVector weightVec;
	int degPoly = InputData->getValue("MLSPolynomialDegree");
	for (int i = 0; i < neighbours.size(); i++) {

		// Weight function
		double weight = cubicSplineWgtCalc(iPoint, coords[neighbours[i]],
				radiusVec, logFile);
		weightVec.push_back(weight);

#ifdef _MLSDebugMode_
		printVector(radiusVec, "radiusVec", logFile);
#endif

		//Basis function
		dbVector nCoords(ndim);
		if (activate_SBM == true) {
			for (int j = 0; j < ndim; j++) {
				nCoords[j] = (coords[neighbours[i]][j] - iPoint[j])
						/ radiusVec[j];
			}

#ifdef _MLSDebugMode_
			printVector(nCoords, "Shifted Coords", logFile);
#endif

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
	dbVector iPointBasisVec;
	if (activate_SBM == true) {
		dbVector newiPoint(ndim, 0);
		iPointBasisVec = PascalBasisCalc(newiPoint, degPoly, logFile);
	} else {
		iPointBasisVec = PascalBasisCalc(iPoint, degPoly, logFile);
	}

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

//#ifdef _MLSDebugMode_
	printMatrix(inverseMomentMat, "Inverse Moment Matrix", logFile);
	printVector(iPointBasisVec, "iPointBasisVec", logFile);
	printVector(interpolants, "Interpolants", logFile);

	double intSum = 0;
	for (int i = 0; i < interpolants.size(); i++)
		intSum += interpolants[i];

	logFile << "Sum of interpolants: " << intSum << endl;
//#endif

	return interpolants;

}
;

/*!****************************************************************************/
//!
//dbVector Interpolation::MLSCalcCheck(dbVector& iPoint,dbMatrix& coords,
//		dbVector& radiusVec,InputFileData* InputData,ofstream& logFile){
//
//}
/*!****************************************************************************/
//! Function to setup the Pascal polynomial's basis
dbVector Interpolation::PascalBasisCalc(dbVector& coord, int& orderPoly,
		ofstream& logFile) {

	// Determine the no of dimensions
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
/*!****************************************************************************/
//! Debugging
void Interpolation::MLSUnitTest(ofstream& logFile) {

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

