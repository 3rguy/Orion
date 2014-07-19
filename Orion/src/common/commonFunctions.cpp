/*
 * commonFunctions.cpp
 *
 *  Created on: 03 Jul 2014
 *      Author: ritesh
 */

#include "commonFunctions.h"

/***********************************************************************/
void clearArray(blVector& v) {

	using namespace std;

	v.assign(v.size(), false);

}

//void clearArray(intVector& v) {
//
//	using namespace std;
//
//	v.assign(v.size(), 0);
//
//}

//void clearArray(intMatrix& T) {
//
//	using namespace std;
//
//	for (int i = 0; i < T.size(); i++)
//		T[i].assign(T[i].size(), 0);
//
//}

//void clearArray(intMatrix3& T) {
//
//	using namespace std;
//
//	for (int i = 0; i < T.size(); i++)
//		for (int j = 0; j < T[i].size(); j++)
//			T[i][j].assign(T[i][j].size(), 0);
//
//}

//void clearArray(dbVector& v) {
//
//	using namespace std;
//
//	v.assign(v.size(), 0);
//
//}

//void clearArray(dbMatrix& T) {
//
//	using namespace std;
//
//	for (int i = 0; i < T.size(); i++)
//		T[i].assign(T[i].size(), 0);
//
//}

//void clearArray(dbMatrix3& T) {
//
//	using namespace std;
//
//	for (int i = 0; i < T.size(); i++)
//		for (int j = 0; j < T[i].size(); j++)
//			T[i][j].assign(T[i][j].size(), 0);
//
//}

void clearArray(dbMatrix4& T) {

	using namespace std;

	for (int i = 0; i < T.size(); i++)
		for (int j = 0; j < T[i].size(); j++)
			for (int k = 0; k < T[i][j].size(); k++)
				T[i][j][k].assign(T[i][j][k].size(), 0);

}

void clearArray(dbMatrix5& T) {

	using namespace std;

	for (int i = 0; i < T.size(); i++)
		for (int j = 0; j < T[i].size(); j++)
			for (int k = 0; k < T[i][j].size(); k++)
				for (int l = 0; l < T[i][j][k].size(); l++)

					T[i][j][k][l].assign(T[i][j][k][l].size(), 0);

}

// set all entries of a map<string,bool> to false
void clearMap(std::map<std::string, bool>& data, std::ofstream& logFile) {

	using namespace std;

	// Loop over all map entries.
	for (map<string, bool>::iterator q = data.begin(); q != data.end(); ++q)

		q->second = false;

}

// set all entries of a map<string,bool> to false
void clearMap(std::map<std::string, bool>& data) {

	using namespace std;

	// Loop over all map entries.
	for (map<string, bool>::iterator q = data.begin(); q != data.end(); ++q)

		q->second = false;

}

/***********************************************************************/
// Quicksort a integer vector.
// left ... first element
// right ... last element
void sortIntVector(intVector& vec, int left, int right) {

	int i, j;
	int l, m;

	if (left < right) {
		m = vec[left];
		i = left;
		j = right;

		while (i <= j) {

			while (vec[i] < m)
				i++;
			while (vec[j] > m)
				j--;

			if (i <= j) {
				l = vec[i];
				vec[i] = vec[j];
				vec[j] = l;

				i++;
				j--;
			}
		}

		sortIntVector(vec, left, j);
		sortIntVector(vec, i, right);
	}
}

/***********************************************************************/
// Quicksort a double vector (increasing order).
// left = startIdx
// right = endIdx (note: not endIdx+1 !)
void sortDoubleVector(dbVector& vec, int left, int right) {

	int i, j;
	double l, m;

	if (left < right) {
		m = vec[left];
		i = left;
		j = right;

		while (i <= j) {

			while (vec[i] < m)
				i++;
			while (vec[j] > m)
				j--;

			if (i <= j) {
				l = vec[i];
				vec[i] = vec[j];
				vec[j] = l;

				i++;
				j--;
			}
		}

		sortDoubleVector(vec, left, j);
		sortDoubleVector(vec, i, right);
	}
}

/***********************************************************************/
// Quicksort a integer array together with its indices (increasing order).
void sortValuesIdx(intVector& values, intVector& idx, int left, int right) {

	int i, j, k;
	int l, m;

	if (left < right) {
		m = values[left];
		i = left;
		j = right;

		while (i <= j) {

			while (values[i] < m)
				i++;
			while (values[j] > m)
				j--;

			if (i <= j) {
				l = values[i];
				values[i] = values[j];
				values[j] = l;

				k = idx[i];
				idx[i] = idx[j];
				idx[j] = k;

				i++;
				j--;
			}
		}

		sortValuesIdx(values, idx, left, j);
		sortValuesIdx(values, idx, i, right);
	}
}

/***********************************************************************/
// Quicksort a double array together with its indices.
void sortValuesIdx(dbVector& values, intVector& idx, int left, int right) {

	int i, j, k;
	double l, m;

	if (left < right) {
		m = values[left];
		i = left;
		j = right;

		while (i <= j) {

			while (values[i] < m)
				i++;
			while (values[j] > m)
				j--;

			if (i <= j) {
				l = values[i];
				values[i] = values[j];
				values[j] = l;

				k = idx[i];
				idx[i] = idx[j];
				idx[j] = k;

				i++;
				j--;
			}
		}

		sortValuesIdx(values, idx, left, j);
		sortValuesIdx(values, idx, i, right);
	}
}

/***********************************************************************/
// Sort the particles in 3 dimensions.
void sortParticles(dbMatrix& ptcleCoords, intVector& idx, int sortingCoord,
		int numOfPlanes, std::ofstream& logFile) {

	using namespace std;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int sortingCoord1, sortingCoord2, sortingCoord3;

	if (sortingCoord == 1) {
		sortingCoord1 = 0;
		sortingCoord2 = 1;
		sortingCoord3 = 2;
	} else if (sortingCoord == 2) {
		sortingCoord1 = 1;
		sortingCoord2 = 0;
		sortingCoord3 = 2;
	} else if (sortingCoord == 3) {
		sortingCoord1 = 2;
		sortingCoord2 = 0;
		sortingCoord3 = 1;
	} else {
		logFile << "Particles can not be sorted with respect to coordinate\n "
				<< "direction '" << sortingCoord << "' in function "
				<< "'sortParticles'!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int particlesNum = ptcleCoords[0].size();

	// Determination of sorting parameters.
	int numOfPlaneElems, numOfRowElems;

	int cubicRoot = (int) ceil(pow((double) particlesNum, 1.0 / 3.0));
	int square = (int) pow((double) cubicRoot, 2);

	if (numOfPlanes > 0) {
		numOfRowElems = (int) ceil(sqrt((double) particlesNum / numOfPlanes));
		numOfPlaneElems = (int) pow((double) numOfRowElems, 2);
	} else {
		numOfPlanes = cubicRoot;
		numOfRowElems = cubicRoot;
		numOfPlaneElems = square;
	}

	intVector& ptclsIdx = idx;
	dbVector ptclsVec(particlesNum);

	if (idx.size() < particlesNum)
		ptclsIdx.resize(particlesNum);

	logFile << "######################################################" << endl;
	logFile << "************* particle sorting algorithm *************" << endl;
	logFile << "sorting coordinate 1: " << sortingCoord1 << endl;
	logFile << "sorting coordinate 2: " << sortingCoord2 << endl;
	logFile << "sorting coordinate 3: " << sortingCoord3 << endl;
	logFile << "number of particles: " << particlesNum << endl;
	logFile << "row elements = " << numOfRowElems << "(" << cubicRoot << ")"
			<< endl;
	logFile << "plane elements = " << numOfPlaneElems << "(" << square << ")"
			<< endl;
	logFile << "planes = " << numOfPlanes << "(" << cubicRoot << ")" << endl;

	// Loop over all particles and store their coordinates 1.
	for (int i = 0; i < particlesNum; i++) {
		ptclsIdx[i] = i;
		ptclsVec[i] = ptcleCoords[sortingCoord1][i];
	}

	// Sort the nodes according their coordinate 1.
	sortValuesIdx(ptclsVec, ptclsIdx, 0, particlesNum - 1);

#ifdef _commonDebugMode_
	logFile<<"************* sorted particles according "<<sortingCoord1
	<<"-direction ***********"<<endl;
	for(int i=0;i<particlesNum;i++)
	logFile<<"new: "<<i<<"; old: "<<ptclsIdx[i]
	<<" coord = "<<ptclsVec[i]<<endl;
#endif

	int planeStartElem = 0;
	int planeEndElem = -1;

	int rowStartElem, rowEndElem;

	// Loop over all planes.
	for (int i = 0; i < numOfPlanes; i++) {

		planeEndElem += numOfPlaneElems;

		// Store 2-direction coordinates of a portion of particles of the
		// current plane.
		if (planeEndElem > particlesNum - 1)
			planeEndElem = particlesNum - 1;

		// Store 2-direction coordinates of a portion of particles of the
		// current plane.
		for (int k = planeStartElem; k <= planeEndElem; k++)
			ptclsVec[k] = ptcleCoords[sortingCoord2][ptclsIdx[k]];

		// Sort the nodes according their coordinate 2.
		sortValuesIdx(ptclsVec, ptclsIdx, planeStartElem, planeEndElem);

#ifdef _commonDebugMode_
		logFile<<"************* sorted particles according 2"
		<<"-direction ***********"<<endl;
		for(int k=planeStartElem;k<=planeEndElem;k++)
		logFile<<"new: "<<k<<"; old: "<<ptclsIdx[k]
		<<" coord = "<<ptclsVec[k]<<endl;
#endif

		rowStartElem = planeStartElem;
		rowEndElem = planeStartElem - 1;

		// Loop over all rows within the current plane.
		for (int j = 0; j < numOfRowElems; j++) {

			rowEndElem += numOfRowElems;

			if (rowEndElem > planeEndElem)
				rowEndElem = planeEndElem;

			// Store 3-direction coordinates of a portion of particles of the
			// current row.
			for (int k = rowStartElem; k <= rowEndElem; k++)
				ptclsVec[k] = ptcleCoords[sortingCoord3][ptclsIdx[k]];

			// Sort the nodes according their coordinate 3.
			sortValuesIdx(ptclsVec, ptclsIdx, rowStartElem, rowEndElem);

#ifdef _commonDebugMode_
			logFile<<"************* sorted particles according 3"
			<<"-direction ***********"<<endl;
			for(int k=rowStartElem;k<=rowEndElem;k++)
			logFile<<"new: "<<k<<"; old: "<<ptclsIdx[k]
			<<" coord = "<<ptclsVec[k]<<endl;
#endif

			rowStartElem += numOfRowElems;
		}

		planeStartElem += numOfPlaneElems;
	}

	//*********************************************************************
	// Store the particles vector in a new proper order.
	dbMatrix oldPtcleCoords = ptcleCoords;

	for (int i = 0; i < particlesNum; i++) {
		ptcleCoords[0][i] = oldPtcleCoords[0][ptclsIdx[i]];
		ptcleCoords[1][i] = oldPtcleCoords[1][ptclsIdx[i]];
		ptcleCoords[2][i] = oldPtcleCoords[2][ptclsIdx[i]];
	}

}

/***********************************************************************/
/***********************************************************************/
// Search for a certain integer-vector entry and return the position.
int findIntVecPos(int& entry, int startIdx, int endIdx, intVector& vec) {

	int position = -1;

	for (int i = startIdx; i < endIdx; i++)

		if (entry == vec[i]) {
			position = i;
			break;
		}

	return position;

}

/***********************************************************************/
/***********************************************************************/
// Search for a certain integer-vector entry and return the position.
int findDoubleVecPos(double& entry, int startIdx, int endIdx, dbVector& vec) {

	int position = -1;

	for (int i = startIdx; i < endIdx; i++)

		if (entry == vec[i]) {
			position = i;
			break;
		}

	return position;

}

/***********************************************************************/
// Search in a certain matrix column for a certain integer-matrix entry
// and return the position.
int findIntMatPos(int& entry, int startIdx, int endIdx, int vecPos,
		intMatrix& mat) {

	int position = -1;

	for (int i = startIdx; i < endIdx; i++)

		if (entry == mat[i][vecPos]) {
			position = i;
			break;
		}

	return position;

}

/***********************************************************************/
// Search for a certain row in an integer matrix and return the position.
int findIntMatPos(intVector& entries, int startIdx, int endIdx,
		std::string& mode, intMatrix& mat) {

	using namespace std;

	int position = -1;

	if (mode.compare("arbitrary") == 0) {

		intVector matEntries;
		intVector vecEntries = entries;

		sortIntVector(vecEntries, 0, vecEntries.size() - 1);

		for (int i = startIdx; i < endIdx; i++) {

			matEntries = mat[i];
			sortIntVector(matEntries, 0, matEntries.size() - 1);

			if (vecEntries.size() > matEntries.size())
				break;

			for (int j = 0; j < matEntries.size(); j++) {

				if (vecEntries[j] != matEntries[j])
					break;

				else if (j + 1 == vecEntries.size()) {
					position = i;
					break;
				}

			}

		}

	}

	else if (mode.compare("strict") == 0) {

		for (int i = startIdx; i < endIdx; i++) {

			if (entries.size() > mat[i].size())
				break;

			for (int j = 0; j < mat[i].size(); j++) {

				if (entries[j] != mat[i][j])
					break;

				else if (j + 1 == entries.size()) {
					position = i;
					break;
				}

			}

		}

	}

	else {
		cerr << "In commonFunctions::findIntMatPos mode '" << mode
				<< "' is not " << "supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return position;

}

/***********************************************************************/
// Search for a certain row in double matrix and return the position.
int findDoubleMatPos(dbVector& entries, int startIdx, int endIdx,
		std::string mode, double tol, dbMatrix& mat, std::ofstream& logFile) {

	using namespace std;

	int position = -1;

#ifdef _commonDebugMode_
	logFile<<"****************** findDoubleMatPos ******************"<<endl;
	for(int i=0;i<entries.size();i++)
	logFile<<"vec["<<i<<"] = "<<entries[i]<<endl;
	logFile<<"------------------------------------------------------"<<endl;
	for(int i=0;i<mat.size();i++) {
		for(int j=0;j<mat[i].size();j++)
		logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i][j]<<endl;
		logFile<<"---------"<<endl;
	}
#endif

	if (mode.compare("strict") == 0) {

		// loop over all matrix's rows
		for (int i = startIdx; i < endIdx; i++) {

			if (entries.size() > mat[i].size())
				break;

			// loop over all entries of the matrix current row
			for (int j = 0; j < mat[i].size(); j++) {

				// not matching entries
				if (fabs(entries[j] - mat[i][j]) > fabs(tol))
					break;

				// all values are found
				else if (j + 1 == entries.size()) {
					position = i;
					i = endIdx;
					break;
				}

			}

		}

	}

	else {
		cerr << "In commonFunctions::findDoubleMatPos mode '" << mode
				<< "' is not " << "supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return position;

}

/***********************************************************************/
// Compare two integer vectors.
bool compareIntVecs(std::string& mode, intVector& small, intVector& big,
		std::ofstream& logFile) {

	using namespace std;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	bool flag = false;

	if (mode.compare("arbitrary-subvector") == 0) {

		intVector vecBig = big;
		intVector vecSmall = small;

		if (vecBig.size() < vecSmall.size()) {
			cerr << "In commonFunctions::compareIntVecs not enough entries in "
					<< "vector (rank = " << rank << ")!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		sortIntVector(vecBig, 0, vecBig.size() - 1);
		sortIntVector(vecSmall, 0, vecSmall.size() - 1);

#ifdef _commonDebugMode_
		logFile<<"***************** input vectors ********************"<<endl;
		for(int i=0;i<vecSmall.size();i++)
		logFile<<"vec1["<<i<<"] = "<<vecSmall[i]<<endl;
		logFile<<"----------------------------------------------------"<<endl;
		for(int i=0;i<vecBig.size();i++)
		logFile<<"vec2["<<i<<"] = "<<vecBig[i]<<endl;
		logFile<<"**************** sorted vectors ********************"<<endl;
		for(int i=0;i<vecSmall.size();i++)
		logFile<<"vec1["<<i<<"] = "<<vecSmall[i]<<endl;
		logFile<<"----------------------------------------------------"<<endl;
		for(int i=0;i<vecBig.size();i++)
		logFile<<"vec2["<<i<<"] = "<<vecBig[i]<<endl;
		logFile<<"****************************************************"<<endl;
#endif

		int neededHits = vecSmall.size();

		for (int i = vecBig.size() - 1; i >= 0; i--) {

			for (int j = 0; j < vecSmall.size(); j++) {

#ifdef _commonDebugMode_
				logFile<<"small["<<j<<"] = "<<vecSmall[j]<<"; big["<<i<<"] = "
				<<vecBig[i]<<"; neededHits = "<<neededHits<<endl;
#endif

				// match
				if (vecSmall[j] == vecBig[i]) {
					vecSmall.erase(vecSmall.begin() + j);
					--j;
					--neededHits;

					// all found
					if (neededHits == 0) {
						flag = true;
						j = vecSmall.size();
						i = 0;
					}
					// not enough entries left to satisfy the needed hits
					else if (i < neededHits) {
						j = vecSmall.size();
						i = 0;
					}

					break;
				}

			}

#ifdef _commonDebugMode_
			logFile<<"----------------------------------------------------"<<endl;
#endif

		}

	} else {
		cerr << "In commonFunctions::compareIntVecs mode '" << mode
				<< "' is not " << "supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return flag;
}

/***********************************************************************/
// remove redundant entries of a integer vector
void removeRedundantEntries(std::vector<int>& vec, int startIdx, int& endIdx,
		std::ofstream& logFile) {

	using namespace std;

	vector<int> tmp(endIdx - startIdx);

	if (tmp.size() > 0) {

		for (int i = 0, j = startIdx; j < endIdx; i++, j++)

			tmp[i] = vec[j];

		sortIntVector(tmp, 0, tmp.size() - 1);

		int m = startIdx + 1;
		vec[startIdx] = tmp[0];

		for (int i = 1; i < tmp.size(); i++)

			// skip already stored elements
			if (tmp[i - 1] != tmp[i]) {

				vec[m] = tmp[i];
				m++;

			}

		endIdx = m;
	}

}

/***********************************************************************/
// remove redundant entries in a certain column of a matrix
void removeRedundantEntries(std::vector<std::vector<int> >& mat, int idx) {

	using namespace std;

	int pos;

	// loop over rows
	for (int i = 1; i < mat.size(); i++) {

		pos = 1;

		while (pos != -1) {

			// Search for a certain integer-matrix entry and return the position.
			pos = findIntMatPos(mat[i][idx], 0, i, idx, mat);

			if (pos != -1)
				mat.erase(mat.begin() + i);

		}

	}

}

/**********************************************************************/
// find max value in a vector
int findMaxValue(std::string& mode, std::vector<int>& vec) {

	using namespace std;

	if (vec.size() == 0) {
		cerr << "In commonFunctions::findMaxValue vector is not initialized!"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int idx = 0;

	for (int i = 1; i < vec.size(); i++) {

		if (vec[i] > vec[i - 1])
			idx = i;

	}

	return idx;

}

/**********************************************************************/
// find min value in a vector
int findMinValue(std::vector<int>& vec) {

	using namespace std;

	if (vec.size() == 0) {
		cerr << "In commonFunctions::findMinValue vector is not initialized!"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int idx = 0;

	for (int i = 1; i < vec.size(); i++) {

		if (vec[i] < vec[i - 1])
			idx = i;

	}

	return idx;

}

/***********************************************************************/
/***********************************************************************/
// Calculate scalar product of two vector.
void scalarProduct(dbVector& a, dbVector& b, double& scalar,
		std::ofstream& logFile) {

	using namespace std;

	scalar = 0;

	if (a.size() != b.size()) {
		logFile << "Not compatible vectors for vector scalar product!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	scalar = 0;

	for (int i = 0; i < a.size(); i++)
		scalar += a[i] * b[i];

}

/***********************************************************************/
// Calculate scalar product of two matrices or tensors.
void scalarProduct(dbMatrix& A, dbMatrix& B, double& scalar, bool aTrans,
		bool bTrans, std::ofstream& logFile) {

	using namespace std;

	scalar = 0;

	if ((!aTrans && !bTrans) || (aTrans && bTrans)) {

		if (A.size() != B.size() || A[0].size() != B[0].size()) {
			logFile << "Not compatible matrices for matrix scalar product!"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for (int i = 0; i < A.size(); i++)

			for (int j = 0; j < A[0].size(); j++)

				scalar += A[i][j] * B[i][j];

	} else {

		if (A.size() != B[0].size() || A[0].size() != B.size()) {
			logFile << "Not compatible matrices for matrix scalar product!"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for (int i = 0; i < A.size(); i++)

			for (int j = 0; j < A[0].size(); j++)

				scalar += A[i][j] * B[j][i];

	}

}

/***********************************************************************/
// Calculate inner product of matrix or tensor T and vector u.
void innerTensorProduct(dbMatrix& T, dbVector& u, dbVector& result, bool tTrans,
		std::ofstream& logFile) {

	using namespace std;

	if (result.size() != 0)
		clearArray(result);

	if (T.size() == 0 || u.size() == 0) {
		logFile << "Not compatible matrix for inner product!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (!tTrans) {

		if (T[0].size() != u.size()) {
			logFile << "Not compatible matrix for inner product!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (result.size() < T.size())
			result.resize(T.size());

		for (int i = 0; i < T.size(); i++)

			for (int j = 0; j < u.size(); j++)
				result[i] += T[i][j] * u[j];
	} else if (tTrans) {

		if (T.size() != u.size()) {
			logFile << "Not compatible matrix for inner product!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (result.size() < T[0].size())
			result.resize(T[0].size());

		for (int i = 0; i < T[0].size(); i++)

			for (int j = 0; j < u.size(); j++)
				result[i] += T[j][i] * u[j];

	}

}

/***********************************************************************/
// Calculate inner product of matrix or tensor T and tensor S.
void innerTensorProduct(dbMatrix& T, dbMatrix& S, dbMatrix& Result, bool tTrans,
		bool sTrans, std::ofstream& logFile) {

	using namespace std;

	if (Result.size() != 0)
		clearArray(Result);

	if (T.size() == 0 || S.size() == 0) {
		logFile << "Not compatible matrices for inner product!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (!tTrans && !sTrans) {

		if (T[0].size() != S.size()) {
			logFile << "Not compatible matrices for inner product!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (Result.size() < T.size()) {
			Result.resize(T.size());

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S[0].size());

		}

		if (Result[0].size() < S[0].size())

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S[0].size());

		for (int i = 0; i < T.size(); i++)

			for (int j = 0; j < S[0].size(); j++)

				for (int k = 0; k < T[0].size(); k++)

					Result[i][j] += T[i][k] * S[k][j];
	} else if (!tTrans && sTrans) {

		if (T[0].size() != S[0].size()) {
			logFile << "Not compatible matrices for inner product!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (Result.size() < T.size()) {
			Result.resize(T.size());

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S.size());

		}

		if (Result[0].size() < S.size())

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S.size());

		for (int i = 0; i < T.size(); i++)

			for (int j = 0; j < S.size(); j++)

				for (int k = 0; k < T[0].size(); k++)

					Result[i][j] += T[i][k] * S[j][k];
	} else if (tTrans && !sTrans) {

		if (T.size() != S.size()) {
			logFile << "Not compatible matrices for inner product!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (Result.size() < T[0].size()) {
			Result.resize(T[0].size());

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S[0].size());

		}

		if (Result[0].size() < S[0].size())

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S[0].size());

		for (int i = 0; i < T[0].size(); i++)

			for (int j = 0; j < S[0].size(); j++)

				for (int k = 0; k < T.size(); k++)

					Result[i][j] += T[k][i] * S[k][j];
	} else if (tTrans && sTrans) {

		if (T.size() != S[0].size()) {
			logFile << "Not compatible matrices for inner product!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (Result.size() < T[0].size()) {
			Result.resize(T[0].size());

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S.size());

		}

		if (Result[0].size() < S.size())

			for (int i = 0; i < Result.size(); i++)
				Result[i].resize(S.size());

		for (int i = 0; i < T[0].size(); i++)

			for (int j = 0; j < S.size(); j++)

				for (int k = 0; k < T.size(); k++)

					Result[i][j] += T[k][i] * S[j][k];
	}

}

/***********************************************************************/
// Calculate cross product of two vectors u and v .
void crossProduct(dbVector& u, dbVector& v, dbVector& result, double& det,
		std::ofstream& logFile) {

	using namespace std;

	if (u.size() != v.size()) {
		logFile << "Not usable vectors for vector cross product!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (result.size() != 0)
		clearArray(result);

	if (result.size() < u.size())
		result.resize(u.size());

	// Get indices of permutations.
	intMatrix permutations = getPermutations(u.size());

	// Loop over all permutations unequal zero.
	for (int i = 0; i < 2 * u.size(); i++) {
		result[permutations[i][0]] += det * u[permutations[i][1]]
				* v[permutations[i][2]] * permutations[i][3];
	}

}

/***********************************************************************/
// Calculate cross product of two vectors u and v .
void crossProduct(dbVector& u, dbVector& v, dbVector& result) {

	using namespace std;

	if (u.size() != v.size()) {
		cerr << "Not usable vectors for vector cross product!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (result.size() != 0)
		clearArray(result);

	if (result.size() < u.size())
		result.resize(u.size());

	// Get indices of permutations.
	intMatrix permutations = getPermutations(u.size());

	// Loop over all permutations unequal zero.*/
	for (int i = 0; i < 2 * u.size(); i++) {
		result[permutations[i][0]] += u[permutations[i][1]]
				* v[permutations[i][2]] * permutations[i][3];
	}

}

/***********************************************************************/
// Set the Kronecker symbol.
dbMatrix getKroneckerSymbol(int size) {

	using namespace std;

	dbMatrix delta(size, dbVector(size));

	for (int i = 0; i < size; i++)
		delta[i][i] = 1.0;

	return delta;
}

/***********************************************************************/
// Calculate the permutations unequal zero if the size of used arrays
// is 'size'.
intMatrix getPermutations(int size) {

	using namespace std;

	intMatrix permutations(2 * size, intVector(4));

	// Set all odd permutations.
	int counter;

	for (int i = 0; i < size; i++) {

		counter = i;

		// Set indices i,j,k.
		for (int j = 0; j < 3; j++) {
			permutations[i][j] = counter;

			if (++counter == size)
				counter = 0;
		}

		permutations[i][3] = 1;
	}

	// Set all even permutations.
	for (int i = size; i < 2 * size; i++) {

		counter = 2 * size - 1 - i;

		// Set indices i,j,k.
		for (int j = 0; j < 3; j++) {
			permutations[i][j] = counter;

			if (--counter < 0)
				counter = size - 1;
		}

		permutations[i][3] = -1;
	}

	return permutations;
}

/***********************************************************************/
// Determine the missing indices which build permutations for
// a given permutation index.
intMatrix3 getPermutationMissingIdx(int dim, int pos) {

	using namespace std;

	int i, j, ipos, jpos;
	intMatrix3 m;

	pos -= 1;

	switch (dim) {

	case 3:

		ipos = pos - 1;
		jpos = pos + 1;

		if (ipos < 0)
			ipos = dim - 1;

		if (jpos == dim)
			jpos = 0;

		allocateArray(m, dim, 2, 3);

		for (int p = 0; p < dim; p++) {

			i = p - 1;
			j = p + 1;

			if (i < 0)
				i = dim - 1;

			if (j == dim)
				j = 0;

			// positive permutation

			m[p][0][ipos] = i;
			m[p][0][jpos] = j;
			m[p][0][2] = 1;

			// negative permutation
			m[p][1][ipos] = j;
			m[p][1][jpos] = i;
			m[p][1][2] = -1;

		}

		break;

	default:
		cerr << "In commonFunctions::getPermutationMissingIdx dimension "
				<< "not\n supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return m;
}

double getPi() {

	return ((double) 3.14159265358979323846);
}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a single precision matrix.
int calcDetSingle(dbMatrix& Amat, double& det, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetSingle no valid matrix for\n"
				<< "calculating determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	det = 0;

	/*********************************************************************/
	/* Intialisation of the variables used by the external subroutine.*/
	int info = 0;
	int lda = Amat.size();
	int n = Amat[0].size();
	int* ipvt = new int[Amat.size() + 1];

	float* mat = new float[Amat.size() * Amat.size() + 1];

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

#ifdef _commonDebugMode_
	logFile<<"********** not factorized FORTRAN matrix *********"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine sgefa factorizing Amat by
	// gaussian elimination.
	sgefa_(mat, lda, n, ipvt, info);

#ifdef _commonDebugMode_
	logFile<<"******************** factorized matrix **************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
	logFile<<"info "<<info<<endl;
	logFile<<"*********************** pivots **********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	logFile<<"pivot["<<i<<"] = "<<ipvt[i]<<endl;
#endif

	if (info != 0) {
		logFile << "In commonFunctions::calcDetSingle the factorization of\n"
				<< "input matrix could not be calculated!" << endl;
		cerr << "In commonFunctions::calcDetSingle the factorization of\n"
				<< "input matrix could not be calculated!" << endl;
	}

	/*********************************************************************/
	// Call of external FORTRAN subroutine sgedi calculating the
	// determinant of Amat.
	int job = 10;
	float* work = new float[Amat.size() + 1];
	float* result = new float[3];

	sgedi_(mat, lda, n, ipvt, result, work, job);

#ifdef _commonDebugMode_
	logFile<<"************************ results *********************"<<endl;
	logFile<<result[0]<<" "<<result[1]<<" "<<result[2]<<endl;
#endif

	if (fabs(result[0]) >= 1.0 && fabs(result[0]) < 10.0 && result[1] <= 128)

		det = (double) result[0] * pow(10, (double) result[1]);
	else if (fabs(result[0]) >= 1.0 && fabs(result[0]) < 10.0
			&& result[1] > 128)
		det = 10.0e128;
	else
		det = 0.0;

	delete[] ipvt, mat, work, result;

#ifdef _commonDebugMode_
	logFile<<"********************* determinate ********************"<<endl;
	logFile<<"det = "<<det<<endl;
#endif

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a square double precision matrix.
int calcDetDouble(dbMatrix& Amat, double& det) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		cout << "In commonFunctions::calcDetDouble no valid matrix for\n"
				<< "calculating determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Intialisation of the variables used by the external subroutine.
	int info = 0;
	int lda = Amat.size();
	int n = Amat[0].size();
	double rcond;

	//   double* mat = new double[Amat.size()*Amat.size()+1];
	//   double* work = new double[Amat.size()+1];
	//   int* ipvt = new int[Amat.size()+1];
	//   double* determinant = new double[3];
	dbVector mat(Amat.size() * Amat.size());
	dbVector work(Amat.size());
	intVector ipvt(Amat.size());
	dbVector determinant(3);

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgeco factorizing Amat by
	// gaussian elimination.
	// rcond < dbl_epsilon means matrix is singular
	dgeco_(&mat[0], lda, n, &ipvt[0], rcond, &work[0]);

	if (fabs(rcond) < DBL_EPSILON) {
		cout << "In commonFunctions::calcDetDouble the factorization of\n"
				<< "provides a singular matrix!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Call of external FORTRAN subroutine sgedi calculating the
	// determinant of Amat.
	int job = 10;

	dgedi_(&mat[0], lda, n, &ipvt[0], &determinant[0], &work[0], job);

	det = determinant[0] * pow(10.0, determinant[1]);

	//delete[] ipvt,mat,work,determinant;

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a square double precision matrix.
int calcDetDouble(dbMatrix& Amat, double& det, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetDouble no valid matrix for\n"
				<< "calculating determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Intialisation of the variables used by the external subroutine.
	int info = 0;
	int lda = Amat.size();
	int n = Amat[0].size();
	double rcond;

	//   double* mat = new double[Amat.size()*Amat.size()+1];
	//   double* work = new double[Amat.size()+1];
	//   int* ipvt = new int[Amat.size()+1];
	//   double* determinant = new double[3];
	dbVector mat(Amat.size() * Amat.size());
	dbVector work(Amat.size());
	intVector ipvt(Amat.size());
	dbVector determinant(3);

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

#ifdef _commonDebugMode_
	logFile<<"********** not factorized FORTRAN matrix *********"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgeco factorizing Amat by
	// gaussian elimination.
	// rcond < dbl_epsilon means matrix is singular
	dgeco_(&mat[0], lda, n, &ipvt[0], rcond, &work[0]);

#ifdef _commonDebugMode_
	logFile<<"******************** factorized matrix **************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
	logFile<<"rcond "<<rcond<<endl;
	logFile<<"*********************** pivots **********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	logFile<<"pivot["<<i<<"] = "<<ipvt[i]<<endl;
#endif

	if (fabs(rcond) < DBL_EPSILON) {
		logFile << "In commonFunctions::calcDetDouble the factorization of\n"
				<< "provides a singular matrix!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgedi calculating the
	// determinant of Amat.
	int job = 10;

	dgedi_(&mat[0], lda, n, &ipvt[0], &determinant[0], &work[0], job);

	det = determinant[0] * pow(10.0, determinant[1]);

#ifdef _commonDebugMode_
	double linpackDet = det;
	logFile<<"********************* determinate ********************"<<endl;
	logFile<<"det[0]="<<determinant[0]<<" det[1]="<<determinant[1]<<endl;
	logFile<<"detLinpack = "<<det<<endl;
	//   calcDetDoubleSmall(Amat,det,logFile);
	//   logFile<<"detSmall = "<<det<<endl;
	//   calcDetDoubleTens(Amat,det,logFile);
	//   logFile<<"detTens = "<<det<<endl;
	//   det = linpackDet;
#endif

	//delete[] ipvt,mat,work,determinant;

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a small double precision matrix.
void calcDetDoubleSmall(dbMatrix& Amat, double& det, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetDouble no valid matrix for\n"
				<< "calculating determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int col1, col2;
	double product1, product2;

	det = 0.0;

	for (int i = 0; i < Amat.size(); i++) {

		col1 = i;
		col2 = Amat.size() - 1 - i;
		product1 = product2 = 1;

		for (int r = 0; r < Amat.size(); r++) {

			// even permutations
			if (col1 == Amat.size())
				col1 = 0;

			product1 *= Amat[r][col1];
			col1++;

			// odd permutations
			if (col2 < 0)
				col2 = Amat.size() - 1;

			product2 *= Amat[r][col2];
			col2--;

#ifdef _commonDebugMode_
			logFile<<"product1 = "<<product1<<endl;
			logFile<<"product2 = "<<product2<<endl;
#endif
		}

		det += product1;
		det -= product2;

	}

#ifdef _commonDebugMode_
	logFile<<"********************* determinate ********************"<<endl;
	logFile<<"det = "<<det<<endl;
#endif

}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a double precision tensor.
void calcDetDoubleTens(dbMatrix& Amat, double& det, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetDoubleTens no valid matrix for\n"
				<< "calculating determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	intMatrix e = getPermutations(Amat.size());
	det = 0;

	// det A = 1/6 e_{irl} e_{jsm} A_{ij} A_{rs} A_{lm}

	// Loop over all six possible permutations.
	for (int l = 0; l < e.size(); l++)

		// Loop over all six possible permutations.
		for (int m = 0; m < e.size(); m++)

			det += e[l][3] * e[m][3] * Amat[e[l][0]][e[m][0]]
					* Amat[e[l][1]][e[m][1]] * Amat[e[l][2]][e[m][2]];

	det *= 1.0 / 6.0;

#ifdef _commonDebugMode_
	logFile<<"********************* permutations *******************"<<endl;
	for(int m=0;m<e.size();m++) {
		logFile<<m+1<<".) ";
		for(int n=0;n<e[m].size();n++)
		logFile<<e[m][n]<<" ";
		logFile<<endl;
	}
	logFile<<"********************* determinate ********************"<<endl;
	logFile<<"det = "<<det<<endl;
#endif

}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a dense square double precision matrix
// making use of the Cayley-Hamilton theorem
int calcDetDoubleDense(dbMatrix& Amat, double& det, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetDoubleDense no valid matrix\n"
				<< "to calculate its determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int info = 0;
	int Asize = Amat.size();

	if (Asize == 3) {

		// compute the invariants

		// tr A = I1
		double trA = 0;

		for (int i = 0; i < Asize; i++)
			trA += Amat[i][i];

		// tr A^2 = A_ik A_ki
		double trA2 = 0;

		for (int i = 0; i < Asize; i++)

			for (int k = 0; k < Asize; k++)

				trA2 += Amat[i][k] * Amat[k][i];

		// trA^3 = A_ik*A_kl*A_li
		double trA3 = 0;

		for (int i = 0; i < Asize; i++)

			for (int k = 0; k < Asize; k++)

				for (int l = 0; l < Asize; l++)

					trA3 += Amat[i][k] * Amat[k][l] * Amat[l][i];

		// I3 = 1/3*(trA^3 - 3/2*trA*trA^2 + 1/2*(trA)^3)

		det = 1.0 / 3.0 * (trA3 - 3.0 / 2.0 * trA * trA2 + 0.5 * pow(trA, 3.0));

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		int flag;
		double lapackDet = 0;
		if(computeNorm(Amat,3,logFile) != 0)
		flag = calcDetDoubleSparse(Amat,lapackDet,logFile);
		logFile<<"********************* determinate ********************"<<endl;
		logFile<<"det = "<<det<<" ?= "<<lapackDet<<endl;
#endif

		return info;

	} else

		return (calcDetDoubleSparse(Amat, det, logFile));

}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a dense square double precision matrix
// making use of the Cayley-Hamilton theorem
int calcDetDoubleDense(dbMatrix& Amat, double& det) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		cout << "In commonFunctions::calcDetDoubleDense no valid matrix\n"
				<< "to calculate its determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int info = 0;
	int Asize = Amat.size();

	if (Asize == 3) {

		// compute the invariants

		// tr A = I1
		double trA = 0;

		for (int i = 0; i < Asize; i++)
			trA += Amat[i][i];

		// tr A^2 = A_ik A_ki
		double trA2 = 0;

		for (int i = 0; i < Asize; i++)

			for (int k = 0; k < Asize; k++)

				trA2 += Amat[i][k] * Amat[k][i];

		// trA^3 = A_ik*A_kl*A_li
		double trA3 = 0;

		for (int i = 0; i < Asize; i++)

			for (int k = 0; k < Asize; k++)

				for (int l = 0; l < Asize; l++)

					trA3 += Amat[i][k] * Amat[k][l] * Amat[l][i];

		// I3 = 1/3*(trA^3 - 3/2*trA*trA^2 + 1/2*(trA)^3)

		det = 1.0 / 3.0 * (trA3 - 3.0 / 2.0 * trA * trA2 + 0.5 * pow(trA, 3.0));

		return info;

	} else

		return (calcDetDoubleSparse(Amat, det));

}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a sparse square double precision matrix
// using LAPACK routines.
int calcDetDoubleSparse(dbMatrix& Amat, double& det, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetDoubleSparse no valid matrix\n"
				<< "to calculate its determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int info = 0;

	// Intialisation of the variables used by the external subroutine.
	int lda = Amat.size();
	int n = Amat[0].size();
	double rcond;

	//   double* mat = new double[Amat.size()*Amat.size()+1];
	//   double* work = new double[Amat.size()+1];
	//   int* ipvt = new int[Amat.size()+1];
	//   double* determinant = new double[3];
	dbVector mat(Amat.size() * Amat.size());
	dbVector work(Amat.size());
	intVector ipvt(Amat.size());
	dbVector determinant(3);

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

#ifdef _commonDebugMode_
	logFile<<"********** not factorized FORTRAN matrix *********"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgeco factorizing Amat by
	// gaussian elimination.
	// rcond < dbl_epsilon means matrix is singular
	dgeco_(&mat[0], lda, n, &ipvt[0], rcond, &work[0]);

#ifdef _commonDebugMode_
	logFile<<"******************** factorized matrix **************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
	logFile<<"rcond "<<rcond<<endl;
	logFile<<"*********************** pivots **********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	logFile<<"pivot["<<i<<"] = "<<ipvt[i]<<endl;
#endif

	if (fabs(rcond) < DBL_EPSILON) {
		logFile << "In commonFunctions::calcDetDoubleSparse the factorization\n"
				<< "provides a singular matrix!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgedi calculating the
	// determinant of Amat.
	int job = 10;

	dgedi_(&mat[0], lda, n, &ipvt[0], &determinant[0], &work[0], job);

	det = determinant[0] * pow(10.0, determinant[1]);

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
	logFile<<"********************* determinate ********************"<<endl;
	logFile<<"det = "<<det<<endl;
#endif

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Calculate the determinante of a sparse square double precision matrix
// using LAPACK routines.
int calcDetDoubleSparse(dbMatrix& Amat, double& det) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		cout << "In commonFunctions::calcDetDoubleSparse no valid matrix\n"
				<< "to calculate its determinant!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int info = 0;

	// Intialisation of the variables used by the external subroutine.
	int lda = Amat.size();
	int n = Amat[0].size();
	double rcond;

	//   double* mat = new double[Amat.size()*Amat.size()+1];
	//   double* work = new double[Amat.size()+1];
	//   int* ipvt = new int[Amat.size()+1];
	//   double* determinant = new double[3];
	dbVector mat(Amat.size() * Amat.size());
	dbVector work(Amat.size());
	intVector ipvt(Amat.size());
	dbVector determinant(3);

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgeco factorizing Amat by
	// gaussian elimination.
	// rcond < dbl_epsilon means matrix is singular
	dgeco_(&mat[0], lda, n, &ipvt[0], rcond, &work[0]);

	if (fabs(rcond) < DBL_EPSILON) {
		cout << "In commonFunctions::calcDetDoubleSparse the factorization\n"
				<< "provides a singular matrix!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Call of external FORTRAN subroutine dgedi calculating the
	// determinant of Amat.
	int job = 10;

	dgedi_(&mat[0], lda, n, &ipvt[0], &determinant[0], &work[0], job);

	det = determinant[0] * pow(10.0, determinant[1]);

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Calculate the single precision inverse of a matrix.
void calcInvSingle(dbMatrix& Amat, dbMatrix& inverse, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetSingle no valid matrix for\n"
				<< "calculating inverse!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Intialisation of the variables used by the external subroutine.
	float det = 0;
	int info = 0;
	int lda = Amat.size();
	int n = Amat[0].size();

	int* ipvt = new int[Amat.size()];
	float* mat = new float[Amat.size() * Amat.size()];

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

#ifdef _commonDebugMode_
	logFile<<"********** not factorized FORTRAN matrix *********"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine sgefa factorizing Amat by
	// gaussian elimination.
	sgefa_(mat, lda, n, ipvt, info);

	if (info != 0) {
		logFile << "In commonFunctions::calcInvSingle the factorization of\n"
				<< "input matrix could not be calculated!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _commonDebugMode_
	logFile<<"******************** factorized matrix **************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
	logFile<<"info "<<info<<endl;
	logFile<<"*********************** pivots **********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	logFile<<"pivot["<<i<<"] = "<<ipvt[i]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine sgedi calculating the
	// inverse of Amat.
	int job = 01;
	float* work = new float[Amat.size()];
	float dummy[2];

	sgedi_(mat, lda, n, ipvt, &dummy[0], work, job);

#ifdef _commonDebugMode_
	logFile<<"**************** fortran inverse *********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<" = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Check the pre-allocated storage.
	if (inverse.size() < Amat.size()) {
		inverse.resize(Amat.size());

		for (int i = 0; i < Amat.size(); i++)
			inverse[i].resize(Amat.size());

	}

	// Store the calculated inverse in a C++ array.
	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			inverse[i][j] = mat[i * Amat.size() + j];

	delete[] ipvt, mat, work;
}

/***********************************************************************/
/***********************************************************************/
// Calculate the double precision inverse of a matrix.
void calcInvDouble(dbMatrix& Amat, dbMatrix& inverse, std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcDetSingle no valid matrix for\n"
				<< "calculating inverse!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//   if(Amat.size() == 3)
	//     calcInv3x3(Amat,inverse,logFile);

	//   else {

	/*********************************************************************/
	// Intialisation of the variables used by the external subroutine.
	int info = 0;
	int lda = Amat.size();
	int n = Amat[0].size();

	int lwork = 3 * Amat.size();

	//   double* work = new double[lwork+1];
	//   double* mat = new double[lda*n+1];
	//   int* ipvt = new int[Amat.size()+1];
	dbVector work(lwork);
	dbVector mat(lda * n);
	intVector ipvt(Amat.size());

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

#ifdef _commonDebugMode_
	logFile<<"********** not inverted FORTRAN matrix *************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine 'dgetrf' factorizing Amat by
	// gaussian elimination.
	dgetrf_(lda, n, &mat[0], lda, &ipvt[0], info);

	if (info != 0) {
		logFile << "In commonFunctions::calcInvDouble the factorization of\n"
				<< "input matrix could not be calculated!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _commonDebugMode_
	logFile<<"*********************** pivots **********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	logFile<<"pivot["<<i<<"] = "<<ipvt[i]<<endl;
	logFile<<"**************** factorized matrix ******************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine 'dgetri' to calculate the
	// the inverse of matrix Amat by
	dgetri_(n, &mat[0], lda, &ipvt[0], &work[0], lwork, info);

#ifdef _commonDebugMode_
	logFile<<"**************** fortran inverse *********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
	dbMatrix one(Amat.size(),dbVector(Amat.size()));
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	for(int k=0;k<Amat.size();k++)
	one[i][j] += Amat[i][k]*mat[k*Amat.size()+j];
	logFile<<"************** Kronecker delta(?) ********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"delta["<<i<<"]["<<j<<"] = "<<one[i][j]<<endl;
#endif

	if (info != 0) {
		logFile << "In commonFunctions::calcInvDouble the inverse could not \n"
				<< "be calculated!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Check the pre-allocated storage.
	if (inverse.size() < Amat.size()) {
		inverse.resize(Amat.size());

		for (int i = 0; i < Amat.size(); i++)
			inverse[i].resize(Amat.size());

	}

	// Store the calculated inverse in a C++ array.
	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			inverse[i][j] = mat[i * Amat.size() + j];

	//delete[] ipvt,mat,work;

}

/************************************************************************/
// Invert 3x3 matrix without using fortran.
void calcInv3x3(dbMatrix& A, dbMatrix& inv, std::ofstream& logFile) {

	using namespace std;

	double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
			+ A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2])
			+ A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
	allocateArray(inv, 3, 3);
	inv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
	inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / det;
	inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
	inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / det;
	inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
	inv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / det;
	inv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
	inv[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) / det;
	inv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;

}

/***********************************************************************/
/***********************************************************************/
// Calculate the double precision inverse of a dense square matrix
// making use of the Cayley-Hamilton theorem
int calcInvDoubleDense(dbMatrix& Amat, dbMatrix& inverse,
		std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcInvDoubleDense no valid matrix\n"
				<< "to calculate its inverse!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int info = 0;
	int Asize = Amat.size();

	if (Asize == 3) {

		dbMatrix AA;
		innerTensorProduct(Amat, Amat, AA, false, false, logFile);

		if (inverse.size() != Asize)
			allocateArray(inverse, Asize, Asize);

		clearArray(inverse);

		// compute the invariants

		// I1 = trA
		double I1 = 0;

		for (int i = 0; i < Asize; i++)
			I1 += Amat[i][i];

		// I2 = 0.5 [ (tr A)^2 - tr A^2 ] = 0.5*[ I1^2 - A_ik A_ki ]

		// tr A^2 = A_ik A_ki
		double I2;
		double trA2 = 0;

		for (int i = 0; i < Asize; i++)

			trA2 += AA[i][i];

		I2 = 0.5 * (pow(I1, 2.0) - trA2);

		// I3 = 1/3*(trA^3 - 3/2*trA*trA^2 + 1/2*(trA)^3)
		double I3 = 0;
		double trA3 = 0;

		// trA^3 = A_ik*A_kl*A_li

		for (int i = 0; i < Asize; i++)

			for (int k = 0; k < Asize; k++)

				for (int l = 0; l < Asize; l++)

					trA3 += Amat[i][k] * Amat[k][l] * Amat[l][i];

		I3 = 1.0 / 3.0 * (trA3 - 3.0 / 2.0 * I1 * trA2 + 0.5 * pow(I1, 3.0));

		/**********************************************************************/
		// A^-1 = 1/I3*(A^2 - I1*A + I2*1)
		for (int i = 0; i < Asize; i++) {

			inverse[i][i] += I2;

			for (int j = 0; j < Asize; j++) {

				inverse[i][j] += AA[i][j] - I1 * Amat[i][j];
				inverse[i][j] *= 1.0 / I3;
			}

		}

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		int flag;
		dbMatrix lapackInverse;
		if(computeNorm(Amat,3,logFile) != 0)
		flag = calcInvDoubleSparse(Amat,lapackInverse,logFile);
		logFile<<"********************* inverse ***********************"<<endl;
		for(int i=0;i<Asize;i++)
		for(int j=0;j<Asize;j++)
		logFile<<"Ainv["<<i<<"]["<<j<<"] = "<<inverse[i][j]
		<<" ?= "<<lapackInverse[i][j]<<endl;
#endif

		return info;

	}

	else

		return (calcInvDoubleSparse(Amat, inverse, logFile));

}

/***********************************************************************/
/***********************************************************************/
// Calculate the double precision inverse of a sparse square matrix using
// LAPACK routines.
int calcInvDoubleSparse(dbMatrix& Amat, dbMatrix& inverse,
		std::ofstream& logFile) {

	using namespace std;

	if (Amat.size() != Amat[0].size()) {
		logFile << "In commonFunctions::calcInvDoubleSparse no valid matrix\n"
				<< "to calculate its inverse!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Intialisation of the variables used by the external subroutine.
	int info = 0;
	int lda = Amat.size();
	int n = Amat[0].size();

	int lwork = 3 * Amat.size();

	//   double* work = new double[lwork+1];
	//   double* mat = new double[lda*n+1];
	//   int* ipvt = new int[Amat.size()+1];
	dbVector work(lwork);
	dbVector mat(lda * n);
	intVector ipvt(Amat.size());

	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			mat[i * Amat.size() + j] = Amat[i][j];

#ifdef _commonDebugMode_
	logFile<<"********** not inverted FORTRAN matrix *************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine 'dgetrf' factorizing Amat by
	// gaussian elimination.
	dgetrf_(lda, n, &mat[0], lda, &ipvt[0], info);

	if (info != 0) {
		logFile << "In commonFunctions::calcInvDoubleSparse the factorization\n"
				<< "of input matrix could not be calculated!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

#ifdef _commonDebugMode_
	logFile<<"*********************** pivots **********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	logFile<<"pivot["<<i<<"] = "<<ipvt[i]<<endl;
	logFile<<"**************** factorized matrix ******************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
#endif

	/*********************************************************************/
	// Call of external FORTRAN subroutine 'dgetri' to calculate the
	// the inverse of matrix Amat by
	dgetri_(n, &mat[0], lda, &ipvt[0], &work[0], lwork, info);

#ifdef _commonDebugMode_
	logFile<<"**************** fortran inverse *********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*Amat.size()+j]<<endl;
	dbMatrix one(Amat.size(),dbVector(Amat.size()));
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	for(int k=0;k<Amat.size();k++)
	one[i][j] += Amat[i][k]*mat[k*Amat.size()+j];
	logFile<<"************** Kronecker delta(?) ********************"<<endl;
	for(int i=0;i<Amat.size();i++)
	for(int j=0;j<Amat.size();j++)
	logFile<<"delta["<<i<<"]["<<j<<"] = "<<one[i][j]<<endl;
#endif

	if (info != 0) {
		logFile << "In commonFunctions::calcInvDoubleSparse the inverse could\n"
				<< "not be calculated!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/*********************************************************************/
	// Check the pre-allocated storage.
	int Asize = Amat.size();

	if (inverse.size() != Asize)
		allocateArray(inverse, Asize, Asize);

	clearArray(inverse);

	// Store the calculated inverse in a C++ array.
	for (int i = 0; i < Amat.size(); i++)

		for (int j = 0; j < Amat.size(); j++)

			inverse[i][j] = mat[i * Amat.size() + j];

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Check the Positive Definiteness of a rank 4 Tensor
void posDefCheck4(dbMatrix4& C4, dbVector& det, std::ofstream& logFile) {

	using namespace std;
	intMatrix v2m = matrixToVectorFull(3);
	dbMatrix posDefcTensor;

	if (posDefcTensor.size() < 9) {
		allocateArray(posDefcTensor, 9, 9);
	}

	if (det.size() != posDefcTensor.size()) {
		allocateArray(det, posDefcTensor.size());
	}

	clearArray(posDefcTensor);
	clearArray(det);

	//Converting 4th order Tensor to Matrix
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++) {
					posDefcTensor[v2m[i][j]][v2m[k][l]] += C4[i][j][k][l];
				}

	//Calculating the Positive-Definiteness
	for (int i = 1; i < (1 + posDefcTensor.size()); i++) {
		double det2 = 0;
		dbMatrix cSemi;
		if (cSemi.size() < i)
			allocateArray(cSemi, i, i);
		clearArray(cSemi);

		for (int j = 0; j < i; j++)
			for (int k = 0; k < i; k++) {
				cSemi[j][k] = posDefcTensor[j][k];
			}
		calcDetDoubleSparse(cSemi, det2, logFile);
		det[i - 1] = det2;
	}

}

/***********************************************************************/
/***********************************************************************/
// Check the Positive Definiteness of a Matrix or rank 2 Tensor
void posDefCheck(dbMatrix& C2, dbVector& det, std::ofstream& logFile) {

	using namespace std;

	if (det.size() != C2.size()) {
		allocateArray(det, C2.size());
	}
	clearArray(det);

	//Calculating the Positive-Definiteness
	for (int i = 1; i < (1 + C2.size()); i++) {
		double det2 = 0;
		dbMatrix cSemi;
		if (cSemi.size() < i)
			allocateArray(cSemi, i, i);
		clearArray(cSemi);

		for (int j = 0; j < i; j++)
			for (int k = 0; k < i; k++) {
				cSemi[j][k] = C2[j][k];
			}
		calcDetDoubleSparse(cSemi, det2, logFile);
		det[i - 1] = det2;

	}

}

/************************************************************************/
/************************************************************************/
// Compute a vector norm.
double computeNorm(dbVector& aVector, int type, std::ofstream& logFile) {

	using namespace std;

	double norm = 0;

	switch (type) {
	case 1:

		for (int i = 0; i < aVector.size(); i++)
			norm += aVector[i];

		break;
	case 2:

		for (int i = 0; i < aVector.size(); i++)
			norm += pow(aVector[i], 2);

		norm = sqrt(norm);

		break;
	case 3:

		for (int i = 0; i < aVector.size(); i++)

			if (fabs(aVector[i]) > norm)
				norm = fabs(aVector[i]);

		break;
	default:
		logFile << "Chosen convergence norm type isn't available, "
				<< "check file input.dat!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"NORM = "<<norm<<endl;
#endif

	return norm;

}

/***********************************************************************/
void calcL1Norm(dbVector& oldGlobalForceVec, double& residualNorm,
		std::ofstream& logFile) {

	using namespace std;

	residualNorm = 0;

	for (int i = 0; i < oldGlobalForceVec.size(); i++)
		residualNorm += fabs(oldGlobalForceVec[i]);

}

void calcL2Norm(dbVector& oldGlobalForceVec, double& residualNorm,
		std::ofstream& logFile) {

	using namespace std;

	double summe = 0;

	for (int i = 0; i < oldGlobalForceVec.size(); i++)
		summe += pow(oldGlobalForceVec[i], 2);

	residualNorm = sqrt(summe);

}

void calcLInftyNorm(dbVector& oldGlobalForceVec, double& residualNorm,
		std::ofstream& logFile) {

	using namespace std;

	residualNorm = fabs(oldGlobalForceVec[0]);

	for (int i = 0; i < oldGlobalForceVec.size(); i++)

		if (fabs(oldGlobalForceVec[i]) < residualNorm)
			residualNorm = fabs(oldGlobalForceVec[i]);

}

/************************************************************************/
/************************************************************************/
// Compute a matrix norm.
double computeNorm(dbMatrix& aMatrix, int type, std::ofstream& logFile) {

	using namespace std;

	double maxEigenvalue, norm, sum;
	dbMatrix dummyMatrix;

	switch (type) {

	// maximum absolute column sum norm
	// || A ||_1 = max { sum_{i=1}^n |a_ij| }_j
	case 1:

		norm = 0;

		for (int j = 0; j < aMatrix[0].size(); j++) {
			sum = 0;

			for (int i = 0; i < aMatrix.size(); i++)

				sum += fabs(aMatrix[i][j]);

			if (sum > norm)
				norm = sum;

		}

		break;

		// spectral norm, which is the square root of the maximum eigenvalue
		// of A^H A ( A^H is the conjugate transpose )
		// || A ||_2 = ( maximum eigenvalue of A^H A )^{1/2}
	case 2:

		// Compute A^H A
		dummyMatrix = dbMatrix(aMatrix[0].size(), dbVector(aMatrix[0].size()));

		for (int i = 0; i < aMatrix[0].size(); i++)

			for (int j = 0; j < aMatrix[0].size(); j++)

				for (int k = 0; k < aMatrix.size(); k++)

					dummyMatrix[i][j] = aMatrix[k][i] * aMatrix[k][j];

		if (computeMaxEigenValue(dummyMatrix, maxEigenvalue, logFile) != 0) {
			logFile << "Calulation of the maximal eigenvalue in \n"
					<< "commonFunctions::computeNorm failed!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		// || A ||_2 = ( maximum eigenvalue of A^H A )^{1/2}
		norm = sqrt(maxEigenvalue);

		break;

		// Frobenius or Euclidean norm of A(m,n)
		// || A ||_F = ( sum_{i=1}^m sum_{j=1}^n a_ij^2 )^{1/2}
	case 3:

		norm = 0;

		for (int i = 0; i < aMatrix.size(); i++)

			for (int j = 0; j < aMatrix[0].size(); j++)

				norm += pow(aMatrix[i][j], 2);

		norm = sqrt(norm);

		break;

		// maximum absolute row sum norm
		// || A ||_{infinity} = max { sum_{j=1}^n | |a_ij| }_i
	case 4:

		norm = 0;

		for (int i = 0; i < aMatrix.size(); i++) {
			sum = 0;

			for (int j = 0; j < aMatrix[0].size(); j++)

				sum += fabs(aMatrix[i][j]);

			if (sum > norm)
				norm = sum;

		}

		break;

	default:
		logFile << "Chosen matrix norm type isn't available!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"NORM = "<<norm<<endl;
#endif

	return norm;

}

/***********************************************************************/
void getFEMMeshData(intVector& data, std::map<std::string, double>& params) {

	using namespace std;

	// Select a mesh type.
	switch (data[0]) {

	// tetraedral element
	case 1:

		switch (data[1]) {

		// 4 node tetraedral element
		case 1:
			params["nodesPerVolumeElement"] = 4;
			params["nodesPerSurfaceElement"] = 3;
			params["nodesPerLineElement"] = 2;

			break;
		default:
			cerr << "Finte element type " << data[0] << " of order " << data[1]
					<< " isn't " << "supported!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		break;

		// hexahedral element
	case 2:

		switch (data[1]) {

		// 8 node hexahedral element
		case 1:
			params["nodesPerVolumeElement"] = 8;
			params["nodesPerSurfaceElement"] = 4;
			params["nodesPerLineElement"] = 2;
			break;

			// 27 node hexahedral element
		case 2:
			params["nodesPerVolumeElement"] = 27;
			params["nodesPerSurfaceElement"] = 9;
			params["nodesPerLineElement"] = 3;
			break;

		default:
			cerr << "Finte element type " << data[0] << " of order " << data[1]
					<< " isn't " << "supported!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		break;

	default:
		cerr << "Finte element type " << data[0] << " of order " << data[1]
				<< " isn't " << "supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/************************************************************************/
void getGaussQuadratureData(intVector& data,
		std::map<std::string, double>& params) {

	using namespace std;

	// Select a back ground mesh type.
	switch (data[0]) {

	// tetrahedral element
	case 1:

		// volume Gauss quadrature
		switch (data[1]) {

		case 1:
			params["gaussPointsPerVolumeElement"] = 1;
			break;

		case 2:
			params["gaussPointsPerVolumeElement"] = 4;
			break;

		case 3:
			params["gaussPointsPerVolumeElement"] = 5;
			break;

		default:
			cerr << "Volume Gauss quadrature order " << data[1]
					<< " isn't supported!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		// surface Gauss quadrature
		switch (data[2]) {

		case 1:
			params["gaussPointsPerSurfaceElement"] = 1;
			break;

		case 2:
			params["gaussPointsPerSurfaceElement"] = 3;
			break;

		case 3:
			params["gaussPointsPerSurfaceElement"] = 4;
			break;

		default:
			cerr << "Surface Gauss quadrature order " << data[3]
					<< " isn't supported!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		break;

		/********************************************************************/
		// hexahedral element
	case 2:

		// volume Gauss quadrature
		switch (data[1]) {

		case 1:
			params["gaussPointsPerVolumeElement"] = 1;
			break;

		case 2:
			params["gaussPointsPerVolumeElement"] = 8;
			break;

		case 3:
			params["gaussPointsPerVolumeElement"] = 27;
			break;

		case 4:
			params["gaussPointsPerVolumeElement"] = 64;
			break;

		case 5:
			params["gaussPointsPerVolumeElement"] = 125;
			break;

		default:
			cerr << "Volume Gauss quadrature order " << data[1]
					<< " isn't supported!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		// surface Gauss quadrature
		switch (data[2]) {

		case 1:
			params["gaussPointsPerSurfaceElement"] = 1;
			break;

		case 2:
			params["gaussPointsPerSurfaceElement"] = 4;
			break;

		case 3:
			params["gaussPointsPerSurfaceElement"] = 9;
			break;

		case 4:
			params["gaussPointsPerSurfaceElement"] = 16;
			break;

		case 5:
			params["gaussPointsPerSurfaceElement"] = 25;
			break;

		default:
			cerr << "Surface Gauss quadrature order " << data[2]
					<< " isn't supported!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		break;

	default:
		cerr << "FEM element-type " << data[0] << " isn't supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

	/**********************************************************************/
	// line Gauss quadrature
	switch (data[3]) {

	case 1:
		params["gaussPointsPerLineElement"] = 1;
		break;

	case 2:
		params["gaussPointsPerLineElement"] = 2;
		break;

	case 3:
		params["gaussPointsPerLineElement"] = 3;
		break;

	case 4:
		params["gaussPointsPerLineElement"] = 4;
		break;

	case 5:
		params["gaussPointsPerLineElement"] = 5;
		break;

	default:
		cerr << "Line Gauss quadrature order " << data[3] << " isn't supported!"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/***********************************************************************/
intMatrix matrixToVector(int size) {

	using namespace std;

	intMatrix m(size, intVector(size));

	switch (size) {

	case 1:
		m[0][0] = 0;

		break;

	case 2:
		m[0][0] = 0;
		m[1][1] = 1;
		m[0][1] = 2;
		m[1][0] = 2;

		break;

	case 3:
		m[0][0] = 0;
		m[1][1] = 1;
		m[2][2] = 2;
		m[0][1] = 3;
		m[1][0] = 3;
		m[1][2] = 4;
		m[2][1] = 4;
		m[2][0] = 5;
		m[0][2] = 5;

		break;

	case 4:
		m[0][0] = 0;
		m[1][1] = 1;
		m[2][2] = 2;
		m[0][1] = 3;
		m[1][0] = 3;
		m[1][2] = 4;
		m[2][1] = 4;
		m[2][0] = 5;
		m[0][2] = 5;

		m[3][3] = 6;
		m[0][3] = 7;
		m[3][0] = 7;
		m[1][3] = 8;
		m[3][1] = 8;
		m[2][3] = 9;
		m[3][2] = 9;

		break;

	case 5:
		m[0][0] = 0;
		m[1][1] = 1;
		m[2][2] = 2;
		m[0][1] = 3;
		m[1][0] = 3;
		m[1][2] = 4;
		m[2][1] = 4;
		m[2][0] = 5;
		m[0][2] = 5;

		m[3][3] = 6;
		m[4][4] = 7;
		m[0][3] = 8;
		m[3][0] = 8;
		m[0][4] = 9;
		m[4][0] = 9;
		m[1][3] = 10;
		m[3][1] = 10;
		m[1][4] = 11;
		m[4][1] = 11;
		m[2][3] = 12;
		m[3][2] = 12;
		m[2][4] = 13;
		m[4][2] = 13;
		m[3][4] = 14;
		m[4][3] = 14;

		break;

	case 6:
		m[0][0] = 0;
		m[1][1] = 1;
		m[2][2] = 2;
		m[0][1] = 3;
		m[1][0] = 3;
		m[1][2] = 4;
		m[2][1] = 4;
		m[2][0] = 5;
		m[0][2] = 5;

		m[3][3] = 6;
		m[4][4] = 7;
		m[5][5] = 8;
		m[0][3] = 9;
		m[3][0] = 9;
		m[0][4] = 10;
		m[4][0] = 10;
		m[0][5] = 11;
		m[5][0] = 11;
		m[1][3] = 12;
		m[3][1] = 12;
		m[1][4] = 13;
		m[4][1] = 13;
		m[1][5] = 14;
		m[5][1] = 14;
		m[2][3] = 15;
		m[3][2] = 15;
		m[2][4] = 16;
		m[4][2] = 16;
		m[2][5] = 17;
		m[5][2] = 17;
		m[3][4] = 18;
		m[4][3] = 18;
		m[3][5] = 19;
		m[5][3] = 19;
		m[4][5] = 20;
		m[5][4] = 20;

		break;

	default:
		cerr << "Matrix (" << size << "," << size
				<< ") to vector conversion is \n" << "not supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return m;
}

/***********************************************************************/
intMatrix matrixToVectorFull(int size) {

	using namespace std;

	intMatrix m(size, intVector(size));

	switch (size) {

	case 1:
		m[0][0] = 0;

		break;

	case 2:
		m[0][0] = 0;
		m[1][1] = 1;
		m[0][1] = 2;
		m[1][0] = 3;

		break;

	case 3:
		m[0][0] = 0;
		m[1][1] = 1;
		m[2][2] = 2;
		m[0][1] = 3;
		m[1][0] = 4;
		m[1][2] = 5;
		m[2][1] = 6;
		m[2][0] = 7;
		m[0][2] = 8;

		break;

	case 4:
		m[0][0] = 0;
		m[1][1] = 1;
		m[2][2] = 2;
		m[0][1] = 3;
		m[1][0] = 4;
		m[1][2] = 5;
		m[2][1] = 6;
		m[2][0] = 7;
		m[0][2] = 8;

		m[3][3] = 9;
		m[0][3] = 10;
		m[3][0] = 11;
		m[1][3] = 12;
		m[3][1] = 13;
		m[2][3] = 14;
		m[3][2] = 15;

		break;

	default:
		cerr << "Matrix (" << size << "," << size
				<< ") to vector conversion is \n" << "not supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return m;
}

intMatrix vectorToMatrix(int size) {

	using namespace std;

	int numOfEntries = (int) (size * size + size) / 2;

	intMatrix m(numOfEntries, intVector(2));

	switch (size) {

	case 1:

		m[0][0] = 0;
		m[0][1] = 0;

		break;

	case 2:

		m[0][0] = 0;
		m[0][1] = 0;

		m[1][0] = 1;
		m[1][1] = 1;

		m[2][0] = 0;
		m[2][1] = 1;

		break;

	case 3:

		m[0][0] = 0;
		m[0][1] = 0;

		m[1][0] = 1;
		m[1][1] = 1;

		m[2][0] = 2;
		m[2][1] = 2;

		m[3][0] = 0;
		m[3][1] = 1;

		m[4][0] = 1;
		m[4][1] = 2;

		m[5][0] = 2;
		m[5][1] = 0;

		break;

	case 4:

		m[0][0] = 0;
		m[0][1] = 0;

		m[1][0] = 1;
		m[1][1] = 1;

		m[2][0] = 2;
		m[2][1] = 2;

		m[3][0] = 0;
		m[3][1] = 1;

		m[4][0] = 1;
		m[4][1] = 2;

		m[5][0] = 2;
		m[5][1] = 0;

		m[6][0] = 3;
		m[6][1] = 3;

		m[7][0] = 0;
		m[7][1] = 3;

		m[8][0] = 1;
		m[8][1] = 3;

		m[9][0] = 2;
		m[9][1] = 3;

		break;

	case 5:

		m[0][0] = 0;
		m[0][1] = 0;

		m[1][0] = 1;
		m[1][1] = 1;

		m[2][0] = 2;
		m[2][1] = 2;

		m[3][0] = 0;
		m[3][1] = 1;

		m[4][0] = 1;
		m[4][1] = 2;

		m[5][0] = 2;
		m[5][1] = 0;

		m[6][0] = 3;
		m[6][1] = 3;

		m[7][0] = 4;
		m[7][1] = 4;

		m[8][0] = 0;
		m[8][1] = 3;

		m[9][0] = 0;
		m[9][1] = 4;

		m[10][0] = 1;
		m[10][1] = 3;

		m[11][0] = 1;
		m[11][1] = 4;

		m[12][0] = 2;
		m[12][1] = 3;

		m[13][0] = 2;
		m[13][1] = 4;

		m[14][0] = 3;
		m[14][1] = 4;

		break;

	case 6:

		m[0][0] = 0;
		m[0][1] = 0;

		m[1][0] = 1;
		m[1][1] = 1;

		m[2][0] = 2;
		m[2][1] = 2;

		m[3][0] = 0;
		m[3][1] = 1;

		m[4][0] = 1;
		m[4][1] = 2;

		m[5][0] = 2;
		m[5][1] = 0;

		m[6][0] = 3;
		m[6][1] = 3;

		m[7][0] = 4;
		m[7][1] = 4;

		m[8][0] = 5;
		m[8][1] = 5;

		m[9][0] = 0;
		m[9][1] = 3;

		m[10][0] = 0;
		m[10][1] = 4;

		m[11][0] = 0;
		m[11][1] = 5;

		m[12][0] = 1;
		m[12][1] = 3;

		m[13][0] = 1;
		m[13][1] = 4;

		m[14][0] = 1;
		m[14][1] = 5;

		m[15][0] = 2;
		m[15][1] = 3;

		m[16][0] = 2;
		m[16][1] = 4;

		m[17][0] = 2;
		m[17][1] = 5;

		m[18][0] = 3;
		m[18][1] = 4;

		m[19][0] = 3;
		m[19][1] = 5;

		m[20][0] = 4;
		m[20][1] = 5;

		break;

	default:
		cerr << "Vector (" << numOfEntries << ") to matrix conversion "
				<< "is not supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return m;
}

intMatrix vectorToMatrix(intVector& idx) {

	using namespace std;

	int size = idx.size();
	int numOfEntries = (int) (size * size + size) / 2;

	intMatrix m(numOfEntries, intVector(2));

	switch (size) {

	case 1:

		m[0][0] = idx[0];
		m[0][1] = idx[0];

		break;

	case 2:

		m[0][0] = idx[0];
		m[0][1] = idx[0];

		m[1][0] = idx[1];
		m[1][1] = idx[1];

		m[2][0] = idx[0];
		m[2][1] = idx[1];

		break;

	case 3:

		m[0][0] = idx[0];
		m[0][1] = idx[0];

		m[1][0] = idx[1];
		m[1][1] = idx[1];

		m[2][0] = idx[2];
		m[2][1] = idx[2];

		m[3][0] = idx[0];
		m[3][1] = idx[1];

		m[4][0] = idx[1];
		m[4][1] = idx[2];

		m[5][0] = idx[2];
		m[5][1] = idx[0];

		break;

	default:
		cerr << "Vector (" << numOfEntries << ") to matrix conversion "
				<< "is not supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return m;
}

/***********************************************************************/
// Get the symmetric matrix which is stored as a vector.
dbMatrix convertVectorToMatrix(int size, dbVector& v) {

	using namespace std;

	int numOfEntries = (int) (size * size + size) / 2;
	dbMatrix m(size, dbVector(size));

	if (v.size() != numOfEntries) {
		cerr << "In convertVectorToMatrix vector has not enough entries!"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	switch (size) {

	case 1:
		m[0][0] = v[0];

		break;

	case 2:
		m[0][0] = v[0];
		m[1][1] = v[1];
		m[0][1] = v[2];
		m[1][0] = v[2];

		break;

	case 3:
		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		m[0][1] = v[3];
		m[1][0] = v[3];
		m[1][2] = v[4];
		m[2][1] = v[4];
		m[2][0] = v[5];
		m[0][2] = v[5];

		break;

	case 4:
		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		m[0][1] = v[3];
		m[1][0] = v[3];
		m[1][2] = v[4];
		m[2][1] = v[4];
		m[2][0] = v[5];
		m[0][2] = v[5];

		m[3][3] = v[6];
		m[0][3] = v[7];
		m[3][0] = v[7];
		m[1][3] = v[8];
		m[3][1] = v[8];
		m[2][3] = v[9];
		m[3][2] = v[9];

		break;

	case 5:
		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		m[0][1] = v[3];
		m[1][0] = v[3];
		m[1][2] = v[4];
		m[2][1] = v[4];
		m[2][0] = v[5];
		m[0][2] = v[5];

		m[3][3] = v[6];
		m[4][4] = v[7];
		m[0][3] = v[8];
		m[3][0] = v[8];
		m[0][4] = v[9];
		m[4][0] = v[9];
		m[1][3] = v[10];
		m[3][1] = v[10];
		m[1][4] = v[11];
		m[4][1] = v[11];
		m[2][3] = v[12];
		m[3][2] = v[12];
		m[2][4] = v[13];
		m[4][2] = v[13];
		m[3][4] = v[14];
		m[4][3] = v[14];

		break;

	case 6:
		m[0][0] = v[0];
		m[1][1] = v[1];
		m[2][2] = v[2];
		m[0][1] = v[3];
		m[1][0] = v[3];
		m[1][2] = v[4];
		m[2][1] = v[4];
		m[2][0] = v[5];
		m[0][2] = v[5];

		m[3][3] = v[6];
		m[4][4] = v[7];
		m[5][5] = v[8];
		m[0][3] = v[9];
		m[3][0] = v[9];
		m[0][4] = v[10];
		m[4][0] = v[10];
		m[0][5] = v[11];
		m[5][0] = v[11];
		m[1][3] = v[12];
		m[3][1] = v[12];
		m[1][4] = v[13];
		m[4][1] = v[13];
		m[1][5] = v[14];
		m[5][1] = v[14];
		m[2][3] = v[15];
		m[3][2] = v[15];
		m[2][4] = v[16];
		m[4][2] = v[16];
		m[2][5] = v[17];
		m[5][2] = v[17];
		m[3][4] = v[18];
		m[4][3] = v[18];
		m[3][5] = v[19];
		m[5][3] = v[19];
		m[4][5] = v[20];
		m[5][4] = v[20];

		break;

	default:
		cerr << "In convertVectorToMatrix vector to matrix (" << size << ","
				<< size << ")\n" << "conversion is not supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return m;
}

/***********************************************************************/
// Computation of roots of polynom 3rd degree ax^3 + bx^2 + cx + d
void polynom3rdRoots(dbVector& coefficients, dbVector& roots,
		std::ofstream& logFile) {

	using namespace std;

	if (coefficients.size() < 4 || coefficients[0] == 0) {
		logFile << "In order to calculate the roots of polynom of"
				<< "third degree four valid coefficients are needed!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	double a = coefficients[0];
	double b = coefficients[1];
	double c = coefficients[2];
	double d = coefficients[3];

	// substitution of x = y-b/(3*a) =>
	// r*y^3 + p*y + q = 0
	// where
	// r = a
	// p = c - b^2/(3*a)
	// q = 2*b^3/(27*a^2) - bc/(3a) + d

	double p = c - pow(b, 2) / (3 * a);
	double q = 2 * pow(b, 3) / (27 * pow(a, 2)) - b * c / (3 * a) + d;

	// substitution of y = z - p/(3*a*z) and multiplying through by z^3 =>
	// a z^6 + q*z^3 - p^3/(27a^2) = 0 and
	// z^3 = (-q+-sqrt(q^2+4*p^3/(27*a))/(2*a)

	double z1 = pow((-q + sqrt(pow(q, 2) + 4 * pow(p, 3) / (27 * a))) / (2 * a),
			1.0 / 3.0);
	double z2 = -z1;

	double z3 = pow((-q - sqrt(pow(q, 2) + 4 * pow(p, 3) / (27 * a))) / (2 * a),
			1.0 / 3.0);
	double z4 = -z1;

	// backsubstitution
	roots = dbVector(4);
	roots[0] = z1 - p / (3 * a * z1) - b / (3 * a);
	roots[1] = z2 - p / (3 * a * z2) - b / (3 * a);
	roots[2] = z3 - p / (3 * a * z3) - b / (3 * a);
	roots[3] = z4 - p / (3 * a * z4) - b / (3 * a);

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"*** computation of the roots of polynom 3rd degree **"<<endl;
	logFile<<"r1 = "<<roots[0]<<endl;
	logFile<<"r2 = "<<roots[1]<<endl;
	logFile<<"r3 = "<<roots[2]<<endl;
	logFile<<"r4 = "<<roots[3]<<endl;
	for(int i=0;i<roots.size();i++)
	logFile<<"p(r"<<i<<") = "<<a*pow(roots[i],3)+b*pow(roots[i],2)+
	c*roots[i]+d<<endl;
#endif

}

/***********************************************************************/

dbVector quadraticPolynom(double x) {

	using namespace std;

	vector<double> basis(3);

	basis[0] = 1.0;          // 1
	basis[1] = x;            // x
	basis[2] = x * x;          // x^2

	return (basis);
}

dbVector dQuadraticPolynom(double x) {

	using namespace std;

	vector<double> basis(3);

	basis[0] = 0.0;            // 0
	basis[1] = 1.0;            // 1
	basis[2] = 2.0 * x;          // 2x

	return (basis);
}

dbVector d2QuadraticPolynom(double x) {

	using namespace std;

	vector<double> basis(3);

	basis[0] = 0;            // 0
	basis[1] = 0;            // 0
	basis[2] = 2.0;          // 2

	return (basis);
}

/***********************************************************************/

dbVector cubicPolynom(double x) {

	using namespace std;

	vector<double> basis(4);

	basis[0] = 1.0;            // 1
	basis[1] = x;            // x
	basis[2] = x * x;          // x^2
	basis[3] = x * x * x;        // x^3

	return (basis);
}

dbVector dCubicPolynom(double x) {

	using namespace std;

	vector<double> basis(4);

	basis[0] = 0;            // 0
	basis[1] = 1.0;            // 1
	basis[2] = 2.0 * x;          // 2x
	basis[3] = 3.0 * x * x;        // 3x^2

	return (basis);
}

dbVector d2CubicPolynom(double x) {

	using namespace std;

	vector<double> basis(4);

	basis[0] = 0;            // 0
	basis[1] = 0;            // 0
	basis[2] = 2.0;            // 2
	basis[3] = 6.0 * x;          // 6x

	return (basis);
}

dbVector d3CubicPolynom(double x) {

	using namespace std;

	vector<double> basis(4);

	basis[0] = 0;            // 0
	basis[1] = 0;            // 0
	basis[2] = 0;            // 0
	basis[3] = 6.0;          // 6

	return (basis);
}

/***********************************************************************/

dbVector quarticPolynom(double x) {

	using namespace std;

	vector<double> basis(5);

	basis[0] = 1.0;          // 1
	basis[1] = x;            // x
	basis[2] = x * x;          // x^2
	basis[3] = x * x * x;        // x^3
	basis[4] = x * x * x * x;      // x^4

	return (basis);
}

dbVector dQuarticPolynom(double x) {

	using namespace std;

	vector<double> basis(5);

	basis[0] = 0.0;            // 0
	basis[1] = 1.0;            // 1
	basis[2] = 2.0 * x;          // 2x
	basis[3] = 3.0 * x * x;        // 3x^2
	basis[4] = 4.0 * x * x * x;      // 4x^3

	return (basis);
}

dbVector d2QuarticPolynom(double x) {

	using namespace std;

	vector<double> basis(5);

	basis[0] = 0;            // 0
	basis[1] = 0;            // 0
	basis[2] = 2.0;          // 2
	basis[3] = 6.0 * x;        // 6x
	basis[4] = 12.0 * x * x;     // 12x^2

	return (basis);
}

dbVector d3QuarticPolynom(double x) {

	using namespace std;

	vector<double> basis(5);

	basis[0] = 0;            // 0
	basis[1] = 0;            // 0
	basis[2] = 0;            // 0
	basis[3] = 6.0;          // 6
	basis[4] = 24.0 * x;       // 24x

	return (basis);
}

/***********************************************************************/
// Compute all eigenvalues of a real symmetric matrix.
int calcEigenValues(dbMatrix& aMatrix, dbVector& eigenvalues,
		std::ofstream& logFile) {

	using namespace std;

	int info = 0;

	char jobz = 'N';
	char uplo = 'U';

	int n = aMatrix.size();
	int lda = n;
	int lwork = 1 + 6 * n + 2 * n * n + 1;
	int liwork = 3 + 5 * n + 1;

	lwork *= 2;
	liwork *= 2;

	//   double* mat = new double[n*n+1];
	//   double* w = new double[n+1];
	//   double* work = new double[lwork+1];
	//   int* iwork = new int[liwork+1];
	dbVector mat(n * n);
	dbVector w(n);
	dbVector work(lwork);
	intVector iwork(liwork);

	for (int i = 0; i < n; i++)

		for (int j = 0; j < n; j++) {

			if (aMatrix[i][j] != aMatrix[j][i]) {
				logFile << "In commonFunctions::calcEigenValues only\n"
						<< "symmetric matrices are admissible!" << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}

			mat[i * n + j] = aMatrix[i][j];

		}

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"************** stored fortran matrix ****************"<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*n+j]<<endl;
#endif

	// call of the LAPACK routine
	dsyevd_(jobz, uplo, n, &mat[0], lda, &w[0], &work[0], lwork, &iwork[0],
			liwork, info);

	if (info != 0) {
		logFile << "In commonFunctions::calcEigenValues computation\n "
				<< "eigenvalues failed!" << endl;
	}

#ifdef _commonDebugMode_
	logFile<<"************** computed eigenvalues ****************"<<endl;
	for(int i=0;i<n;i++) {
		logFile<<"e[i]= "<<w[i]<<endl;
	}
#endif

	// Store the real and the imaginary part of all eigenvalues.
	eigenvalues.resize(n);

	for (int i = 0; i < n; i++)
		eigenvalues[i] = w[i];

	//delete[] mat,work,w,iwork;

	return info;
}

/***********************************************************************/
// Compute all eigenvalues of a real general matrix.
int calcEigenValues(dbMatrix& aMatrix, dbMatrix& eigenvalues,
		std::ofstream& logFile) {

	using namespace std;

	char jobvl = 'N';
	char jobvr = 'N';

	int lda = aMatrix.size();
	int n = aMatrix[0].size();

	int ldvl = 1;
	int ldvr = 1;

	int lwork = 4 * aMatrix.size();

	//   double* wr = new double[n+1];
	//   double* wi = new double[n+1];
	//   double* vl;
	//   double* vr;
	//   double* work = new double[lwork+1];
	//   double* mat = new double[lda*n+1];

	dbVector wr(n);
	dbVector wi(n);
	dbVector vl(ldvl);
	dbVector vr(ldvr);
	dbVector work(lwork);
	dbVector mat(lda * n);

	for (int i = 0; i < lda; i++)

		for (int j = 0; j < n; j++)

			mat[j * n + i] = aMatrix[i][j];

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"************** stored fortran matrix ****************"<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*n+j]<<endl;
#endif

	int info = 0;

	// call of the LAPACK routine
	dgeev_(jobvl, jobvr, n, &mat[0], lda, &wr[0], &wi[0], &vl[0], ldvl, &vr[0],
			ldvr, &work[0], lwork, info);

	if (info != 0) {
		logFile << "In commonFunctions::calcEigenvalues computation of \n"
				<< "eigenvalues failed!" << endl;
	}

#ifdef _commonDebugMode_
	logFile<<"************** computed eigenvalues *****************"<<endl;
	for(int i=0;i<n;i++) {
		logFile<<"re["<<i<<"] = "<<wr[i]<<endl;
		logFile<<"im["<<i<<"] = "<<wi[i]<<endl;
	}
#endif

	// Store the real and the imaginary part of all eigenvalues.
	eigenvalues = dbMatrix(n, dbVector(2));

	for (int i = 0; i < n; i++) {
		eigenvalues[i][0] = wr[i];
		eigenvalues[i][1] = wi[i];
	}

	//delete[] mat,work,wr,wi;

	return info;
}

/***********************************************************************/
// Compute all eigenvalues of a generalized eigenvalue problem
// consisting of real general matrices.
int calcEigenValues(dbMatrix& aMatrix, dbMatrix& bMatrix, dbMatrix& eigenvalues,
		std::ofstream& logFile) {

	using namespace std;

	char jobvl = 'N';
	char jobvr = 'N';

	int lda = aMatrix.size();
	int n = aMatrix[0].size();
	double* matA = new double[lda * n + 1];

	for (int i = 0; i < lda; i++)

		for (int j = 0; j < n; j++)

			matA[i * n + j] = aMatrix[j][i];

	int ldb = bMatrix.size();
	double* matB = new double[ldb * n + 1];

	for (int i = 0; i < ldb; i++)

		for (int j = 0; j < n; j++)

			matB[i * n + j] = bMatrix[i][j];

	if (lda != ldb || aMatrix[0].size() != bMatrix[0].size()) {
		logFile << "In function commonFunctions::calcEigenvalues matrix\n"
				<< "dimensions are not allowed to be different!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	double* alphaR = new double[n + 1];
	double* alphaI = new double[n + 1];
	double* beta = new double[n + 1];

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"************** stored fortran matrices **************"<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"matA["<<i<<"]["<<j<<"] = "<<matA[i*n+j]<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"matB["<<i<<"]["<<j<<"] = "<<matB[i*n+j]<<endl;
#endif

	double* vl;
	int ldvl = 1;
	double* vr;
	int ldvr = 1;

	int lwork = 10 * aMatrix.size();
	double* work = new double[lwork + 1];

	int info = 0;

	// call of the LAPACK routine
	dggev_(jobvl, jobvr, n, matA, lda, matB, ldb, alphaR, alphaI, beta, vl,
			ldvl, vr, ldvr, work, lwork, info);

	if (info != 0) {
		logFile << "In commonFunctions::calcEigenvalues computation of \n"
				<< "eigenvalues failed!" << endl;
	}

#ifdef _commonDebugMode_
	logFile<<"************** computed eigenvalues *****************"<<endl;
	for(int i=0;i<n;i++) {
		logFile<<"re["<<i<<"] = "<<alphaR[i]<<endl;
		if(fabs(beta[i]) > DBL_EPSILON)
		logFile<<"im["<<i<<"] = "<<alphaI[i]/beta[i]<<endl;
		else
		logFile<<"im["<<i<<"] = "<<alphaI[i]<<"/"<<beta[i]<<endl;
	}
#endif

	// Store the real and the imaginary part of all eigenvalues.
	eigenvalues = dbMatrix(n, dbVector(2));

	for (int i = 0; i < n; i++) {

		if (fabs(beta[i]) > DBL_EPSILON)
			eigenvalues[i][0] = alphaR[i] / beta[i];

		else if (alphaR[i] == 0)
			eigenvalues[i][0] = 0;

		else
			eigenvalues[i][0] = DBL_MAX;

		if (fabs(beta[i]) > DBL_EPSILON)
			eigenvalues[i][1] = alphaI[i] / beta[i];

		else if (alphaI[i] == 0)
			eigenvalues[i][1] = 0;

		else
			eigenvalues[i][1] = DBL_MAX;

	}

	delete[] alphaR, alphaI, beta, matA, matB, work;

	return info;
}

/***********************************************************************/
// Compute the maximum eigenvalue of a real general matrix.
int computeMaxEigenValue(dbMatrix& aMatrix, double& maxEigenValue,
		std::ofstream& logFile) {

	using namespace std;

	maxEigenValue = 0;
	dbMatrix eigenvalues;

	int info = calcEigenValues(aMatrix, eigenvalues, logFile);

	for (int i = 0; i < eigenvalues.size(); i++) {

		// check the real part of the eigenvalue.
		if (eigenvalues[i][0] > maxEigenValue)
			maxEigenValue = eigenvalues[i][0];

#ifdef _commonDebugMode_
		if(eigenvalues[i][0] < 0)
		logFile<<"negativ real part eigenvalue = "<<eigenvalues[i][0]<<endl;

		if(eigenvalues[i][1] > 0.0000001)
		logFile<<"large imaginary eigenvalue part = "<<eigenvalues[i][1]<<endl;
#endif

	}

#ifdef _commonDebugMode_
	logFile<<"************** computed eigenvalues *****************"<<endl;
	for(int i=0;i<eigenvalues.size();i++) {
		logFile<<"re["<<i<<"] = "<<eigenvalues[i][0]<<endl;
		logFile<<"im["<<i<<"] = "<<eigenvalues[i][1]<<endl;
	}
#endif

	return info;
}

/***********************************************************************/
// Compute all eigenvalues and eigenvectors of a real symmetric matrix.
int calcEigenValuesVectors(dbMatrix& aMatrix, dbVector& eigenvalues,
		dbMatrix& eigenvectors, std::ofstream& logFile) {

	using namespace std;

	int info = 0;

	char jobz = 'V';
	char uplo = 'U';

	int n = aMatrix.size();
	int lda = n;
	int lwork = 1 + 6 * n + 2 * n * n + 1;
	int liwork = 3 + 5 * n + 1;

	lwork *= 2;
	liwork *= 2;

	//   double* mat = new double[n*n+1];
	//   double* w = new double[n+1];
	//   double* work = new double[lwork+1];
	//   int* iwork = new int[liwork+1];

	dbVector mat(n * n);
	dbVector w(n);
	dbVector work(lwork);
	intVector iwork(liwork);

	for (int i = 0; i < n; i++)

		for (int j = 0; j < n; j++) {

			if (aMatrix[i][j] != aMatrix[j][i]) {
				logFile << "In commonFunctions::calcEigenValuesVectors only\n"
						<< "symmetric matrices are admissible!" << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}

			mat[j * n + i] = aMatrix[i][j];

		}

#ifdef _commonDebugMode_
	logFile<<"*****************************************************"<<endl;
	logFile<<"************** stored fortran matrix ****************"<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"mat["<<i<<"]["<<j<<"] = "<<mat[i*n+j]<<endl;
#endif

	// call of the LAPACK routine
	dsyevd_(jobz, uplo, n, &mat[0], lda, &w[0], &work[0], lwork, &iwork[0],
			liwork, info);

	if (info != 0) {
		logFile << "In commonFunctions::calcEigenValuesVectors computation\n "
				<< "eigenvalues failed!" << endl;
	}

#ifdef _commonDebugMode_
	logFile<<"************** computed eigenvalues ****************"<<endl;
	for(int i=0;i<n;i++) {
		logFile<<"e[i]= "<<w[i]<<endl;
	}
	logFile<<"************** computed eigenvectors ****************"<<endl;
	for(int i=0;i<n;i++) {
		logFile<<"v[i]= ";
		for(int j=0;j<n;j++)
		logFile<<mat[i*n+j]<<" ";
		logFile<<endl;
	}
#endif

	// Store the real and the imaginary part of all eigenvalues.
	resizeArray(eigenvalues, n);

	for (int i = 0; i < n; i++)
		eigenvalues[i] = w[i];

	// Store all eigenvectors.
	resizeArray(eigenvectors, n, n);

	for (int i = 0; i < n; i++)

		for (int j = 0; j < n; j++)

			eigenvectors[j][i] = mat[i * n + j];

	//delete[] mat,work,w,iwork;

	return info;
}

/***********************************************************************/
// Compute all eigenvalues and eigenvectors of a real general matrix.
int calcEigenValuesVectors(dbMatrix& A, dbMatrix& eigenvalues,
		dbMatrix& eigenvectors, std::ofstream& logFile) {

	using namespace std;

	int info = 0;

	char jobvl = 'N';
	char jobvr = 'V';

	int lda = A.size();
	int n = A[0].size();

	int ldvl = n;
	int ldvr = n;

	int lwork = 10 * A.size();

	//   double* matA = new double[lda*n+1];
	//   double* WR = new double[n+1];
	//   double* WI = new double[n+1];
	//   double* vl = new double[ldvl*n+1];
	//   double* vr = new double[ldvr*n+1];
	//   double* work = new double[lwork+1];

	dbVector work(lwork);
	dbVector matA(lda * n);
	dbVector WR(n);
	dbVector WI(n);
	dbVector vl(ldvl * n);
	dbVector vr(ldvr * n);

	for (int i = 0; i < lda; i++)

		for (int j = 0; j < n; j++)

			matA[j * n + i] = A[i][j];

#ifdef _commonDebugMode_
	logFile<<"****************************************************"<<endl;
	logFile<<"************** stored fortran matrices *************"<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"matA["<<i<<"]["<<j<<"] = "<<matA[i*n+j]<<endl;
#endif

	// call of the LAPACK routine
	dgeev_(jobvl, jobvr, n, &matA[0], lda, &WR[0], &WI[0], &vl[0], ldvl, &vr[0],
			ldvr, &work[0], lwork, info);

	if (info != 0) {
		logFile << "In commonFunctions::calcEigenValuesVectors computation\n "
				<< "eigenvalues failed!" << endl;
	}

#ifdef _commonDebugMode_
	logFile<<"************** computed eigenvalues *****************"<<endl;
	for(int i=0;i<n;i++) {
		logFile<<"E_re["<<i<<"] = "<<WR[i]<<endl;
	}
	for(int i=0;i<n;i++) {
		logFile<<"E_im["<<i<<"] = "<<WI[i]<<endl;
	}
#endif

	// Store the real and the imaginary part of all eigenvalues.
	resizeArray(eigenvalues, n, 2);

	for (int i = 0; i < n; i++) {

		eigenvalues[i][0] = WR[i];
		eigenvalues[i][1] = WI[i];

		if (fabs(WI[i]) > 1.0e-04) {
			logFile
					<< "In commonFunctions::calcEigenValuesVectors computation\n "
					<< "eigenvalues failed due to complex eigenvalue " << WI[i]
					<< " !" << endl;

		}

	}

	// Store all eigenvectors.
	resizeArray(eigenvectors, n, n);

	for (int i = 0; i < n; i++)

		for (int j = 0; j < n; j++)

			eigenvectors[j][i] = vr[i * n + j];

	//delete[] WR,WI,matA,work,vl,vr;

	return info;
}

/***********************************************************************/
// Compute the singular value decomposition of a real general matrix
// A = USV^T.
int computeGeneralMatrixSVD(dbMatrix& A, dbVector& S, dbMatrix& U, dbMatrix& V,
		std::ofstream& logFile) {

	using namespace std;

	int size = A.size(); // size=4

	char jobu = 'A';
	char jobvt = 'A';

	int m = A.size();			// m=5
	int n = A[0].size();		// n=4
	int lda = A.size();		// lda=5
	int ldu = A.size();		// ldu=5
	int ldvt = A[0].size();	// ldvt=4
	int lwork = 10 * m;			// lwork=50

#ifdef _commonDebugMode_SVD_
	logFile << "m = " << m << endl;
	logFile << "n = " << n << endl;
	logFile << "lda = " << lda << endl;
	logFile << "ldu = " << ldu << endl;
	logFile << "ldvt = " << ldvt << endl;
	logFile << "lwork = " << lwork << endl;
#endif

	int lds;

	if (m < n)
		lds = m;		//lds = n = 4
	else
		lds = n;

	dbVector matA(lda * n); 	//matA.size()=20
	dbVector matS(lds);		//matS.size()=4
	unsigned long long ldu_ = (unsigned long long) (ldu);
	unsigned long long matU_size = (ldu_ * ldu_);

	//double* matU = new double[matU_size];

	//  logFile << "creating matU of size: " << 1e+2 << endl;
	//  dbVector matU_2(1e+2);
	//  dbVector().swap(matU_2);

	//  logFile << "creating matU of size: " << 1e+4 << endl;
	//  dbVector matU_4(1e+4);
	//  dbVector().swap(matU_4);

	//  logFile << "creating matU of size: " << 1e+6 << endl;
	//  dbVector matU_6(1e+6);
	//  dbVector().swap(matU_6);

	//  logFile << "creating matU of size: " << 1e+9 << endl;
	//  dbVector matU_9(1e+9);
	//  dbVector().swap(matU_9);

	double* matU = new double[matU_size];
	dbVector matVt(ldvt * n);	//matVt.size()=16
	dbVector work(lwork);		//work.size()=50

	// store the original matrix
	for (int i = 0; i < lda; i++)
		for (int j = 0; j < n; j++)
			matA[j * lda + i] = A[i][j];

#ifdef _commonDebugMode_SVD_
	logFile<<"******************************************************"<<endl;
	logFile<<"******************* original matrix ******************"<<endl;
	for(int i=0;i<A.size();i++)
	for(int j=0;j<A[i].size();j++)
	logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]<<endl;
	logFile<<"************** stored fortran matrix ****************"<<endl;
	for(int i=0;i<lda;i++)
	for(int j=0;j<n;j++)
	logFile<<"matA["<<i<<"]["<<j<<"] = "<<matA[i*n+j]<<endl;
#endif

	// compute the singular value decomposition
	int info = 0;
	dgesvd_(jobu, jobvt, m, n, &matA[0], lda, &matS[0], &matU[0], ldu,
			&matVt[0], ldvt, &work[0], lwork, info);

	if (info != 0) {
		cout << "In commonFunctions::computeGeneralSymmetricSVD singular\n"
				<< "value decomposition failed." << endl;
	}

	// store the singular values
	allocateArray(S, lds);

	bool flag = true;

	for (int i = 0; i < lds; i++) {
		S[i] = matS[i];
	}

	// unique and nonzero singular values
	if (flag) {

		// store the left and right singular vectors of A,
		// ie. U and V, respectively with
		// A = U diagA V^T

		resizeArray(U, ldu, ldu);
		resizeArray(V, ldvt, n);

		for (int i = 0; i < ldu; i++)

			for (int j = 0; j < ldu; j++)

				U[j][i] = matU[i * ldu + j];

		for (int i = 0; i < n; i++)

			for (int j = 0; j < ldvt; j++)

				V[i][j] = matVt[i * n + j]; // transpose it returned!

	}

	else {

		for (int i = 0; i < lds; i++)
			S[i] = A[i][i];

		U = getKroneckerSymbol(size);
		V = getKroneckerSymbol(size);

	}

	//delete[] matA,matS,matU,matVt,work;

#ifdef _commonDebugMode_SVD_
	logFile<<"******************************************************"<<endl;
	logFile<<"****************** singular values *******************"<<endl;
	for(int i=0;i<S.size();i++)
	logFile<<"S["<<i<<"] = "<<S[i]<<endl;
	logFile<<"*************** left singular vectors ****************"<<endl;
	for(int i=0;i<U.size();i++)
	for(int j=0;j<U[i].size();j++)
	logFile<<"U["<<i<<"]["<<j<<"] = "<<U[i][j]<<endl;
	logFile<<"*************** right singular vectors ****************"<<endl;
	for(int i=0;i<V.size();i++)
	for(int j=0;j<V[i].size();j++)
	logFile<<"V["<<i<<"]["<<j<<"] = "<<V[i][j]<<endl;
	logFile<<"********************** A = UDV^T **********************"<<endl;
	double value;
	for(int i=0;i<U.size();i++)
	for(int j=0;j<V.size();j++) {
		value = 0;
		for(int k=0;k<S.size();k++)
		value += U[i][k]*S[k]*V[j][k];
		logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]<<" ?= "<<value<<endl;
	}
#endif

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Compute the natural logarithm of a diagonalizable symmetric matrix
// using the spectral decomposition.
int spectralCalcMatrixLn(dbMatrix& A, dbMatrix& lnA, std::ofstream& logFile) {

	using namespace std;

	int size = A.size();
	dbVector E;
	dbMatrix V, invV;

	int info = calcEigenValuesVectors(A, E, V, logFile);

#ifdef _commonDebugMode_
	logFile<<"******************************************************"<<endl;
	logFile<<"****************** eigenvalues ***********************"<<endl;
	for(int i=0;i<size;i++) {
		logFile<<"E[i]= "<<E[i]<<endl;
	}
	logFile<<"******************** eigenvectors ********************"<<endl;
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	logFile<<"V["<<i<<"]["<<j<<"] = "<<V[i][j]<<endl;
#endif

	//   if(size == 3)
	//     calcInv3x3(V,invV,logFile);
	//   else
	calcInvDoubleSparse(V, invV, logFile);

	// A'_ij = V^{-1}_ik A_kl V_lj --> A' ... diagonalized matrix

	dbMatrix diagA;
	allocateArray(diagA, size, size);

	for (int i = 0; i < size; i++)

		for (int j = 0; j < size; j++)

			for (int k = 0; k < size; k++)

				for (int l = 0; l < size; l++)

					diagA[i][j] += invV[i][k] * A[k][l] * V[l][j];

#ifdef _commonDebugMode_
	logFile<<"******************************************************"<<endl;
	logFile<<"******************** eigenvectors ********************"<<endl;
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	logFile<<"V["<<i<<"]["<<j<<"] = "<<V[i][j]<<endl;
	logFile<<"*************** inverted eigenvectors ****************"<<endl;
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	logFile<<"invV["<<i<<"]["<<j<<"] = "<<invV[i][j]<<endl;
	logFile<<"***************** diagonalized matrix ****************"<<endl;
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	logFile<<"diagA["<<i<<"]["<<j<<"] = "<<diagA[i][j]<<endl;
	logFile<<"******************* D = e(ln(D)) *********************"<<endl;
	for(int i=0;i<size;i++)
	logFile<<"diagA["<<i<<"]["<<i<<"] = "<<diagA[i][i]<<" ?= "
	<<exp(log(diagA[i][i]))<<endl;
	logFile<<"******************** A = VDV^-1 **********************"<<endl;
	double value;
	dbMatrix dummy;
	allocateArray(dummy,size,size);
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++) {
		value = 0;
		for(int k=0;k<size;k++)
		for(int l=0;l<size;l++)
		value += V[i][k]*diagA[k][l]*invV[l][j];
		logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]<<" ?= "<<value<<endl;
	}
#endif

	// ln A = V B V^{-1}  --> B_{ii} = ln (A'_{ii})

	resizeArray(lnA, size, size);
	clearArray(lnA);

	for (int i = 0; i < size; i++)

		diagA[i][i] = log(diagA[i][i]);

	for (int i = 0; i < size; i++)

		for (int j = 0; j < size; j++)

			for (int k = 0; k < size; k++)

				for (int l = 0; l < size; l++)

					lnA[i][j] += V[i][k] * diagA[k][l] * invV[l][j];

#ifdef _commonDebugMode_
	logFile<<"************** ln of diagonalized matrix *************"<<endl;
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	logFile<<"ln(diagA)["<<i<<"]["<<j<<"] = "<<diagA[i][j]<<endl;
	logFile<<"******************* ln of given matrix ***************"<<endl;
	for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	logFile<<"lnA["<<i<<"]["<<j<<"] = "<<lnA[i][j]<<endl;
#endif

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Compute the natural logarithm of a diagonalizable non-symmetric matrix.
int calcNonsymMatrixLn(dbMatrix& A, dbMatrix& lnA, std::ofstream& logFile) {

	using namespace std;

	int size = A.size();
	int info = 0;

	bool usedSVD = false;

	// use singular value decomposition
	if (usedSVD) {

		dbVector S;
		dbMatrix U, V;
		info = computeGeneralMatrixSVD(A, S, U, V, logFile);

		if (info != 0) {
			logFile
					<< "In commonFunctions::calcNonsymMatrixExp singular value\n"
					<< "decomposition failed." << endl;
		}

		// ln(A) = U ln(S) V^T

		resizeArray(lnA, size, size);
		clearArray(lnA);

		for (int i = 0; i < size; i++)

			S[i] = log(S[i]);

		for (int i = 0; i < size; i++)

			for (int j = 0; j < size; j++)

				for (int k = 0; k < size; k++)

					lnA[i][j] += U[i][k] * S[k] * V[j][k];

#ifdef _commonDebugMode_
		logFile<<"**************** ln of given matrix ****************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"lnA["<<i<<"]["<<j<<"] = "<<lnA[i][j]<<endl;
		logFile<<"******************** exp(ln A) *********************"<<endl;
		dbMatrix expA,dexpA;
		dbMatrix3 d2expA;
		calcExponential(lnA,expA,dexpA,d2expA,0,logFile);
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]
		<<" ?= exp(ln(A))["<<i<<"]["<<j<<"] = "<<expA[i][j]<<endl;
#endif

	}

	/**********************************************************************/
	// use eigenvalues and vectors
	else {

		dbMatrix E, V, invV;
		info = calcEigenValuesVectors(A, E, V, logFile);

#ifdef _commonDebugMode_
		logFile<<"****************************************************"<<endl;
		logFile<<"****************** eigenvalues *********************"<<endl;
		for(int i=0;i<size;i++) {
			logFile<<"E_re[i]= "<<E[i][0]<<endl;
			logFile<<"E_im[i]= "<<E[i][1]<<endl;
		}
		logFile<<"****************** eigenvectors ********************"<<endl;
		for(int i=0;i<size;i++)
		for(int k=0;k<size;k++)
		logFile<<"V["<<i<<"]["<<k<<"] = "<<V[i][k]<<endl;
#endif

		calcInvDoubleSparse(V, invV, logFile);

		// A'_ij = V^{-1}_ik A_kl V_lj --> A' ... diagonalized matrix

		dbMatrix diagA;
		allocateArray(diagA, size, size);

		for (int i = 0; i < size; i++)

			for (int j = 0; j < size; j++)

				for (int k = 0; k < size; k++)

					for (int l = 0; l < size; l++)

						diagA[i][j] += invV[i][k] * A[k][l] * V[l][j];

#ifdef _commonDebugMode_
		logFile<<"******************************************************"<<endl;
		logFile<<"******************** eigenvectors ********************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"V["<<i<<"]["<<j<<"] = "<<V[i][j]<<endl;
		logFile<<"*************** inverted eigenvectors ****************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"invV["<<i<<"]["<<j<<"] = "<<invV[i][j]<<endl;
		logFile<<"***************** diagonalized matrix ****************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"diagA["<<i<<"]["<<j<<"] = "<<diagA[i][j]<<endl;
		logFile<<"******************* D = e(ln(D)) *********************"<<endl;
		for(int i=0;i<size;i++)
		logFile<<"diagA["<<i<<"]["<<i<<"] = "<<diagA[i][i]<<" ?= "
		<<exp(log(diagA[i][i]))<<endl;
		logFile<<"******************** A = VDV^-1 **********************"<<endl;
		double value;
		dbMatrix dummy;
		allocateArray(dummy,size,size);
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++) {
			value = 0;
			for(int k=0;k<size;k++)
			for(int l=0;l<size;l++)
			value += V[i][k]*diagA[k][l]*invV[l][j];
			logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]<<" ?= "<<value<<endl;
		}
#endif

		// ln A = V B V^{-1}  --> B_{ii} = ln (A'_{ii})

		resizeArray(lnA, size, size);
		clearArray(lnA);

		for (int i = 0; i < size; i++)

			diagA[i][i] = log(diagA[i][i]);

		for (int i = 0; i < size; i++)

			for (int j = 0; j < size; j++)

				for (int k = 0; k < size; k++)

					for (int l = 0; l < size; l++)

						lnA[i][j] += V[i][k] * diagA[k][l] * invV[l][j];

#ifdef _commonDebugMode_
		logFile<<"************ ln of diagonalized matrix *************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"ln(diagA)["<<i<<"]["<<j<<"] = "<<diagA[i][j]<<endl;
		logFile<<"***************** ln of given matrix ***************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"lnA["<<i<<"]["<<j<<"] = "<<lnA[i][j]<<endl;
		logFile<<"******************** exp(ln A) *********************"<<endl;
		dbMatrix expA,dexpA;
		dbMatrix3 d2expA;
		calcExponential(lnA,expA,dexpA,d2expA,0,logFile);
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]
		<<" ?= exp(ln(A))["<<i<<"]["<<j<<"] = "<<expA[i][j]<<endl;
#endif

	}

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Compute the natural logarithm of a diagonalizable non-symmetric
// matrix using a Taylor series with an initial value A0 = 1
int taylorApproxDenseMatrixLn(dbMatrix& A, dbMatrix& lnA,
		std::ofstream& logFile) {

	using namespace std;

	int info = 0;

	if (A.size() == 0) {

		logFile << "In commonFunctions::taylorApproxDenseMatrixLn no valid\n"
				<< "matrix to calculate logarithm!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (A.size() != A[0].size()) {

		logFile << "In commonFunctions::taylorApproxDenseMatrixLn no valid\n"
				<< "matrix to calculate logarithm!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (lnA.size() != A.size())
		allocateArray(lnA, A.size(), A.size());

	clearArray(lnA);

	// ln (A) = ln (1 - W) =
	//
	//           = -W -1/2 W^2 -1/3 W^3 - ....
	//
	// with W = 1 - A

	dbMatrix W = getKroneckerSymbol(A.size());

	// lnA -= W
	for (int i = 0; i < A.size(); i++)

		for (int j = 0; j < A[i].size(); j++)

			W[i][j] -= A[i][j];

	if (fabs(computeNorm(W, 1, logFile)) >= 1.0) {
		logFile << "In commonFunctions::taylorApproxMatrixLn norm of input\n"
				<< "matrix || 1 - A || must be smaller than 2." << endl;
		for (int i = 0; i < A.size(); i++)
			for (int j = 0; j < A[i].size(); j++)
				logFile << "A[" << i << "][" << j << "]=" << A[i][j] << endl;
		for (int i = 0; i < W.size(); i++)
			for (int j = 0; j < W[i].size(); j++)
				logFile << "W[" << i << "][" << j << "]=" << W[i][j] << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// apply Cayley-Hamilton theorem to the Taylor series
	//
	// ln(1-W) = alpha0(I1,I2,I3) 1 + alpha1(I1,I2,I3) W + alpha2(I1,I2,I3) W^2

	if (W.size() == 3) {

		int Wsize = W.size();
		dbMatrix delta = getKroneckerSymbol(Wsize);
		dbMatrix WW;
		innerTensorProduct(W, W, WW, false, false, logFile);

		// compute the invariants

		// I1 = tr W
		double I1 = 0;

		for (int i = 0; i < Wsize; i++)
			I1 += W[i][i];

		// I2 = 0.5 [ (tr W)^2 - tr W^2 ] = 0.5*[ I1^2 - W_ik W_ki ]

		double I2;
		double trW2 = 0;

		for (int i = 0; i < Wsize; i++)

			trW2 += WW[i][i];

		I2 = 0.5 * (pow(I1, 2.0) - trW2);

		// I3 = det W
		double I3;
		calcDetDoubleDense(W, I3, logFile);

		// -------------------------------------------------------------------
		// compute initial alpha0,alpha1,alpha2

		double n = 3;
		double tol = 1.0;
		double oldGamma0, oldGamma1, oldGamma2;
		double gamma0, gamma1, gamma2;

		gamma0 = I3;
		gamma1 = -I2;
		gamma2 = I1;

		double alpha0 = -1.0 / 3.0 * gamma0;
		double alpha1 = -1.0 - 1.0 / 3.0 * gamma1;
		double alpha2 = -1.0 / 2.0 - 1.0 / 3.0 * gamma2;

		oldGamma0 = gamma0;
		oldGamma1 = gamma1;
		oldGamma2 = gamma2;

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		logFile<<"I1="<<I1<<endl;
		logFile<<"I2="<<I2<<endl;
		logFile<<"I3="<<I3<<endl;
		logFile<<"============================"<<endl;
		logFile<<"n="<<n<<endl;
		logFile<<"tol="<<fabs(oldGamma0-gamma0)+fabs(oldGamma1-gamma1)
		+fabs(oldGamma2-gamma2)<<endl;
		logFile<<"----------------------------"<<endl;
		logFile<<"oldGamma0="<<oldGamma0<<endl;
		logFile<<"oldGamma1="<<oldGamma1<<endl;
		logFile<<"oldGamma2="<<oldGamma2<<endl;
		logFile<<"----------------------------"<<endl;
		logFile<<"gamma0="<<gamma0<<endl;
		logFile<<"gamma1="<<gamma1<<endl;
		logFile<<"gamma2="<<gamma2<<endl;
		logFile<<"----------------------------"<<endl;
		logFile<<"alpha0="<<alpha0<<endl;
		logFile<<"alpha1="<<alpha1<<endl;
		logFile<<"alpha2="<<alpha2<<endl;
#endif

		n++;

		// -------------------------------------------------------------------
		// correct alpha0,alpha1,alpha2 iteratively
		while (tol > 1.0e-8 && n < 100) {

			gamma0 = I3 * oldGamma2;
			gamma1 = oldGamma0 - I2 * oldGamma2;
			gamma2 = oldGamma1 + I1 * oldGamma2;

			// alpha0 += -\sum\limits_{n=4}^N\,\frac{1}{n}\,\gamma_0^{(n)};
			alpha0 -= 1.0 / n * gamma0;

			// alpha1 += -\sum\limits_{n=4}^N\,\frac{1}{n}\,\gamma_1^{(n)};
			alpha1 -= 1.0 / n * gamma1;

			// alpha2 += -\sum\limits_{n=4}^N\,\frac{1}{n}\,\gamma_2^{(n)}
			alpha2 -= 1.0 / n * gamma2;

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
			logFile<<"============================"<<endl;
			logFile<<"n="<<n<<endl;
			logFile<<"tol="<<fabs(oldGamma0-gamma0)+fabs(oldGamma1-gamma1)
			+fabs(oldGamma2-gamma2)<<endl;
			logFile<<"----------------------------"<<endl;
			logFile<<"oldGamma0="<<oldGamma0<<endl;
			logFile<<"oldGamma1="<<oldGamma1<<endl;
			logFile<<"oldGamma2="<<oldGamma2<<endl;
			logFile<<"----------------------------"<<endl;
			logFile<<"gamma0="<<gamma0<<endl;
			logFile<<"gamma1="<<gamma1<<endl;
			logFile<<"gamma2="<<gamma2<<endl;
			logFile<<"----------------------------"<<endl;
			logFile<<"alpha0="<<alpha0<<endl;
			logFile<<"alpha1="<<alpha1<<endl;
			logFile<<"alpha2="<<alpha2<<endl;
#endif

			oldGamma0 = gamma0;
			oldGamma1 = gamma1;
			oldGamma2 = gamma2;

			tol = fabs(oldGamma0 - gamma0) + fabs(oldGamma1 - gamma1)
					+ fabs(oldGamma2 - gamma2);
			n++;

		}

		if (n >= 100) {
			info = -1;
			logFile
					<< "In commonFunctions::taylorApproxMatrixLn approximation of "
					<< "matrix logarithm has not converged after\n" << n
					<< " iterations: error = " << tol << " !" << endl;
		}

		// ln A = alpha0*1 + alpha1*W + alpha2*W^2

		for (int i = 0; i < Wsize; i++)

			for (int j = 0; j < Wsize; j++)

				lnA[i][j] = alpha0 * delta[i][j] + alpha1 * W[i][j]
						+ alpha2 * WW[i][j];

	}

	/**********************************************************************/
	// just the Taylor series
	else {

		// lnA -= W
		for (int i = 0; i < A.size(); i++)

			for (int j = 0; j < A[i].size(); j++)

				lnA[i][j] -= W[i][j];

		// -----------------------------------------------------------
		// add higher order terms until Frobenius norm small enough:

		int n = 2;

		dbMatrix WW;
		dbMatrix WWW = W;
		dbMatrix E = W;

		while (fabs(computeNorm(E, 1, logFile)) > 1.0e-08 && n < 100) {

			WW = WWW;
			innerTensorProduct(WW, W, WWW, false, false, logFile);

			// lnA -= 1/n W^n
			for (int i = 0; i < A.size(); i++) {
				for (int j = 0; j < A.size(); j++) {

					E[i][j] = 1.0 / n * WWW[i][j];
					lnA[i][j] -= E[i][j];
				}
			}

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
			logFile<<"----------------------------------------------------"<<endl;
			logFile<<"n="<<n<<endl;
			logFile<<"error="<<fabs(computeNorm(E,1,logFile))<<endl;
			for(int i=0;i<E.size();i++)
			for(int j=0;j<E[i].size();j++)
			logFile<<"E["<<i<<"]["<<j<<"]="<<E[i][j]<<endl;
			for(int i=0;i<lnA.size();i++)
			for(int j=0;j<lnA[i].size();j++)
			logFile<<"lnA["<<i<<"]["<<j<<"]="<<lnA[i][j]<<endl;
#endif

			n++;

		}

		if (n >= 100) {
			info = -1;
			logFile
					<< "In commonFunctions::taylorApproxMatrixLn approximation of "
					<< "matrix logarithm has not converged after\n" << n
					<< " iterations: error = "
					<< fabs(computeNorm(E, 1, logFile)) << " !" << endl;
		}

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		dbMatrix lnAtest;
		calcNonsymMatrixLn(A,lnAtest,logFile);
		logFile<<"----------------------------------------------------"<<endl;
		for(int i=0;i<lnA.size();i++)
		for(int j=0;j<lnA[i].size();j++)
		logFile<<"lnA["<<i<<"]["<<j<<"]="<<lnA[i][j]<<" ?= "
		<<lnAtest[i][j]<<endl;
#endif

	}

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Compute the natural logarithm of a diagonalizable non-symmetric
// matrix (Cheng, Higham, Kenney and Laub 2001, p. 1118).
int highamApproxDenseMatrixLn(dbMatrix& A, dbMatrix& lnA,
		std::ofstream& logFile) {

	using namespace std;

	int info = 0;

	if (A.size() == 0) {

		logFile << "In commonFunctions::highamApproxMatrixLn no valid matrix\n"
				<< "to calculate logarithm!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (A.size() != A[0].size()) {

		logFile << "In commonFunctions::highamApproxMatrixLn no valid matrix\n"
				<< "to calculate logarithm!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/**********************************************************************/

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
	dbMatrix E,V,invV;
	info = calcEigenValuesVectors(A,E,V,logFile);
	logFile<<"****************************************************"<<endl;
	logFile<<"****************** eigenvalues *********************"<<endl;
	for(int i=0;i<E.size();i++) {
		logFile<<"E_re[i]= "<<E[i][0]<<endl;
		logFile<<"E_im[i]= "<<E[i][1]<<endl;
	}
	for(int i=0;i<E.size();i++) {

		if(E[i][0] < 0 && E[i][1] != 0) {
			logFile <<"In commonFunctions::highamApproxMatrixLn no valid matrix\n"
			<<"to calculate logarithm, because real part of eigenvalue\n"
			<<"smaller than 0!"<<endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
#endif

	/**********************************************************************/

	if (lnA.size() != A.size())
		allocateArray(lnA, A.size(), A[0].size());

	clearArray(lnA);

	dbMatrix delta = getKroneckerSymbol(A.size());

	int iStep = 0;
	int kStep = 0;
	double Y1norm, M1norm, gamma, gamma2, detM;
	dbMatrix Y, M, Yold, Mold, Y1, M1, MoldInv, lnY, lnM;

	allocateArray(Y, A.size(), A[0].size());
	allocateArray(M, A.size(), A[0].size());
	allocateArray(Y1, A.size(), A[0].size());
	allocateArray(M1, A.size(), A[0].size());

	double n = A.size();
	double dtol = 1.0e-08 * fabs(computeNorm(A, 1, logFile)) / 4.0;

	Y = A;

	// Y1 = || 1 - A ||
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < A[i].size(); j++)
			Y1[i][j] = delta[i][j] - Y[i][j];

	Y1norm = fabs(computeNorm(Y1, 1, logFile));

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
	logFile<<"******************************************************"<<endl;
	logFile<<iStep<<".iStep: Y1norm="<<Y1norm<<" dtol="<<dtol<<endl;
#endif

	// iterate until Y is small enough
	while (Y1norm >= 1.0 && iStep < 100) {

		iStep += 1;

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		logFile<<"****************************************************"<<endl;
		logFile<<iStep<<".iStep: Y1norm="<<Y1norm<<endl;
#endif

		if (iStep == 1) {
			Yold = A;
			Mold = A;
		} else {
			Yold = Y;
			Mold = Y;
		}

		for (int i = 0; i < A.size(); i++)
			for (int j = 0; j < A[i].size(); j++)
				M1[i][j] = delta[i][j] - Mold[i][j];

		M1norm = fabs(computeNorm(M1, 1, logFile));

		if (M1norm < 1.0)
			M1norm = fabs(M1norm + log(1.0 - M1norm));

		kStep = 0;

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		logFile<<"-------------------------------------------------"<<endl;
		logFile<<kStep<<".kStep: M1norm="<<M1norm<<endl;
#endif

		// -------------------------------------------------------------------

		while (M1norm >= dtol / pow(4.0, iStep - 1.0) && kStep < 100) {

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
			logFile<<"-------------------------------------------------"<<endl;
			logFile<<kStep<<".kStep: M1norm="<<M1norm<<endl;
#endif

			calcDetDoubleDense(Mold, detM, logFile);
			calcInvDoubleDense(Mold, MoldInv, logFile);

			gamma = fabs(pow(detM, -1.0 / (2.0 * n)));
			gamma2 = pow(gamma, 2.0);

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
			logFile<<"-------------------------------------------------"<<endl;
			for(int i=0;i<A.size();i++)
			for(int j=0;j<A[i].size();j++)
			logFile<<"Mold["<<i<<"]["<<j<<"]="<<Mold[i][j]<<endl;
			for(int i=0;i<A.size();i++)
			for(int j=0;j<A[i].size();j++)
			logFile<<"MoldInv["<<i<<"]["<<j<<"]="<<MoldInv[i][j]<<endl;
			logFile<<"det(Mold)="<<detM<<endl;
			logFile<<"gamma="<<gamma<<endl;
#endif

			for (int i = 0; i < A.size(); i++)
				for (int j = 0; j < A[i].size(); j++)

					M[i][j] = 0.5
							* (delta[i][j]
									+ 0.5
											* (gamma2 * Mold[i][j]
													+ MoldInv[i][j] / gamma2));

			for (int i = 0; i < A.size(); i++)
				for (int j = 0; j < A[i].size(); j++)

					Y[i][j] = 0.5 * gamma * Yold[i][j]
							* (delta[i][j] + MoldInv[i][j] / gamma2);

			// ----
			// prepare for next kStep-iteration
			Yold = Y;
			Mold = M;

			for (int i = 0; i < A.size(); i++)
				for (int j = 0; j < A[i].size(); j++)
					M1[i][j] = delta[i][j] - Mold[i][j];

			M1norm = fabs(computeNorm(M1, 1, logFile));

			if (M1norm < 1.0)
				M1norm = fabs(M1norm + log(1.0 - M1norm));

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
			logFile<<"-------------------------------------------------"<<endl;
			logFile<<"M1norm="<<M1norm<<" <? "<<dtol/pow(4.0,iStep-1.0)<<endl;
			for(int i=0;i<A.size();i++)
			for(int j=0;j<A[i].size();j++)
			logFile<<"M["<<i<<"]["<<j<<"]="<<M[i][j]<<endl;
			for(int i=0;i<A.size();i++)
			for(int j=0;j<A[i].size();j++)
			logFile<<"Y["<<i<<"]["<<j<<"]="<<Y[i][j]<<endl;
#endif

			kStep++;

			if (kStep >= 100) {
				info = -1;
				break;
			}

		}

		// -------------------------------------------------------------------
		// compute termination norm of iStep-loop
		for (int i = 0; i < A.size(); i++)
			for (int j = 0; j < A[i].size(); j++)
				Y1[i][j] = delta[i][j] - Y[i][j];

		Y1norm = fabs(computeNorm(Y1, 1, logFile));

		// iStep-loop not converged
		if (info != 0 || iStep >= 100 || kStep == 0) {
			info = -1;
			logFile
					<< "In commonFunctions::highamApproxMatrixLn approximation of "
					<< "matrix logarithm has not converged after\n" << iStep
					<< " outer and " << kStep << " inner iterations: error = "
					<< Y1norm << " and " << M1norm << ", respectively." << endl;
			break;
		}
		// lnA += - ln(M_1) - 2ln(M_2) - ... - 2^{s-1}*ln(M^s)
		else {

			//       approxMatrixLn(M,lnM,logFile);

			//       for(int i=0;i<A.size();i++)
			// 	for(int j=0;j<A[i].size();j++)

			// 	  lnA[i][j] -= pow(2.0,iStep-1.0)*lnM[i][j];

			for (int i = 0; i < A.size(); i++)
				for (int j = 0; j < A[i].size(); j++)

					lnA[i][j] -= pow(2.0, iStep - 1.0)
							* (M[i][j] - delta[i][j]);

		}

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		logFile<<"---------"<<endl;
		logFile<<"Y1norm="<<Y1norm<<" <? 0.99 ("<<iStep<<".iStep)"<<endl;
		for(int i=0;i<A.size();i++)
		for(int j=0;j<A[i].size();j++)
		logFile<<"Y["<<i<<"]["<<j<<"]="<<Y[i][j]<<endl;
		logFile<<"---------"<<endl;
		for(int i=0;i<A.size();i++)
		for(int j=0;j<A[i].size();j++)
		logFile<<"Y1["<<i<<"]["<<j<<"]="<<Y1[i][j]<<endl;
		logFile<<"---------"<<endl;
		for(int i=0;i<A.size();i++)
		for(int j=0;j<A[i].size();j++)
		logFile<<"lnA["<<i<<"]["<<j<<"]="<<lnA[i][j]<<endl;
#endif

	}

	// lnA = 2^s*ln(Y_s) - ln(M_1) - 2ln(M_2) - ... - 2^{s-1}*ln(M^s)

	if (iStep > 0 && info == 0) {

		info = taylorApproxDenseMatrixLn(Y, lnY, logFile);

		for (int i = 0; i < A.size(); i++)
			for (int j = 0; j < A[i].size(); j++)

				lnA[i][j] += pow(2.0, iStep) * lnY[i][j];

	}

	else if (info == 0)

		info = taylorApproxDenseMatrixLn(A, lnA, logFile);

	if (info != 0) {
		dbMatrix E, V, invV;
		int flag = calcEigenValuesVectors(A, E, V, logFile);
		logFile << "****************************************************"
				<< endl;
		logFile << "****************** eigenvalues *********************"
				<< endl;
		for (int i = 0; i < E.size(); i++) {
			logFile << "E_re[i]= " << E[i][0] << endl;
			logFile << "E_im[i]= " << E[i][1] << endl;
		}
		for (int i = 0; i < E.size(); i++) {

			if (E[i][0] < 0 && E[i][1] != 0) {
				logFile
						<< "In commonFunctions::highamApproxMatrixLn no valid matrix to\n"
						<< "calculate logarithm, because real part of eigenvalue\n"
						<< "smaller than 0!" << endl;
			}
		}
	}

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
	logFile<<"iStep"<<iStep<<endl;
	logFile<<"kStep"<<kStep<<endl;
	logFile<<"******************************************************"<<endl;
	dbMatrix lnB;
	calcNonsymMatrixLn(A,lnB,logFile);
	for(int i=0;i<A.size();i++)
	for(int j=0;j<A[i].size();j++)
	logFile<<"lnA["<<i<<"]["<<j<<"]="<<lnA[i][j]<<" ?= "
	<<lnB[i][j]<<endl;
#endif

	return info;
}

/***********************************************************************/
/***********************************************************************/
// Compute the exponential of a diagonalizable non-symmetric
// dense matrix using a Taylor series with an initial value A0 = 0
int taylorApproxDenseMatrixExp(dbMatrix& A, dbMatrix& expA,
		std::ofstream& logFile) {

	using namespace std;

	int info = 0;

	if (A.size() == 0) {

		logFile << "In commonFunctions::approxMatrixExp no valid matrix to\n"
				<< "calculate exponential!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (A.size() != A[0].size()) {

		logFile << "In commonFunctions::approxMatrixExp no valid matrix to\n"
				<< "calculate exponential!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (expA.size() < A.size())
		allocateArray(expA, A.size(), A.size());

	clearArray(expA);

	double n, tol;

	// apply Cayley-Hamilton theorem to the Taylor series
	//
	// exp A = alpha0(I1,I2,I3) 1 + alpha1(I1,I2,I3) A + alpha2(I1,I2,I3) A^2

	if (A.size() == 3) {

		int Asize = A.size();
		dbMatrix delta = getKroneckerSymbol(Asize);
		dbMatrix AA;
		innerTensorProduct(A, A, AA, false, false, logFile);

		// compute the invariants

		// I1 = tr A
		double I1 = 0;

		for (int i = 0; i < Asize; i++)
			I1 += A[i][i];

		// I2 = 0.5 [ (tr A)^2 - tr A^2 ] = 0.5*[ I1^2 - A_ik A_ki ]

		double I2;
		double trA2 = 0;

		for (int i = 0; i < Asize; i++)

			trA2 += AA[i][i];

		I2 = 0.5 * (pow(I1, 2) - trA2);

		// I3 = 1/3*(trA^3 - 3/2*trA*trA^2 + 1/2*(trA)^3)
		double I3 = 0;
		double trA3 = 0;

		// trA^3 = A_ik*A_kl*A_li

		for (int i = 0; i < Asize; i++)

			for (int k = 0; k < Asize; k++)

				for (int l = 0; l < Asize; l++)

					trA3 += A[i][k] * A[k][l] * A[l][i];

		I3 = 1.0 / 3.0 * (trA3 - 3.0 / 2.0 * I1 * trA2 + 0.5 * pow(I1, 3.0));

		// -------------------------------------------------------------------
		// compute initial alpha0,alpha1,alpha2

		n = 3;
		double fac = 1.0 * 2.0 * 3.0;
		tol = 1.0;
		double oldGamma0, oldGamma1, oldGamma2;
		double gamma0, gamma1, gamma2;

		gamma0 = I3;
		gamma1 = -I2;
		gamma2 = I1;

		// alpha0 = 1 + (1/3!)*I3   + sum_{n=4}^N (1/n!)*gamma0^(n)
		// alpha1 = 1 - (1/3!)*I2   + sum_{n=4}^N (1/n!)*gamma1^(n)
		// alpha2 = 1/2 + (1/3!)*I1 + sum_{n=4}^N (1/n!)*gamma2^(n)

		double alpha0 = 1.0 + 1.0 / fac * I3;
		double alpha1 = 1.0 - 1.0 / fac * I2;
		double alpha2 = 0.5 + 1.0 / fac * I1;

		oldGamma0 = gamma0;
		oldGamma1 = gamma1;
		oldGamma2 = gamma2;

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
		logFile<<"I1="<<I1<<endl;
		logFile<<"I2="<<I2<<endl;
		logFile<<"I3="<<I3<<endl;
		logFile<<"============================"<<endl;
		logFile<<"n="<<n<<endl;
		logFile<<"fac="<<fac<<endl;
		logFile<<"tol="<<fabs(oldGamma0-gamma0)+fabs(oldGamma1-gamma1)
		+fabs(oldGamma2-gamma2)<<endl;
		logFile<<"----------------------------"<<endl;
		logFile<<"oldGamma0="<<oldGamma0<<endl;
		logFile<<"oldGamma1="<<oldGamma1<<endl;
		logFile<<"oldGamma2="<<oldGamma2<<endl;
		logFile<<"----------------------------"<<endl;
		logFile<<"gamma0="<<gamma0<<endl;
		logFile<<"gamma1="<<gamma1<<endl;
		logFile<<"gamma2="<<gamma2<<endl;
		logFile<<"----------------------------"<<endl;
		logFile<<"alpha0="<<alpha0<<endl;
		logFile<<"alpha1="<<alpha1<<endl;
		logFile<<"alpha2="<<alpha2<<endl;
#endif

		n++;

		// -------------------------------------------------------------------
		// correct alpha0,alpha1,alpha2 iteratively
		while (tol > 1.0e-8 && n < 100) {

			fac *= n;
			gamma0 = I3 * oldGamma2;
			gamma1 = oldGamma0 - I2 * oldGamma2;
			gamma2 = oldGamma1 + I1 * oldGamma2;

			// alpha0 += sum_{n=4}^N (1/n!)*gamma0^(n);
			alpha0 += 1.0 / fac * gamma0;

			// alpha1 += sum_{n=4}^N (1/n!)*gamma1^(n);
			alpha1 += 1.0 / fac * gamma1;

			// alpha2 += sum_{n=4}^N (1/n!)*gamma2^(n)
			alpha2 += 1.0 / fac * gamma2;

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
			logFile<<"============================"<<endl;
			logFile<<"n="<<n<<endl;
			logFile<<"fac="<<fac<<endl;
			logFile<<"tol="<<fabs(oldGamma0-gamma0)+fabs(oldGamma1-gamma1)
			+fabs(oldGamma2-gamma2)<<endl;
			logFile<<"----------------------------"<<endl;
			logFile<<"oldGamma0="<<oldGamma0<<endl;
			logFile<<"oldGamma1="<<oldGamma1<<endl;
			logFile<<"oldGamma2="<<oldGamma2<<endl;
			logFile<<"----------------------------"<<endl;
			logFile<<"gamma0="<<gamma0<<endl;
			logFile<<"gamma1="<<gamma1<<endl;
			logFile<<"gamma2="<<gamma2<<endl;
			logFile<<"----------------------------"<<endl;
			logFile<<"alpha0="<<alpha0<<endl;
			logFile<<"alpha1="<<alpha1<<endl;
			logFile<<"alpha2="<<alpha2<<endl;
#endif

			oldGamma0 = gamma0;
			oldGamma1 = gamma1;
			oldGamma2 = gamma2;

			tol = fabs(oldGamma0 - gamma0) + fabs(oldGamma1 - gamma1)
					+ fabs(oldGamma2 - gamma2);
			n++;

		}

		if (n >= 100) {
			info = -1;
			logFile
					<< "In commonFunctions::taylorApproxMatrixLn approximation of "
					<< "matrix logarithm has not converged after\n" << n
					<< " iterations: error = " << tol << " !" << endl;
		}

		// exp A = alpha0*1 + alpha1*A + alpha2*A^2

		for (int i = 0; i < Asize; i++)

			for (int j = 0; j < Asize; j++)

				expA[i][j] = alpha0 * delta[i][j] + alpha1 * A[i][j]
						+ alpha2 * AA[i][j];

	}

	/**********************************************************************/
	// just the Taylor series
	else {

		// exp(A) = sum_{n=0}^\infty 1/n! A^n = 1 + A + A^2/(2!) + A^3/(3!) ...

		// expA = 1 + A
		expA = A;

		for (int i = 0; i < A.size(); i++)
			expA[i][i] += 1.0;

#ifdef _commonDebugMode_
		logFile<<"----------------------------------------------------"<<endl;
		logFile<<"n=1"<<endl;
		logFile<<"factorial=1"<<endl;
		for(int i=0;i<expA.size();i++)
		for(int j=0;j<expA[i].size();j++)
		logFile<<"expA["<<i<<"]["<<j<<"]="<<expA[i][j]<<endl;
#endif

		// expA += sum_{n=2}^\infty 1/n! A^n = 1 + A + A^2/(2!) + A^3/(3!) ...
		dbMatrix Aold, Anew, E;

		int iterationStep = 1;
		n = 2.0;
		double factorial = 1.0;
		Anew = A;
		E = getKroneckerSymbol(A.size());

		while (tol > 1.0e-13 && n < 100) {

			factorial *= n;
			Aold = Anew;

			clearArray(Anew);

			for (int i = 0; i < A.size(); i++) {
				for (int j = 0; j < A.size(); j++) {
					for (int k = 0; k < A.size(); k++)
						Anew[i][j] += Aold[i][k] * A[k][j];

					E[i][j] = 1.0 / factorial * Anew[i][j];
					expA[i][j] += E[i][j];
				}
			}

#ifdef _commonDebugMode_
			logFile<<"----------------------------------------------------"<<endl;
			logFile<<"n="<<n<<endl;
			logFile<<"factorial="<<factorial<<endl;
			for(int i=0;i<E.size();i++)
			for(int j=0;j<E[i].size();j++)
			logFile<<"E["<<i<<"]["<<j<<"]="<<E[i][j]<<endl;
			for(int i=0;i<expA.size();i++)
			for(int j=0;j<expA[i].size();j++)
			logFile<<"expA["<<i<<"]["<<j<<"]="<<expA[i][j]<<endl;
#endif

			tol = fabs(computeNorm(E, 3, logFile));
			n++;

		}

	}

	if (n >= 100) {
		info = -1;
		logFile << "In commonFunctions::taylorApproxMatrixExp approximation "
				<< "of matrix expontial\n" << "has not converged after " << n
				<< " iterations" << ": error = " << tol << " !" << endl;
	}

#if defined _commonDebugMode_ || defined _kinematicsDebugMode_
	dbMatrix expAtest;
	int flag =
	calcNonsymMatrixExp(A,expAtest,logFile);
	logFile<<"******************************************************"<<endl;
	logFile<<"*********************** exp(A) ***********************"<<endl;
	for(int r=0;r<expA.size();r++)
	for(int s=0;s<expA[r].size();s++)
	logFile<<"exp(A)["<<r<<"]["<<s<<"] = "
	<<expA[r][s]<<" ?= "<<expAtest[r][s]<<endl;
#endif

	return info;

}

/***********************************************************************/
/***********************************************************************/
// Compute the exponential of a diagonalizable non-symmetric matrix.
int calcNonsymMatrixExp(dbMatrix& A, dbMatrix& expA, std::ofstream& logFile) {

	using namespace std;

	int info = 0;
	int size = A.size();

	bool usedSVD = false;

	// use singular value decomposition
	if (usedSVD) {

		dbVector S;
		dbMatrix U, V;
		info = computeGeneralMatrixSVD(A, S, U, V, logFile);

		if (info != 0) {
			logFile
					<< "In commonFunctions::calcNonsymMatrixExp singular value\n"
					<< "decomposition failed." << endl;
		}

		// exp(A) = U exp(S) V^T

		resizeArray(expA, size, size);
		clearArray(expA);

		for (int i = 0; i < size; i++)

			S[i] = exp(S[i]);

		for (int i = 0; i < size; i++)

			for (int j = 0; j < size; j++)

				for (int k = 0; k < size; k++)

					expA[i][j] += U[i][k] * S[k] * V[j][k];

#ifdef _commonDebugMode_
		logFile<<"******************* exp of given matrix ***************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"expA["<<i<<"]["<<j<<"] = "<<expA[i][j]<<endl;
		logFile<<"********************** ln(exp A) *********************"<<endl;
		dbMatrix lnA;
		calcNonsymMatrixLn(expA,lnA,logFile);
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]
		<<" ?= ln(exp(A))["<<i<<"]["<<j<<"] = "<<lnA[i][j]<<endl;
#endif

	}

	/**********************************************************************/
	// use eigenvectors
	else {

		dbMatrix E, V, invV;
		info = calcEigenValuesVectors(A, E, V, logFile);

#ifdef _commonDebugMode_
		logFile<<"****************************************************"<<endl;
		logFile<<"****************** eigenvalues *********************"<<endl;
		for(int i=0;i<size;i++) {
			logFile<<"E_re[i]= "<<E[i][0]<<endl;
			logFile<<"E_im[i]= "<<E[i][1]<<endl;
		}
		logFile<<"******************** eigenvectors ******************"<<endl;
		for(int i=0;i<size;i++)
		for(int k=0;k<size;k++)
		logFile<<"V["<<i<<"]["<<k<<"] = "<<V[i][k]<<endl;
#endif

		calcInvDoubleSparse(V, invV, logFile);

		// A'_ij = V^{-1}_ik A_kl V_lj --> A' ... diagonalized matrix

		dbMatrix diagA;
		allocateArray(diagA, size, size);

		for (int i = 0; i < size; i++)

			for (int j = 0; j < size; j++)

				for (int k = 0; k < size; k++)

					for (int l = 0; l < size; l++)

						diagA[i][j] += invV[i][k] * A[k][l] * V[l][j];

#ifdef _commonDebugMode_
		logFile<<"****************************************************"<<endl;
		logFile<<"******************** eigenvectors ******************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"V["<<i<<"]["<<j<<"] = "<<V[i][j]<<endl;
		logFile<<"*************** inverted eigenvectors **************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"invV["<<i<<"]["<<j<<"] = "<<invV[i][j]<<endl;
		logFile<<"***************** diagonalized matrix **************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"diagA["<<i<<"]["<<j<<"] = "<<diagA[i][j]<<endl;
		logFile<<"******************* D = e(ln(D)) *******************"<<endl;
		for(int i=0;i<size;i++)
		logFile<<"diagA["<<i<<"]["<<i<<"] = "<<diagA[i][i]<<" ?= "
		<<exp(log(diagA[i][i]))<<endl;
		logFile<<"******************** A = VDV^-1 ********************"<<endl;
		double value;
		dbMatrix dummy;
		allocateArray(dummy,size,size);
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++) {
			value = 0;
			for(int k=0;k<size;k++)
			for(int l=0;l<size;l++)
			value += V[i][k]*diagA[k][l]*invV[l][j];
			logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]<<" ?= "<<value<<endl;
		}
#endif

		// ln A = V B V^{-1}  --> B_{ii} = ln (A'_{ii})

		resizeArray(expA, size, size);
		clearArray(expA);

		for (int i = 0; i < size; i++)

			diagA[i][i] = exp(diagA[i][i]);

		for (int i = 0; i < size; i++)

			for (int j = 0; j < size; j++)

				for (int k = 0; k < size; k++)

					for (int l = 0; l < size; l++)

						expA[i][j] += V[i][k] * diagA[k][l] * invV[l][j];

#ifdef _commonDebugMode_
		logFile<<"*********** exp of diagonalized matrix *************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"exp(diagA)["<<i<<"]["<<j<<"] = "<<diagA[i][j]<<endl;
		logFile<<"***************** exp of given matrix **************"<<endl;
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"expA["<<i<<"]["<<j<<"] = "<<expA[i][j]<<endl;
		logFile<<"******************** ln(exp A) *********************"<<endl;
		dbMatrix lnA;
		calcNonsymMatrixLn(expA,lnA,logFile);
		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		logFile<<"A["<<i<<"]["<<j<<"] = "<<A[i][j]
		<<" ?= ln(exp(A))["<<i<<"]["<<j<<"] = "<<lnA[i][j]<<endl;
#endif

	}

	return info;
}

/************************************************************************/
/************************************************************************/
// Compute the exponential of a matrix.
//void calcExponential(dbMatrix& A, dbMatrix& expA, dbMatrix& dexpA,
//		dbMatrix3& d2expA, int derivOrder, std::ofstream& logFile) {
//
//	using namespace std;
//
//	int rank;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	int usedDims = 3;
//	int fullSize = A.size() * A.size();
//
//	intMatrix vFull = matrixToVectorFull(A.size());
//	dbMatrix delta = getKroneckerSymbol(A.size());
//
//	// tensor^2,tensor^3
//	dbMatrix AA, AAA;
//	innerTensorProduct(A, A, AA, false, false, logFile);
//	innerTensorProduct(AA, A, AAA, false, false, logFile);
//
//	/**********************************************************************/
//
//	bool dexpo = true;
//
//	if (!dexpo) {
//
//		if (expA.size() < A.size())
//			allocateArray(expA, A.size(), A.size());
//
//		clearArray(expA);
//
//		for (int i = 0; i < A.size(); i++) {
//
//			expA[i][i] += 1.0;
//
//			for (int j = 0; j < A.size(); j++) {
//
//				expA[i][j] += A[i][j];
//				expA[i][j] += 1.0 / 2.0 * AA[i][j];
//				expA[i][j] += 1.0 / 6.0 * AAA[i][j];
//
//				for (int k = 0; k < A.size(); k++)
//
//					expA[i][j] += 1.0 / 24.0 * AA[i][k] * AA[k][j];
//
//			}
//
//		}
//
//#ifdef _commonDebugMode_
//		logFile<<"******************************************************"<<endl;
//		logFile<<"*********************** exp(A) ***********************"<<endl;
//		for(int r=0;r<expA.size();r++)
//		for(int s=0;s<expA[r].size();s++)
//		logFile<<"exp(A)["<<r<<"]["<<s<<"] = "
//		<<expA[r][s]<<endl;
//#endif
//
//		// -------------------------------------------------------------------
//
//		if (derivOrder > 0) {
//
//			if (dexpA.size() < fullSize)
//				allocateArray(dexpA, fullSize, fullSize);
//
//			clearArray(dexpA);
//
//			// dexp(A)_{ij} / dA_{rs} = \d_{ir}\,\d_{js}
//			//
//			// + 1/2 (d_{ir} \A_{sj} + d_{js} A_{ir})
//			//
//			// + 1/6 (d_ir A_{sm} A_{mj} + A_{ir} A_{sj}
//			//        + d_{js} A_{im} A_{mr})
//			// + 1/24 ( d_{ir} A_{sp} A_{pq} A_{qj}
//			//        + A_{ir} A_{sp} A_{pj}
//			//        + A_{ip} A_{pr} A_{sj}
//			//        + d_{js} A_{ip} A_{pq} A_{qr} ) + ...
//
//			for (int i = 0; i < A.size(); i++)
//
//				for (int j = 0; j < A.size(); j++)
//
//					for (int r = 0; r < A.size(); r++)
//
//						for (int s = 0; s < A.size(); s++)
//
//							dexpA[vFull[i][j]][vFull[r][s]] = delta[i][r]
//									* delta[j][s]
//
//									+ 1.0 / 2.0
//											* (delta[i][r] * A[s][j]
//													+ delta[j][s] * A[i][r])
//
//									+ 1.0 / 6.0
//											* (delta[i][r] * AA[s][j]
//													+ A[i][r] * A[s][j]
//													+ delta[j][s] * AA[i][r])
//
//									+ 1.0 / 24.0
//											* (delta[i][r] * AAA[s][j]
//													+ A[i][r] * AA[s][j]
//													+ AA[i][r] * A[s][j]
//													+ delta[j][s] * AAA[i][r]);
//
//		}
//
//#ifdef _commonDebugMode_
//		logFile<<"******************************************************"<<endl;
//		logFile<<"*********************** dexp(A) **********************"<<endl;
//		for(int r=0;r<dexpA.size();r++)
//		for(int s=0;s<dexpA[r].size();s++)
//		logFile<<"d(exp(A))/dA["<<r<<"]["<<s<<"] = "
//		<<dexpA[r][s]<<endl;
//#endif
//
//	}
//
//	/**********************************************************************/
//	// compute expA and dexp(A)/dA ( Carlo's dexpo FORTRAN routine is used )
//	else {
//
//		int fullSize = A.size() * A.size();
//
//		if (expA.size() < A.size())
//			allocateArray(expA, A.size(), A.size());
//
//		clearArray(expA);
//
//		if (dexpA.size() < fullSize)
//			allocateArray(dexpA, fullSize, fullSize);
//
//		clearArray(dexpA);
//
//		// fortran arrays need to be vectors
//
//		dbVector fortranStrain, fortranExp, fortranDexp;
//
//		if (fortranStrain.size() == 0)
//			fortranStrain.resize(fullSize);
//
//		for (int r = 0; r < A.size(); r++)
//
//			for (int s = 0; s < A[r].size(); s++)
//
//				fortranStrain[r * A[r].size() + s] = A[r][s];
//
//#ifdef _commonDebugMode_
//		logFile<<"######################################################"<<endl;
//		logFile<<"******************* fortran strain *******************"<<endl;
//		for(int r=0;r<fortranStrain.size();r++)
//		logFile<<"A["<<r<<"] = "<<fortranStrain[r]<<endl;
//#endif
//
//		/**********************************************************************/
//		// compute exp(A) and dexp(A)/dA
//		if (fortranExp.size() == 0)
//			fortranExp.resize(fullSize);
//
//		if (fortranDexp.size() == 0)
//			fortranDexp.resize((int) pow(A.size(), 4.0));
//
//		int ierr = 0;
//		double tol = 1.0e-12;
//
//		dexpo_(&fortranStrain[0], &fortranExp[0], &fortranDexp[0], tol, rank,
//				ierr);
//
//		if (ierr != 0) {
//			logFile << "In commonFunctions::calcExponential failed!" << endl;
//			MPI_Abort(MPI_COMM_WORLD, 1);
//		}
//
//		// store expA
//
//		for (int k = 0; k < usedDims; k++)
//
//			for (int l = 0; l < usedDims; l++)
//
//				expA[k][l] = fortranExp[k * usedDims + l];
//
//		// store dexpA
//
//		for (int k = 0; k < fullSize; k++)
//
//			for (int l = 0; l < fullSize; l++)
//
//				dexpA[k][l] = fortranDexp[k * fullSize + l];
//
//#ifdef _commonDebugMode_
//		logFile<<"******************************************************"<<endl;
//		logFile<<"*********************** exp(A) ***********************"<<endl;
//		for(int r=0;r<expA.size();r++)
//		for(int s=0;s<expA[r].size();s++)
//		logFile<<"exp(A)["<<r<<"]["<<s<<"] = "
//		<<expA[r][s]<<endl;
//		logFile<<"*********************** dexp(A) **********************"<<endl;
//		for(int r=0;r<dexpA.size();r++)
//		for(int s=0;s<dexpA[r].size();s++)
//		logFile<<"d(exp(A))/dA["<<r<<"]["<<s<<"] = "
//		<<dexpA[r][s]<<endl;
//#endif
//
//	}
//
//	/**********************************************************************/
//
//	if (derivOrder > 1) {
//
//		allocateArray(d2expA, fullSize, fullSize, fullSize);
//		clearArray(d2expA);
//
//		// dexp(D)_{ij}/ dD_{pq}dD_{rs}) =
//		//
//		// 1.0/2.0*(delta_{ir}delta_{sp}delta_{jq} + delta_{js}delta_{ip}delta_{rq})
//		//
//		//+ 1.0/6.0*(delta_{ir}delta_{sp}D_{qj} + delta_{ir}delta_{jq}D_{sp}
//		//
//		//     + delta_{ip}delta_{rq}D_{sj} + delta_{sp}delta_{jq}D_{ir}
//		//
//		//     + delta_{js}delta_{ip}D_{qr} + delta_{js}\delta_{rq}D_{ip})
//		//
//		//+ 1.0/24.0*(delta_{ir}delta_{sp}D^2_{qj}
//		//     + delta_{ir}D_{sp}D_{qj} + delta_{ir}delta_{jq}D^2_{sp}
//		//     + delta_{ip}delta_{rq}D^2_{sj} + delta_{sp}D_{ir}D_{qj}
//		//     + delta_{jq}D_{ir}D_{sp}
//		//     + delta_{ip}D_{qr}D_{sj} + delta_{rq}D_{ip}D_{sj}
//		//     + delta_{sp}delta_{jq}D^2_{ir}
//		//     + delta_{js}delta_{ip}D^2_{qr} + delta_{js}D_{ip}D_{qr}
//		//     + delta_{js}delta_{rq}D^2_{ip}) + ...
//
//		for (int i = 0; i < A.size(); i++)
//
//			for (int j = 0; j < A.size(); j++)
//
//				for (int r = 0; r < A.size(); r++)
//
//					for (int s = 0; s < A.size(); s++)
//
//						for (int p = 0; p < A.size(); p++)
//
//							for (int q = 0; q < A.size(); q++)
//
//								d2expA[vFull[i][j]][vFull[r][s]][vFull[p][q]] =
//
//								1.0 / 2.0
//										* (delta[i][r] * delta[s][p]
//												* delta[j][q]
//												+ delta[j][s] * delta[i][p]
//														* delta[r][q])
//
//										+ 1.0 / 6.0
//												* (delta[i][r] * delta[s][p]
//														* A[q][j]
//														+ delta[i][r]
//																* delta[j][q]
//																* A[s][p]
//
//														+ delta[i][p]
//																* delta[r][q]
//																* A[s][j]
//														+ delta[s][p]
//																* delta[j][q]
//																* A[i][r]
//
//														+ delta[j][s]
//																* delta[i][p]
//																* A[q][r]
//														+ delta[j][s]
//																* delta[r][q]
//																* A[i][p])
//
//										+ 1.0 / 24.0
//												* (delta[i][r] * delta[s][p]
//														* AA[q][j]
//														+ delta[i][r] * A[s][p]
//																* A[q][j]
//														+ delta[i][r]
//																* delta[j][q]
//																* AA[s][p]
//														+ delta[i][p]
//																* delta[r][q]
//																* AA[s][j]
//														+ delta[s][p] * A[i][r]
//																* A[q][j]
//														+ delta[j][q] * A[i][r]
//																* A[s][p]
//														+ delta[i][p] * A[q][r]
//																* A[s][j]
//														+ delta[r][q] * A[i][p]
//																* A[s][j]
//														+ delta[s][p]
//																* delta[j][q]
//																* AA[i][r]
//														+ delta[j][s]
//																* delta[i][p]
//																* AA[q][r]
//														+ delta[j][s] * A[i][p]
//																* A[q][r]
//														+ delta[j][s]
//																* delta[r][q]
//																* AA[i][p]);
//
//#ifdef _constitutiveDebugMode_
//		logFile<<"******************************************************"<<endl;
//		logFile<<"************************* d2exp(A) *******************"<<endl;
//		for(int r=0;r<d2expA.size();r++)
//		for(int s=0;s<d2expA[r].size();s++)
//		for(int t=0;t<d2expA[r][s].size();t++)
//		logFile<<"d2(exp(A))/dAdA["<<r<<"]["<<s<<"]["<<t<<"] = "
//		<<d2expA[r][s][t]<<endl;
//#endif
//
//	}
//
//}

/***********************************************************************/
/***********************************************************************/
// Vector and matrix allocation routines
//void allocateArray(intVector& vec,int idx) {
//
//  using namespace std;
//
//  vec.resize(idx);
//
//}
//
//void allocateArray(dbVector& vec,int idx) {
//
//  using namespace std;
//
//  vec.resize(idx);
//
//}
//void allocateArray(blMatrix& mat,int idx1) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//}
//
//void allocateArray(dbMatrix& mat,int idx1) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//}
//void allocateArray(blMatrix& mat,int idx1,int idx2) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//  for(int i=0;i<idx1;i++)
//    mat[i].resize(idx2);
//
//}
//void allocateArray(intMatrix& mat,int idx1,int idx2) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//  for(int i=0;i<idx1;i++)
//    mat[i].resize(idx2);
//
//}
//void allocateArray(intMatrix3& mat,int idx1,int idx2,int idx3) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//  for(int i=0;i<idx1;i++) {
//    mat[i].resize(idx2);
//
//    for(int j=0;j<idx2;j++) {
//      mat[i][j].resize(idx3);
//
//    }
//
//  }
//
//}
void allocateArray(intMatrix3& outMat, intMatrix3& inMat) {

	using namespace std;

	outMat.resize(inMat.size());

	for (int i = 0; i < inMat.size(); i++) {
		outMat[i].resize(inMat[i].size());

		for (int j = 0; j < inMat[i].size(); j++)
			outMat[i][j].resize(inMat[i][j].size());

	}

}

//void allocateArray(dbMatrix& mat,int idx1,int idx2) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//  for(int i=0;i<idx1;i++)
//    mat[i].resize(idx2);
//
//}

//void allocateArray(dbMatrix3& mat,int idx) {
//
//  using namespace std;
//
//  mat.resize(idx);
//}

//void allocateArray(dbMatrix3& mat,int idx1,int idx2,int idx3) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//  for(int i=0;i<idx1;i++) {
//    mat[i].resize(idx2);
//
//    for(int j=0;j<idx2;j++) {
//      mat[i][j].resize(idx3);
//
//    }
//
//  }
//
//}

//void allocateArray(dbMatrix4& mat,int idx) {
//
//  using namespace std;
//
//  mat.resize(idx);
//}

void allocateArray(dbMatrix4& mat, int idx1, int idx2, int idx3, int idx4) {

	using namespace std;

	mat.resize(idx1);

	for (int i = 0; i < idx1; i++) {
		mat[i].resize(idx2);

		for (int j = 0; j < idx2; j++) {
			mat[i][j].resize(idx3);

			for (int k = 0; k < idx3; k++) {
				mat[i][j][k].resize(idx4);

			}

		}

	}

}

void allocateArray(dbMatrix5& mat, int idx1, int idx2, int idx3, int idx4,
		int idx5) {

	using namespace std;

	mat.resize(idx1);

	for (int i = 0; i < idx1; i++) {
		mat[i].resize(idx2);

		for (int j = 0; j < idx2; j++) {
			mat[i][j].resize(idx3);

			for (int k = 0; k < idx3; k++) {
				mat[i][j][k].resize(idx4);

				for (int l = 0; l < idx4; l++) {
					mat[i][j][k][l].resize(idx5);

				}

			}

		}

	}

}

//void allocateArray(ldbMatrix& mat,int idx1,int idx2) {
//
//  using namespace std;
//
//  mat.resize(idx1);
//
//  for(int i=0;i<idx1;i++)
//    mat[i].resize(idx2);
//
//}

//void pushBackVector(intVector& vec, int n) {
//
//	using namespace std;
//
//	vec.resize(vec.size() + 1);
//	vec[vec.size() - 1] = n;
//}

//void pushBackVector(dbVector& vec, double n) {
//
//	using namespace std;
//
//	vec.resize(vec.size() + 1);
//	vec[vec.size() - 1] = n;
//}

//void pushBackVector(intMatrix& mat, intVector& vec) {
//
//	using namespace std;
//
//	mat.resize(mat.size() + 1);
//	mat[mat.size() - 1] = vec;
//}

//void pushBackVector(dbMatrix& mat, dbVector& vec) {
//
//	using namespace std;
//
//	mat.resize(mat.size() + 1);
//	mat[mat.size() - 1] = vec;
//}

void pushBackVector(intMatrix& oldMat, intMatrix& deltaMat) {

	using namespace std;

	int oldSize = oldMat.size();

	oldMat.resize(oldSize + deltaMat.size());

	for (int i = oldSize, j = 0; i < oldSize + deltaMat.size(); i++, j++)
		oldMat[i] = deltaMat[j];
}

/************************************************************************/
// Vector resize routines which shrink and fit also its reserved capacity
//void resizeArray(intVector& vec, int size) {
//
//	vec.resize(size);
//	intVector(vec).swap(vec);
//}

//void resizeArray(dbVector& vec, int size) {
//
//	vec.resize(size);
//	dbVector(vec).swap(vec);
//}

//void resizeArray(intMatrix& mat, int size) {
//
//	mat.resize(size);
//	intMatrix(mat).swap(mat);
//}

//void resizeArray(dbMatrix& mat, int size) {
//
//	mat.resize(size);
//	dbMatrix(mat).swap(mat);
//
//}

void resizeArray(dbMatrix& mat, int size1, int size2) {

	mat.resize(size1);

	for (int i = 0; i < mat.size(); i++)
		mat[i].resize(size2);

	dbMatrix(mat).swap(mat);

}

void resizeArray(dbMatrix3& mat, int size) {

	mat.resize(size);
	dbMatrix3(mat).swap(mat);

}

void resizeArray(dbMatrix4& mat, int size) {

	mat.resize(size);
	dbMatrix4(mat).swap(mat);

}

/***********************************************************************/
// convert a integer number to a string
std::string convertIntToString(int number) {

	using namespace std;

	int digit, digit1, digit2;
	double floatPart, intPart;

	string convertedNumber;

	double value = (double) number;

	bool flag = true;

	while (flag) {

		value /= 10.0;
		floatPart = modf(value, &intPart);
		value = floatPart * 10.0;

		digit1 = (int) ceil(value);
		digit2 = (int) floor(value);

		if (fabs(digit1 - value) > 0.1)

			digit = digit2 + 48;

		else

			digit = digit1 + 48;

		convertedNumber.insert(convertedNumber.begin(), (char) digit);

		if (intPart == 0)
			break;

		value = intPart;
	}

	return convertedNumber;
}

/***********************************************************************/
// Check if a double precision number equals to zero within numerical
// precision limits.
bool checkZero(double& value) {

	if (fabs(value) < DBL_EPSILON * 1.0e+03)
		return true;

	else
		return false;

}

/***********************************************************************/
// Determine the twist of a point given by its reference configuration
// and its current configuration around a given axis going through the
// origin.
double calcPlacementTwist(dbVector& reference, dbVector& current,
		dbVector& rotAxis, std::ofstream& logFile) {

	using namespace std;

	if (reference.size() != 3 || current.size() != 3 || rotAxis.size() != 3) {
		logFile << "In commonFunctions::calcPlacementTwist input vectors\n"
				<< "are not admissible!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	double axisAbs = 0;

	dbVector xref = reference;
	dbVector xcur = current;
	dbVector axis = rotAxis;

	for (int j = 0; j < axis.size(); j++)
		axisAbs += axis[j] * axis[j];

	axisAbs = sqrt(axisAbs);

	// normalize the rotation axis
	if (axisAbs != 1.0)

		for (int j = 0; j < axis.size(); j++)
			axis[j] /= axisAbs;

#ifdef _commonDebugMode_
	logFile<<"******************* calculate twist ******************"<<endl;
	logFile<<"axis: "<<axis[0]<<" "<<axis[1]<<" "<<axis[2]<<endl;
	logFile<<"xref: "<<xref[0]<<" "<<xref[1]<<" "<<xref[2]<<endl;
	logFile<<"xcur: "<<xcur[0]<<" "<<xcur[1]<<" "<<xcur[2]<<endl;
#endif

	// calculate vectors projected onto the rotation plane

	dbVector xrefp(xref.size());
	dbVector xcurp(xcur.size());
	intMatrix e = getPermutations(xref.size());

	// a'_m = n x a x n = e_ijk e_klm n_i a_j n_l

	// a'_m = n_i a_m n_i - n_m a_i n_i

	for (int m = 0; m < xrefp.size(); m++)

		for (int i = 0; i < xcurp.size(); i++) {

			xrefp[m] += axis[i] * xref[m] * axis[i]
					- axis[m] * xref[i] * axis[i];
			xcurp[m] += axis[i] * xcur[m] * axis[i]
					- axis[m] * xcur[i] * axis[i];

		}

#ifdef _commonDebugMode_
	logFile<<"------------------------------------------------------"<<endl;
	logFile<<"xrefp: "<<xrefp[0]<<" "<<xrefp[1]<<" "<<xrefp[2]<<endl;
	logFile<<"xcurp: "<<xcurp[0]<<" "<<xcurp[1]<<" "<<xcurp[2]<<endl;
#endif

	// xrefp * xcurp = |xrefp||xcurp| cos(twist)

	double xcurpAbs = 0;
	double xrefpAbs = 0;
	double twist = 0;

	for (int i = 0; i < xcurp.size(); i++) {

		xcurpAbs += xcurp[i] * xcurp[i];
		xrefpAbs += xrefp[i] * xrefp[i];
		twist += xrefp[i] * xcurp[i];

	}

	xcurpAbs = sqrt(xcurpAbs);
	xrefpAbs = sqrt(xrefpAbs);

	twist = twist / (xcurpAbs * xrefpAbs);

#ifdef _commonDebugMode_
	logFile<<"------------------------------------------------------"<<endl;
	logFile<<"cos(twist)="<<twist<<endl;
#endif

	twist = acos(twist);

	return twist;

}

/***********************************************************************/
// Might be buggy ?!
// Determine the twist of a point given by its origin 'reference' and its
// image 'current' around a given axis going through a certain point
// (fix point).
double calcPlacementTwist(dbVector& reference, dbVector& current,
		dbVector& rotAxis, dbVector& fixPoint, std::ofstream& logFile) {

	using namespace std;

	if (reference.size() != 3 || current.size() != 3 || rotAxis.size() != 3
			|| fixPoint.size() != 3) {
		logFile << "In commonFunctions::calcPlacementTwist input vectors\n"
				<< "are not admissible!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	double axisAbs = 0;

	dbVector xref = reference;
	dbVector xcur = current;
	dbVector axis = rotAxis;

	for (int j = 0; j < xref.size(); j++) {
		xref[j] -= fixPoint[j];
		xcur[j] -= fixPoint[j];
		axisAbs += axis[j] * axis[j];
	}

	axisAbs = sqrt(axisAbs);

	// normalize the rotation axis
	if (axisAbs != 1.0)

		for (int j = 0; j < axis.size(); j++)
			axis[j] /= axisAbs;

	// calculate vectors projected onto the rotation plane

	double scalar = 0;
	for (int j = 0; j < xref.size(); j++)
		scalar += xref[j] * axis[j];

	for (int j = 0; j < xref.size(); j++)
		xref[j] = xref[j] - scalar * axis[j];

	scalar = 0;
	for (int j = 0; j < xcur.size(); j++)
		scalar += xcur[j] * axis[j];

	for (int j = 0; j < xcur.size(); j++)
		xcur[j] = xcur[j] - scalar * axis[j];

	scalar = 0;
	double scalar1 = 0;
	double scalar2 = 0;
	for (int j = 0; j < xcur.size(); j++) {

		scalar += xref[j] * xcur[j];
		scalar1 += xref[j] * xref[j];
		scalar2 += xcur[j] * xcur[j];

	}

	scalar1 = sqrt(scalar1);
	scalar2 = sqrt(scalar2);

	double twist = 0;

	if (fabs(scalar1 * scalar2) > 0)
		twist = acos(scalar / (scalar1 * scalar2));

	return twist;

}

/***********************************************************************/
// Project a vector onto a plane defined by its normal vector.
void projectVecToPlane(dbVector& a, dbVector& normal, dbVector& b,
		std::ofstream& logFile) {

	using namespace std;

	if (a.size() != normal.size()) {
		logFile << "In commonFunctions::projectVecToPlane input vectors\n"
				<< "are not admissible!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int size = a.size();

	if (b.size() != size)
		b.resize(size);

	clearArray(b);

	intMatrix e = getPermutations(size);

	// b_m = n x a x n = e_ijk e_klm n_i a_j n_l =

	//      = n_i a_m n_i - n_m a_i n_i

	for (int m = 0; m < size; m++)

		for (int i = 0; i < size; i++)

			b[m] += normal[i] * a[m] * normal[i] - normal[m] * a[i] * normal[i];

#ifdef _commonDebugMode_
	logFile<<"*************** project vector to plane **************"<<endl;
	logFile<<"normal: ";
	for(int i=0;i<size;i++)
	logFile<<normal[i]<<" ";
	logFile<<endl;
	logFile<<"a: ";
	for(int i=0;i<size;i++)
	logFile<<a[i]<<" ";
	logFile<<endl;
	logFile<<"b: ";
	for(int i=0;i<size;i++)
	logFile<<b[i]<<" ";
	logFile<<endl;
#endif

}

/***********************************************************************/
/***********************************************************************/
// Calculate an orthogonal basis to a given vector v
void calcOrthogonalBasis(dbVector& vec, dbMatrix& basis,
		std::ofstream& logFile) {

	using namespace std;

	if (vec.size() != 3) {
		logFile << "In commonFunctions::calcOrthogonalBasis invalid input\n"
				<< " vector!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	allocateArray(basis, 3, 3);
	clearArray(basis);

	// normalized v
	double absVec;
	scalarProduct(vec, vec, absVec, logFile);
	absVec = sqrt(absVec);

	if (absVec < DBL_EPSILON) {
		logFile << "In commonFunctions::calcOrthogonalBasis given vector is 0"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	for (int i = 0; i < vec.size(); i++)
		vec[i] /= absVec;

	basis[0] = vec;
	dbVector& a = basis[1];
	dbVector& b = basis[2];

	// ---------------------------------------------------------------------
	// compute second basis vector a
	//
	// v1*a1 + v2*a2 + v3*a3 = 0 (from the vector dot product with
	//                            v \perp a)
	// a3 = (-v1 - v2)/v3 (choose a1 = a2 = 1)

	if (vec[0] != 0) {

		a[1] = 1.0;
		a[2] = 1.0;
		a[0] = (-vec[1] - vec[2]) / vec[0];
	} else if (vec[1] != 0) {

		a[0] = 1.0;
		a[2] = 1.0;
		a[1] = (-vec[0] - vec[2]) / vec[1];
	} else if (vec[2] != 0) {

		a[0] = 1.0;
		a[1] = 1.0;
		a[2] = (-vec[0] - vec[1]) / vec[2];
	}

	// normalized a
	scalarProduct(a, a, absVec, logFile);
	absVec = sqrt(absVec);

	for (int i = 0; i < a.size(); i++)
		a[i] /= absVec;

	// ---------------------------------------------------------------------
	// compute third basis vector a

	crossProduct(vec, a, b);

}

/***********************************************************************/
/***********************************************************************/
// Calculate rotation tensor.
void calcRotationMatrix(dbVector& rotVec, dbMatrix& rTens,
		std::ofstream& logFile) {

	using namespace std;

	if (rTens.size() < 3) {
		rTens.resize(3);

		for (int i = 0; i < 3; i++)
			rTens[i].resize(3);

	}

	clearArray(rTens);

	//---------------------------------------------------------------------

	ldbVector rots(3);

	for (int k = 0; k < rots.size(); k++)
		rots[k] = (long double) rotVec[k];

	long double absRotVec = (long double) sqrt(
			rots[0] * rots[0] + rots[1] * rots[1] + rots[2] * rots[2]);

	long double cosinusTerm, sinusTerm;

	if (absRotVec > DBL_EPSILON) {

		sinusTerm = sin(absRotVec) / absRotVec;
		cosinusTerm = ((long double) 1.0 - cos(absRotVec))
				/ (absRotVec * absRotVec);
	} else {
		sinusTerm = 1.0;
		cosinusTerm = 0.5;
	}

#ifdef  _kinematicsDebugMode_
	logFile<<"****************************************************"<<endl;
	logFile<<"calculating rotation tensor"<<endl;
	logFile<<"rot[0] = "<<rots[0]<<"rot[1] = "<<rots[1]
	<<"rot[2] = "<<rots[2]<<endl;
	logFile<<"absRotVec "<<absRotVec<<endl;
	logFile<<"sinusTerm "<<sinusTerm<<endl;
	logFile<<"cosinusTerm "<<cosinusTerm<<endl;
#endif

	dbMatrix delta = getKroneckerSymbol(3);

	rTens[0][0] = delta[0][0]
			- cosinusTerm * (rots[1] * rots[1] + rots[2] * rots[2]);

	rTens[0][1] = -sinusTerm * rots[2] + cosinusTerm * rots[0] * rots[1];

	rTens[0][2] = sinusTerm * rots[1] + cosinusTerm * rots[0] * rots[2];

	rTens[1][0] = sinusTerm * rots[2] + cosinusTerm * rots[1] * rots[0];

	rTens[1][1] = delta[1][1]
			- cosinusTerm * (rots[0] * rots[0] + rots[2] * rots[2]);

	rTens[1][2] = -sinusTerm * rots[0] + cosinusTerm * rots[1] * rots[2];

	rTens[2][0] = -sinusTerm * rots[1] + cosinusTerm * rots[2] * rots[0];

	rTens[2][1] = sinusTerm * rots[0] + cosinusTerm * rots[2] * rots[1];

	rTens[2][2] = delta[2][2]
			- cosinusTerm * (rots[0] * rots[0] + rots[1] * rots[1]);

}

/************************************************************************/
/************************************************************************/
// Determine the rotation matrix corresponding the rotation of the
// circumferential tangent vector around the cylinder's axis
// so that the tangent goes through the given position
//
// axis    ... cylinder axis (x1,x2,x3 Cartesian axis)
// X       ... position vector
// R       ... rotation matrix (output)
// dR      ... spatial derivatives of rotation matrix (output)
//
void setCylinderMetrics(dbVector& axis, dbVector& X, dbMatrix& R, dbMatrix3& dR,
		std::ofstream& logFile) {

	using namespace std;

	double theta, rho;

	double pi = getPi();
	int usedDims = axis.size();

	if (usedDims < 3) {

		cerr << "In commonFunctions::setRotationMatrix only two or three\n"
				<< "dimensional space is allowed!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (X.size() < usedDims) {

		cerr << "In commonFunctions::setRotationMatrix not enough entries in\n"
				<< "position vector!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// project position vector onto the rotation plane
	dbVector XP;
	projectVecToPlane(X, axis, XP, logFile);

	dbMatrix delta = getKroneckerSymbol(usedDims);
	R = delta;

	if (dR.size() < 3)
		allocateArray(dR, 3, 3, 3);

	clearArray(dR);

	// ---------------------------------------------------------------------
	// cylinder axis = X1 => tangent = X3
	if (axis[0] != 0 && axis[1] == 0 && axis[2] == 0) {

		// radius
		rho = sqrt(pow(XP[1], 2) + pow(XP[2], 2));

		// angle
		if (XP[2] >= 0)
			theta = acos(XP[1] / rho);

		else
			theta = 2.0 * pi - acos(XP[1] / rho);

		if (!checkZero(rho)) {

			R[1][1] = cos(theta);
			R[1][2] = sin(theta);
			R[2][1] = -sin(theta);
			R[2][2] = cos(theta);

			// ---------------------
			// spatial derivatives with respect to x2 and x3

			// d(cos(theta))/dx_i =  d(cos(theta))/dtheta dtheta/dx_i = -sin(theta) dtheta/dy dy/dx_i

			// d(sin(theta))/dx_i = d(sin(theta))/dtheta dtheta/dx_i = cos(theta) dtheta/dy dtheta/dx_i

			// with y = x2/(x2^2 = x3^2)^0.5

			double DcosDth = -sin(theta);
			double DsinDth = cos(theta);

			double y = XP[1] / rho;

			// dtheta/dy = -1/(1-y^2)^0.5 x3 >=  0
			// dtheta/dy = +1/(1-y^2)^0.5 x3 <  0
			double DthDy;

			if (XP[2] >= 0)
				DthDy = -1.0 / pow(1.0 - y * y, 0.5);

			else
				DthDy = 1.0 / pow(1.0 - y * y, 0.5);

			// dy/dx2 = (rho-x2^2 rho^{-1}) / rho^2
			//
			// dy/dx3 = -x2 x3 / rho^3

			double DyDx2 = (rho - pow(XP[1], 2.0) * 1.0 / rho) / pow(rho, 2.0);
			double DyDx3 = -XP[1] * XP[2] / pow(rho, 3.0);

			dR[1][1][1] = DcosDth * DthDy * DyDx2;
			dR[1][1][2] = DsinDth * DthDy * DyDx2;
			dR[1][2][1] = -DsinDth * DthDy * DyDx2;
			dR[1][2][2] = DcosDth * DthDy * DyDx2;

			dR[2][1][1] = DcosDth * DthDy * DyDx3;
			dR[2][1][2] = DsinDth * DthDy * DyDx3;
			dR[2][2][1] = -DsinDth * DthDy * DyDx3;
			dR[2][2][2] = DcosDth * DthDy * DyDx3;

		}

#ifdef _commonDebugMode_
		logFile<<"rho="<<rho<<endl;
		if(XP[2] >= 0)
		theta = acos(XP[1]/rho);
		else
		theta = 2.0*pi - acos(XP[1]/rho);
		logFile<<"theta="<<theta<<endl;
		logFile<<"R[1][1]="<<R[1][1]<<" ?= "<<cos(theta)<<endl;
		logFile<<"R[1][2]="<<R[1][2]<<" ?= "<<sin(theta)<<endl;
		logFile<<"R[2][1]="<<R[2][1]<<" ?= "<<-sin(theta)<<endl;
		logFile<<"R[2][2]="<<R[2][2]<<" ?= "<<cos(theta)<<endl;
#endif

	}

	// ---------------------------------------------------------------------
	// cylinder axis = X2 => tangent = X1
	else if (axis[0] == 0 && axis[1] != 0 && axis[2] == 0) {

		// radius
		rho = sqrt(pow(XP[2], 2) + pow(XP[0], 0));

		// angle
		if (XP[0] >= 0)
			theta = acos(XP[2] / rho);

		else
			theta = 2.0 * pi - acos(XP[2] / rho);

		if (!checkZero(rho)) {

			R[0][0] = cos(theta);
			R[0][2] = -sin(theta);
			R[2][0] = sin(theta);
			R[2][2] = cos(theta);

			// ---------------------
			// spatial derivatives with respect to x2 and x3

			// d(cos(theta))/dx_i =  d(cos(theta))/dtheta dtheta/dx_i = -sin(theta) dtheta/dy dy/dx_i

			// d(sin(theta))/dx_i = d(sin(theta))/dtheta dtheta/dx_i = cos(theta) dtheta/dy dtheta/dx_i

			// with y = x3/(x3^2 = x0^2)^0.5

			double DcosDth = -sin(theta);
			double DsinDth = cos(theta);

			double y = XP[2] / rho;

			// dtheta/dy = -1/(1-y^2)^0.5; x1 >=  0
			// dtheta/dy = +1/(1-y^2)^0.5; x1 <  0
			double DthDy;

			if (XP[0] >= 0)
				DthDy = -1.0 / pow(1.0 - y * y, 0.5);

			else
				DthDy = 1.0 / pow(1.0 - y * y, 0.5);

			// dy/dx3 = (rho-x3^2 rho^{-1}) / rho^2
			//
			// dy/dx1 = -x3 x1 / rho^3

			double DyDx3 = (rho - pow(XP[2], 2.0) * 1.0 / rho) / pow(rho, 2.0);
			double DyDx1 = -XP[2] * XP[0] / pow(rho, 3.0);

			dR[2][0][0] = DcosDth * DthDy * DyDx3;
			dR[2][0][2] = -DsinDth * DthDy * DyDx3;
			dR[2][2][0] = DsinDth * DthDy * DyDx3;
			dR[2][2][2] = DcosDth * DthDy * DyDx3;

			dR[0][0][0] = DcosDth * DthDy * DyDx1;
			dR[0][0][2] = -DsinDth * DthDy * DyDx1;
			dR[0][2][0] = DsinDth * DthDy * DyDx1;
			dR[0][2][2] = DcosDth * DthDy * DyDx1;

		}

#ifdef _commonDebugMode_
		logFile<<"rho="<<rho<<endl;
		if(XP[0] >= 0)
		theta = acos(XP[2]/rho);
		else
		theta = 2.0*pi - acos(XP[2]/rho);
		logFile<<"theta="<<theta<<endl;
		logFile<<"R[0][0]="<<R[0][0]<<" ?= "<<cos(theta)<<endl;
		logFile<<"R[0][2]="<<R[0][2]<<" ?= "<<-sin(theta)<<endl;
		logFile<<"R[2][0]="<<R[2][0]<<" ?= "<<sin(theta)<<endl;
		logFile<<"R[2][2]="<<R[2][2]<<" ?= "<<cos(theta)<<endl;
#endif

		logFile
				<< "In commonFunctions::setCylinderMetrics chosen rotation direction "
				<< "is not tested yet!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// ---------------------------------------------------------------------
	// cylinder axis = X3 => tangent = X2
	else if (axis[0] == 0 && axis[1] == 0 && axis[2] != 0) {

		// radius
		rho = sqrt(pow(XP[0], 2) + pow(XP[1], 2));

		// angle
		if (XP[1] >= 0)
			theta = acos(XP[0] / rho);

		else
			theta = 2.0 * pi - acos(XP[0] / rho);

		if (!checkZero(rho)) {

			R[0][0] = cos(theta);
			R[0][1] = sin(theta);
			R[1][0] = -sin(theta);
			R[1][1] = cos(theta);

			// ---------------------
			// spatial derivatives with respect to x2 and x3

			// d(cos(theta))/dx_i =  d(cos(theta))/dtheta dtheta/dx_i = -sin(theta) dtheta/dy dy/dx_i

			// d(sin(theta))/dx_i = d(sin(theta))/dtheta dtheta/dx_i = cos(theta) dtheta/dy dtheta/dx_i

			// with y = x1/(x1^2 = x2^2)^0.5

			double DcosDth = -sin(theta);
			double DsinDth = cos(theta);

			double y = XP[0] / rho;

			// dtheta/dy = -1/(1-y^2)^0.5 x2 >=  0
			// dtheta/dy = +1/(1-y^2)^0.5 x2 <  0
			double DthDy;

			if (XP[1] >= 0)
				DthDy = -1.0 / pow(1.0 - y * y, 0.5);

			else
				DthDy = 1.0 / pow(1.0 - y * y, 0.5);

			// dy/dx1 = (rho-x1^2 rho^{-1}) / rho^2
			//
			// dy/dx2 = -x1 x2 / rho^3

			double DyDx1 = (rho - pow(XP[0], 2.0) * 1.0 / rho) / pow(rho, 2.0);
			double DyDx2 = -XP[0] * XP[1] / pow(rho, 3.0);

			dR[0][0][0] = DcosDth * DthDy * DyDx1;
			dR[0][0][1] = DsinDth * DthDy * DyDx1;
			dR[0][1][0] = -DsinDth * DthDy * DyDx1;
			dR[0][1][1] = DcosDth * DthDy * DyDx1;

			dR[1][0][0] = DcosDth * DthDy * DyDx2;
			dR[1][0][1] = DsinDth * DthDy * DyDx2;
			dR[1][1][0] = -DsinDth * DthDy * DyDx2;
			dR[1][1][1] = DcosDth * DthDy * DyDx2;

		}

#ifdef _commonDebugMode_
		logFile<<"rho="<<rho<<endl;
		if(XP[1] >= 0)
		theta = acos(XP[0]/rho);
		else
		theta = 2.0*pi - acos(XP[0]/rho);
		logFile<<"theta="<<theta<<endl;
		logFile<<"R[0][0]="<<R[0][0]<<" ?= "<<cos(theta)<<endl;
		logFile<<"R[0][1]="<<R[0][1]<<" ?= "<<sin(theta)<<endl;
		logFile<<"R[1][0]="<<R[1][0]<<" ?= "<<-sin(theta)<<endl;
		logFile<<"R[1][1]="<<R[1][1]<<" ?= "<<cos(theta)<<endl;
#endif

		logFile
				<< "In commonFunctions::setCylinderMetrics chosen rotation direction "
				<< "is not tested yet!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	else {

		cerr
				<< "In commonFunctions::setCylinderMetrics chosen rotation direction "
				<< "is not supported!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/***********************************************************************/
/***********************************************************************/
// compute the distance between two points
double calcDistance(dbVector& P, dbVector& Q) {

	using namespace std;

	if (P.size() != Q.size()) {
		cerr << "In commonFunctions::calcDistance coordinate dimensions don't\n"
				<< "match!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	double distance = 0;

	for (int i = 0; i < P.size(); i++)

		distance += pow(P[i] - Q[i], 2);

	return (sqrt(distance));

}

// populate a symmetric matrix that has upper triangle filled.
void populateSymmetricMatrix(dbMatrix& sym, int size) {
	for (int i = 1; i < size; i++)
		for (int j = 0; j < i; j++)
			sym[i][j] = sym[j][i];
}
/***********************************************************************/
/***********************************************************************/
// Determine an equal distribution of an array over the processors.
void getParallelArrayDistribution(int globalNum, int& startIdx, int& endIdx,
		std::ofstream& logFile) {

	using namespace std;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int defaultPortion = (int) ceil((double) globalNum / size);

	if (defaultPortion * rank < globalNum
			&& defaultPortion * (rank + 1) <= globalNum) {

		startIdx = defaultPortion * rank;
		endIdx = defaultPortion * (rank + 1);
	} else if (defaultPortion * rank <= globalNum
			&& defaultPortion * (rank + 1) >= globalNum) {

		startIdx = defaultPortion * rank;
		endIdx = globalNum;
	} else {
		startIdx = endIdx = globalNum;
	}

}

/************************************************************************/
// returns the process' physical memory usage in Mbyte.
int process_RAM_usage() {

	using namespace std;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat", ios_base::in);

	if (!stat_stream) {
		cerr << "In function commonFunctions::process_RAM_usage can't open"
				<< "file '/proc/self/stat'!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime
			>> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue
			>> starttime >> vsize >> rss; // don't care about the rest

	long page_size = sysconf(_SC_PAGE_SIZE);
	int mem = convertByteToMByte(rss * page_size);

	return (mem);

}

/************************************************************************/
// returns the process' memory usage in Mbyte.
int process_mem_usage() {

	using namespace std;

	// 'file' stat seems to give the most reliable results

	ifstream stat_stream("/proc/self/stat", ios_base::in);

	if (!stat_stream) {
		cerr << "In function commonFunctions::process_mem_usage can't open"
				<< "file '/proc/self/stat'!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// dummy vars for leading entries in stat that we don't care about

	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want

	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime
			>> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue
			>> starttime >> vsize >> rss; // don't care about the rest

	long page_size = sysconf(_SC_PAGE_SIZE);
	long int virtualMem = vsize; // virtual memory
	long int physicalMem = rss * page_size; // physical memory allocated

	int mem = convertByteToMByte(virtualMem) + convertByteToMByte(physicalMem);

	return (mem);

}

/************************************************************************/
// returns the system's total physical memory in Byte
/*int getTotalRAM() {

 FILE* stream = popen( "head -n1 /proc/meminfo", "r" );
 std::ostringstream output;
 int bufsize = 128;

 while( !feof( stream ) && !ferror( stream )) {
 char buf[bufsize];
 int bytesRead = fread( buf, 1, bufsize, stream );
 output.write( buf, bytesRead );
 }

 std::string result = output.str();

 std::string label, ram;
 std::istringstream iss(result);
 iss >> label;
 iss >> ram;

 return (atoi(ram.c_str())*1024);
 }*/

/************************************************************************/
// returns the total physical memory of one node in Mbyte
int getTotalRAM() {

	using namespace std;

	string key1, key2;
	long int value;

	ifstream inputFile;
	inputFile.open("/proc/meminfo");

	if (!inputFile) {
		cerr << "In function commonFunctions::getTotalRAM can't open"
				<< "file '/proc/meminfo'!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	while (!inputFile.eof()) {

		inputFile >> key1 >> value >> key2;
		//cout<<"total-key "<<key1<<" value "<<value<<" "<<key2<<endl;

		if (key1 == "MemTotal:")
			break;

	}

	int mem = convertkByteToMByte(value);

	inputFile.close();
	return (mem);
}

/************************************************************************/
// returns the free physical memory of one node in Mbyte
int getFreeRAM() {

	using namespace std;

	string key1, key2;
	long int value;

	ifstream inputFile;
	inputFile.open("/proc/meminfo");

	if (!inputFile) {
		cerr << "In function commonFunctions::getFreeRAM can't open"
				<< "file '/proc/meminfo'!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int mem = 0;
	bool flagFree = false;
	bool flagCached = false;

	while (!inputFile.eof()) {

		inputFile >> key1 >> value >> key2;
		//cout<<"free-key "<<key1<<" value "<<value<<" "<<key2<<endl;

		if (flagFree && flagCached)
			break;

		else if (key1 == "MemFree:") {

			mem += convertkByteToMByte(value);
			flagFree = true;

		}

		else if (key1 == "Cached:") {

			mem += convertkByteToMByte(value);
			flagCached = true;

		}

	}

	inputFile.close();
	return (mem);
}

/************************************************************************/
// check the RAM usage
void checkRAMusage(std::ofstream& logFile) {

	using namespace std;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int usedRAM = process_RAM_usage();
	int freeRAM = getFreeRAM();
	int totalRAM = getTotalRAM();

	if ((double) freeRAM / totalRAM < 0.15) {

		logFile << "In function checkRAMusage free memory insufficient:\n"
				<< "used RAM (by process) = " << usedRAM << " MB\n"
				<< "free RAM (on mainboard) = " << freeRAM << " MB\n"
				<< "total RAM (on mainboard) = " << totalRAM << " MB" << endl;

		cout << "used RAM (by process " << rank << ") = " << usedRAM << " MB\n"
				<< "free RAM (on mainboard) = "
				<< (double) freeRAM / totalRAM * 100.0 << " %" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	else {

		logFile << "used RAM (by process) = " << usedRAM << " MB\n"
				<< "free RAM (on mainboard) = " << freeRAM << " MB\n"
				<< "total RAM (on mainboard) = " << totalRAM << " MB" << endl;

		int proc = 0;
		int usedGlobalRAM = 0;
		intVector recvBuf(size);
		MPI_Allgather(&usedRAM, 1, MPI_INT, &recvBuf[0], 1, MPI_INT,
				MPI_COMM_WORLD);

		for (int i = 0; i < size; i++)

			if (recvBuf[i] > usedGlobalRAM) {

				usedGlobalRAM = recvBuf[i];
				proc = i;
			}

		if (rank == proc)
			cout << "max used RAM (by process " << proc << ") = "
					<< usedGlobalRAM << " MB\n" << "free RAM (on mainboard) = "
					<< (double) freeRAM / totalRAM * 100.0 << " %" << endl;

	}

}

/************************************************************************/
// check the current local RAM usage in Mbyte
void checkRAMusage(int numberOfProcsPerNode, std::ofstream& logFile) {

	using namespace std;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int usedRAM = process_RAM_usage();

	// ---------------------------------------------------------------------
	// available physical memory
	int totalRAM = getTotalRAM();

	// compute locally available RAM
	totalRAM /= numberOfProcsPerNode;

	// ---------------------------------------------------------------------
	// currently free memory

	int freeRAM = getFreeRAM(); // free memory in Mbyte

	// compute locally free RAM
	freeRAM /= numberOfProcsPerNode;

	logFile << "used RAM (by process) = " << usedRAM << " MB ("
			<< (double) usedRAM / totalRAM * 100 << " %)\n"
			<< "free RAM (for process) = " << freeRAM << " MB ("
			<< (double) freeRAM / totalRAM * 100 << " %)" << endl;

}

int convertByteToMByte(long int value) {

	return (value / (int) pow(2.0, 20));
}

int convertkByteToMByte(long int value) {

	return (value / (int) pow(2.0, 10));
}

/************************************************************************/
// Geometric information
//convert from cartesian to prolate spheroidal coordinates.
void cartToProSpheroid(dbVector& cartCoords, dbVector& proSpheroidCoords,
		std::map<std::string, double>& modelData, std::ofstream& logFile) {
	using namespace std;

	double C = modelData["prolateC"];
	int usedDims = (int) modelData["usedDimensions"];
	double pi = getPi();

	double x = cartCoords[0];
	double y = cartCoords[1];
	double z = cartCoords[2];

	allocateArray(proSpheroidCoords, usedDims);

	if (x == 0 && y == 0)
		proSpheroidCoords[2] = 0;
	else if (x >= 0)
		proSpheroidCoords[2] = atan(y / x);
	else
		proSpheroidCoords[2] = pi + atan(y / x);

	double b = -(pow(x, 2) + pow(y, 2) + pow(z, 2) + pow(C, 2)) / pow(C, 2);
	double c = pow(z, 2) / pow(C, 2);

	proSpheroidCoords[0] = acosh(sqrt((-b + sqrt(pow(b, 2) - 4 * c)) / 2));

	if (z >= 0)
		proSpheroidCoords[1] = acos(sqrt((-b - sqrt(pow(b, 2) - 4 * c)) / 2));
	else {
		//TEMP: problem when argument approx -1
		if (-sqrt((-b - sqrt(pow(b, 2) - 4 * c)) / 2) <= -0.999999
				&& -sqrt((-b - sqrt(pow(b, 2) - 4 * c)) / 2) > -1.01)
			proSpheroidCoords[1] = acos(-0.99999999999999);
		else
			proSpheroidCoords[1] = acos(
					-sqrt((-b - sqrt(pow(b, 2) - 4 * c)) / 2));
	}

}

//convert from prolate spheroidal to cartesian coordinates.
void proSpheroidToCart(dbVector& cartCoords, dbVector& proSpheroidCoords,
		std::map<std::string, double>& modelData) {
	using namespace std;

	double C = modelData["prolateC"];
	int usedDims = (int) modelData["usedDimensions"];

	double mu = proSpheroidCoords[0];
	double nu = proSpheroidCoords[1];
	double phi = proSpheroidCoords[2];

	allocateArray(cartCoords, usedDims);

	cartCoords[0] = C * sinh(mu) * sin(nu) * cos(phi);
	cartCoords[1] = C * sinh(mu) * sin(nu) * sin(phi);
	cartCoords[2] = C * cosh(mu) * cos(nu);

}

//get prolate spheroidal unit basis vectors at a point.
void getProSpheroidVecs(dbVector& proSpheroidCoords, dbMatrix& unitVectors,
		std::map<std::string, double>& modelData) {
	using namespace std;

	double C = modelData["prolateC"];
	int usedDims = (int) modelData["usedDimensions"];

	double mu = proSpheroidCoords[0];
	double nu = proSpheroidCoords[1];
	double phi = proSpheroidCoords[2];

	allocateArray(unitVectors, usedDims, usedDims);

	double scale = 1
			/ sqrt(
					pow(cosh(mu), 2) * pow(sin(nu), 2)
							+ pow(sinh(mu), 2) * pow(cos(nu), 2));

	//radial unit vector
	unitVectors[0][0] = scale * cosh(mu) * sin(nu) * cos(phi);
	unitVectors[0][1] = scale * cosh(mu) * sin(nu) * sin(phi);
	unitVectors[0][2] = scale * sinh(mu) * cos(nu);

	//Longitudinal unit vector
	unitVectors[1][0] = -scale * sinh(mu) * cos(nu) * cos(phi);
	unitVectors[1][1] = -scale * sinh(mu) * cos(nu) * sin(phi);
	unitVectors[1][2] = +scale * cosh(mu) * sin(nu);

	scale = 1 / (sinh(mu) * sin(nu));
	//circumfrential unit vector
	unitVectors[2][0] = -scale * sinh(mu) * sin(nu) * sin(phi);
	unitVectors[2][1] = scale * sinh(mu) * sin(nu) * cos(phi);
	unitVectors[2][2] = 0;

}

/************************************************************************/

//convert from cartesian to prolate spheroidal coordinates more accurately
void cart2ProSpheroid(dbVector& cartCoords, dbVector& proSpheroidCoords,
		std::map<std::string, double>& modelData, std::ofstream& logFile) {
	using namespace std;

	double C = modelData["prolateC"];
	int usedDims = (int) modelData["usedDimensions"];
	double pi = getPi();

	double x = cartCoords[0];
	double y = cartCoords[1];
	double z = cartCoords[2];

	allocateArray(proSpheroidCoords, usedDims);

	proSpheroidCoords[2] = atan2(y, x);
	if (proSpheroidCoords[2] < 0)
		proSpheroidCoords[2] = 2 * pi + proSpheroidCoords[2];

	double d1 = sqrt(pow((z + C), 2) + pow(x, 2) + pow(y, 2));
	double d2 = sqrt(pow((z - C), 2) + pow(x, 2) + pow(y, 2));

	double a = (d1 - d2) / (2 * C);
	double b = (d1 + d2) / (2 * C);

	proSpheroidCoords[0] = acosh(b);
	proSpheroidCoords[1] = acos(a);

	//  if(z>=0)
	//    proSpheroidCoords[1] = acos(sqrt((-b-sqrt(pow(b,2)-4*c))/2));
	//  else{
	//    //TEMP: problem when argument approx -1
	//    if(-sqrt((-b-sqrt(pow(b,2)-4*c))/2) <= -0.999999 && -sqrt((-b-sqrt(pow(b,2)-4*c))/2) > -1.01)
	//      proSpheroidCoords[1] = acos(-0.99999999999999);
	//    else
	//      proSpheroidCoords[1] = acos(-sqrt((-b-sqrt(pow(b,2)-4*c))/2));
	//  }

}

/***********************************************************************/

//combine n sorted vectors of differing lengths into a single vector
//recursive function
//rows represent each vector;
void combineSortedVecs(intMatrix& inputMat) {
	using namespace std;

	//sizes could be a problem, hopefully they will be the same or problems don't occur...
	//perhaps quite costly copying...
	int numVecs = inputMat.size();
	if (inputMat.size() != 1) {
		for (int i = 0; i < numVecs / 2; i++) {
			//reduce the two vectors into 1 and save to position of first
			// vectors 2*i and 2*i+1 are combined and then fill vector i
			combine2SortedVecs(inputMat[2 * i], inputMat[2 * i + 1]);
			inputMat[i] = inputMat[2 * i];
		}
		//the additional vector is moved to numVecs/2 + 1 and resize matrix
		if (numVecs % 2 == 1) {
			inputMat[numVecs / 2] = inputMat[numVecs - 1];
			inputMat.resize(numVecs / 2 + 1);
		} else
			inputMat.resize(numVecs / 2);

		combineSortedVecs(inputMat);
	}

}

// combines sorted vectors inputVec1 and inputVec2 and puts the sorted vector into inputVec1
void combine2SortedVecs(intVector& inputVec1, intVector& inputVec2) {
	using namespace std;

	int startIdx = 0;
	for (int i = 0; i < inputVec2.size(); i++) {
		for (int j = startIdx; j < inputVec1.size(); j++) {
			if (inputVec2[i] == inputVec1[j]) {
				startIdx = j;
				break;
			} else if (inputVec2[i] < inputVec1[j]) {
				inputVec1.insert(inputVec1.begin() + j, inputVec2[i]);
				startIdx = j;
				break;
			} else if (j + 1 == inputVec1.size()) {
				inputVec1.insert(inputVec1.end(), inputVec2[i]);
				startIdx = j + 1;
				break;
			}
		}
	}
}

// ========================================================================
// ========================================================================
// Ritesh common functions

// Print matrix
//void printMatrix(dbMatrix& matrix) {
//
//	using namespace std;
//
//	cout << "Displaying matrix of size [" << matrix.size() << " x "
//			<< matrix[0].size() << " ]:" << endl;
//	for (int i = 0; i < matrix.size(); i++) {
//			for (int j = 0; j < matrix[i].size(); j++) {
//			cout << matrix[i][j] << "\t";
//		}
//		cout << endl;
//	}
//	cout << endl;
//}
//
//void printMatrix(dbMatrix& matrix, std::ofstream& logFile) {
//
//	using namespace std;
//
//	logFile << "Size:" << matrix.size() << " x " << matrix[0].size() << endl;
//	for (int i = 0; i < matrix.size(); i++) {
//			for (int j = 0; j < matrix[i].size(); j++) {
//			logFile << matrix[i][j] << "\t";
//		}
//		logFile << endl;
//	}
//}
//
//void printMatrix(dbMatrix& matrix, std::string& msg) {
//
//	using namespace std;
//
//	cout << msg << " [ " << matrix.size() << " x " << matrix[0].size()
//			<< " ]:" << endl;
//	for (int i = 0; i < matrix.size(); i++) {
//		for (int j = 0; j < matrix[i].size(); j++) {
//			cout << matrix[i][j] << "\t";
//		}
//		cout << endl;
//	}
//	cout << endl;
//}

void printMatrix(dbMatrix& matrix, const char* msg, std::ofstream& logFile) {

	using namespace std;

	int nRow = matrix.size();
	int nCol = 0;

	if (nRow > 0)
		nCol = matrix[0].size();

	logFile << msg << " [ " << nRow << " x " << nCol << " ]:" << endl;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			logFile << matrix[i][j] << "\t";
		}
		logFile << endl;
	}
	logFile << endl;
}

void printMatrix(intMatrix& matrix, const char* msg, std::ofstream& logFile) {

	using namespace std;

	int nRow = matrix.size();
	int nCol = 0;

	if (nRow > 0)
		nCol = matrix[0].size();

	logFile << msg << " [ " << nRow << " x " << nCol << " ]:" << endl;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			logFile << matrix[i][j] << "\t";
		}
		logFile << endl;
	}
	logFile << endl;
}

//void printVector(dbVector& A) {
//
//	using namespace std;
//
//	cout << "Displaying vector of size:" << A.size() << endl;
//	for (int i = 0; i < A.size(); i++) {
//		cout << A[i] << "\t";
//	}
//	cout << endl << endl;
//}
//
//void printVector(dbVector& A, std::ofstream& logFile) {
//
//	using namespace std;
//
//	logFile << "Displaying vector of size:" << A.size() << endl;
//	for (int i = 0; i < A.size(); i++) {
//		logFile << A[i] << "\t";
//	}
//	logFile << endl << endl;
//}
//
//void printVector(dbVector& A, std::string& msg) {
//
//	using namespace std;
//
//	cout << msg << " [" << A.size() << "]" << endl;
//	for (int i = 0; i < A.size(); i++) {
//		cout << A[i] << "\t";
//	}
//	cout << endl << endl;
//}

void printVector(dbVector& A, const char* msg, std::ofstream& logFile) {

	using namespace std;

	logFile << msg << " [" << A.size() << "]" << endl;
	for (int i = 0; i < A.size(); i++) {
		logFile << A[i] << "\t";
	}
	logFile << endl << endl;
}

//void printVector(intVector& A) {
//
//	using namespace std;
//
//	cout << "Displaying vector of size:" << A.size() << endl;
//	for (int i = 0; i < A.size(); i++) {
//		cout << A[i] << "\t";
//	}
//	cout << endl << endl;
//}
//
//void printVector(intVector& A, std::ofstream& logFile) {
//
//	using namespace std;
//
//	logFile << "Displaying vector of size:" << A.size() << endl;
//	for (int i = 0; i < A.size(); i++) {
//		logFile << A[i] << "\t";
//	}
//	logFile << endl << endl;
//}
//
//void printVector(intVector& A, std::string& msg) {
//
//	using namespace std;
//
//	cout << msg << " [" << A.size() << "]" << endl;
//	for (int i = 0; i < A.size(); i++) {
//		cout << A[i] << "\t";
//	}
//	cout << endl << endl;
//}

void printVector(intVector& A, const char* msg, std::ofstream& logFile) {

	using namespace std;

	logFile << msg << " [" << A.size() << "]" << endl;
	for (int i = 0; i < A.size(); i++) {
		logFile << A[i] << "\t";
	}
	logFile << endl << endl;
}




