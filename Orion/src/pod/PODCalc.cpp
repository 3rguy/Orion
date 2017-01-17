#include "PODCalc.h"

using namespace std;

/*!****************************************************************************/
/*!****************************************************************************/
//! PODCalc constructor
PODCalc::PODCalc(dbMatrix& fullMatrix, InputFileData* InputData,
		ofstream& logFile){

	// POD Preprocessing
	setDataMatrix(fullMatrix,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
//! PODCalc constructor
PODCalc::PODCalc(dbMatrix& fullMatrix,dbMatrix& reducedMatrix, double enLev,
		InputFileData* InputData, ofstream& logFile){

	// POD Preprocessing
	setDataMatrix(fullMatrix,InputData,logFile);

	// Define the energy level to be conserved
	setEnergyLevel(enLev);

	// Calculate the POMs and POVs
	setPOVsandPOMs(InputData,logFile);

	// Determine the number of POMs needed to meet the minimum energy level
	// requirement
	energyConservCalc(logFile);

	// Reduce Matrix's size
	compressMatrix(reducedMatrix,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void PODCalc::setDataMatrix(dbMatrix& mat,InputFileData* InputData,
		ofstream& logFile){

	dataMat = mat;

	if(InputData->getValue("PODMeanCalculation") == 1){

		meanVec = dbVector(mat.size(),0);

		for(int r = 0; r <mat.size(); r++){
			for(int c = 0; c <mat[r].size(); c++){
				meanVec[r] += mat[r][c];
			}
			meanVec[r] = meanVec[r] / mat[r].size();
		}

		for(int r = 0; r <mat.size(); r++){
			for(int c = 0; c <mat[r].size(); c++){
				dataMat[r][c] = mat[r][c] - meanVec[r];
			}
		}
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the Singular Value Decomposition and set the POVs and POMs
void PODCalc::SVDCalc(ofstream& logFile) {

	// Carry out the SVD of matrix "dataMat"
	dbMatrix V;
	dbVector singularValues;

	SVD(dataMat,POMs,V,singularValues,logFile);

	POVs.resize(singularValues.size());

	int n = singularValues.size();
	logFile << "number of singular values, n: " << n << endl;
	for(int i=0; i<singularValues.size();i++){
		POVs[i] = (singularValues[i]*singularValues[i])/n;
		logFile << "POVs[i] = (" << singularValues[i] <<"*" << singularValues[i] << ")/"<< n <<" = " << POVs[i] << endl;
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the Singular Value Decomposition of the A matrix
void PODCalc::SVD(dbMatrix& A, dbMatrix& U, dbMatrix& V, dbVector& S,
				  ofstream& logFile) {

	int flag = computeGeneralMatrixSVD(A, S, U, V, logFile);

	if (flag != 0){
		logFile << "SVD calculation fail" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	else
		logFile << "SVD calculation successful" << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the Karhunen-Loevre Decomposition and set the POVs and POMs
void PODCalc::KLDCalc(ofstream& logFile) {

	KLD(logFile);
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the POVs and POMs using the method of snapshots
void PODCalc::snapshotCalc(ofstream& logFile) {

	/* Algorithm adapted from:
		N.J. Falkiewicz and E.S. Cesnik, Proper Orthogonal Decomposition for
		Reduction-Order Thermal Solution in Hypersonic Aerothermoelastic
		Simulations. AIAA Journal, 994:1009(Vol.49,No.5), May 2011
	*/

	// Calculate (A^T)A
	dbMatrix snapshotMat;
//	printMatrix(dataMat,"dataMat",logFile);
	innerTensorProduct(dataMat, dataMat, snapshotMat, true, false, logFile);
//	printMatrix(snapshotMat,"snapshotMat",logFile);

	// Compute 1/n * (A^T)A
	int n = dataMat[0].size();
	for (int i = 0; i < snapshotMat.size(); i++) {
		for (int j = 0; j < snapshotMat[i].size(); j++) {
			snapshotMat[i][j] = snapshotMat[i][j] / (n);
		}
	}

//	printMatrix(snapshotMat,"snapshotMat",logFile);

	// Carry out the SVD of the snapshot matrix
	dbMatrix U, V;
	dbVector S;
	SVD(snapshotMat, U, V, S, logFile);

//	printMatrix(U,"U",logFile);
//	printVector(S,"S",logFile);
//	printMatrix(V,"V",logFile);

	// Calculate the POM
	POMs.resize(dataMat.size(), dbVector(S.size(), 0));
	dbVector eigVector(U.size(),0);
	dbVector POMs_i;
	for (int i = 0; i < U[0].size(); i++) {

		// Select one vector
		for (int j = 0; j < U.size(); j++)
			eigVector[j] = U[j][i];

		innerTensorProduct(dataMat,eigVector,POMs_i,false,logFile);
//		printVector(POMs_i,"POMs_i",logFile);

		for (int k = 0; k < POMs_i.size(); k++){
//			logFile << "sqrt(S[k]*n): " << sqrt(S[i]*n) <<endl;
			POMs_i[k] = POMs_i[k]/sqrt(S[i]*n);
		}

		for (int l = 0; l < POMs_i.size(); l++)
			POMs[l][i] = POMs_i[l];
	}

	// Assign the POVs
	POVs = S;
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the POVs and POMs using the method of snapshots
void PODCalc::snapshotCalc_two(ofstream& logFile) {

	/* Algorithm adapted from:
		N.J. Falkiewicz and E.S. Cesnik, Proper Orthogonal Decomposition for
		Reduction-Order Thermal Solution in Hypersonic Aerothermoelastic
		Simulations. AIAA Journal, 994:1009(Vol.49,No.5), May 2011
	*/

	// Calculate (A^T)A
	dbMatrix snapshotMat;
//	printMatrix(dataMat,"dataMat",logFile);
	innerTensorProduct(dataMat, dataMat, snapshotMat, true, false, logFile);
//	printMatrix(snapshotMat,"snapshotMat",logFile);

	// Compute 1/n * (A^T)A
	int n = dataMat[0].size();
	for (int i = 0; i < snapshotMat.size(); i++) {
		for (int j = 0; j < snapshotMat[i].size(); j++) {
			snapshotMat[i][j] = snapshotMat[i][j] / (n);
		}
	}

	// Carry out eigenvalue decomposition
	dbVector eigenvalues;
	dbMatrix eigenvectors;
	calcEigenValuesVectors(snapshotMat,eigenvalues,eigenvectors,logFile);

	reverse(eigenvalues.begin(),eigenvalues.end());
	POVs = eigenvalues;

	for(int i = 0 ; i < eigenvectors.size(); i++){
		reverse(eigenvectors[i].begin(),eigenvectors[i].end());
	}

	POMs = dbMatrix(dataMat.size(),dbVector(dataMat[0].size(),0));

	for (int i = 0; i < eigenvectors[0].size(); i++) {

		dbVector POMs_i;
		dbVector eigVector(eigenvectors.size(),0);

		// Select one vector
		for (int j = 0; j < eigenvectors.size(); j++)
			eigVector[j] = eigenvectors[j][i];

		innerTensorProduct(dataMat,eigVector,POMs_i,false,logFile);

		for (int k = 0; k < POMs_i.size(); k++){
			POMs_i[k] = POMs_i[k]/sqrt(POVs[i]*n);
		}

		for (int l = 0; l < POMs_i.size(); l++)
			POMs[l][i] = POMs_i[l];
	}

#ifdef _PODCalcDebugMode_
	printVector(POVs,"POVs(before clean)",logFile);
	printMatrix(POMs,"POMs(before clean)",logFile);
#endif

	cleanPOVsAndPOMs(POVs,POMs,logFile);

#ifdef _PODCalcDebugMode_
	printVector(POVs,"POVs(after clean)",logFile);
	printMatrix(POMs,"POMs(after clean)",logFile);
#endif
}

/*!****************************************************************************/
/*!****************************************************************************/
void PODCalc::cleanPOVsAndPOMs(dbVector& vec,dbMatrix& mat,ofstream& logFile){

	intVector idList;

	for(int i=0; i<vec.size(); i++){
		if(vec[i] < 0)
			idList.push_back(i);
	}

	for(int i = idList.size()-1; i == 0 ; i++){
		vec.erase(vec.begin() + idList[i]);

		for(int j = 0 ; j < mat.size(); j++){
			mat[j].erase(mat[j].begin() + idList[i]);
		}
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the Karhunen-Loevre Decomposition of the A matrix
void PODCalc::KLD(ofstream& logFile) {

	/* Algorithm adapted from:
		N.J. Falkiewicz and E.S. Cesnik, Proper Orthogonal Decomposition for
		Reduction-Order Thermal Solution in Hypersonic Aerothermoelastic
		Simulations. AIAA Journal, 994:1009(Vol.49,No.5), May 2011
	*/

	// Calculate A(A^T)

	printMatrix(dataMat,"dataMat",logFile);

	dbMatrix KLDMat;
//	printMatrix(dataMat,"dataMat",logFile);
	innerTensorProduct(dataMat, dataMat, KLDMat, false, true, logFile);
//	printMatrix(KLDMat,"KLDMat",logFile);

	// Compute 1/n * A(A^T)
	int n = dataMat[0].size();
	for (int i = 0; i < KLDMat.size(); i++) {
		for (int j = 0; j < KLDMat[i].size(); j++) {
			KLDMat[i][j] = KLDMat[i][j] / (n);
		}
	}

	dbVector eigenvalues;
	dbMatrix eigenvectors;
	calcEigenValuesVectors(KLDMat,eigenvalues,eigenvectors,logFile);

	printVector(eigenvalues,"eigenvalues",logFile);
	printMatrix(eigenvectors,"eigenvectors",logFile);

	reverse(eigenvalues.begin(),eigenvalues.end());
	POVs = eigenvalues;

	for(int i = 0 ; i < eigenvectors.size(); i++){
		reverse(eigenvectors[i].begin(),eigenvectors[i].end());
	}
	POMs = eigenvectors;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the Proper Orthogonal Values and Modes using the selected method
void PODCalc::setPOVsandPOMs(InputFileData* InputData,ofstream& logFile){

	int choice = InputData->getValue("PODCalculationType") ;

	switch (choice) {
	case 1:
		SVDCalc(logFile);
		break;

	case 2:
		KLDCalc(logFile);
		break;

	case 3:
		snapshotCalc_two(logFile);
		break;

	default:
		cout << "Specified PODCalculationType["<< choice <<"] is invalid" << endl;
		break;
	}
}

/*!****************************************************************************/
/*!****************************************************************************/
//! Retrieve the POMs conserved based on the energy level
dbMatrix& PODCalc::getPOMsConserved(ofstream& logFile){

	dbMatrix POMsSelected;
	resizeArray(POMsSelected,numOfPOVsConserved,POMs[0].size());

	for(int i=0;i<numOfPOVsConserved;i++){
		for(int j=0;j<POMs[i].size();j++){
			POMsSelected[i][j] = POMs[i][j];
		}
	}

	return POMsSelected;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Retrieve the POVs conserved based on the energy level
dbVector& PODCalc::getPOVsConserved(ofstream& logFile){

	dbVector POVsSelected;
	resizeArray(POVsSelected,numOfPOVsConserved);

	for(int i=0;i<numOfPOVsConserved;i++){
		POVsSelected[i] = POVs[i];
	}

	return POVsSelected;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Set the number of POVs to be conserved, depending on the energy level
void PODCalc::energyConservCalc(ofstream& logFile){

	double totEnergy = 0;

	printVector(POVs,"POVs",logFile);

	//Calculate the total amount of energy
	for(int i=0;i<POVs.size();i++){
		totEnergy += POVs[i];
	}

	double percEnergy=0;
	double cumEnergy=0;

	for(int j=0;j < POVs.size();j++){

		// calculate the cumulative energy
		cumEnergy += POVs[j];

		// Calculate the percentage of energy
		energyConserved = (cumEnergy/totEnergy)*100;

		numOfPOVsConserved = j+1;

		// Check if energy level has been reached
		if(energyConserved >= energyLevel){
			break;
		}
	}

#ifdef _PODCalcDebugMode_
	logFile << "Set-up number of POVs to be conserved" << endl;
	logFile << "-------------------------------------" << endl;
	logFile << "Energy conserved: " << energyConserved <<endl;
	logFile << "Number of POVs conserved: " << numOfPOVsConserved <<endl;
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Reduce the size of the data matrix using the generated POMs and POVs
void PODCalc::compressMatrix(dbMatrix& reducedMatrix,InputFileData* InputData,
		ofstream& logFile){

#ifdef _PODCalcDebugMode_
	logFile << "******* PODCalc::compressMatrix -> Before Compression *******" << endl;

	logFile << "Data Matrix" << endl;
	logFile << "-----------" << endl;
	printMatrix(dataMat,"",logFile);
#endif

	reducedMatrix.resize(numOfPOVsConserved, dbVector(dataMat[0].size(),0));

	// Assemble selected POMs in Matrix form
	for(int i = 0; i < numOfPOVsConserved; i++){
		for(int j = 0; j < dataMat.size(); j++){
			for(int k = 0; k < dataMat[j].size(); k++){
				reducedMatrix[i][k] += POMs[j][i] * dataMat[j][k];
			}
		}
	}

#ifdef _PODCalcDebugMode_
	logFile << "******* PODCalc::compressMatrix -> After Compression *******" << endl;

	logFile << "Reduced Matrix" << endl;
	logFile << "--------------" << endl;
	printMatrix(reducedMatrix,"",logFile);
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Reduce the size of a full matrix using the generated POMs and POVs
void PODCalc::compressMatrix(dbMatrix& fullMatrix,dbMatrix& reducedMatrix,
		InputFileData* InputData,ofstream& logFile){


	reducedMatrix.resize(numOfPOVsConserved, dbVector(fullMatrix[0].size(),0));

	// Assemble selected POMs in Matrix form
	for(int i = 0; i < numOfPOVsConserved; i++){
		for(int j = 0; j < fullMatrix.size(); j++){
			for(int k = 0; k < fullMatrix[j].size(); k++){
				reducedMatrix[i][k] += POMs[j][i] * fullMatrix[j][k];
			}
		}
	}

#ifdef _PODCalcDebugMode_
	logFile << "******* PODCalc::compressMatrix *******" << endl;

	logFile << "Full Matrix" << endl;
	logFile << "-----------" << endl;
	printMatrix(fullMatrix,"",logFile);

	logFile << "Reduced Matrix" << endl;
	logFile << "--------------" << endl;
	printMatrix(reducedMatrix,"",logFile);
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Reduce the size of a full vector using the generated POMs and POVs
void PODCalc::compressVector(dbVector& fullVector,dbVector& reducedVector,
		InputFileData* InputData,ofstream& logFile){

	reducedVector.resize(numOfPOVsConserved,0);

	// Assemble selected POMs in Matrix form
	for(int i = 0; i < numOfPOVsConserved; i++){
		for(int j = 0; j < fullVector.size(); j++){
			reducedVector[i] += POMs[j][i] * fullVector[j];
		}
	}

#ifdef _PODCalcDebugMode_
	logFile << "******* PODCalc::compressVector *******" << endl;

	logFile << "Full Vector" << endl;
	logFile << "-----------" << endl;
	printVector(fullVector,"",logFile);

	logFile << "Reduced Vector" << endl;
	logFile << "--------------" << endl;
	printVector(reducedVector,"",logFile);
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! expand any matrix using the Proper Orthogonal Modes
void PODCalc::expandMatrix(dbMatrix& reducedMatrix,dbMatrix& fullMatrix,
		InputFileData* InputData,ofstream& logFile){

#ifdef _checkMode_
	if(reducedMatrix.size() !=  numOfPOVsConserv){
		cout << "ERROR: In PODCalc::expandMatrix, the number of POMs "
				"preserved and the row number in reducedMatrix are not "
				"the same" << endl;
		logFile << "ERROR: In PODCalc::expandMatrix, the number of POMs "
				"preserved and the row number in reducedMatrix are not "
				"the same" << endl;

		MPI_Abort(MPI_COMM_WORLD,1);
	}
#endif

	fullMatrix.resize(reducedMatrix.size(),dbVector(reducedMatrix[0].size(),0));

	for(int i = 0; i < POMs.size(); i++){
		for(int k = 0; k < reducedMatrix[0].size(); k++){
			for(int j = 0; j < numOfPOVsConserved; j++){
				fullMatrix[i][k] += POMs[i][j] * reducedMatrix[j][k];
			}
		}
	}

#ifdef _PODCalcDebugMode_
	logFile << "******* PODCalc::expandMatrix *******" << endl;

	logFile << "Reduced Matrix" << endl;
	logFile << "--------------" << endl;
	printMatrix(reducedMatrix,"",logFile);

	logFile << "Full Matrix" << endl;
	logFile << "-----------" << endl;
	printMatrix(fullMatrix,"",logFile);
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! expand any matrix using the Proper Orthogonal Modes
void PODCalc::expandVector(dbVector& reducedVector,dbVector& fullVector,
		InputFileData* InputData,ofstream& logFile){

#ifdef _checkMode_
	if(reducedVector.size() !=  numOfPOVsConserv){
		cout << "ERROR: In PODCalc::expandVector, the number of POMs "
				"preserved and the row number in reducedVector are not "
				"the same" << endl;
		logFile << "ERROR: In PODCalc::expandVector, the number of POMs "
				"preserved and the row number in reducedVector are not "
				"the same" << endl;

		MPI_Abort(MPI_COMM_WORLD,1);
	}
#endif

	fullVector.resize(POMs.size(),0);

	for(int i = 0; i < POMs.size(); i++){
		for(int j = 0; j < numOfPOVsConserved; j++){
			fullVector[i] += POMs[i][j] * reducedVector[j];
		}
	}

#ifdef _PODCalcDebugMode_
	logFile << "******* PODCalc::expandVector *******" << endl;

	logFile << "Reduced Vector" << endl;
	logFile << "--------------" << endl;
	printVector(reducedVector,"",logFile);

	logFile << "Full Vector" << endl;
	logFile << "-----------" << endl;
	for(int i=0;i < fullVector.size(); i++){
		logFile << fullVector[i] << endl;
	}
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! expand any matrix using the Proper Orthogonal Modes
void PODCalc::addMean(dbVector& vec, InputFileData* InputData,
		ofstream& logFile) {

	if(vec.size() == meanVec.size()){
		for (int i = 0; i < vec.size(); i++) {
			vec[i] = vec[i] + meanVec[i];
		}
	}else{

		cout << "In PODCalc::addMean, size of vec is not equal to size of meanVec" << endl;
		logFile << "In PODCalc::addMean, size of vec is not equal to size of meanVec" << endl;
		MPI_Abort(MPI_COMM_WORLD,1);

	}

}
