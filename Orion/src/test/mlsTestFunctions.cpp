/*
 * mlsTestFunctions.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: rama
 */

#include "mlsTestFunctions.h"

void mlsTest(InputFileData* InputData,ofstream& logFile){

	int choice = InputData->getValue("mlsTest_testChoice");

	switch(choice){

	case 1:
		mlsTest_one(InputData,logFile);
		break;

	case 2:
		mlsTest_two(InputData,logFile);
		break;

	default:
		logFile << "ERROR: mlsTest_testChoice value is incorrect. Check input file." << endl;
		cout << "ERROR: mlsTest_testChoice value is incorrect. Check input file." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

}

// *****************************************************************************
// *****************************************************************************
void mlsTest_one(InputFileData* InputData,ofstream& logFile){

	logFile << "******* mlsTest_one *******" << endl;

	DataContainer* MLSData = new DataContainer();

	vector<Particle> ptclList;
	create_domain(ptclList,MLSData,InputData,logFile);

	vector<Particle> ptcls;
	particlesSelection(ptcls,ptclList,MLSData,InputData,logFile);

	setupFunction(ptcls,ptclList,MLSData,InputData,logFile);

	// *************************************************************************
	// Calculating the approximate solutions
	findSupports(ptcls,ptclList,MLSData,InputData,logFile);

	calcInterpolants(ptcls,ptclList,MLSData,InputData,logFile);

	calcApproxPolynomialFunction(ptcls,ptclList,MLSData,InputData,logFile);

	// *************************************************************************

	errorCalculation(ptcls,ptclList,MLSData,InputData,logFile);

	plotShapeFunctions(ptcls,ptclList,MLSData,InputData,logFile);

}

// *****************************************************************************
// *****************************************************************************
void mlsTest_two(InputFileData* InputData,ofstream& logFile){

	logFile << "******* mlsTest_two *******" << endl;

	DataContainer* MLSData = new DataContainer();

	vector<Particle> ptclList;
	create_domain(ptclList,MLSData,InputData,logFile);


	// create sampling points
	vector<Particle> ptcls;
	create_samplingPoints(ptcls,MLSData,InputData,logFile);

	setupFunction(ptcls,ptclList,MLSData,InputData,logFile);

	// *************************************************************************
	// Calculating the approximate solutions
	findSupports(ptcls,ptclList,MLSData,InputData,logFile);

	calcInterpolants(ptcls,ptclList,MLSData,InputData,logFile);

	calcApproxPolynomialFunction(ptcls,ptclList,MLSData,InputData,logFile);

	errorCalculation(ptcls,ptclList,MLSData,InputData,logFile);

	int nDim = InputData->getValue("mlsTest_nDim");
	if(nDim == 1){
		plot1DShapeFunctionsOfSamplePoints(ptcls,ptclList,MLSData,InputData,logFile);
		plot1DErrorOfSamplePoints(ptcls,ptclList,MLSData,InputData,logFile);
	}
	else if(nDim == 2){
		plot2DShapeFunctionsOfSamplePoints(ptcls,ptclList,MLSData,InputData,logFile);
		plot2DErrorOfSamplePoints(ptcls,ptclList,MLSData,InputData,logFile);
	}

}

// *****************************************************************************
// *****************************************************************************
void create_domain(vector<Particle>& ptclList,DataContainer* MLSData,InputFileData* InputData,
		ofstream& logFile){

	logFile << "******* create_domain *******" << endl;

	int nDim = InputData->getValue("mlsTest_nDim");
	int numCoords;
	dbMatrix coords;

	if(InputData->getValue("mlsTest_readDomainFromFile") == 1){
		coords = read_nodes(MLSData,InputData,logFile);
	}else{
		if(InputData->getValue("mlsTest_radomNodes") == 1){
			coords = create_randomNodes(MLSData,InputData,logFile);
		}else{
			coords = create_nodes(MLSData,InputData,logFile);
		}
	}

	numCoords = coords.size();

	// =========================================================================
	if (InputData->getValue("mlsTest_saveDomainToFile") == 1) {

		string fileName = "mlsTest_domainNodes.grf";
		ofstream writeGraphRes(fileName);

		for (int i = 0; i < coords.size(); i++) {
			for (int j = 0; j < nDim; j++) {
				writeGraphRes << coords[i][j] << " ";
			}
			writeGraphRes << endl;
		}
	}
	// =========================================================================

	MLSData->setValue("numCoords",numCoords);

	ptclList = vector<Particle>(numCoords,Particle(nDim));
	for(int i=0; i<numCoords; i++){
		ptclList[i].getCoords() = coords[i];
	}

	// Printing Coordinates
	logFile << " ******* Domain *******" << endl;
	logFile << "Number of points: " << ptclList.size() << endl;
	for(int i = 0 ; i < ptclList.size(); i++){
		logFile << "Point[" << i << "] ";
		dbVector& coords = ptclList[i].getCoords();
		for(int j=0; j < coords.size(); j++){
			logFile << coords[j] << "\t" ;
		}
		logFile << endl;
	}
}

// *****************************************************************************
// *****************************************************************************
dbMatrix read_nodes(DataContainer* MLSData, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "******* read_nodes *******" << endl;

	cout << "Reading domain points from file" << endl;

	string fileName = "mlsTest_domainNodes.grf";
	ifstream readDomainNodes(fileName);
	string line;

	dbMatrix coords;

	if(!readDomainNodes){
		logFile << "Domain file: " << fileName << " has not been found" << endl;
		cout << "Domain file: " << fileName << " has not been found" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	getline(readDomainNodes, line);
	while (readDomainNodes.good()) {

		stringstream numLine(line);
		double num;
		string numStr;

		dbVector coordVec;

		while (getline(numLine,numStr, ' ')) {
			num = stod(numStr);
//			logFile << num << endl;
			coordVec.push_back(num);
		}

		coords.push_back(coordVec);

		getline(readDomainNodes, line);
	}

	printMatrix(coords,"read coords",logFile);

	readDomainNodes.close();

	if(coords.size() == 0){
		logFile << "No coordinates were read from: " << fileName << endl;
		cout << "No coordinates were read from: " << fileName << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	return coords;
}

// *****************************************************************************
// *****************************************************************************
dbMatrix create_nodes(DataContainer* MLSData,InputFileData* InputData,
		ofstream& logFile){

	logFile << "******* create_nodes *******" << endl;

	cout << "Generating domain points" << endl;

	int nDim = InputData->getValue("mlsTest_nDim");

	// Input parameters
	double startCoord = InputData->getValue("mlsTest_DomainStart");
	double endCoord = InputData->getValue("mlsTest_DomainEnd");

	double nodal_dist = InputData->getValue("mlsTest_nodalSpacing");

	// Output parameters
	dbMatrix coords;

	intVector SBM_dim(1);
	SBM_dim[0] = nDim;
	dbVector SBM_coeff(1);
	SBM_coeff[0] = 1;

//	Interpolation MLSTestInt.MLSUnitTest_meshMultiDim(length, nodal_dist,
//			nDim, dbMatrix& coords, int& numCoords, intVector SBM_dim,
//			dbVector SBM_coeff, ofstream& logFile);

	// =========================================================================
	// Define the coordinate list in a specific dimension
	dbVector coordList;
	for (int i = 0; startCoord + (i * nodal_dist) <= endCoord; i++)
		coordList.push_back(startCoord + (i * nodal_dist));

	printVector(coordList, "coordList", logFile);

	dbMatrix x;
	x.push_back(coordList);
	if (nDim > 1) {
		for (int i = 1; i < nDim; i++) {
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
		coords[i].resize(nDim);
		for (int j = 0; j < nDim; j++)
			coords[i][j] = x[j][i];
	}

	return coords;

}

// *****************************************************************************
// *****************************************************************************
dbMatrix create_randomNodes(DataContainer* MLSData,InputFileData* InputData,
		ofstream& logFile){

	logFile << "******* create_randomNodes *******" << endl;

	cout << "Generating domain points" << endl;

	int nDim = InputData->getValue("mlsTest_nDim");

	// Input parameters
	double startCoord = InputData->getValue("mlsTest_DomainStart");
	double endCoord = InputData->getValue("mlsTest_DomainEnd");

	double nodal_dist = InputData->getValue("mlsTest_nodalSpacing");

	// Output parameters
	dbMatrix coords;

	intVector SBM_dim(1);
	SBM_dim[0] = nDim;
	dbVector SBM_coeff(1);
	SBM_coeff[0] = 1;

//	Interpolation MLSTestInt.MLSUnitTest_meshMultiDim(length, nodal_dist,
//			nDim, dbMatrix& coords, int& numCoords, intVector SBM_dim,
//			dbVector SBM_coeff, ofstream& logFile);

	// =========================================================================
	// Define the coordinate list in a specific dimension
	dbVector coordList;
	for (int i = 0; startCoord + (i * nodal_dist) <= endCoord; i++)
		coordList.push_back(startCoord + (i * nodal_dist));

	printVector(coordList, "coordList", logFile);

	dbMatrix x;
	x.push_back(coordList);
	if (nDim > 1) {
		for (int i = 1; i < nDim; i++) {
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
		coords[i].resize(nDim);
		for (int j = 0; j < nDim; j++)
			coords[i][j] = x[j][i];
	}

	printMatrix(coords,"coords[before]",logFile);

	srand(time(0));

	// Make points random
	for (int i = 0; i < coords.size(); i++) {
		for (int j = 0; j < nDim; j++){
			if(coords[i][j] != startCoord && coords[i][j] != endCoord){

				double rNum = ((double)rand()/RAND_MAX);

				double coef = 0;
				if(rNum > 0.5)
					coef = 1;
				else
					coef = -1;

				double noise = coef * rNum * nodal_dist * 0.5;
				logFile << "coef: " << coef << endl;
				logFile << "rNum: " << rNum << endl;
				logFile << "nodal_dist: " << nodal_dist << endl;
				logFile << "Noise: " << noise << endl;

				coords[i][j] = coords[i][j] + noise;
			}
		}
	}

	printMatrix(coords,"coords[after]",logFile);

//	ofstream writeDomainNodes("mlsTest_DomainRandomNodes.grf");
//	for (int i = 0; i < coords.size(); i++) {
//			for (int j = 0; j < nDim; j++){
//				writeDomainNodes << coords[i][j] << "\t";
//			}
//			writeDomainNodes << endl;
//	}


//	cout << "Controlled ending" << endl;
//	MPI_Abort(MPI_COMM_WORLD, 1);


	return coords;

}

// *****************************************************************************
// *****************************************************************************
void create_samplingPoints(vector<Particle>& ptcls,DataContainer* MLSData,
		InputFileData* InputData,ofstream& logFile){

	logFile << "******* create_samplingPoints *******" << endl;

	double old_nodal_dist = InputData->getValue("mlsTest_nodalSpacing");
	double new_nodal_dist = InputData->getValue("mlsTest_nodalSampleSpacing");

	InputData->setValue("mlsTest_nodalSpacing",new_nodal_dist);

	int domainSource = InputData->getValue("mlsTest_readDomainFromFile");
	InputData->setValue("mlsTest_readDomainFromFile",0);

	int saveDomain = InputData->getValue("mlsTest_saveDomainToFile");
	InputData->setValue("mlsTest_saveDomainToFile",0);

	int randomPoint = InputData->getValue("mlsTest_radomNodes");
	InputData->setValue("mlsTest_radomNodes",0);

	DataContainer* temp = new DataContainer();
	create_domain(ptcls,temp,InputData,logFile);

	InputData->setValue("mlsTest_radomNodes",randomPoint);

	InputData->setValue("mlsTest_saveDomainToFile",saveDomain);

	InputData->setValue("mlsTest_readDomainFromFile",domainSource);

	InputData->setValue("mlsTest_nodalSpacing",old_nodal_dist);

}

// *****************************************************************************
// *****************************************************************************
void particlesSelection(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,	InputFileData* InputData, ofstream& logFile){

	logFile << "******* particlesSelection *******" << endl;

	int pSelectType = InputData->getValue("mlsTest_particleSelectionType");

	switch(pSelectType){
	case 1:
	{
		int nRandPoints = InputData->getValue("mlsTest_numOfRandomPoints");
		int nCoords = MLSData->getInt("numCoords");

		intVector rPointList(nRandPoints,0);

		for(int i = 0 ; i < nRandPoints; i++){
			srand (time(NULL));
			rPointList[i] = rand() % (nCoords-1) + 0;
			ptcls.push_back(ptclList[rPointList[i]]);
			ptcls[i].getID() = rPointList[i];
		}

		logFile << "Random point chosen: " ;
		for(int i=0; i< nRandPoints; i++){
			logFile << rPointList[i] << " " ;
		}
		logFile << endl;

		break;

	}
	case 2:
	{
		int nSpecPoints = InputData->getValue("mlsTest_numOfSpecificPoints");
		int nCoords = MLSData->getInt("numCoords");

		intVector sPointList(nSpecPoints,0);

		string temp = "mlsTest_specificPointID";
		for(int i=0; i < nSpecPoints; i++){
			string inputName = temp + std::to_string(i);
			sPointList[i] = InputData->getValue(inputName.c_str());

			if(sPointList[i] < nCoords){
				ptcls.push_back(ptclList[sPointList[i]]);
				ptcls[i].getID() = sPointList[i];
			}
			else{
				logFile << "'" << inputName
						<<"' is greater than the number of nodes: "
						<< nCoords << endl;
				cout << "'" << inputName
						<<"' is greater than the number of nodes: "
						<< nCoords << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}

		break;
	}

	case 3:
	{
		int nDim = InputData->getValue("mlsTest_nDim");

		string coord = "mlsTest_coord";
		dbVector coords(nDim,0);

		for(int i=0; i<nDim; i++){
			string lookupName = coord + std::to_string(i);
			coords[i] = InputData->getValue(lookupName.c_str());
		}

		ptcls.resize(1,Particle(nDim));
		ptcls[0].getCoords() = coords;

		break;
	}

	case 4:

		for(int i = 0 ; i < ptclList.size(); i++){
			ptcls.push_back(ptclList[i]);
			ptcls[i].getID() = i;
		}

		break;

	default:
		logFile << "ERROR: mlsTest_particleSelectionType valid choices are:"
				" 1, 2 and 3" << endl;
		cout << "ERROR: mlsTest_particleSelectionType valid choices are:"
				" 1, 2 and 3" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

}

// *****************************************************************************
// *****************************************************************************
void setupFunction(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,	InputFileData* InputData, ofstream& logFile){

	logFile << "******* setupPolynomialFunction *******" << endl;

	int nDim = InputData->getValue("mlsTest_nDim");

	int nPoints = ptclList.size();

	dbVector exactPolyValuePtclList(nPoints);
	dbVector exactPolyValuePtcls(ptcls.size());

	int polyFuncCoefType = InputData->getValue("mlsTest_polyFuncCoefType");
	double polyFuncCoef = InputData->getValue("mlsTest_polyFuncCoef");

	// Set the coefficients of the polynomial
	dbVector polyCoefVec(nDim);
	if (polyFuncCoefType == 0) {
		for (int i = 0; i < nDim; i++) {
			srand(time(NULL));
			polyCoefVec[i] = rand() % 10 + 1;
		}
	} else if (polyFuncCoefType == 1) {
		for (int i = 0; i < nDim; i++) {

			polyCoefVec[i] = polyFuncCoef;
		}

	} else {
		logFile << "'mlsTest_polyFuncType' is invalid" << endl;
		cout << "'mlsTest_polyFuncType' is invalid" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	printVector(polyCoefVec, "Polynomial Function Coefficients", logFile);

	MLSData->setValue("polynomialCoefficientVector",polyCoefVec);

	// =========================================================================
	for (int i = 0; i < nPoints; i++) {
			dbVector& coords = ptclList[i].getCoords();
			exactPolyValuePtclList[i] = calcFunction(coords,MLSData,InputData,
																	   logFile);
	}

	for (int i = 0; i < ptcls.size(); i++) {
		dbVector& coords = ptcls[i].getCoords();
		exactPolyValuePtcls[i] = calcFunction(coords,MLSData,InputData,logFile);
	}
	// =========================================================================
	/*
	int polyFuncDeg = InputData->getValue("mlsTest_polyFuncDeg");
	// Calculate the exact polynomial function values
	for (int i = 0; i < nPoints; i++) {
		dbVector& coords = ptclList[i].getCoords();
		for (int j = 0; j < nDim; j++) {
			exactPolyValuePtclList[i] += (polyCoefVec[j]
					* pow(coords[j], polyFuncDeg) );
//			exactPolyValuePtclList[i] += sin(coords[j]);
		}
		exactPolyValuePtclList[i] += 100000;
//		exactPolyValuePtclList[i] += 100;
	}

	for (int i = 0; i < ptcls.size(); i++) {
		dbVector& coords = ptcls[i].getCoords();
		for (int j = 0; j < nDim; j++) {
			exactPolyValuePtcls[i] += (polyCoefVec[j]
					* pow(coords[j], polyFuncDeg) );
//			exactPolyValuePtcls[i] += sin(coords[j]);
		}
		exactPolyValuePtcls[i] += 100000;
//		exactPolyValuePtcls[i] += 100 ;
	}
	*/
	// =========================================================================

	// Save the polynomial function value vector
	MLSData->setValue("exactPolyValuePtclList", exactPolyValuePtclList);
	MLSData->setValue("exactPolyValuePtcls", exactPolyValuePtcls);
}

// *****************************************************************************
// *****************************************************************************
double calcFunction(dbVector& coords,DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile){

	int nDim = InputData->getValue("mlsTest_nDim");
	int funcType = InputData->getValue("mlsTest_funcType");

	double value = 0;

	switch (funcType) {

	// -------------------------------------------------------------------------
	// Polynomial Function: f(x,y) = 100000 + (x^polyFuncDeg) + (y^polyFuncDeg)
	case 1: {
		int polyFuncDeg = InputData->getValue("mlsTest_polyFuncDeg");
		dbVector polyCoefVec = MLSData->getDbVector(
				"polynomialCoefficientVector");

		for (int j = 0; j < nDim; j++) {
			value += (polyCoefVec[j] * pow(coords[j], polyFuncDeg));
		}
		value += 100000;

		break;
	}

	// -------------------------------------------------------------------------
	// Polynomial Function: f(x,y) = 100 + sin(x) + sin(y)
	case 2: {

		for (int j = 0; j < nDim; j++) {
			value += sin(coords[j]);
		}
		value += 100;

		break;

	}
	default:
		logFile << "mlsTest_funcType is not valid." << endl;
		cout << "mlsTest_funcType is not valid." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);

	}

	return value;
}

// *****************************************************************************
// *****************************************************************************
bool isInInfluenceZone(Particle& ptcl_A, Particle& ptcl_B, double& influenceRadius,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile){


	double dist = 0;

	dbVector coordVec_A = ptcl_A.getCoords();
	dbVector coordVec_B = ptcl_B.getCoords();

	if(coordVec_A.size() == coordVec_B.size()){
		dist = calcDistance(coordVec_A,coordVec_B);
	}
	else{
		logFile << "In isInInfluenceZone, coordVec_A.size() != coordVec_B.size()"
				<< endl;
		cout << "In isInInfluenceZone, coordVec_A.size() != coordVec_B.size()"
		     << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if(dist <= influenceRadius)
		return true;
	else
		return false;


}

// *****************************************************************************
// *****************************************************************************
void findSupports(vector<Particle>& ptcls,vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData, ofstream& logFile){

	logFile << "******* findSupports *******" << endl;

	double influenceRadius = InputData->getValue("mlsTest_influenceRadius");
	double numPtclsSupport = InputData->getValue("mlsTest_numPtclsSupport");


	for(int i=0; i<ptcls.size(); i++){

		logFile << "**************************************************" << endl;
		logFile << "Point of interest [" << i << "]: ";
		for(int j=0; j<ptcls[i].getCoords().size(); j++){
			logFile << ptcls[i].getCoords()[j] << " ";
		}
		logFile << endl;

		dbVector ptclDistance(ptclList.size(),0);

		for(int j=0; j<ptclList.size(); j++){
			ptclDistance[j] = calcDistance(ptcls[i].getCoords(),
														ptclList[j].getCoords());
			logFile << "Point[" << j << "]: ";
			for(int k=0; k<ptclList[j].getCoords().size(); k++){
				logFile << ptclList[j].getCoords()[k] << " ";
			}
			logFile << "= " << ptclDistance[j] << endl;
		}

//		dbVector sPtclDistance(numPtclsSupport,numeric_limits<double>::max());
//		intVector sPtcls(numPtclsSupport,0);
//
//		for(int j=0; j<ptclList.size(); j++){
//			for(int k=0; k<sPtcls.size(); k++){
//				if(i != j && sPtclDistance[k] > ptclDistance[j]){
//					sPtclDistance[k] = ptclDistance[j];
//					sPtcls[k] = j;
//					break;
//				}
//			}
//		}

		dbVector sPtclDistance;
		intVector sPtcls;

		for(int j=0; j<ptclList.size(); j++){

//			if(ptcls[i].getCoords() != ptclList[j].getCoords() &&
//					influenceRadius >= ptclDistance[j]){
			if(influenceRadius >= ptclDistance[j]){
				sPtcls.push_back(j);
			}
		}

		ptcls[i].getSupportPtcls() = sPtcls;

		logFile << "Particle[" << i << "] Supports: ";
		for(int j = 0 ; j < ptcls[i].getSupportPtcls().size(); j++){
			logFile << ptcls[i].getSupportPtcls()[j] << ", ";
		}
		logFile << endl;
	}

}

// *****************************************************************************
// *****************************************************************************
void calcInterpolants(vector<Particle>& ptcls, vector<Particle>& ptclList,
		DataContainer* MLSData, InputFileData* InputData, ofstream& logFile) {

	logFile << "******* calcInterpolants *******" << endl;

	int nDim = InputData->getValue("mlsTest_nDim");
	double radFactor = InputData->getValue("mlsTest_radiusFactor");

	double influenceRadius = InputData->getValue("mlsTest_influenceRadius");
	int radType = InputData->getValue("mlsTest_radiusType");

	for (int i = 0; i < ptcls.size(); i++) {

		// Determine the support radius of each particle
		dbVector radVec(nDim, 0);

		intVector& sPtcls = ptcls[i].getSupportPtcls();

		dbVector& pCoords = ptcls[i].getCoords();

		dbMatrix sCoordList(sPtcls.size(), dbVector());

		if (radType == 0) {

			for (int j = 0; j < sPtcls.size(); j++) {

				dbVector& sCoords = ptclList[sPtcls[j]].getCoords();
				sCoordList[j] = sCoords;

				for (int k = 0; k < sCoords.size(); k++) {
					double dist = abs(pCoords[k] - sCoords[k]);
					if (dist > radVec[k])
						radVec[k] = dist * radFactor;
				}

			}

		} else if (radType == 1) {

			for (int j = 0; j < sPtcls.size(); j++) {
				sCoordList[j] = ptclList[sPtcls[j]].getCoords();
			}

			for (int k = 0; k < nDim; k++) {
				radVec[k] = influenceRadius;
			}
		}

		printVector(radVec, "final chosen radius vector", logFile);

		dbVector interpolantVec;
		Interpolation* MLSInt = new Interpolation(pCoords, sCoordList, radVec,
				interpolantVec, InputData, logFile);

		ptcls[i].getShapeFuncs() = interpolantVec;

	}
}

// *****************************************************************************
// *****************************************************************************
void calcApproxPolynomialFunction(vector<Particle>& ptcls,
		vector<Particle>& ptclList,DataContainer* MLSData,
		InputFileData* InputData,ofstream& logFile){

	logFile << "******* calcApproxPolynomialFunction *******" << endl;

	dbVector exactPolyValuePtclList = MLSData->getDbVector("exactPolyValuePtclList");
	dbVector approxPolyValuesPtcls(ptcls.size(),0);

	for(int i=0; i<ptcls.size(); i++){

		logFile << "--------------" << endl;
		logFile << "Particle[" << i << "]" << endl;

		intVector& sPtcls = ptcls[i].getSupportPtcls();
		dbVector& interpolants = ptcls[i].getShapeFuncs();

		for(int j=0; j<sPtcls.size(); j++){
			approxPolyValuesPtcls[i] += exactPolyValuePtclList[sPtcls[j]] * interpolants[j];
			logFile << "exactPolyValueList[sPtcls[j]]: " << exactPolyValuePtclList[sPtcls[j]] << endl;
			logFile << "interpolants[j]: " << interpolants[j] << endl;
			logFile << "approxPolyValues: " << approxPolyValuesPtcls[i] << endl;

		}
	}

	MLSData->setValue("approxPolyValuesPtcls",approxPolyValuesPtcls);

}


// *****************************************************************************
// *****************************************************************************
void errorCalculation(vector<Particle>& ptcls, vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile){

	logFile << "******* errorCalculation *******" << endl;

	int numCoords = ptcls.size();
	dbVector exactPolyValuePtclList = MLSData->getDbVector("exactPolyValuePtclList");
	dbVector exactPolyValuePtcls = MLSData->getDbVector("exactPolyValuePtcls");


	dbVector approxPolyValuesPtcls = MLSData->getDbVector("approxPolyValuesPtcls");

	dbMatrix exactPolyMat(numCoords,dbVector(1));
	dbMatrix approxPolyMat(numCoords,dbVector(1));

	for(int i=0; i<numCoords; i++){
		exactPolyMat[i][0] = exactPolyValuePtcls[i];
		approxPolyMat[i][0] = approxPolyValuesPtcls[i];
	}


	printMatrix(exactPolyMat,"exactPolyMat",logFile);
	printMatrix(approxPolyMat,"approxPolyMat",logFile);

	double L2Error = ErrorCalc::relativeL2norm(exactPolyMat,approxPolyMat,
															InputData, logFile);

	logFile << "L2 Error norm = " << L2Error << endl;
	cout << "L2 Error norm = " << L2Error << endl;

}

// *****************************************************************************
// *****************************************************************************
void plotShapeFunctions(vector<Particle>& ptcls, vector<Particle>& ptclList,
		DataContainer* MLSData,InputFileData* InputData,ofstream& logFile){

	logFile << "******* plotShapeFunctions *******" << endl;

	int nDim = InputData->getValue("mlsTest_nDim");

	string headerline = "#";
	for(int j=0; j < nDim; j++){
		headerline.append("x" + std::to_string(j) + "\t");
	}
	headerline.append("shapeFunctions");


	for(int i = 0 ; i < ptcls.size(); i++){

		string fileName = "particleShapeFunction" + std::to_string(i) + ".grf";
		ofstream writeGraphRes(fileName);

		// write header line first.
		writeGraphRes << headerline << endl;

		// assume the domain have zero shape function overall.
		dbVector shapeFuncVector(ptclList.size(),0);

		// extract the shape functions and fill them in the shapeFuncVector.
		intVector& sPtcls = ptcls[i].getSupportPtcls();
		dbVector& sPtclsShapeValues = ptcls[i].getShapeFuncs();

		for(int j=0; j < sPtcls.size(); j++){
			shapeFuncVector[sPtcls[j]] = sPtclsShapeValues[j];
		}


		// write coordinates and shape functions to file
		for(int j=0; j < ptclList.size(); j++){
			dbVector& coords = ptclList[j].getCoords();

			for(int k=0; k < coords.size(); k++){
				writeGraphRes << coords[k] << "\t";
			}
			writeGraphRes << shapeFuncVector[j] << endl;
		}

		writeGraphRes.close();

	}
}

// *****************************************************************************
// *****************************************************************************
void plot1DShapeFunctionsOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "******* plot1DShapeFunctionsOfSamplePoints *******" << endl;

	int nSpecPoints = InputData->getValue("mlsTest_numOfSpecificPoints");
	int nCoords = MLSData->getInt("numCoords");

	int nDim = InputData->getValue("mlsTest_nDim");
	string headerline = "#";
	for (int j = 0; j < nDim; j++) {
		headerline.append("x" + std::to_string(j) + "\t");
	}
	headerline.append("shapeFunctions");

	intVector selectedPointList(nSpecPoints, 0);

	string temp = "mlsTest_specificPointID";
	for (int i = 0; i < nSpecPoints; i++) {
		string inputName = temp + std::to_string(i);
		selectedPointList[i] = InputData->getValue(inputName.c_str());

		intVector selectedSamplePointList;
		dbVector selectedSamplePointListShapeFunction;
		if (selectedPointList[i] < nCoords && selectedPointList[i] > -1) {

			for (int j = 0; j < ptcls.size(); j++) {
				intVector& supportSamplePoints = ptcls[j].getSupportPtcls();
				dbVector& supportSamplePointsShapeFunc =
						ptcls[j].getShapeFuncs();
				for (int k = 0; k < supportSamplePoints.size(); k++) {
					if (supportSamplePoints[k] == selectedPointList[i]) {

						selectedSamplePointList.push_back(j);
						selectedSamplePointListShapeFunction.push_back(
								supportSamplePointsShapeFunc[k]);

						break;
					}
				}
			}

			printVector(selectedSamplePointList,"selectedSamplePointList",logFile);
			printVector(selectedSamplePointListShapeFunction,"selectedSamplePointListShapeFunction",logFile);

			string fileName = "particleShapeFunction" + std::to_string(i)
					+ ".grf";
			ofstream writeGraphRes(fileName);

			// write header line first.
			writeGraphRes << headerline << endl;

			// assume the domain have zero shape function overall.
			dbVector shapeFuncVector(ptcls.size(), 0);

			for (int j = 0; j < selectedSamplePointList.size(); j++) {
				shapeFuncVector[selectedSamplePointList[j]]
				                = selectedSamplePointListShapeFunction[j];
			}

			// write coordinates and shape functions to file
			for (int j = 0; j < ptcls.size(); j++) {
				dbVector& coords = ptcls[j].getCoords();

				for (int k = 0; k < coords.size(); k++) {
					writeGraphRes << coords[k] << "\t";
				}
				writeGraphRes << shapeFuncVector[j] << endl;
			}

			writeGraphRes.close();

		} else {
			logFile << "'" << inputName
					<< "' is greater than the number of nodes: " << nCoords
					<< endl;
			cout << "'" << inputName
					<< "' is greater than the number of nodes: " << nCoords
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

}

// *****************************************************************************
// *****************************************************************************
void plot1DErrorOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "******* plot1DErrorOfSamplePoints *******" << endl;

	int numCoords = ptcls.size();

	dbVector exactPolyValuePtcls = MLSData->getDbVector("exactPolyValuePtcls");
	dbVector approxPolyValuesPtcls = MLSData->getDbVector("approxPolyValuesPtcls");

	printVector(exactPolyValuePtcls,"exactPolyValuePtcls",logFile);
	printVector(approxPolyValuesPtcls,"approxPolyValuesPtcls",logFile);

	dbVector errorPolyValuesPtcls(numCoords,0);

	for(int i=0; i<numCoords; i++){
		errorPolyValuesPtcls[i]
		     = fabs((approxPolyValuesPtcls[i]-exactPolyValuePtcls[i])
		                     	 	 	 	 	 	 /exactPolyValuePtcls[i]);

		logFile << "errorPolyValuesPtcls[i]: " << errorPolyValuesPtcls[i] << endl;
	}

	printVector(errorPolyValuesPtcls,"errorPolyValuesPtcls",logFile);

	// =========================================================================

	int nSpecPoints = InputData->getValue("mlsTest_numOfSpecificPoints");
	int nCoords = MLSData->getInt("numCoords");

	int nDim = InputData->getValue("mlsTest_nDim");
	string headerline = "#";
	for (int j = 0; j < nDim; j++) {
		headerline.append("x" + std::to_string(j) + "\t");
	}
	headerline.append("shapeFunctions");

	intVector selectedPointList(nSpecPoints, 0);

	string fileName_err = "particleErrorValues.grf";
	ofstream writeGraphErr(fileName_err);
	writeGraphErr.precision(15);

	string fileName_ext = "particleExactValues.grf";
	ofstream writeGraphExt(fileName_ext);
	writeGraphExt.precision(15);

	string fileName_app = "particleApproxValues.grf";
	ofstream writeGraphApp(fileName_app);
	writeGraphApp.precision(15);

	// write header line first.
	writeGraphErr << headerline << endl;
	writeGraphExt << headerline << endl;

	// write coordinates and shape functions to file
	printVector(errorPolyValuesPtcls,"errorPolyValuesPtcls",logFile);
	for (int j = 0; j < ptcls.size(); j++) {

		dbVector& coords = ptcls[j].getCoords();

		for (int k = 0; k < coords.size(); k++) {
			writeGraphErr << coords[k] << "\t";
			writeGraphExt << coords[k] << "\t";
			writeGraphApp << coords[k] << "\t";
		}

		logFile << "errorPolyValuesPtcls[j]: " << errorPolyValuesPtcls[j] << endl;
		writeGraphErr << errorPolyValuesPtcls[j] << endl;
		writeGraphExt << exactPolyValuePtcls[j] << endl;
		writeGraphApp << approxPolyValuesPtcls[j] << endl;
	}

	writeGraphErr.close();
	writeGraphExt.close();
	writeGraphApp.close();
}

// *****************************************************************************
// *****************************************************************************
void plot2DShapeFunctionsOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile){

	logFile << "******* plot2DShapeFunctionsOfSamplePoints *******" << endl;

	int nSpecPoints = InputData->getValue("mlsTest_numOfSpecificPoints");
	int nCoords = MLSData->getInt("numCoords");

	int nDim = InputData->getValue("mlsTest_nDim");
	string headerline = "#";
	for (int j = 0; j < nDim; j++) {
		headerline.append("x" + std::to_string(j) + "\t");
	}
	headerline.append("shapeFunctions");

	intVector selectedPointList(nSpecPoints, 0);

	string temp = "mlsTest_specificPointID";
	for (int i = 0; i < nSpecPoints; i++) {
		string inputName = temp + std::to_string(i);
		selectedPointList[i] = InputData->getValue(inputName.c_str());

		intVector selectedSamplePointList;
		dbVector selectedSamplePointListShapeFunction;
		if (selectedPointList[i] < nCoords && selectedPointList[i] > -1) {

			for (int j = 0; j < ptcls.size(); j++) {
				intVector& supportSamplePoints = ptcls[j].getSupportPtcls();
				dbVector& supportSamplePointsShapeFunc =
						ptcls[j].getShapeFuncs();
				for (int k = 0; k < supportSamplePoints.size(); k++) {
					if (supportSamplePoints[k] == selectedPointList[i]) {

						selectedSamplePointList.push_back(j);
						selectedSamplePointListShapeFunction.push_back(
								supportSamplePointsShapeFunc[k]);

						break;
					}
				}
			}

			printVector(selectedSamplePointList,"selectedSamplePointList",logFile);
			printVector(selectedSamplePointListShapeFunction,"selectedSamplePointListShapeFunction",logFile);

			string fileName = "particleShapeFunction" + std::to_string(i)
					+ ".grf";
			ofstream writeGraphRes(fileName);

			// write header line first.
			writeGraphRes << headerline << endl;

			// assume the domain have zero shape function overall.
			dbVector shapeFuncVector(ptcls.size(), 0);

			for (int j = 0; j < selectedSamplePointList.size(); j++) {
				shapeFuncVector[selectedSamplePointList[j]]
				                = selectedSamplePointListShapeFunction[j];
			}

			// write coordinates and shape functions to file
			double current, previous;
			for (int j = 0; j < ptcls.size(); j++) {

				dbVector& coords = ptcls[j].getCoords();
				current = coords[coords.size()-1];

				if(j != 0 && previous != current){
					previous = current;
					writeGraphRes << endl;
				}


				for (int k = 0; k < coords.size(); k++) {
					writeGraphRes << coords[k] << "\t";
				}
				writeGraphRes << shapeFuncVector[j] << endl;

				previous = current;
			}

			writeGraphRes.close();

		} else {
			logFile << "'" << inputName
					<< "' is greater than the number of nodes: " << nCoords
					<< endl;
			cout << "'" << inputName
					<< "' is greater than the number of nodes: " << nCoords
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
}

// *****************************************************************************
// *****************************************************************************
void plot2DErrorOfSamplePoints(vector<Particle>& ptcls,
		vector<Particle>& ptclList, DataContainer* MLSData,
		InputFileData* InputData, ofstream& logFile) {

	logFile << "******* plot2DErrorOfSamplePoints *******" << endl;

	int numCoords = ptcls.size();

	dbVector exactPolyValuePtcls = MLSData->getDbVector("exactPolyValuePtcls");
	dbVector approxPolyValuesPtcls = MLSData->getDbVector("approxPolyValuesPtcls");

	dbVector errorPolyValuesPtcls(numCoords,0);

	for(int i=0; i<numCoords; i++){
		errorPolyValuesPtcls[i]
		     = abs((approxPolyValuesPtcls[i]-exactPolyValuePtcls[i])
		                      	 	 	 	 	 	 /exactPolyValuePtcls[i]);
	}

	// =========================================================================

	int nSpecPoints = InputData->getValue("mlsTest_numOfSpecificPoints");
	int nCoords = MLSData->getInt("numCoords");

	int nDim = InputData->getValue("mlsTest_nDim");
	string headerline = "#";
	for (int j = 0; j < nDim; j++) {
		headerline.append("x" + std::to_string(j) + "\t");
	}
	headerline.append("shapeFunctions");

	intVector selectedPointList(nSpecPoints, 0);

	string fileName_err = "particleErrorValues.grf";
	ofstream writeGraphErr(fileName_err);
	writeGraphErr.precision(15);

	string fileName_ext = "particleExactValues.grf";
	ofstream writeGraphExt(fileName_ext);
	writeGraphExt.precision(15);

	string fileName_app = "particleApproxValues.grf";
	ofstream writeGraphApp(fileName_app);
	writeGraphApp.precision(15);

	// write header line first.
	writeGraphErr << headerline << endl;
	writeGraphExt << headerline << endl;
	writeGraphApp << headerline << endl;

	// write coordinates and shape functions to file
	double current, previous;
	for (int j = 0; j < ptcls.size(); j++) {

		dbVector& coords = ptcls[j].getCoords();
		current = coords[coords.size() - 1];

		if (j != 0 && previous != current) {
			previous = current;
			writeGraphErr << endl;
			writeGraphExt << endl;
			writeGraphApp << endl;
		}

		for (int k = 0; k < coords.size(); k++) {
			writeGraphErr << coords[k] << "\t";
			writeGraphExt << coords[k] << "\t";
			writeGraphApp << coords[k] << "\t";
		}
		writeGraphErr << errorPolyValuesPtcls[j] << endl;
		writeGraphExt << exactPolyValuePtcls[j] << endl;
		writeGraphApp << approxPolyValuesPtcls[j] << endl;

		previous = current;
	}

	writeGraphErr.close();
	writeGraphExt.close();
	writeGraphApp.close();
}
