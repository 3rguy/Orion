#include "CardiacFibreGenerator.h"

CardiacFibreGenerator::CardiacFibreGenerator(InputFileData* InputData,
		  std::ofstream& logFile){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  pi = getPi();

  if(rank == 0)
	  cout << "generating fibre distribution" << endl;
	logFile << "generating fibre distribution" << endl;

	// define the vertical axis
	z.resize(3,0);
	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = 1.0;

	vector<Particle> ptcls;

	vector<FEMElement> surfaceElems;
	dbVector surfaceAngles;

	readMeshFile(ptcls,surfaceElems,surfaceAngles,InputData,logFile);

	/**************************************************************************/
	// Calculate all surface normals
	calcAllSurfaceNormals(surfaceElems,InputData,logFile);

	/**************************************************************************/
	// Compute the averaged surface normal at each node
	intVector sNodes;
	calcNodalNormal(ptcls,surfaceElems,sNodes,InputData,logFile);

	/**************************************************************************/
	// Calculate the circumferential direction at each node
	calcNodalCircumDirections(ptcls,sNodes,InputData,logFile);

	/**************************************************************************/
	// Calculate the outward pointing normal at each node
	calcOutwardNormal(ptcls,sNodes,InputData,logFile);

	/**************************************************************************/
	// Calculate the sheet normal at each node
	calcSheetNormal(ptcls,sNodes,InputData,logFile);

	/**************************************************************************/
	// Calculate the fibre direction projection vector at each node
	calcFibreDirectProjection(ptcls,sNodes,InputData,logFile);

	/**************************************************************************/
	// Calculate the fibre direction vector at each node
	calcFibreDirection(ptcls,sNodes,InputData,logFile);

	/**************************************************************************/
	// Calculate the fibre direction vector at each node
	calcOrthogonalFibreDirection(ptcls,sNodes,InputData,logFile);

	/**************************************************************************/
	// Generate result for seska
	if(rank == 0)
	  writingResultToFile(ptcls,sNodes,InputData,logFile);

	// Generate result for GiD
	generateFEMresFile(InputData,logFile);

	if(rank == 0)
	  saveResultsToFile_flexible(ptcls,InputData,logFile);

  if(rank == 0)
    cout << "finshed generating fibre distribution" << endl;
  logFile << "finished generating fibre distribution" << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Normalise any given vector
dbVector CardiacFibreGenerator::normaliseVec(dbVector& vec, std::ofstream& logFile){

	using namespace std;

	double normVal = computeNorm(vec,2,logFile);

	dbVector aVec(vec.size());
	for(int i=0; i<vec.size();i++){
		aVec[i] = vec[i]/normVal;
	}

	return aVec;
}

/******************************************************************************/
/******************************************************************************/
void CardiacFibreGenerator::readMeshFile(vector<Particle>& ptcls, vector<FEMElement>& surfaceElems,
				  dbVector& surfaceAngles, InputFileData* InputData,
				  std::ofstream& logFile) {

	int ptcleNum, ID;
	string name, token;
	double value;

	vector<FEMElement> epiElems;
	vector<FEMElement> endoElems;
	vector<FEMElement>* surfElems = NULL;

	dbVector epiAngleVec;
	dbVector endoAngleVec;
	dbVector* angleVec = NULL;

	int usedDims = 3;
	string filename = "mesh.dat";
	double minCoord = 0;
	double maxCoord = 0;

	int elemType = 0;

	// Open file with particle coordinates.
	ifstream meshFile(filename.c_str());

	if (meshFile) {

		/**********************************************************************/
		// read sample particle coordinates
		meshFile >> name;
		meshFile >> ptcleNum;

		ptcls = vector<Particle>(ptcleNum, Particle(0));

		for (int i = 0; i < ptcls.size(); i++) {
			Particle& ptcle = ptcls[i];
			vector<double>& coords = ptcle.getCoords();
			clearArray(coords);

			if (usedDims == 3)
				meshFile >> ID >> name >> coords[0] >> coords[1] >> coords[2];
			else {
				cerr << "Only 3D supported.\n";
				MPI_Abort(MPI_COMM_WORLD, 1);
			}

		}

		/**********************************************************************/
		// read sample surface helix angle and nodes
		string line;
		string resultOne = "Epicardial_Helix_Angle";
		string resultTwo = "Endocardial_Helix_Angle";
		int surfNum = 0;
		while (meshFile.good()) {
			getline(meshFile, line);

			if (line == "Epicardial_Helix_Angle") {

				meshFile >> surfNum;

				getline(meshFile, line); // get rid of blank line

				epiElems.resize(surfNum, FEMElement(3));
				surfElems = &epiElems;

				elemType = 1;

				epiAngleVec.resize(surfNum);
				angleVec = &epiAngleVec;

			} else if (line == "Endocardial_Helix_Angle") {

				meshFile >> surfNum;

				getline(meshFile, line); // get rid of blank line

				endoElems.resize(surfNum, FEMElement(3));
				surfElems = &endoElems;

				elemType = 2;

				endoAngleVec.resize(surfNum);
				angleVec = &endoAngleVec;
			}

			if (surfElems != NULL) {
				intVector surfNodes;
				int nodeID;
				for (int i = 0; i < surfNum; i++) {

					getline(meshFile, line);

					istringstream iss_line(line);
					string nodes_str, surfID_str, angle_str;

					getline(iss_line, surfID_str, ':');
					getline(iss_line, angle_str, ':');
					getline(iss_line, nodes_str, ':');

					istringstream nodes_line(nodes_str);
					while (nodes_line >> nodeID) {
						surfNodes.push_back(nodeID);
					}

					surfElems->at(i).setGlobalID(atoi(surfID_str.c_str()));

					surfElems->at(i).getNodes() = surfNodes;

					surfElems->at(i).getElemType() = elemType;

					angleVec->at(i) = atof(angle_str.c_str());

					surfNodes.clear();
				}
				surfNum = 0;
				surfElems = NULL;
				angleVec = NULL;
			}
		}

	} else {
		cerr << "File 'mesh.dat' does not exist!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (epiAngleVec.size() == 0) {
		cerr << "In 'mesh.dat', epicardial_helix_angle missing!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	} else if (endoAngleVec.size() == 0) {
		cerr << "In 'mesh.dat', endocardial_helix_angle missing!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Assign coordinate list to each epi-cardial surface and angle to each
	// particle
	for (int i = 0; i < epiElems.size(); i++) {
		intVector& nodesVec = epiElems[i].getNodes();
		dbMatrix coordMat(nodesVec.size());
		for (int j = 0; j < nodesVec.size(); j++) {
			coordMat[j] = ptcls[nodesVec[j] - 1].getCoords();
			ptcls[nodesVec[j] - 1].setBeta(epiAngleVec[i]);
		}
		epiElems[i].getNodalCoords() = coordMat;
	}

	// Assign coordinate list to each endo-cardial surface and angle to each
	// particle
	for (int i = 0; i < endoElems.size(); i++) {
		intVector& nodesVec = endoElems[i].getNodes();
		dbMatrix coordMat(nodesVec.size());
		for (int j = 0; j < nodesVec.size(); j++) {
			coordMat[j] = ptcls[nodesVec[j] - 1].getCoords();
			ptcls[nodesVec[j] - 1].setBeta(endoAngleVec[i]);
		}
		endoElems[i].getNodalCoords() = coordMat;
	}

	generateFEMresFile(InputData, logFile);

	// Merge surface elements list
	surfaceElems = epiElems;
	surfaceElems.insert(surfaceElems.end(), endoElems.begin(), endoElems.end());

	// Merge surface angles list
	surfaceAngles = epiAngleVec;
	surfaceAngles.insert(surfaceAngles.end(), endoAngleVec.begin(),
															endoAngleVec.end());


	// ------------------------------------------------------------------------
	// Debugging stuff

#ifdef _anisotropyInterpolationDebugMode_
	logFile << "****************************************************" << endl;
	logFile << "**************** particle coordinates **************" << endl;
	for (int i = 0; i < ptcls.size(); i++) {
		Particle& ptcle = ptcls[i];
		vector<double>& coords = ptcle.getCoords();
		logFile << i + 1 << ": ";
		for (int j = 0; j < coords.size(); j++) {
			logFile << coords[j] << " ";
		}
		logFile << endl;
	}

	logFile << "****************************************************" << endl;
	logFile << "****************** surface elements ****************" << endl;
	logFile << "epicardial_helix_angle" << endl;
	for (int i = 0; i < epiElems.size(); i++) {
		logFile << "(" << i << ") SurfaceID[" << epiElems[i].getGlobalID()
				<< "|" << epiAngleVec[i] << "]: ";
		intVector& nodesList = epiElems[i].getNodes();
		for (int j = 0; j < nodesList.size(); j++) {
			logFile << nodesList[j] << "(";
			dbVector& ptclCoords = ptcls[nodesList[j] - 1].getCoords();
			for (int k = 0; k < ptclCoords.size(); k++) {
				logFile << ptclCoords[k] << ",";
			}
			logFile << ") | " ;
		}
		logFile << endl;
	}
	logFile << endl;

	logFile << "endocardial_helix_angle" << endl;
	for (int i = 0; i < endoElems.size(); i++) {
		logFile << "(" << i << ") SurfaceID[" << endoElems[i].getGlobalID()
				<< "|" << endoAngleVec[i] << "]: ";
		intVector& nodesList = endoElems[i].getNodes();
		for (int j = 0; j < nodesList.size(); j++) {
			logFile << nodesList[j] << "(";
			dbVector& ptclCoords = ptcls[nodesList[j] - 1].getCoords();
			for (int k = 0; k < ptclCoords.size(); k++) {
				logFile << ptclCoords[k] << ",";
			}
			logFile << ") | " ;
		}
		logFile << endl;
	}
	logFile << endl;
#endif

}

/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate all surface normal
void CardiacFibreGenerator::generateFEMresFile(InputFileData* InputData,
		std::ofstream& logFile){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0)
	system("cp fem.msh fem_fibre.msh");

}


/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate all surface normal
void CardiacFibreGenerator::calcAllSurfaceNormals(vector<FEMElement>& sElems,InputFileData* InputData,
		std::ofstream& logFile){

	// Compute the surface normal of each epi-cardial surface

#ifdef _anisotropyInterpolationDebugMode_
	logFile << "***** Surface Normal of elements *****" << endl;
#endif

	for (int i = 0; i < sElems.size(); i++) {

#ifdef _anisotropyInterpolationDebugMode_
		logFile << "Elem[" << sElems[i].getGlobalID() << "|" << i << "]"
				<< endl;
#endif

		sElems[i].getSurfaceNormal() = calcSurfaceNormal(
				sElems[i].getNodalCoords(), InputData, logFile);

#ifdef _anisotropyInterpolationDebugMode_
		logFile << "--------------------------------------------------" << endl;
#endif

	}
}


/*!****************************************************************************/
/*!****************************************************************************/
//! Calculate the surface normal from a set of nodal coordinates
dbVector CardiacFibreGenerator::calcSurfaceNormal(dbMatrix& surfacePointsCoords,
		InputFileData* InputData,
		std::ofstream& logFile){

	using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//	logFile << "surfacePointsCoords:" << endl;
//	for(int i = 0; i < surfacePointsCoords.size(); i++){
//		for(int j = 0; j < surfacePointsCoords[i].size(); j++){
//			logFile << surfacePointsCoords[i][j] << " ";
//		}
//		logFile << endl;
//	}

	if(surfacePointsCoords.size() < 3){
		logFile <<"In FEMGeometryExt::caclSurfaceNormal: Too few coordinates"
				" to compute surface normal"<<endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	dbVector vecOne(surfacePointsCoords[0].size());
	dbVector vecTwo(surfacePointsCoords[0].size());

	for(int i=0; i<surfacePointsCoords[0].size();i++){

		if(abs(surfacePointsCoords[1][i] - surfacePointsCoords[0][i])<1e-14)
			vecOne[i] = 0;
		else
			vecOne[i] = surfacePointsCoords[1][i] - surfacePointsCoords[0][i];

		if(abs(surfacePointsCoords[2][i] - surfacePointsCoords[0][i])<1e-14)
			vecTwo[i] = 0;
		else
			vecTwo[i] = surfacePointsCoords[2][i] - surfacePointsCoords[0][i];
	}

	dbVector surfaceNormalVec;
	crossProduct(vecOne,vecTwo,surfaceNormalVec);

	// Invert the direction of the surface normal
	for(int i = 0 ; i < surfaceNormalVec.size(); i++){
		surfaceNormalVec[i] = surfaceNormalVec[i]*-1;
	}

//	printVector(vecOne,"vecOne",logFile);
//	printVector(vecTwo,"vecTwo",logFile);
//	printVector(surfaceNormalVec,"surfaceNormalVec",logFile);

	return surfaceNormalVec;
}

/*!****************************************************************************/
/*!****************************************************************************/
// Compute the averaged surface normal at each node
void CardiacFibreGenerator::calcNodalNormal(vector<Particle>& ptcls, vector<FEMElement>& sElems,
		intVector& sNodes, InputFileData* InputData, std::ofstream& logFile) {

//	logFile << "***** Compiling all surface nodes *****" << endl;
	sNodes = compileSurfaceNodes(ptcls, sElems, InputData, logFile);
//	printVector(sNodes,"sNodes",logFile);

//	logFile << "***** Average nodal normals of all surface nodes *****" << endl;
	calcAveragedNodalNormals(ptcls, sElems, sNodes, InputData, logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
// Compile list of surface nodes
intVector CardiacFibreGenerator::compileSurfaceNodes(vector<Particle>& ptcls,
		vector<FEMElement>& elemVec, InputFileData* InputData,
		std::ofstream& logFile) {

	intVector nodalVec;

	int pos;
	for (int i = 0; i < elemVec.size(); i++) {
		intVector& nodes = elemVec[i].getNodes();
		for (int j = 0; j < nodes.size(); j++) {
			ptcls[nodes[j] - 1].getElems().push_back(i);
			ptcls[nodes[j] - 1].getMaterialID() = elemVec[i].getElemType();
			pos = findIntVecPos(nodes[j], 0, nodalVec.size(), nodalVec);
			if (pos == -1) {
				nodalVec.resize(nodalVec.size() + 1);
				nodalVec[nodalVec.size() - 1] = nodes[j];
			}
		}
	}

	sortIntVector(nodalVec, 0, nodalVec.size()-1);

	// Enforce epicardial nodes
	for (int i = 0; i < elemVec.size(); i++) {
		if (elemVec[i].getElemType() == 1) {
			intVector& nodes = elemVec[i].getNodes();
			for (int j = 0; j < nodes.size(); j++) {
				ptcls[nodes[j] - 1].getMaterialID() = elemVec[i].getElemType();
			}
		}
	}

	return nodalVec;
}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the nodal normals by averaging the normal of its
// surrounding/neighbouring surfaces.
void CardiacFibreGenerator::calcAveragedNodalNormals(vector<Particle>& ptcls,
		vector<FEMElement>& elemVec, intVector& sNodes,
		InputFileData* InputData, std::ofstream& logFile) {

	for (int i = 0; i < sNodes.size(); i++) {
		intVector surfElems = ptcls[sNodes[i] - 1].getElems();

		eliminateSpecificSurface(surfElems, int(1), elemVec, InputData,
				logFile);

		dbVector avgNormals(3, 0);
		int numOfSurfaces = surfElems.size();

//		logFile << "Node[" << sNodes[i] << "] with elements(" << numOfSurfaces
//				<< "): ";
//		for (int j = 0; j < surfElems.size(); j++) {
//			logFile << surfElems[j] << " ";
//		}
//		logFile << endl;

		// Calculate the average normal from the selected surfaces
		for (int j = 0; j < numOfSurfaces; j++) {

			dbVector& sNormal = elemVec[surfElems[j]].getSurfaceNormal();

//			logFile << "[" << surfElems[j] << "] ";
//			for (int l = 0; l < sNormal.size(); l++) {
//				logFile << sNormal[l] << ", ";
//			}
//			logFile << endl;

			for (int k = 0; k < sNormal.size(); k++) {
				avgNormals[k] += sNormal[k] / numOfSurfaces;
			}
		}

//		logFile << "[Avg]: ";
//		for (int l = 0; l < avgNormals.size(); l++) {
//			logFile << avgNormals[l] << ", ";
//		}
//		logFile << endl;

		avgNormals = normaliseVec(avgNormals, logFile);

		// store the averaged surface normal at each node
		ptcls[sNodes[i] - 1].getAllSurfaceNormals().resize(1);
		ptcls[sNodes[i] - 1].getAllSurfaceNormals()[0] = avgNormals;
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the nodal normals by averaging the normal of its
// surrounding/neighbouring surfaces.
void CardiacFibreGenerator::calcNodalCircumDirections(vector<Particle>& ptcls, intVector& sNodes,
	InputFileData* InputData, std::ofstream& logFile){

//	logFile << "##### calcNodalCircumDirections #####" << endl;

	double orientationIndexEpi = 1.0;
	double orientationIndexEndo = -1.0;

	// For all epi-cardial surface node
	for(int i = 0; i < sNodes.size(); i++){
//		logFile << "----- Node[" << sNodes[i]-1 << "] -----" <<endl;

		Particle& ptcle = ptcls[sNodes[i]-1];
		dbVector& n = ptcle.getAllSurfaceNormals()[0];
//		printVector(n,"n",logFile);

		dbVector c;
		crossProduct(z,n,c);
//		printVector(c,"c",logFile);

		if(ptcle.getMaterialID() == 1){
			std::transform(c.begin(), c.end(), c.begin(),
						   std::bind1st(std::multiplies<double>(),orientationIndexEpi));
		}
		else if (ptcle.getMaterialID() == 2){
			std::transform(c.begin(), c.end(), c.begin(),
					        std::bind1st(std::multiplies<double>(),orientationIndexEndo));
		}
		else if (ptcle.getMaterialID() > 2){
			logFile << "In calcNodalCircumDirections, materialID["
					<< ptcle.getMaterialID() << "] is not supported" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		ptcle.getAllSurfaceNormals().resize(2);
		ptcle.getAllSurfaceNormals()[1] = c;
//		printVector(ptcle.getAllSurfaceNormals()[1],"c+-",logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the outward pointing normal at each node
void CardiacFibreGenerator::calcOutwardNormal(vector<Particle>& ptcls, intVector& sNodes,
		InputFileData* InputData, std::ofstream& logFile){

//	logFile << "##### calcOutwardNormal #####" << endl;

	// For all surface nodes
	for(int i = 0; i < sNodes.size(); i++){
		Particle& ptcle = ptcls[sNodes[i]-1];
		dbVector n_cz;

//		printVector(ptcle.getAllSurfaceNormals()[1],"c",logFile);
		crossProduct(ptcle.getAllSurfaceNormals()[1],z,n_cz);

		ptcle.getAllSurfaceNormals().resize(3);
		ptcle.getAllSurfaceNormals()[2] = n_cz;
//		printVector(ptcle.getAllSurfaceNormals()[2],"n_cz",logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the sheet normal at each node
void CardiacFibreGenerator::calcSheetNormal(vector<Particle>& ptcls, intVector& sNodes,
		InputFileData* InputData, std::ofstream& logFile){

//	logFile << "##### calcSheetNormal #####" << endl;

	for(int i = 0; i < sNodes.size(); i++){
		Particle& ptcle = ptcls[sNodes[i]-1];
		dbVector& n = ptcle.getAllSurfaceNormals()[0];
		dbVector& n_cz = ptcle.getAllSurfaceNormals()[2];
		double val;

//		printVector(n_cz,"n_cz",logFile);
//		printVector(n,"n",logFile);
		scalarProduct(n_cz,n,val,logFile);
//		logFile << "dot(n_cz,n): " << val << endl;

		double sign_val = copysign(1.0,val);
//		logFile << "sign of dot(n_cz,n): " << sign_val << endl;

		dbVector s(n.size());
		std::transform(n.begin(), n.end(), s.begin(),
				std::bind1st(std::multiplies<double>(),sign_val));

		s = normaliseVec(s,logFile);

		ptcle.getAllSurfaceNormals().resize(4);
		ptcle.getAllSurfaceNormals()[3] = s;
//		printVector(ptcle.getAllSurfaceNormals()[3],"s",logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the outward pointing normal at each node
void CardiacFibreGenerator::calcFibreDirectProjection(vector<Particle>& ptcls, intVector& sNodes,
		InputFileData* InputData, std::ofstream& logFile){

//	logFile << "##### calcFibreDirectProjection #####" << endl;

	// For all surface nodes
	for(int i = 0; i < sNodes.size(); i++){
		Particle& ptcle = ptcls[sNodes[i]-1];

		dbVector& c = ptcle.getAllSurfaceNormals()[1];
		dbVector p(c.size());

		double cos_alpha = cos ( ptcle.getBeta() * pi / 180.0 );
		double sin_alpha = sin ( ptcle.getBeta() * pi / 180.0 );

//		printVector(c,"c",logFile);
//		logFile <<"cos_alpha: " << cos_alpha << endl;
//		logFile <<"sin_alpha: " << sin_alpha << endl;

		for(int j = 0; j < c.size(); j++){
			p[j] = (cos_alpha * c[j]) + (sin_alpha * z[j]);
		}

		ptcle.getAllSurfaceNormals().resize(5);
		ptcle.getAllSurfaceNormals()[4] = p;
//		printVector(ptcle.getAllSurfaceNormals()[4],"p",logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the outward pointing normal at each node
void CardiacFibreGenerator::calcFibreDirection(vector<Particle>& ptcls, intVector& sNodes,
		InputFileData* InputData, std::ofstream& logFile){

//	logFile << "##### calcFibreDirection #####" << endl;

	for (int i = 0; i < sNodes.size(); i++) {
		Particle& ptcle = ptcls[sNodes[i] - 1];

		dbVector n_cz = ptcle.getAllSurfaceNormals()[2];
		dbVector s = ptcle.getAllSurfaceNormals()[3];
		dbVector p = ptcle.getAllSurfaceNormals()[4];
//		printVector(n_cz,"n_cz",logFile);
//		printVector(s,"s",logFile);
//		printVector(p,"p",logFile);

		double p_dot_s;
		scalarProduct(p, s, p_dot_s, logFile);
//		logFile << "p_dot_s: " << p_dot_s << endl;

		double n_cz_dot_s;
		scalarProduct(n_cz, s, n_cz_dot_s, logFile);
//		logFile << "n_cz_dot_s: " << n_cz_dot_s << endl;

		dbVector f(n_cz.size());
		for (int j = 0; j < f.size(); j++) {
			f[j] = (-1 * p_dot_s * n_cz[j]) + (n_cz_dot_s * p[j]);
		}

		f = normaliseVec(f,logFile);

		ptcle.getAllSurfaceNormals().resize(6);
		ptcle.getAllSurfaceNormals()[5] = f;
//		printVector(f,"f",logFile);
	}

}

void CardiacFibreGenerator::calcOrthogonalFibreDirection(vector<Particle>& ptcls, intVector& sNodes,
		InputFileData* InputData, std::ofstream& logFile){

//	logFile << "##### calcFibreDirection #####" << endl;

	for (int i = 0; i < sNodes.size(); i++) {
		Particle& ptcle = ptcls[sNodes[i] - 1];

		dbVector s = ptcle.getAllSurfaceNormals()[3];
		dbVector f = ptcle.getAllSurfaceNormals()[5];
//		printVector(s,"s",logFile);
//		printVector(f,"f",logFile);

		dbVector m;
		crossProduct(f,s,m);
//		printVector(m,"m",logFile);

		m = normaliseVec(m,logFile);

		ptcle.getAllSurfaceNormals().resize(7);
		ptcle.getAllSurfaceNormals()[6] = m;
//		printVector(m,"m",logFile);
	}

}

/*!****************************************************************************/
/*!****************************************************************************/
// calculate the outward pointing normal at each node
void CardiacFibreGenerator::writingResultToFile(vector<Particle>& ptcls, intVector& sNodes,
		InputFileData* InputData,
		std::ofstream& logFile) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//	sortIntVector(surfaceNodes, 0, surfaceNodes.size()-1);

	/**************************************************************************/
	// write anisotropy direction in file
	string filename = "anisotropy.dat";
	ofstream anisotropyFile(filename.c_str());
	anisotropyFile.precision(15);

	// -------------------------------------------------------------------------
	// write coordinates
	anisotropyFile << "COORDINATES" << endl;
	anisotropyFile << sNodes.size() << endl;

	for (int i = 0; i < sNodes.size(); i++) {

		Particle& ptcle = ptcls[sNodes[i] - 1];
		vector<double>& coords = ptcle.getCoords();

		anisotropyFile << i + 1 << " : ";

		for (int j = 0; j < coords.size(); j++){
			anisotropyFile << coords[j] << " ";
		}
		anisotropyFile << endl;

	}

	// -------------------------------------------------------------------------
	// write directions
	anisotropyFile << "DIRECTIONS" << endl;
	anisotropyFile << sNodes.size() << endl;

	// -------------------------------------------------------------------------
	// loop over all ptcls
	for (int i = 0; i < sNodes.size(); i++) {

		Particle& ptcle = ptcls[sNodes[i] - 1];
		vector<double>& resultVec = ptcle.getAllSurfaceNormals()[5];

		anisotropyFile << i + 1 << " : transverse-isotropy" << endl;

		// loop over all directions
		for (int j = 0; j < resultVec.size(); j++) {
			anisotropyFile << resultVec[j] << " ";
		}
		anisotropyFile << endl;
	}

	// -------------------------------------------------------------------------
//	dbMatrix resultMat(3,dbVector());
//	// loop over all ptcls
//	for (int i = 0; i < sNodes.size(); i++) {
//
//		Particle& ptcle = ptcls[sNodes[i] - 1];
//		resultMat[0] = ptcle.getAllSurfaceNormals()[3];
//		resultMat[1] = ptcle.getAllSurfaceNormals()[5];
//		resultMat[2] = ptcle.getAllSurfaceNormals()[6];
//
//		anisotropyFile << i + 1 << " : orthotropy" << endl;
//
//		// loop over all directions
//		for (int j = 0; j < resultMat.size(); j++) {
//			for(int k = 0; k < resultMat[j].size(); k++){
//				anisotropyFile << resultMat[j][k] << " ";
//			}
//			anisotropyFile << endl;
//		}
//	}
//	// -------------------------------------------------------------------------

	anisotropyFile.flush();
	anisotropyFile.close();

}

/*!****************************************************************************/
/*!****************************************************************************/
void CardiacFibreGenerator::printVector(intVector& A, const char* msg, std::ofstream& logFile) {

	using namespace std;

	logFile << msg << " [" << A.size() << "]" << endl;
	for (int i = 0; i < A.size(); i++) {
		logFile << A[i] << "\t";
	}
	logFile << endl << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
void CardiacFibreGenerator::printVector(dbVector& A, const char* msg, std::ofstream& logFile) {

	using namespace std;

	logFile << msg << " [" << A.size() << "]" << endl;
	for (int i = 0; i < A.size(); i++) {
		logFile << A[i] << "\t";
	}
	logFile << endl << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
// Outputting calculated matrix to fem.res file
void CardiacFibreGenerator::saveResultsToFile_res_format(std::ofstream& femResFile,
		std::string& resultName, dbMatrix& resultMatrix, double& stepValue,
		InputFileData* InputData, ofstream& logFile) {


	string headerLine = "GiD Post Results File 1.0";

	string resultType = "Vector";
	string my_location = "OnNodes";
	string analysis_name = " ";

	if (femResFile.is_open()) {

		femResFile << headerLine << endl;

		// Result details
		int numDofs = resultMatrix[0].size();

		// Setting up the result type header line
		string my_type_header = "ComponentNames \"DOF_1\",";
		if (numDofs > 1) {
			for (int i = 2; i < numDofs + 1; i++) {
				stringstream i_ss;
				i_ss << i;
				my_type_header += " \"DOF_" + i_ss.str() + "\",";
			}
		}

		// Setting up the Result line
		femResFile << "Result " << "\"" << resultName << "\" \""
				<< analysis_name << "\" " << stepValue << " " << resultType
				<< " " << my_location << endl;

		femResFile << my_type_header << endl;

		femResFile << "Values" << endl;

		int nodeNum = 1;
		for (int r = 0; r < resultMatrix.size(); r++) {

			femResFile << nodeNum;

			for(int s = 0; s < numDofs; s++){
				femResFile << " " << resultMatrix[r][s];
			}

			femResFile << endl;
			nodeNum++;

		}

		femResFile << "End Values" << endl;
	}
	femResFile.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
// Outputting calculated matrix to fem.res file
void CardiacFibreGenerator::saveResultsToFile(vector<Particle>& ptcls,
		InputFileData* InputData, ofstream& logFile){

//	int resultType = 0;
//	std::string resultName = "averaged_nodal_normal_n";
	// *********************************************************************
//	int resultType = 1;
//	std::string resultName = "circumferential_direction_c";

//	int resultType = 2;
//	std::string resultName = "outward_pointing_normal_n_cz";

	int resultType = 3;
	std::string resultName = "sheet_normal_s";

//	int resultType = 4;
//	std::string resultName = "fibre direction projection_p";

//	int resultType = 5;
//	std::string resultName = "fibre_direction_f";


	dbMatrix resultMatrix(ptcls.size(),dbVector());
	dbVector dummyVec(3,0);
	for(int i = 0; i < ptcls.size(); i++){

		if(ptcls[i].getAllSurfaceNormals().size()>0){
			dbVector& resultVec = ptcls[i].getAllSurfaceNormals()[resultType];
			resultMatrix[i] = resultVec;
		}
		else
			resultMatrix[i] = dummyVec;
	}

	string filename = "fem_fibre.res";
	ofstream femResFile(filename.c_str());

	double stepValue = 0.001;

	saveResultsToFile_res_format(femResFile,resultName,resultMatrix,stepValue,
			InputData, logFile);

	femResFile.close();

//	cout << "Result saved: " << resultName << endl;
//	logFile << "Result saved: " << resultName << endl;

}

/*!****************************************************************************/
/*!****************************************************************************/
// Outputting calculated matrix to fem.res file
void CardiacFibreGenerator::saveResultsToFile_res_format_flexible(std::string& fileName,
		std::vector<std::string>& resultName,
		std::vector<dbMatrix>& resultMatrix, dbVector& stepValue,
		InputFileData* InputData, ofstream& logFile) {

	ofstream femResFile(fileName.c_str());

	string headerLine = "GiD Post Results File 1.0";

	string resultType = "Vector";
	string my_location = "OnNodes";
	string analysis_name = " ";

	if (femResFile.is_open()) {

		femResFile << headerLine << endl;

		for (int i = 0; i < resultName.size(); i++) {

			// Setting up the Result line
			femResFile << "Result " << "\"" << resultName[i] << "\" \""
					<< analysis_name << "\" " << stepValue[i] << " "
					<< resultType << " " << my_location << endl;


			// Setting up the result type header line
			int numDofs = resultMatrix[i][0].size();
			string my_type_header = "ComponentNames \"DOF_1\",";
			if (numDofs > 1) {
				for (int i = 2; i < numDofs + 1; i++) {
					stringstream i_ss;
					i_ss << i;
					my_type_header += " \"DOF_" + i_ss.str() + "\",";
				}
			}
			femResFile << my_type_header << endl;

			// Outputting results
			femResFile << "Values" << endl;
			int nodeNum = 1;
			for (int r = 0; r < resultMatrix[i].size(); r++) {

				femResFile << nodeNum;

				for (int s = 0; s < numDofs; s++) {
					femResFile << " " << resultMatrix[i][r][s];
				}

				femResFile << endl;
				nodeNum++;

			}
			femResFile << "End Values" << endl;
		}
	}
	femResFile.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
// Outputting calculated matrix to fem.res file
void CardiacFibreGenerator::saveResultsToFile_res_format_flexible_resultTypes(std::string& fileName,
		std::vector<std::string>& vectorResultName,
		std::vector<dbMatrix>& vectorResultMatrix,
		std::vector<std::string>& scalarResultName,
		std::vector<dbVector>& scalarResultVector, dbVector& stepValue,
		InputFileData* InputData, ofstream& logFile) {

	ofstream femResFile(fileName.c_str());

	string headerLine = "GiD Post Results File 1.0";

	string my_location = "OnNodes";
	string analysis_name = " ";

	if (femResFile.is_open()) {

		femResFile << headerLine << endl;

		string resultType = "Vector";
		for (int i = 0; i < vectorResultName.size(); i++) {

			// Setting up the Result line
			femResFile << "Result " << "\"" << vectorResultName[i] << "\" \""
					<< analysis_name << "\" " << stepValue[i] << " "
					<< resultType << " " << my_location << endl;

			// Setting up the result type header line
			int numDofs = vectorResultMatrix[i][0].size();
			string my_type_header = "ComponentNames \"DOF_1\",";
			if (numDofs > 1) {
				for (int i = 2; i < numDofs + 1; i++) {
					stringstream i_ss;
					i_ss << i;
					my_type_header += " \"DOF_" + i_ss.str() + "\",";
				}
			}
			femResFile << my_type_header << endl;

			// Outputting results
			femResFile << "Values" << endl;
			int nodeNum = 1;
			for (int r = 0; r < vectorResultMatrix[i].size(); r++) {

				femResFile << nodeNum;

				for (int s = 0; s < numDofs; s++) {
					femResFile << " " << vectorResultMatrix[i][r][s];
				}

				femResFile << endl;
				nodeNum++;

			}
			femResFile << "End Values" << endl;
		}

		// --------------------------------------------------------------------
		// Output scalar values
		resultType = "Scalar";
		for (int i = 0; i < scalarResultName.size(); i++) {

			// Setting up the Result line
			femResFile << "Result " << "\"" << scalarResultName[i] << "\" \""
					<< analysis_name << "\" " << stepValue[i] << " "
					<< resultType << " " << my_location << endl;

			// Setting up the result type header line
			string my_type_header = "ComponentNames \"DOF_1\",";
			femResFile << my_type_header << endl;

			// Outputting results
			femResFile << "Values" << endl;
			int nodeNum = 1;
			for (int r = 0; r < scalarResultVector[i].size(); r++) {

				femResFile << nodeNum << " " << scalarResultVector[i][r] << endl;

				nodeNum++;

			}
			femResFile << "End Values" << endl;
		}

	}

	femResFile.close();
}

/*!****************************************************************************/
/*!****************************************************************************/
// Outputting calculated matrix to fem.res file
void CardiacFibreGenerator::saveResultsToFile_flexible(vector<Particle>& ptcls,
		InputFileData* InputData, ofstream& logFile) {

	intVector resultsToSave(7);
	resultsToSave[0] = 0;
	resultsToSave[1] = 1;
	resultsToSave[2] = 2;
	resultsToSave[3] = 3;
	resultsToSave[4] = 4;
	resultsToSave[5] = 5;
	resultsToSave[6] = 6;

	map<int, string> resultInfo;
	resultInfo[0] = "averaged_nodal_normal_n";
	resultInfo[1] = "circumferential_direction_c";
	resultInfo[2] = "outward_pointing_normal_n_cz";
	resultInfo[3] = "sheet_normal_s";
	resultInfo[4] = "fibre_direction_projection_p";
	resultInfo[5] = "fibre_direction_f";
	resultInfo[6] = "fibre_direction_m";

//	vector<string> resultName(resultsToSave.size());
//	vector<dbMatrix> resultMatrix(resultsToSave.size(),dbMatrix()) ;
//	vector<double> stepValue(resultsToSave.size(),0.001);
//	for(int i = 0; i < resultsToSave.size(); i++){
//		resultName[i] = resultInfo.find(resultsToSave[i])->second;
//
//		dbMatrix resultMat(ptcls.size(),dbVector());
//		dbVector dummyVec(3,0);
//		for(int j = 0; j < ptcls.size(); j++){
//			if(ptcls[j].getAllSurfaceNormals().size()>0){
//				dbVector& resultVec = ptcls[j].getAllSurfaceNormals()[resultsToSave[i]];
//				resultMat[j] = resultVec;
//			}
//			else
//				resultMat[j] = dummyVec;
//		}
//
//		resultMatrix[i] = resultMat;
//	}

//	string fileName = "fem_fibre.res";
//	saveResultsToFile_res_format_flexible(fileName,resultName,resultMatrix,stepValue,
//			InputData, logFile);

//	for (int i = 0; i < resultName.size(); i++) {
//		cout << "Result saved: " << resultName[i] << endl;
//		logFile << "Result saved: " << resultName[i] << endl;
//	}

	// -------------------------------------------------------------------------

	vector<string> vectorResultName(resultsToSave.size());
	vector<dbMatrix> vectorResultMatrix(resultsToSave.size(), dbMatrix());
	vector<double> stepValue(resultsToSave.size(), 0.001);
	for (int i = 0; i < resultsToSave.size(); i++) {
		vectorResultName[i] = resultInfo.find(resultsToSave[i])->second;

		dbMatrix resultMat(ptcls.size(), dbVector());
		dbVector dummyVec(3, 0);
		for (int j = 0; j < ptcls.size(); j++) {
			if (ptcls[j].getAllSurfaceNormals().size() > 0) {
				dbVector& resultVec =
						ptcls[j].getAllSurfaceNormals()[resultsToSave[i]];
				resultMat[j] = resultVec;
			} else
				resultMat[j] = dummyVec;
		}

		vectorResultMatrix[i] = resultMat;
	}

	string scalarResultToSave = "fibre_angle";
	vector<string> scalarResultName(1);
	scalarResultName[0] = scalarResultToSave;

	vector<dbVector> scalarResultMatrix(scalarResultName.size());
	for (int i = 0; i < scalarResultMatrix.size(); i++) {
		for (int j = 0; j < ptcls.size(); j++) {
			scalarResultMatrix[i].push_back(ptcls[j].getBeta());
		}
	}

	string fileName = "fem_fibre.res";
	saveResultsToFile_res_format_flexible_resultTypes(fileName,
			vectorResultName, vectorResultMatrix, scalarResultName,
			scalarResultMatrix, stepValue, InputData, logFile);

//	for(int i = 0; i < vectorResultName.size(); i++){
//		cout << "Result saved: " << vectorResultName[i] << endl;
//		logFile << "Result saved: " << vectorResultName[i] << endl;
//	}

	// -------------------------------------------------------------------------

}

/*!****************************************************************************/
/*!****************************************************************************/
//	enforce epicardium properties on boundary nodes (boundary between
//	epicardium and endocardium mesh)
void CardiacFibreGenerator::eliminateSpecificSurface(intVector& surfElems, int surfID,
		vector<FEMElement>& elemVec,InputFileData* InputData,
		std::ofstream& logFile){

		bool notSame = false;
		for(int j=1; j < surfElems.size(); j++){
//			logFile << "elem(" << j <<"):"
//					<< elemVec[surfElems[j-1]].getElemType() << "!=(?)"
//					<< elemVec[surfElems[j]].getElemType() << endl;
			if(elemVec[surfElems[j-1]].getElemType() !=
										elemVec[surfElems[j]].getElemType() ){
//				logFile << "notSame = true" << endl;
				notSame = true;
			}
		}
		if(notSame == true){
//			printVector(surfElems,"surfElems(before)",logFile);
			intVector dummyVec;
			for(int k = 0; k < surfElems.size(); k++){
				if(elemVec[surfElems[k]].getElemType() == surfID){
					dummyVec.push_back(surfElems[k]);
				}
			}
			surfElems = dummyVec;
//			printVector(surfElems,"surfElems(after)",logFile);
		}
}
