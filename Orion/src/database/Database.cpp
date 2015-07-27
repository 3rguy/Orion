#include "Database.h"

Database::Database(string& databaseName, ofstream& logFile) {

	cout << "Loading-up database" << endl;
	logFile << "Loading-up database" << endl;

	//! Read database file (e.g myDatabase.csv)
	dbMatrix paramMatrix;
	vector<string> fileNameVec;
	intVector anchorPointList;
	readDatabase(databaseName, paramNamesVec, paramMatrix, anchorPointList, fileNameVec,
			logFile);

	//! Create Data and store them
	for (int i = 0; i < paramMatrix.size(); i++) {
		Data a;

		// Define Data
		a.setId(i);
		a.setParamValuesVec(paramMatrix[i]);
		a.setFolderName(fileNameVec[i]);
		a.setAnchorPoint(anchorPointList[i]);

		// Add Data to list
		dataList.resize(dataList.size() + 1);
		dataList[dataList.size() - 1] = a;
	}

#ifdef _DatabaseDebugMode_
	printDatabase(logFile);
#endif
}

void Database::addEntry(Data entry) {
	cout << "Database::addEntry" << endl;

	entry.setId((dataList[dataList.size() - 1].getId()) + 1);
	dataList.resize(dataList.size() + 1);
	dataList[dataList.size() - 1] = entry;
}

/*!***************************************************************************/
/*!***************************************************************************/
//! Delete a specific entry in the database
void Database::deleteEntry(int id) {
	cout << "Database::deleteEntry" << endl;

	for (int i = 0; i < dataList.size(); i++) {
		if (dataList[i].getId() == id) {
			dataList.erase(dataList.begin() + i);
			break;
		}
	}
}

/*!***************************************************************************/
/*!***************************************************************************/
//! Read database file
void Database::readDatabase(string& databaseName, vector<string>& paramNameVec,
		dbMatrix& paramMatrix, intVector& anchorPointList,
		vector<string> &fileNameVec, ofstream& logFile) {

	ifstream myfile(databaseName.c_str());

	if (myfile.is_open()) {
		// Read header line
		string line;
		getline(myfile, line);
		istringstream ss(line);
		while (ss) {
			string s;
			if (!getline(ss, s, ','))
				break;
			paramNameVec.push_back(s);
		}

		// Read header vector to hold parameters' name but not
		// filenames and anchorPoints
		paramNameVec.pop_back();
		paramNameVec.pop_back();

#ifdef _DatabaseDebugMode_
		logFile << "Headers[" << paramNameVec.size() << "]: ";
		for (int i = 0; i < paramNameVec.size(); i++)
			logFile << paramNameVec[i] << " ";
		logFile << endl;
#endif

		double nParameters = paramNameVec.size();

		while (myfile.good()) {

			if (!getline(myfile, line))
				break;

			if(line[0] == '#')
				continue;

			resizeArray(paramMatrix, paramMatrix.size() + 1, nParameters);
			fileNameVec.resize(fileNameVec.size() + 1);
			anchorPointList.resize(anchorPointList.size() + 1);

			istringstream ss(line);
			int counter = 0;
			while (ss) {
				string s;
				if (!getline(ss, s, ',')) // Break line using a specific delimiter
					break;

				if (counter < nParameters) {
					paramMatrix[paramMatrix.size() - 1][counter] = atof(
							s.c_str());
					counter++;
				}
				else if(counter == nParameters){	// Record anchor points
					anchorPointList[anchorPointList.size() - 1] = atof(s.c_str());
					counter++;
				}
				else if(counter == nParameters+1){	// Record file name
					fileNameVec[fileNameVec.size() - 1] = s;
				}
			}
		}
//#ifdef _DatabaseDebugMode_
		logFile << "***** Database[ "<< databaseName <<"] *****" << endl;
		for (int i = 0; i < paramNameVec.size(); i++){
			logFile << paramNameVec[i] << "\t";
		}
		logFile << "Filename \t anchorPoint";
		logFile << endl;
		for(int j=0; j <paramMatrix.size(); j++){
			for(int k=0; k <paramMatrix[j].size(); k++){
				logFile << paramMatrix[j][k] << std::setw(6);
			}
			logFile << anchorPointList[j] << std::setw(6);
			logFile << fileNameVec[j] << endl;
		}
//#endif

		myfile.close();
	} else {
		logFile << "ERROR: Unable to open database file" << endl;
		cout << "ERROR: Unable to open database file" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!***************************************************************************/
/*!***************************************************************************/
//! Print the content of the database
void Database::printDatabase(ofstream& logFile) {

	logFile << "********** Printing database **********" << endl;
	logFile << "---------------------------------------" << endl;
	logFile << "Size of Database: " << dataList.size() << endl;
	for (int i = 0; i < dataList.size(); i++) {

		dbVector pValues = dataList[i].getParamValuesVec();

		logFile << "Id:" << dataList[i].getId() << "\t";
		for (int j = 0; j < pValues.size(); j++) {
			logFile << paramNamesVec[j] << ": " << pValues[j] << "\t";
		}
		logFile << "File Name: " << dataList[i].getFolderName() << endl;
	}
	logFile << endl;
}

/*!***************************************************************************/
/*!***************************************************************************/
int Database::findData(int val) {

	int location = -1;

	for (int i = 0; i < dataList.size(); i++) {
		if (dataList[i].getId() == val) {
			location = i;
			break;
		}
	}

	return location;
}

/*!***************************************************************************/
/*!***************************************************************************/
//! Import the displacement for a specific Data
void Database::importDisp(InputFileData* InputData,int val, ofstream& logFile) {

	int location = findData(val);
	//double start_time = clock();
	dataList[location].readResultFile(InputData,logFile);
	//double end_time = clock();
	//logFile << "Total time taken to load file:" << end_time - start_time
	//		  << "ms" << endl;

}

/*!***************************************************************************/
/*!***************************************************************************/
//! Print the parameter info of a specific Data
void Database::printDataInfo(int val, ofstream& logFile) {
	logFile << "" << endl;

	int location = findData(val);

	dbVector pValues = dataList[location].getParamValuesVec();

	logFile << "Id:" << dataList[location].getId() << "\t";
	for (int j = 0; j < pValues.size(); j++) {
		logFile << paramNamesVec[j] << ": " << pValues[j] << "\t";
	}
	logFile << "File Name: " << dataList[location].getFolderName() << endl;

//	dataList[location].printDispMatrix(logFile);
}
