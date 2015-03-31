//============================================================================
// Name        : Orion.cpp
// Author      : Ritesh Rao Rama
// Version     :
// Copyright   : 
// Description : POD-based real-time calculation
//============================================================================
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include <time.h>
#include <chrono>
#include <ctime>

//#include "clapack.h"
#define PI 3.14159265

//#include "commonTypedefs.h"
//#include "commonFunctions.h"
//#include "fortranFunctions.h"

#include <limits>
#include "InputFileData.h"

//#include "Database.h"

//#include "PODCalc.h"
//#include "FEMGeometryX.h"
//#include "ParticleX.h"
//#include "FEMElementX.h"

#include "defs.h"
#include "ROMCalc.h"
#include "ErrorCalc.h"
#include "Interpolation.h"
#include "testFunctions.h"

#include "Data.h"

using namespace std;

void readAdditionalFiles(InputFileData* InputData, ofstream& logFile);
void readInputFile(string filename, InputFileData* InputData, ofstream& logFile);
void readMaterialFile(string filename, InputFileData* InputData,
		ofstream& logFile);

int main(int argc, char* argv[]) {

	PetscInitialize(&argc, &argv, (char*) 0, 0);

	// *************************************************************************
	// Log writer
	string fileOutput = "oLog.dat";
	ofstream logFile(fileOutput.c_str(), std::ofstream::out | std::ofstream::app);
	logFile.precision(15);

	time_t tim;  //create variable of time_t
	time(&tim); //pass variable tim to time function

	cout << "================================" << endl;
	cout << "ORION: " << ctime(&tim);
	cout << "--------------------------------" << endl;

	logFile << "================================" << endl;
	logFile << "ORION: " << ctime(&tim);
	logFile << "--------------------------------" << endl;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	//*************************************************************************
	// Read input.dat file
	string filename = "input.dat";
	logFile << "Reading: " << filename << endl;
	cout << "Reading: " << filename << endl;
	InputFileData* InputData = new InputFileData(filename,logFile);

	readAdditionalFiles(InputData,logFile);

	int choice = InputData->getValue("orionMode");

	switch (choice) {
	case 0: {
		//*************************************************************************
		// Reduced Order Method Calculation Calculation
		ROMCalc* reducedOrderCalculation = new ROMCalc(InputData, logFile);
		delete reducedOrderCalculation;
		break;
	}

	case 1: {
		//*************************************************************************
		// Error Calculation
		string fileOutput2 = "oLogError.dat";
		ofstream elogFile(fileOutput2.c_str(), std::ofstream::out);
		logFile.precision(15);

		ErrorCalc myError(InputData, elogFile);
		break;
	}

	case 2: {
		//*************************************************************************
		// Testing specific algorithm in Orion
		testFunctions(InputData, logFile);
		break;
	}

	default:
		cout << "ERROR: In main(), choice is not valid" << endl;
		logFile << "ERROR: In main(), choice is not valid" << endl;
	}

	delete InputData;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	cout << "Total calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;
	logFile << "Total calculation time: " << elapsed_seconds.count() << " sec"
			<< endl;

	cout << "--------------------------------" << endl;
	logFile << "--------------------------------" << endl;

	return 0;

}

// *****************************************************************************
// read additional files
void readAdditionalFiles(InputFileData* InputData, ofstream& logFile){

	string inputFile = "oInput.dat";
	readInputFile(inputFile,InputData,logFile);

	string materialFile = "oMaterial.dat";
	readMaterialFile(materialFile,InputData,logFile);

}


// *****************************************************************************
// read input file and add/overwrite parameters in InputData
void readInputFile(string filename, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "###################################################" << endl;
	logFile << "############ " << filename << " ###################" << endl;

	ifstream inputFile(filename.c_str());

	if (inputFile.eof()) {
		logFile << "Input file: " << filename << " contains no data!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (inputFile) {

		logFile << "Reading: " << filename << endl;
		cout << "Reading: " << filename << endl;

		string name, line;
		double value;

		while (inputFile.good()) {

			inputFile >> name >> value;

			InputData->setValue(name.c_str(), value);

			logFile << name << ": " << value << endl;
		}

	}

	inputFile.close();

}

// *****************************************************************************
// read material file and overwrite parameters in InputData->getMaterials()
void readMaterialFile(string filename, InputFileData* InputData,
		ofstream& logFile) {

	logFile << "###################################################" << endl;
	logFile << "############ " << filename << " ###################" << endl;

	ifstream materialFile(filename.c_str());

	if (materialFile.eof()) {
		logFile << "Material file " << filename << " contains no data!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (materialFile) {

		logFile << "Reading: " << filename << endl;
		cout << "Reading: " << filename << endl;

		string name, line;
		double value;

		vector<string> paramNamesVec;
		dbVector paramValuesVec;

		while (materialFile.good()) {

			materialFile >> name >> value;

			paramNamesVec.push_back(name);
			paramValuesVec.push_back(value);

			logFile << paramNamesVec[paramNamesVec.size() - 1] << ": "
					<< paramValuesVec[paramValuesVec.size() - 1] << endl;
		}

		materialFile.close();

		//! Extract materials
		vector<map<string, double> >& materialsSet = InputData->getMaterials();
		if (materialsSet.size() > 1) {
			logFile
					<< "Geometry with multiple material types is not yet supported"
					<< endl;
			cout << "Geometry with multiple material types is not yet supported"
					<< endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! Extract the material's parameters
		map<string, double>& material = materialsSet[0];

		//! Extract the parameters' value
		vector<string> missingParameters;
		for (int i = 0; i < paramNamesVec.size(); i++) {

			//! Search for specific parameter in the material's parameters list
			map<string, double>::iterator it_material = material.find(
					paramNamesVec[i]);

			if (it_material == material.end()) {
				// If specific parameter is not found, record it in
				// missing parameters vector
				missingParameters.resize(missingParameters.size() + 1,
						paramNamesVec[i]);
			} else {
				// If specific parameter is found, record it in parameters vector
				it_material->second = paramValuesVec[i];
			}
		}

		if (missingParameters.size() > 0) {
			logFile << "ERROR: The following parameters are missing from the "
					"input file: " << endl;
			cout << "ERROR: The following parameters are missing from the "
					"input files: " << endl;

			for (int i = 0; i < missingParameters.size(); i++) {
				logFile << missingParameters[i] << " ";
				cout << missingParameters[i] << " ";
			}
			logFile << endl;
			cout << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		//! Check
		logFile << "******************************************************"
				<< endl;
		logFile << "New material parameters" << endl;
		vector<map<string, double> > materialsSetCheck =
				InputData->getMaterials();
		map<string, double> materialCheck = materialsSetCheck[0];
		map<string, double>::iterator it_materialCheck_beg =
				materialCheck.begin();
		map<string, double>::iterator it_materialCheck_end =
				materialCheck.end();
		map<string, double>::iterator it_materialCheck = materialCheck.begin();
		for (it_materialCheck = it_materialCheck_beg;
				it_materialCheck != it_materialCheck_end; it_materialCheck++) {
			logFile << it_materialCheck->first << ": "
					<< it_materialCheck->second << endl;
		}
	}

}
