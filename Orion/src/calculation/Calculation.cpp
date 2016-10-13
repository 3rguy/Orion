/*
 * Calculation.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <Calculation.h>

Calculation::Calculation(InputFileData* InputData, ofstream& logFile) {

	DataContainer* problemData = new DataContainer();

	// Model for reduced order calculation
	int calcMode = InputData->getValue("calculationMode");

	switch(calcMode){
	case 1:
	{
		ModelSet* ModelCalc = new ModelSet(problemData,InputData,logFile);
		delete ModelCalc;
		break;
	}

	case 2:
		// error calculation
		logFile << "ERROR: Error class has not been integrated with the Calculation class yet." << endl;
		cout << "ERROR: Error class has not been integrated with the Calculation class yet" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;

	default:
		logFile << "ERROR: Incorrect model selected. Valid options: 1, 2." << endl;
		cout << "ERROR: Incorrect model selected. Valid options: 1, 2." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

