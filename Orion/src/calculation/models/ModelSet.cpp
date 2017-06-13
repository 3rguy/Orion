/*
 * Model.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <models/ModelSet.h>

ModelSet::ModelSet(DataContainer* problemData, InputFileData* InputData,
		ofstream& logFile) {

	string databaseName = "myDatabase.csv";
	Database myDatabase(databaseName, logFile);

	int modelType = InputData->getValue("modelType");

	switch (modelType) {

	case 1: {
		Standard* StandardCalc = new Standard(myDatabase, problemData,
				InputData, logFile);

		StandardCalc->preProcessing(problemData, InputData, logFile);

		StandardCalc->reducedCalculation(problemData, InputData, logFile);

		StandardCalc->postProcessing(problemData, InputData, logFile);

		delete StandardCalc;

		break;
	}

	case 2: {
		Cardiac* CardiacCalc = new Cardiac(myDatabase, problemData, InputData,
				logFile);

		CardiacCalc->preProcessing(problemData, InputData, logFile);

		CardiacCalc->reducedCalculation(problemData, InputData, logFile);

		CardiacCalc->postProcessing(problemData, InputData, logFile);

		delete CardiacCalc;

		break;
	}

	case 3: {
		Plasticity* PlasticityCalc = new Plasticity(myDatabase, problemData,
				InputData, logFile);

		PlasticityCalc->preProcessing(problemData, InputData, logFile);

		PlasticityCalc->reducedCalculation(problemData, InputData, logFile);

		PlasticityCalc->postProcessing(problemData, InputData, logFile);

		delete PlasticityCalc;

		break;
	}

	case 4: {
		MicroRVE* MicroRVECalc = new MicroRVE(myDatabase, problemData,
				InputData, logFile);

		MicroRVECalc->preProcessing(problemData, InputData, logFile);

		MicroRVECalc->reducedCalculation(problemData, InputData, logFile);

		MicroRVECalc->postProcessing(problemData, InputData, logFile);

		delete MicroRVECalc;

		break;
	}

	default:
		cout << "ERROR: In ModelSet(), modelType choice is invalid." << endl;
		logFile << "ERROR: In ModelSet(), modelType choice is invalid." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		break;
	}

}

