/*
 * Standard.cpp
 *
 *  Created on: Jul 28, 2015
 *      Author: rama
 */

#include <models/Standard.h>

/*!****************************************************************************/
/*!****************************************************************************/
Standard::Standard(Database& myDB, DataContainer* problemData,
			 InputFileData* InputData, ofstream& logFile){

	cout << "Model Calculation: Standard" << endl;
	logFile << "Model Calculation: Standard" << endl;

	myDatabase = myDB;

}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::preProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	initPreProcessing(myDatabase,problemData,InputData,logFile);

	defaultPreProcessing(myDatabase,problemData,InputData,logFile);
}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::reducedCalculation(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	calculation(myDatabase,problemData,InputData,logFile);

}

/*!****************************************************************************/
/*!****************************************************************************/
void Standard::postProcessing(DataContainer* problemData,
		InputFileData* InputData, ofstream& logFile){

	defaultPostProcessing(myDatabase,problemData,InputData,logFile);

}

