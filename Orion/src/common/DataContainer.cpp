#include "DataContainer.h"

DataContainer::DataContainer(){

}

/*!
 * Constructor of the struct 'DataContainer'
 * @param logFile logFile output
 */
DataContainer::DataContainer(std::ofstream& logFile){

	using namespace std;

	logFile << "DataContainer::DataContainer" << endl;
}

/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves the name and integer vector in the intVectorList
 * @param name	name of Vector
 * @param vec	vector of integers
 * @param logFile	logFile output
 */
void DataContainer::setValue(const char* name, intVector vec) {

	string str(name);

	if (findIntVectorList(name) == false)
		intVectorList[str] = vec;
	else {
		cout << "'" << name << "' already exists in DataContainer intVectorList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!
 * A member function to find a specific integer vector in the intVectorList
 * @param name	name of the vector to be looked for
 * @param logFile	logFile output
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findIntVectorList(const char* name) {

	string str(name);
	bool isFound;
	map<string, intVector>::iterator it = intVectorList.find(str);

	if (it == intVectorList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;
}

intVector& DataContainer::getIntVector(const char* name) {

		string str(name);
		bool isFound;
		map<string, intVector>::iterator it = intVectorList.find(str);

		if (it == intVectorList.end()){
			cout << "Cannot find '" << name << "' in intVectorList" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		else
			return it->second;

}


/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves save the name and double vector in the dbVectorList
 * @param name	name of Vector
 * @param vec	vector of doubles
 * @param logFile	logFile output
 */
void DataContainer::setValue(const char* name,dbVector vec){

	string str(name);

	if (findIntVectorList(name) == false)
		dbVectorList[str] = vec;
	else {
		cout << "'" << name << "' already exists in DataContainer dbVectorList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!
 * A member function to find a specific vector of doubles in the dbVectorList
 * @param name	name of the vector to be looked for
 * @param logFile	logFile output
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findDbVectorList(const char* name){

	string str(name);
	bool isFound;

	map<string,dbVector>::iterator it = dbVectorList.find(str);

	if(it == dbVectorList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;
}


dbVector& DataContainer::getDbVector(const char* name) {

		string str(name);
		bool isFound;
		map<string, dbVector>::iterator it = dbVectorList.find(str);

		if (it == dbVectorList.end()){
			cout << "Cannot find '" << name << "' in dbVectorList" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		else
			return it->second;

}

/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves save the name and matrix of integers in the intMatrixList
 * @param name	name of matrix
 * @param mat	matrix of integers
 * @param logFile	logFile output
 */
void DataContainer::setValue(const char* name,intMatrix mat){

	string str(name);

	if (findIntVectorList(name) == false)
		intMatrixList[str] = mat;
	else {
		cout << "'" << name << "' already exists in DataContainer intMatrix"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!
 * A member function to find a specific matrix of integers in the intMatrixList
 * @param name	name of the vector to be looked for
 * @param logFile	logFile output
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findIntMatrixList(const char* name){

	string str(name);
	bool isFound;

	map<string,intMatrix>::iterator it = intMatrixList.find(str);

	if(it == intMatrixList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;

}

intMatrix& DataContainer::getIntMatrix(const char* name) {

		string str(name);
		bool isFound;
		map<string, intMatrix>::iterator it = intMatrixList.find(str);

		if (it == intMatrixList.end()){
			cout << "Cannot find '" << name << "' in intMatrixList" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		else
			return it->second;

}


/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves save the name and matrix of doubles in the dbMatrixList
 * @param name name of matrix
 * @param mat	matrix of doubles
 * @param logFile logFile output
 */
void DataContainer::setValue(const char* name,dbMatrix mat){

	string str(name);
	if (findIntVectorList(name) == false)
		dbMatrixList[str] = mat;
	else {
		cout << "'" << name << "' already exists in DataContainer dbMatrixList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!
 * A member function to find a specific matrix of doubles in the dbMatrixList
 * @param name	name of the vector to be looked for
 * @param logFile	logFile output
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findDbMatrixList(const char* name){

	string str(name);
	bool isFound;

	map<string,dbMatrix>::iterator it = dbMatrixList.find(str);

	if(it == dbMatrixList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;

}

dbMatrix& DataContainer::getDbMatrix(const char* name) {

		string str(name);
		bool isFound;
		map<string, dbMatrix>::iterator it = dbMatrixList.find(str);

		if (it == dbMatrixList.end()){
			cout << "Cannot find '" << name << "' in dbMatrixList" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		else
			return it->second;

}


/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves save the name and double vector in the dbMatrixVecList
 * @param name name of vector list
 * @param matVec vector list with Matrix of doubles
 * @param logFile logFile output
 */
void DataContainer::setValue(const char* name,vector<dbMatrix> matVec){

	string str(name);
	if (findDbMatrixVecList(name) == false)
		dbMatrixVecList[str] = matVec;
	else {
		cout << "'" << name << "' already exists in DataContainer dbMatrixVecList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

/*!
 * A member function to find a specific vector of matrix double in the dbMatrixVecList
 * @param name	name of the vector to be looked for
 * @param logFile	logFile output
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findDbMatrixVecList(const char* name){

	string str(name);
	bool isFound;

	map<string,vector<dbMatrix> >::iterator it = dbMatrixVecList.find(str);

	if(it == dbMatrixVecList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;

}


vector<dbMatrix>& DataContainer::getDbMatrixVec(const char* name) {

		string str(name);
		bool isFound;
		map<string, vector<dbMatrix> >::iterator it = dbMatrixVecList.find(str);

		if (it == dbMatrixVecList.end()){
			cout << "Cannot find '" << name << "' in dbMatrixVecList" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		else
			return it->second;

}

