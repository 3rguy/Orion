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
 * A member function that saves the name and an integer in the intList
 * @param name	name of integer
 * @param val	integer
 */
void DataContainer::setValue(const char* name, int val) {

	string str(name);

	if (findInt(name) == false)
		intList[str] = val;
	else {
		cout << "'" << name << "' already exists in DataContainer intList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!
 * A member function to find a specific integer in the intList
 * @param name	name of the integer to be looked for
 * @param logFile	logFile output
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findInt(const char* name){

	string str(name);
	bool isFound;
	map<string, int>::iterator it = intList.find(str);

	if (it == intList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;

}

/*!
 * A member function to find a specific integer in the intList
 * @param name	name of the integer to be looked for
 * @return	return 'True' if found and 'false' if not.
 */
int& DataContainer::getInt(const char* name){

	string str(name);
	bool isFound;
	map<string, int>::iterator it = intList.find(str);

	if (it == intList.end()){
		cout << "Cannot find '" << name << "' in intList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		return it->second;


}

/*!
 * A member function to delete a specific integer in the intList
 * @param name	name of the integer to be looked for
 */
void DataContainer::deleteInt(const char* name){

	string str(name);
	bool isFound;
	map<string, int>::iterator it = intList.find(str);

	if (it == intList.end()){
		cout << "In deleteInt, cannot find '" << name << "' in intList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		intList.erase(it);

}

/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves the name and integer vector in the doubleList
 * @param name	name of double
 * @param val	double
 */
void DataContainer::setValue(const char* name,double val){

	string str(name);

	if (findDouble(name) == false)
		doubleList[str] = val;
	else {
		cout << "'" << name << "' already exists in DataContainer doubleList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}
/*!
 * A member function that saves the name and a double in the doubleList
 * @param name	name of double
 * @return bool	return 'True' if found and 'false' if not.
 */
bool DataContainer::findDouble(const char* name){

	string str(name);
	bool isFound;
	map<string, double>::iterator it = doubleList.find(str);

	if (it == doubleList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;
}
/*!
 * A member function that saves the name and a double in the doubleList
 * @param name	name of double
 * @param vec	vector of integers
 * @return logFile	return 'True' if found and 'false' if not.
 */
double& DataContainer::getDouble(const char* name){

	string str(name);
	bool isFound;
	map<string, double>::iterator it = doubleList.find(str);

	if (it == doubleList.end()){
		cout << "Cannot find '" << name << "' in doubleList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		return it->second;
}

/*!
 * A member function to delete a specific double in the doubleList
 * @param name	name of double to be looked for
 */
void DataContainer::deleteDouble(const char* name){

	string str(name);
	bool isFound;
	map<string, double>::iterator it = doubleList.find(str);

	if (it == doubleList.end()){
		cout << "In deleteDouble, cannot find '" << name << "' in doubleList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		doubleList.erase(it);

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

	if (findIntVector(name) == false)
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
bool DataContainer::findIntVector(const char* name) {

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

/*!
 * A member function to delete a specific integer in the IntVector
 * @param name	name of the vector of integers to be looked for
 */
void DataContainer::deleteIntVector(const char* name){

	string str(name);
	bool isFound;
	map<string, intVector>::iterator it = intVectorList.find(str);

	if (it == intVectorList.end()){
		cout << "In deleteIntVector, cannot find '" << name << "' in intVectorList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		intVectorList.erase(it);

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

	if (findIntVector(name) == false)
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
bool DataContainer::findDbVector(const char* name){

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

/*!
 * A member function to delete a specific double vector in the dbVectorList
 * @param name	name of the vector of double to be looked for
 */
void DataContainer::deleteDbVector(const char* name){

	string str(name);
	bool isFound;
	map<string, dbVector>::iterator it = dbVectorList.find(str);

	if (it == dbVectorList.end()){
		cout << "In deleteIntVector, cannot find '" << name << "' in dbVectorList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		dbVectorList.erase(it);

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

	if (findIntVector(name) == false)
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
bool DataContainer::findIntMatrix(const char* name){

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

/*!
 * A member function to delete a specific integer matrix in the IntMatrixList
 * @param name	name of the vector of double to be looked for
 */
void DataContainer::deleteIntMatrix(const char* name){

	string str(name);
	bool isFound;
	map<string, intMatrix>::iterator it = intMatrixList.find(str);

	if (it == intMatrixList.end()){
		cout << "In deleteIntMatrix, cannot find '" << name << "' in intMatrixList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		intMatrixList.erase(it);

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
	if (findIntVector(name) == false)
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
bool DataContainer::findDbMatrix(const char* name){

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

/*!
 * A member function to delete a specific double matrix in the dbMatrixList
 * @param name	name of the vector of double to be looked for
 */
void DataContainer::deleteDbMatrix(const char* name){

	string str(name);
	bool isFound;
	map<string, dbMatrix>::iterator it = dbMatrixList.find(str);

	if (it == dbMatrixList.end()){
		cout << "In deleteDbMatrix, cannot find '" << name << "' in dbMatrixList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		dbMatrixList.erase(it);

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
	if (findDbMatrixVec(name) == false)
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
bool DataContainer::findDbMatrixVec(const char* name){

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

/*!
 * A member function to delete a specific double matrix in the dbMatrixList
 * @param name	name of the vector of double to be looked for
 */
void DataContainer::deleteDbMatrixVec(const char* name){

	string str(name);
	bool isFound;
	map<string, vector<dbMatrix> >::iterator it = dbMatrixVecList.find(str);

	if (it == dbMatrixVecList.end()){
		cout << "In deleteDbMatrixVec, cannot find '" << name << "' in dbMatrixVecList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		dbMatrixVecList.erase(it);

}

/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves the name and a string in the intList
 * @param name	name of string
 * @param val	string
 */
void DataContainer::setValue(const char* name, string val) {

	string str(name);

	if (findString(name) == false)
		stringList[str] = val;
	else {
		cout << "'" << name << "' already exists in DataContainer stringList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!
 * A member function to find a specific string in the stringList
 * @param name	name of the string to be looked for
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findString(const char* name){

	string str(name);
	bool isFound;
	map<string, string>::iterator it = stringList.find(str);

	if (it == stringList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;

}

/*!
 * A member function to find a specific string in the stringList
 * @param name	name of the string to be looked for
 * @return	return 'True' if found and 'false' if not.
 */
string& DataContainer::getString(const char* name){

	string str(name);
	bool isFound;
	map<string, string>::iterator it = stringList.find(str);

	if (it == stringList.end()){
		cout << "Cannot find '" << name << "' in stringList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		return it->second;

}

/*!
 * A member function to delete a specific string in the stringList
 * @param name	name of the vector of double to be looked for
 */
void DataContainer::deleteString(const char* name){

	string str(name);
	bool isFound;
	map<string, string >::iterator it = stringList.find(str);

	if (it == stringList.end()){
		cout << "In deletestring, cannot find '" << name << "' in stringList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		stringList.erase(it);

}

/*!****************************************************************************/
/*!****************************************************************************/
/*!
 * A member function that saves the name and a string vector in the stringVectorList
 * @param name	name of string
 * @param val	string vector
 */
void DataContainer::setValue(const char* name, vector<string> val) {

	string str(name);

	if (findString(name) == false)
		stringVectorList[str] = val;
	else {
		cout << "'" << name << "' already exists in DataContainer stringList"
				<< endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

}

/*!
 * A member function to find a specific string vector in the stringVectorList
 * @param name	name of the string vector to be looked for
 * @return	return 'True' if found and 'false' if not.
 */
bool DataContainer::findStringVec(const char* name){

	string str(name);
	bool isFound;
	map<string, vector<string> >::iterator it = stringVectorList.find(str);

	if (it == stringVectorList.end())
		isFound = false;
	else
		isFound = true;

	return isFound;

}

/*!
 * A member function to find a specific string vector in the stringVectorList
 * @param name	name of the string to be looked for
 * @return	return 'True' if found and 'false' if not.
 */
vector<string>& DataContainer::getStringVec(const char* name){

	string str(name);
	bool isFound;
	map<string, vector<string> >::iterator it = stringVectorList.find(str);

	if (it == stringVectorList.end()){
		cout << "Cannot find '" << name << "' in stringVectorList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		return it->second;

}

/*!
 * A member function to delete a specific string vector in the stringVectorList
 * @param name	name of the vector of string vector to be looked for
 */
void DataContainer::deleteStringVec(const char* name){

	string str(name);
	bool isFound;
	map<string, vector<string> >::iterator it = stringVectorList.find(str);

	if (it == stringVectorList.end()){
		cout << "In deletestring, cannot find '" << name << "' in stringVectorList" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		stringVectorList.erase(it);

}
