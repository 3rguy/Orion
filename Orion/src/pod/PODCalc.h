/*
 * PODCalc.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef PODCALC_H_
#define PODCALC_H_


#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "defs.h"
#include "Database.h"
#include "Data.h"
#include "ParticleX.h"
//#include "ParticleDistributionX.h"
#include "InputFileData.h"
#include "GridNodes.h"

class PODCalc
{
public:

    PODCalc(dbMatrix& fullMatrix,dbMatrix& reducedMatrix, double enLev,
    		InputFileData* InputData, ofstream& logFile);

    ~PODCalc(){};

    // Set the energy level
    void setEnergyLevel(double enLev){energyLevel=enLev;};

    // Set data matrix
    void setDataMatrix(dbMatrix& mat,InputFileData* InputData, ofstream& logFile);

    void setPOVsandPOMs(InputFileData* InputData,ofstream& logFile);
    dbMatrix& getPOMs(){return POMs;};
    dbVector& getPOVs(){return POVs;};

    dbVector& getMeanVec(){return meanVec;};

    void SVDCalc(ofstream& logFile);
    void SVD(dbMatrix& A,dbMatrix& U,dbMatrix& V,dbVector& S,
    		 ofstream& logFile);
    void KLDCalc(ofstream& logFile);
    void KLD(ofstream& logFile); //Not yet implemented

    void snapshotCalc(ofstream& logFile);
    void snapshotCalc_two(ofstream& logFile);

    void cleanPOVsAndPOMs(dbVector& vec,dbMatrix& mat,ofstream& logFile);

    dbMatrix& getPOMsConserved(ofstream& logFile);
    dbVector& getPOVsConserved(ofstream& logFile);

    int& getNumPOVConserved(){return numOfPOVsConserved;}

    void energyConservCalc(ofstream& logFile);
    double& getEnergyConserved(){return energyConserved;};

    // Compression functions
    void compressMatrix(dbMatrix& reducedMatrix,InputFileData* InputData,
    		ofstream& logFile);
    void compressMatrix(dbMatrix& fullMatrix,dbMatrix& reducedMatrix,
    		InputFileData* InputData,ofstream& logFile);

    void compressVector(dbVector& fullVector,dbVector& reducedVector,
    		InputFileData* InputData,ofstream& logFile);

    // Expansion functions
    void expandMatrix(dbMatrix& reducedMatrix,dbMatrix& fullMatrix,
    		InputFileData* InputData,ofstream& logFile);

    void expandVector(dbVector& reducedVector,dbVector& fullVector,
    		InputFileData* InputData,ofstream& logFile);

private:

    dbMatrix dataMat;
    dbVector meanVec;

    double energyLevel;
    double energyConserved;

    dbVector POVs;
    dbMatrix POMs;

    int numOfPOVsConserved;
};

#endif /* PODCALC_H_ */
