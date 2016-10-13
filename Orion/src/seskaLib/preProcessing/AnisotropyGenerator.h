/// Project the preferred a material direction onto two opposite planes
/// of a structure, each with a user-specified normal vector 'z' and
/// rotated by a used-specified angle

#ifndef AnisotropyGenerator_h_
#define AnisotropyGenerator_h_

#include "float.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cmath>

#include <functional>   // std::plus
#include <algorithm>	// std::transform

#include "petsc.h"

#include "AnisotropyCondition.h"
#include "commonTypedefs.h"
#include "commonFunctions.h"
#include "ConditionElement.h"
#include "ConditionParticle.h"
#include "MeshlessApproximation.h"

#include "FEMGeometry.h"


using namespace std;

///
/// Paper: Jonathan Wong and Ellen Kuhl,"Generating fibre orientation maps in "
/// "human heart models using Poisson Interpolation", Computer Methods in
/// Biomechanics and Biomedical Engineering, 2012, 1-10.

/// Objective: To compute the anisotropy directions of any surface defined using
/// the epi-cardial and endo-cardial helix angle.

/// Note : The "surfaceNormal" dbMatrix in the Particle class is defined by a
/// set of vectors (in a row-wise manner) as follows:
/// [0] n, averaged nodal normal
/// [1] c, circumferential direction
/// [2] n_cz, outward pointing normal
/// [3] s, sheet normal
/// [4] p, fibre direction projection
/// [5] f, fibre direction
/// [6] m, orthogonal fibre direction to f and s

class AnisotropyGenerator {

  public:

    AnisotropyGenerator(MeshlessApproximation* MeshlessData,
                        InputFileData* InputData,
                        std::map<std::string,double>& calcData,
                        std::map<std::string,double>& modelData,
                        std::map<std::string,double>& constitutiveData,
                        std::ofstream& logFile);
    ~AnisotropyGenerator() {};

    // pre-setup
    void initPreSetup(MeshlessApproximation* MeshlessData,
                      InputFileData* InputData,
                      std::map<std::string,double>& calcData,
                      std::map<std::string,double>& modelData,
                      std::ofstream& logFile);

    void calcAllSurfaceNormals(MeshlessApproximation* MeshlessData,
                               InputFileData* InputData,
                               AnisotropyCondition& condition,
                               std::map<std::string,double>& calcData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    dbVector calcSurfaceNormal(MeshlessApproximation* MeshlessData,
                               InputFileData* InputData,
                               dbMatrix& surfacePointsCoords,
                               std::map<std::string,double>& calcData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);

    void calcNodalNormal(MeshlessApproximation* MeshlessData,
                         InputFileData* InputData,
                         AnisotropyCondition& condition,
                         std::map<std::string,double>& calcData,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile);

    void calcAveragedNodalNormals(MeshlessApproximation* MeshlessData,
                                  InputFileData* InputData,
                                  AnisotropyCondition& condition,
                                  std::map<std::string,double>& calcData,
                                  std::map<std::string,double>& modelData,
                                  std::ofstream& logFile);

    void calcNodalCircumDirections(MeshlessApproximation* MeshlessData,
                                   InputFileData* InputData,
                                   AnisotropyCondition& condition,
                                   std::map<std::string,double>& calcData,
                                   std::map<std::string,double>& modelData,
                                   std::ofstream& logFile);

    void calcOutwardNormal(MeshlessApproximation* MeshlessData,
                           InputFileData* InputData,
                           AnisotropyCondition& condition,
                           std::map<std::string,double>& calcData,
                           std::map<std::string,double>& modelData,
                           std::ofstream& logFile);

    void calcSheetNormal(MeshlessApproximation* MeshlessData,
                         InputFileData* InputData,
                         AnisotropyCondition& condition,
                         std::map<std::string,double>& calcData,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile);

    void calcFibreDirectProjection(MeshlessApproximation* MeshlessData,
                                   InputFileData* InputData,
                                   AnisotropyCondition& condition,
                                   std::map<std::string,double>& calcData,
                                   std::map<std::string,double>& modelData,
                                   std::ofstream& logFile);

    void calcFibreDirection(MeshlessApproximation* MeshlessData,
                            InputFileData* InputData,
                            AnisotropyCondition& condition,
                            std::map<std::string,double>& calcData,
                            std::map<std::string,double>& modelData,
                            std::ofstream& logFile);

    void calcOrthogonalFibreDirection(MeshlessApproximation* MeshlessData,
                                      InputFileData* InputData,
                                      AnisotropyCondition& condition,
                                      std::map<std::string,double>& calcData,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile);

    void storeMaterialDirections(MeshlessApproximation* MeshlessData,
                                 InputFileData* InputData,
                                 AnisotropyCondition& condition,
                                 intVector& ptcleRootList,
                                 std::map<std::string,double>& calcData,
                                 std::map<std::string,double>& modelData,
                                 std::map<std::string,double>& constitutiveData,
                                 std::ofstream& logFile);

    void saveResultsToFile_flexible(MeshlessApproximation* MeshlessData,
                                    InputFileData* InputData,
                                    AnisotropyCondition& condition,
                                    intVector& seskaCondPtcleID,
                                    std::map<std::string,double>& calcData,
                                    std::map<std::string,double>& modelData,
                                    std::ofstream& logFile);

    void rotateVecAboutAxis(dbVector& v, dbVector& k, double& angle,
                                             dbVector& rVec, ofstream& logFile);

    void saveResultsToFile_res_format_flexible_resultTypes(
        std::vector<std::string>& vectorResultName,
        std::vector<dbMatrix>& vectorResultMatrix,
        std::vector<std::string>& scalarResultName,
        std::vector<dbVector>& scalarResultVector,dbVector& stepValue,
        InputFileData* InputData,ofstream& logFile);

    void eliminateSpecificSurface(InputFileData* InputData,AnisotropyCondition& condition,
                                  intVector& surfElems,int surfID,
                                  std::ofstream& logFile);


//    void printVector(intVector& A,const char* msg,std::ofstream& logFile);
//    void printVector(dbVector& A,const char* msg,std::ofstream& logFile);

    double pi;
    intMatrix condPtcleRootLists;
    intMatrix seskaCondPtcleIDLists;

    std::ofstream femResFile;

};

#endif
