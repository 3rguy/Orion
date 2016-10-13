/// Stores all properties of a microspace.
///
/// General: Each point in macro-space can have multiple superimposed
///          micro-spaces where the sum of the associated volume ratios
///          needs to be one.
///
/// Anisotropy use: each micro-space represents one or more material
///                 directions (uniaxial anisotropy (1),
///                             transverse isotropy (1),
///                             orthotropy (3)

#ifndef MicroSpace_h_
#define MicroSpace_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"   
#include "commonTypedefs.h"
#include "InputFileData.h"
#include "MicroGaussPoint.h"


class MicroSpace {


 public:

    MicroSpace() : sarcomereLength(0),restSarcomereLength(0) {};
    ~MicroSpace() {};

    std::vector<MicroGaussPoint>& getMicroIntegrationPoints()
    { return microIntPts; };
  
    /// retrieve the micro-length scale parameter
    dbMatrix& getMicroLengths(InputFileData* InputData,std::ofstream& logFile);

    /// retrieve the micro-rotation vector
    dbVector& getMicroRotationVec(InputFileData* InputData,
                                  std::ofstream& logFile);

    /// retrieve the volume ratio of each microSpace
    void setMicroVolRatio(InputFileData* InputData,std::ofstream& logFile);

    /// set the initial micro-directions (Anisotropy)
    void setInitialMicroDirections(InputFileData* InputData,dbVector& coords,
                                   dbMatrix& macroDirectors,int& macroMatID,
                                   std::map<std::string,double>& calcData,
                                   std::map<std::string,double>& modelData,
                                   std::ofstream& logFile);

    /// set the structural tensors (Anisotropy)
    void setStructuralTensors(InputFileData* InputData,
                              std::map<std::string,double>& calcData,
                              std::map<std::string,double>& modelData,
                              std::ofstream& logFile);

    /// set the micrometrics
    void setMicroMetrics(InputFileData* InputData,dbVector& coords,
                         std::map<std::string,double>& calcData,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile);

    /// set the material properties
    void setMaterialProperties(InputFileData* InputData,
                               std::map<std::string,double>& calcData,
                               std::map<std::string,double>& modelData,
                               std::ofstream& logFile);
  
  int& getMaterialID() { return materialID; };
  int& getConstitutiveLawID() { return constitutiveLawID; };
  int& getMicroSpaceID() { return microSpaceID; };
  bool& getMicroDependentMaterial() { return microDependentMaterial; };

  dbMatrix& getMicroLengths() { return microLengths; };
  dbMatrix& getMicroMetrics() { return microMetrics; };
  dbMatrix3& getMicroMetricsDerivs() { return microMetricsDerivs; };
  dbMatrix& getInitialMicroDirectors() { return initialMicroDirectors; };
  dbVector& getMicroRotationVec() { return microRotVec; };
  double& getMicroVolRatio() { return microVolRatio; };
 
  /// return all structural tensors, one for each direction
  ///(Anisotropy)

  dbMatrix& getStructuralTensor(int directionID) 
    { return structuralTensors[directionID]; };
  dbMatrix3& getStructuralTensors() { return structuralTensors; };

  double& getActiveTension(){ return activeTension; };
  double& getDepolarisationTime() { return depolarisationTime; };
  double& getRestSarcomereLength() { return restSarcomereLength; };
  double& getSarcomereLength() { return sarcomereLength; };

 private:

  int materialID; /// position in InputFileData material-data vector
  int constitutiveLawID;
  int microSpaceID;
  bool microDependentMaterial;
    
  dbMatrix microLengths;
  dbMatrix microMetrics;
  dbMatrix3 microMetricsDerivs;
  dbMatrix initialMicroDirectors;
  dbVector microRotVec;
  double microVolRatio;

  double activeTension;
  double depolarisationTime;
  double restSarcomereLength,sarcomereLength;

  dbMatrix3 structuralTensors;

  std::vector<MicroGaussPoint> microIntPts;


};

#endif
