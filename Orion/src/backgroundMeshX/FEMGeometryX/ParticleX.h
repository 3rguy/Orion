// Stores all properties of a single particle.

#ifndef ParticleX_h_
#define ParticleX_h_

#include <iostream>
#include <vector>

//#include "CustomSpline.h"
#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "InputFileData.h"
#include "IntegrationPointX.h"
//#include "MicroSpace.h"
#include "mpi.h"


class ParticleX {

  public:

    ParticleX(int usedDOF);
    ~ParticleX();

    // Particle identity number.
    void setID(int idx);
    int& getID() { return ID; };

    int& getMaterialID() {return materialID;};
    void setMaterialID(int matID){materialID=matID;};

    // Weight data manipulating 
    dbVector& getAllIntWeights() { return weights; };
    double& getIntWeight(int ID);
    double& getWeight() { return getIntWeight(0); };
    void setWeight(int ID,double value);
    void setWeight(double value) { setWeight(0,value); };

    // Elements to that this particle belongs.
    void setElems(int elem);
    intVector& getElems() { return elems; };

    // Coordinates manipulating 
    void setCoords(double coord1,double coord2,double coord3);
    double& getCoord(int idx) { return coords[idx]; };    
    dbVector& getCoords() { return coords; };

    // Influence radius manipulating
    void setRadii(double& rx,double& ry,double& rz);
    double& getRadius(int dof) { return influenceRadii[dof]; };
    dbVector& getRadii() { return influenceRadii; };

    // Degrees of freedem manipulating 
    void setDOF(int idx,double& value);
    void setDOFs(dbVector& dofs);
    double& getDOF(int idx) { return degreesOfFreedom[idx]; };
    dbVector& getDOF() { return degreesOfFreedom; };
/*
    void setOldDOF(int idx,double& value);
    void setOldDOFs(dbVector& dofs);
    double& getOldDOF(int idx) { return oldDegreesOfFreedom[idx]; };
    dbVector& getOldDOF() { return oldDegreesOfFreedom; };
*/


    double& getStepDOF(int idx) { return stepDegreesOfFreedom[idx]; };
    dbVector& getStepDOF() { return stepDegreesOfFreedom; };

    void setStepDOFMat(dbMatrix mat){stepDegreesOfFreedomMat = mat;};
    dbMatrix& getStepDOFMat(){return stepDegreesOfFreedomMat;};

    double& getDeltaDOF(int idx) { return deltaDOF[idx]; };
    dbVector& getDeltaDOF() { return deltaDOF; };

    // Supporting neighbour particles(neigbour spheres which include this 
    // particle).
    int& getSupportPtcle(int idx) { return supportingPtcls[idx]; };
    intVector& getSupportPtcls() { return supportingPtcls; };

    // Minimum distance between the particle and supporting particles
    void setLocalMinPtcleDistance(double localMinPtcleDist);
    double& getLocalMinPtcleDistance(){ return localMinPtcleDistance;};

    void setBeta(double betaVal);
    double& getBeta(){return beta;};

    void setBetaDerivs(dbVector betaDer);
    dbVector& getBetaDerivs(){return betaDerivs;};
    double& getBetaDerivs(int idx){return betaDerivs[idx];};

    void setShapeFuncsSize(int& size,int& order);
    void setShapeFuncs(int& size,dbVector& sFuncs);
    void setFirstDerivShapes(int& size,dbMatrix& firstDerivShapes);
    void setSecondDerivShapes(int& size,dbMatrix& secondDerivShapes);

    double& getShapeFunc(int idx) { return shapeFuncs[idx]; };
    double& getXDerivShape(int idx) { return firstDerivShapeFuncs[0][idx]; };
    double& getYDerivShape(int idx) { return firstDerivShapeFuncs[1][idx]; };
    double& getZDerivShape(int idx) { return firstDerivShapeFuncs[2][idx]; };
    dbVector& getShapeFuncs() { return shapeFuncs; };
    dbMatrix& getFirstDerivShapes() { return firstDerivShapeFuncs; };
    dbMatrix& getSecondDerivShapes() { return secondDerivShapeFuncs; };

    // Influencing neighbour spheres (sphere cutting if a volume or boundary 
    // integration point is also supported by another neighbouring particle).
    void setInflSpheres(int& size,intVector& iSpheres);
    intVector& getInflSpheres() { return influencingSpheres; };

    // Check if a point is supported by this particle.
    bool querySupported(InputFileData* InputData,dbVector& pointCoords,
			std::map<std::string,double>& modelData,
			std::ofstream& logFile);

    bool querySupported(InputFileData* InputData,dbVector& pointCoords,
            double& ptcleDist,std::map<std::string,double>& modelData,
            std::ofstream& logFile);

    // return customized spline

    //CustomSpline* getSpline() { return &spline; };

    /*******************************************************************/
    // Only used for particles boundary conditions are applied.

    void setSupportNBPtcle(int& DOF,int& ptcle)
      { suppNBoundPtcls[DOF].push_back(ptcle); };
    intVector& getSupportNBPtcls(int& DOF) 
      { return suppNBoundPtcls[DOF]; };

    // Influencing boundary neighbour particles.
    void setSupportBPtcle(int& DOF,int& ptcle)
      { suppBoundPtcls[DOF].push_back(ptcle); };
    intVector& getSupportBPtcls(int& DOF) 
      { return suppBoundPtcls[DOF]; };

    // Shape function values of influencing nonboundary neighbour 
    // particles at its own coordinates.
    void setNBShapeFunc(int& DOF,int& idx,double& value);
    void setNBShapeFuncs(int& DOF,int& suppportSize,
			 dbVector& shapes);
    double& getNBShapeFunc(int& DOF,int idx) { 
      return nBoundShapeFuncs[DOF][idx]; };    
    dbVector& getNBShapeFuncs(int& DOF)  { 
      return nBoundShapeFuncs[DOF]; };

    // Shape function values of influencing boundary neighbour 
    // particles at its own coordinates.
    void setBShapeFunc(int& DOF,int& position,double& value);
    double& getBShapeFunc(int& DOF,int idx) 
      { return boundShapeFuncs[DOF][idx]; };
    dbVector& getBShapeFuncs(int& DOF)  
      { return boundShapeFuncs[DOF]; };

    // Allocate several vectors.
    void setUsedDOF(int& usedDOF);

    /*******************************************************************/
    // Needed stuff for particle integration.

    // surface integration 
    dbVector& getSurfaceNormal(int ID);
    dbMatrix& getAllSurfaceNormals() { return surfaceNormals; };

    // force data manipulating.
    blMatrix& getAllBodyForceDOF() { return bodyForceDOF; };
    dbMatrix& getAllBodyForceLoads() { return bodyForceLoads; };
    blVector& getBodyForceDOF(int ID);
    dbVector& getBodyForce(int ID);

    blMatrix& getAllTractionDOF() { return tractionDOF; };
    dbMatrix& getAllTractionLoads() { return tractionLoads; };
    blVector& getTractionDOF(int ID);
    dbVector& getTraction(int ID);

    blMatrix& getAllLineForceDOF() { return lineForceDOF; };
    dbMatrix& getAllLineForceLoads() { return lineForceLoads; };
    blVector& getLineForceDOF(int ID);
    dbVector& getLineForce(int ID);

    blMatrix& getAllPointForceDOF() { return pointForceDOF; };
    dbMatrix& getAllPointForceLoads() { return pointForceLoads; };
    blVector& getPointForceDOF(int ID);
    dbVector& getPointForce(int ID);

    dbVector& getAllSurfacePressureLoads() { return surfacePressureLoads; };
    double& getSurfacePressure(int ID);

    // deformation boundary conditions
    blMatrix& getAllDeformationBoundDOF() 
      { return deformationBoundDOF; };
    dbMatrix& getAllDeformationBoundConds() 
      { return deformationBoundConds; };
    blVector& getDeformationBoundDOF(int ID);
    dbVector& getDeformationBoundConds(int ID);

    dbVector& getDeltaDeformationBoundConds(int ID);
    dbVector& getInitialDeformationBoundConds(int ID);

    // electric boundary conditions
    blMatrix& getAllElectricBoundDOF() 
      { return electricBoundDOF; };
    dbMatrix& getAllElectricBoundConds() 
      { return electricBoundConds; };
    blVector& getElectricBoundDOF(int ID);
    dbVector& getElectricBoundConds(int ID);

    dbVector& getDeltaElectricBoundConds(int ID);
    dbVector& getInitialElectricBoundConds(int ID);

    // depolarisation boundary conditions
    blMatrix& getAllDepolarisationBoundDOF() 
      { return depolarisationBoundDOF; };
    dbMatrix& getAllDepolarisationBoundConds() 
      { return depolarisationBoundConds; };
    blVector& getDepolarisationBoundDOF(int ID);
    dbVector& getDepolarisationBoundConds(int ID);

    dbVector& getDeltaDepolarisationBoundConds(int ID);
    dbVector& getInitialDepolarisationBoundConds(int ID);

    // micro boundary conditions
    blMatrix& getAllMicroBoundDOF() 
      { return electricBoundDOF; };
    dbMatrix& getAllMicroBoundConds() 
      { return electricBoundConds; };
    blVector& getMicroBoundDOF(int ID);
    dbVector& getMicroBoundConds(int ID);

    dbVector& getDeltaMicroBoundConds(int ID);
    dbVector& getInitialMicroBoundConds(int ID);

    /*******************************************************************/
    // history variables

    dbMatrix& getDeformationGradient();

    // Return the intermediate deformation gradient which is a history
    // variable needed when undeformed configuration is unknown; transfers
    // undeformed to deformed configuration
    dbMatrix& getIntermediateDefGradient();
    double& getIntermediateJacobian();

    dbMatrix& getRotationTens() { return rotationTensor; };

    // Second tensor manipulating
    dbMatrix& getSecondStrainTens() { return secondStrainTens; };

    // plasticity history variables
    dbMatrix& getPlasticityHistory() { return plasticityHistory; };

    dbVector& getInternalTraction() { return internalTraction; };

    /*******************************************************************/
    // Ghost boundary
    int& getMotherPtcle() { return motherPtcle; };
    intVector& getSupportBoundGhostPtcls() 
      { return supportBoundGhostPtcls; };

    /********************************************************************/
    // generalized continuum

    dbVector& getSpinor() { return spinor; };
    dbVector& getDeltaSpinor() { return deltaSpinor; };

    /********************************************************************/
    // generalized continuum
    //std::vector<MicroSpace>& getMicroSpaces() { return microSpaces; };

    /*******************************************************************/
    // anisotropy
    dbMatrix& getMechanicalMaterialDirections()
      { return mechanicalMaterialDirections; };




    /*******************************************************************/
    // active tension (kerckhoffs)

    double& getLc0()
      { return lc0;};
    double& getSarcomereLength()
      { return sarcomereLength;};

    dbVector& getSarcomereLengthHistory()
    { return sarcomereLengthHistory;};


    dbVector& getContractileLengthHistory()
    { return contractileLengthHistory;};

    double& getActiveTension()
      { return activeTension; };
    double& getActiveTensionVariation()
      { return activeTensionVariation; };

    double& getLR()
      { return lR; };

    /********************************************************************/
    //electrophysiology stuff

    double& getDepolarisationTime()
      { return depolarisationTime;};

    /********************************************************************/
    // Clear all class arrays which are not locally needed.
    void clearArrays();


    private:

    int ID,materialID;
    intVector elems;
    dbVector coords;
    dbVector deformations;
    dbVector influenceRadii;
    dbVector degreesOfFreedom;
    dbMatrix stepDegreesOfFreedomMat;
    dbVector stepDegreesOfFreedom; // DOF obtained during one time step
    dbVector deltaDOF;
    intVector supportingPtcls;
    intVector influencingSpheres;
    dbVector shapeFuncs;
    dbMatrix firstDerivShapeFuncs;
    dbMatrix secondDerivShapeFuncs;

    double localMinPtcleDistance;
    double beta;
    dbVector betaDerivs;

    //CustomSpline spline;

    // Only used for particles boundary conditions applied.
    //intVector ownSuppNBCounts;
    intMatrix suppBoundPtcls;
    intMatrix suppNBoundPtcls;
    dbMatrix boundShapeFuncs;
    dbMatrix nBoundShapeFuncs;

    //-------------------------------------------------------------------
    // Needed stuff for particle integration.
    dbVector weights;
    dbMatrix surfaceNormals;

    // loading
    blMatrix bodyForceDOF;
    dbMatrix bodyForceLoads;

    blMatrix tractionDOF;
    dbMatrix tractionLoads;

    blMatrix lineForceDOF;
    dbMatrix lineForceLoads;

    blMatrix pointForceDOF;
    dbMatrix pointForceLoads;

    dbVector surfacePressureLoads;

    // deformation boundary conditions
    blMatrix deformationBoundDOF;
    dbMatrix deformationBoundConds;

    dbMatrix deltaDeformationBoundConds;
    dbMatrix initialDeformationBoundConds;

    // electric boundary conditions
    blMatrix electricBoundDOF;
    dbMatrix electricBoundConds;

    dbMatrix deltaElectricBoundConds;
    dbMatrix initialElectricBoundConds;

    // depolarisation boundary conditions
    blMatrix depolarisationBoundDOF;
    dbMatrix depolarisationBoundConds;

    dbMatrix deltaDepolarisationBoundConds;
    dbMatrix initialDepolarisationBoundConds;

    // micro boundary conditions
    blMatrix microBoundDOF;
    dbMatrix microBoundConds;

    dbMatrix deltaMicroBoundConds;
    dbMatrix initialMicroBoundConds;

    /********************************************************************/
    // history variables
    dbMatrix deformationGradient;
    dbMatrix rotationTensor;
    dbVector internalTraction;
    dbMatrix secondStrainTens;

    // needed when undeformed configuration is unknown; transfers
    // undeformed to deformed configuration
    dbMatrix intermediateDefGradient;
    double intermediateJacobian;

    // history variables for Carlo's viscoplastic constitutive law
    dbMatrix plasticityHistory;

    // ghost boundary
    int motherPtcle;
    intVector supportBoundGhostPtcls;

    /********************************************************************/
    // generalized continuum

    //std::vector<MicroSpace> microSpaces;

    /********************************************************************/
    // anisotropy
    dbMatrix mechanicalMaterialDirections;
    dbMatrix pechanicalMaterialDirections;

    dbVector spinor,deltaSpinor;

    /********************************************************************/
    // active tension (kerckhoffs)

    double lc0;
    double sarcomereLength;
    dbVector sarcomereLengthHistory;
    dbVector contractileLengthHistory;
    double activeTension;
    double activeTensionVariation;
    double lR; //sarcomere length in unloaded config
               // should vary transmurally

    /********************************************************************/
    //electrophysiology stuff

    double depolarisationTime;


};

#endif
