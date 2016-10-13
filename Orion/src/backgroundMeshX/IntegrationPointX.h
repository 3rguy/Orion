/*
 * IntegrationPointX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef INTEGRATIONPOINTX_H_
#define INTEGRATIONPOINTX_H_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"

class IntegrationPointX {

  public:

    IntegrationPointX();
    ~IntegrationPointX() {};

    // Global index datum manipulating
    void setGlobalID(int idx);
    int& getGlobalID() { return globalID; };

    int& getMaterialID() { return materialID; };

    // info relating the element this Gauss point belongs to
    intVector& getElementInfo() { return elementInfo; };

    // Weight data manipulating
    dbVector& getAllIntWeights() { return weights; };
    double& getIntWeight(int ID);
    double& getWeight() { return getIntWeight(0); };
    void setWeight(int ID,double value);
    void setWeight(double value) { setWeight(0,value); };

    // Coordinates data manipulating
    void setCoords(double coord1,double coord2,double coord3);
    void setCoords(dbVector& vec) { coords = vec; };
    double& getCoord(int idx) { return coords[idx]; };
    dbVector& getCoords() { return coords; };

    // Supporting particle data manipulating
    void setSupportPtcle(int particle);
    void checkSuppPtclsSize(int size);
    int& getSupportPtcle(int idx) { return supportingPtcls[idx]; };
    intVector& getSupportPtcls() { return supportingPtcls; };
    int getSupportCounts() { return supportingPtcls.size(); };

    // MaxEnt stuff

    // Local Minimum Distance between particles
    void setLocalMinPtcleDistance(double localMinPtcleDist);
    double& getLocalMinPtcleDistance(){ return localMinPtcleDistance;};

    dbMatrix& getFEShapeFuncDerivs() { return feShapeFuncDerivs; };

    // Shape functions data manipulating
    void setShapeFuncsSize(int& size);
    void setShapeFuncs(int& size,dbVector& sFuncs);
    double& getShapeFunc(int idx) { return shapeFuncs[idx]; };
    dbVector& getShapeFuncs() { return shapeFuncs; };

    // Derivated shape functions data manipulating
    void setFirstDerivShapesSize(int& size);
    void setSecondDerivShapesSize(int& size);
    void setFirstDerivShapes(int& size,dbMatrix& firstDerivShapes);
    void setSecondDerivShapes(int& size,dbMatrix& secondDerivShapes);

    double& getXDerivShape(int idx) { return firstDerivShapeFuncs[0][idx]; };
    double& getYDerivShape(int idx) { return firstDerivShapeFuncs[1][idx]; };
    double& getZDerivShape(int idx) { return firstDerivShapeFuncs[2][idx]; };

    dbMatrix& getFirstDerivShapes() { return firstDerivShapeFuncs; };
    dbMatrix& getSecondDerivShapes() { return secondDerivShapeFuncs; };

    // Idenfification from local to global degrees of freedom.
    void setGlobalDOF(int idx,int value);
    int& getGlobalDOF(int idx) { return globalDOF[idx]; };
    intVector& getAllGlobalDOF() { return globalDOF; };
    void setGlobalDOFSize(int size);

    // Load data manipulating.
    double& getSurfacePressure(int ID);
    double& getSurfacePressure() { return getSurfacePressure(0); };

    blVector& getLineForceDOF(int ID);
    dbVector& getLineForce(int ID);
    blVector& getLineForceDOF() { return getLineForceDOF(0); };
    dbVector& getLineForce() { return getLineForce(0); };

    blVector& getTractionDOF(int ID);
    dbVector& getTraction(int ID);
    blVector& getTractionDOF() { return getTractionDOF(0); };
    dbVector& getTraction() { return getTraction(0); };

    blVector& getBodyForceDOF(int ID);
    dbVector& getBodyForce(int ID);
    blVector& getBodyForceDOF() { return getBodyForceDOF(0); };
    dbVector& getBodyForce() { return getBodyForce(0); };

    blVector& getSurfaceMomentDOF(int ID);
    dbVector& getSurfaceMoment(int ID);
    blVector& getSurfaceMomentDOF() { return getSurfaceMomentDOF(0); };
    dbVector& getSurfaceMoment() { return getSurfaceMoment(0); };

    blVector& getBodyMomentDOF(int ID);
    dbVector& getBodyMoment(int ID);
    blVector& getBodyMomentDOF() { return getBodyMomentDOF(0); };
    dbVector& getBodyMoment() { return getBodyMoment(0); };

    blVector& getSurfaceElectricChargeDOF(int ID);
    dbVector& getSurfaceElectricCharge(int ID);
    blVector& getSurfaceElectricChargeDOF() { return getSurfaceElectricChargeDOF(0); };
    dbVector& getSurfaceElectricCharge() { return getSurfaceElectricCharge(0); };

    blVector& getBodyElectricChargeDOF(int ID);
    dbVector& getBodyElectricCharge(int ID);
    blVector& getBodyElectricChargeDOF() { return getBodyElectricChargeDOF(0); };
    dbVector& getBodyElectricCharge() { return getBodyElectricCharge(0); };

    // Dirichlet boundary conditions
    blVector& getDeformationBoundDOF(int ID);
    dbVector& getDeformationBoundConds(int ID);
    blVector& getDeformationBoundDOF() { return getDeformationBoundDOF(0); };
    dbVector& getDeformationBoundConds() { return getDeformationBoundConds(0); };
    dbVector& getDeltaDeformationBoundConds(int ID);
    dbVector& getInitialDeformationBoundConds(int ID);

    // electric boundary conditions
    blVector& getElectricBoundDOF(int ID);
    dbVector& getElectricBoundConds(int ID);
    blVector& getElectricBoundDOF() { return getElectricBoundDOF(0); };
    dbVector& getElectricBoundConds() { return getElectricBoundConds(0); };
    dbVector& getDeltaElectricBoundConds(int ID);
    dbVector& getInitialElectricBoundConds(int ID);

    // depolarization time boundary conditions
    blVector& getDepolarisationBoundDOF(int ID);
    dbVector& getDepolarisationBoundConds(int ID);
    blVector& getDepolarisationBoundDOF() { return getElectricBoundDOF(0); };
    dbVector& getDepolarisationBoundConds() { return getElectricBoundConds(0); };
    dbVector& getDeltaDepolarisationBoundConds(int ID);
    dbVector& getInitialDepolarisationBoundConds(int ID);
    // Surface normal for boundary Gauss points only.
    void setSurfaceNormal(dbVector& normal);
    dbVector& getSurfaceNormal() { return surfaceNormal; };


    /*******************************************************************/
    // history variables

    dbVector& getStepDisplacement() { return stepDisplacement; };
    dbVector& getVelocity() { return velocity; };

    // needed for energy conserving dynamics (Sansour)
    dbMatrix& getDeformationGradient();

    // Return the intermediate deformation gradient which is a history
    // variable needed when undeformed configuration is unknown; transfers
    // undeformed to deformed configuration
    dbMatrix& getIntermediateDefGradient();
    double& getIntermediateJacobian();

    // Rotation tensor manipulating
    dbMatrix& getRotationTens();

    // Second tensor manipulating
    dbMatrix& getSecondStrainTens() { return secondStrainTens; };

    // plasticity history variables
    dbMatrix& getPlasticityHistory() { return plasticityHistory; };


    /*******************************************************************/
    // Ghost boundary
    intVector& getSupportBoundGhostPtcls()
      { return supportBoundGhostPtcls; };

    // get the penalty parameters for this Gauss point for all three
    // coordinate directions
    dbVector& getPenaltyParameters() { return penaltyParameters; };

    /*******************************************************************/
    // anisotropy
    dbMatrix& getMechanicalMaterialDirections()
      { return mechanicalMaterialDirections; };

    /*******************************************************************/
    // cardiac mechanics

    // active tension (kerckhoffs)
    double& getLc0() { return lc0; };
    double& getSarcomereLength() { return sarcomereLength; };

    dbVector& getSarcomereLengthHistory()
    { return sarcomereLengthHistory;};

    dbVector& getContractileLengthHistory()
    { return contractileLengthHistory;};

    //active tension for the current calculation step.
    //set at step 0 and not changed throughout.
    double& getActiveTension() { return activeTension;};
    double& getActiveTensionVariation()
      { return activeTensionVariation;};

    double& getLR() { return lR;};

    double& getDepolarisationTime() { return depolarisationTime;};


    /*******************************************************************/
    /*******************************************************************/

    int globalID,materialID;
    intVector elementInfo;

    dbVector weights;

    dbVector coords;

    dbVector deformations;
    dbMatrix deformationDerivs;

    intVector supportingPtcls;

    // MaxEnt stuff
    double localMinPtcleDistance;
    dbMatrix feShapeFuncDerivs;

    // meshfree shape functions
    dbVector shapeFuncs;
    dbMatrix firstDerivShapeFuncs;
    dbMatrix secondDerivShapeFuncs;

    intVector globalDOF;

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



    /*******************************************************************/
    // loading
    dbVector surfacePressureLoads;
    blMatrix lineForceDOF;
    dbMatrix lineForceLoads;
    blMatrix tractionDOF;
    dbMatrix tractionLoads;
    blMatrix bodyForceDOF;
    dbMatrix bodyForceLoads;

    blMatrix surfaceMomentDOF;
    dbMatrix surfaceMomentLoads;
    blMatrix bodyMomentDOF;
    dbMatrix bodyMomentLoads;

    blMatrix surfaceElectricChargeDOF;
    dbMatrix surfaceElectricChargeLoads;
    blMatrix bodyElectricChargeDOF;
    dbMatrix bodyElectricChargeLoads;

    dbVector surfaceNormal;

    /*******************************************************************/
    // history variables
    dbVector stepDisplacement;
    dbVector velocity;
    dbMatrix deformationGradient;

    // needed when undeformed configuration is unknown; transfers
    // undeformed to deformed configuration
    dbMatrix intermediateDefGradient;
    double intermediateJacobian;

    dbMatrix rotationTensor;
    dbMatrix secondStrainTens;


    // plasticity history variables
    dbMatrix plasticityHistory;

    /*******************************************************************/
    // Ghost boundary
    intVector supportBoundGhostPtcls;

    dbVector penaltyParameters;

    /********************************************************************/
    // anisotropy
    dbMatrix mechanicalMaterialDirections;

    /********************************************************************/
    // active tension (kerckhoffs)
    double lc0;  //length of contractile element
    double sarcomereLength; //sarcomere length
    dbVector sarcomereLengthHistory;
    dbVector contractileLengthHistory;
    double activeTension; // active tension for current calculation step
    double activeTensionVariation;
    double lR;


    /********************************************************************/
    // electrophysiology

    double depolarisationTime;


};

#endif /* INTEGRATIONPOINTX_H_ */
