/*
 * FEMElementX.h
 *
 *  Created on: 16 Jul 2014
 *      Author: ritesh
 */

#ifndef FEMELEMENTX_H_
#define FEMELEMENTX_H_

#include <vector>

#include "commonTypedefs.h"
#include "ElementTemplateX.h"
#include "GaussPointSetsX.h"

class FEMElementX : public ElementTemplateX {

  public:
    FEMElementX(int usedDOF);
    ~FEMElementX() {};

    // FEM element type.
    void setElemType(int value) {elementType = value;};
    int& getElemType() { return elementType; };
    void setElemOrder(int value) {elementOrder = value;};
    int& getElemOrder() { return elementOrder; };

    // Global identifier manipulating.
    void setGlobalID(int globalID);
    int& getGlobalID() { return globalID; };

    // Return the global ID of the corresponding mother element
    int& getMotherElementID() { return motherElementID; };

    // Material 'membership' manipulating.
    void setMaterialID(int value) {materialID = value;};
    int& getMaterialID() { return materialID; };

    // Manipulating of its connected nodes.
    intVector& getNodes() { return nodes; };
    int& getNode(int idx) { return nodes[idx]; };
    void setNodes(intVector nodesVec){nodes=nodesVec;};

    // Manipulating surface elements
    intVector& getSurfaceElems(){return surfaceElems;};
    void setSurfaceElems(intVector surfElemVec){surfaceElems = surfElemVec;};

    // shape functions
    ElementTemplateX* getVolumeElementTemplate() { return FEVolumeSet; };
    ElementTemplateX* getSurfaceElementTemplate() { return FESurfaceSet; };
    ElementTemplateX* getLineElementTemplate() { return FELineSet; };

    void setVolumeElementTemplate(ElementTemplateX* pt)
      { FEVolumeSet = pt; };

    void setSurfaceElementTemplate(ElementTemplateX* pt)
      { FESurfaceSet = pt; };

    void setLineElementTemplate(ElementTemplateX* pt)
      { FELineSet = pt; };

    dbMatrix& getVolumeShapeFuncOrds() { return *volumeShapeFuncOrds; };
    dbMatrix& getSurfaceShapeFuncOrds() { return *surfaceShapeFuncOrds; };
    dbMatrix& getLineShapeFuncOrds() { return *lineShapeFuncOrds; };
    void setVolumeShapeFuncOrds(dbMatrix& ords) { volumeShapeFuncOrds = &ords; };
    void setSurfaceShapeFuncOrds(dbMatrix& ords) { surfaceShapeFuncOrds = &ords; };
    void setLineShapeFuncOrds(dbMatrix& ords) { lineShapeFuncOrds = &ords; };

    dbMatrix3& getVolumeShapeFuncDerivOrds() { return *volumeShapeFuncDerivOrds; };
    void setVolumeShapeFuncDerivOrds(dbMatrix3& ords) { volumeShapeFuncDerivOrds = &ords; };

//    dbMatrix3& getJacobian(std::vector<ParticleX>& particles,std::ofstream& logFile);

    // integration points
    intVector& getVolumeIntegrationPts() { return volumeIntegrationPts; };
    intVector& getSurfaceIntegrationPts() { return surfaceIntegrationPts; };
    intVector& getLineIntegrationPts() { return lineIntegrationPts; };

    GaussPointSetX* getVolumeGaussSet() { return VolumeGaussSet; };
    GaussPointSetX* getSurfaceGaussSet() { return SurfaceGaussSet; };
    GaussPointSetX* getLineGaussSet() { return LineGaussSet; };
    void setVolumeGaussSet(GaussPointSetX* pt) { VolumeGaussSet = pt; };
    void setSurfaceGaussSet(GaussPointSetX* pt) { SurfaceGaussSet = pt; };
    void setLineGaussSet(GaussPointSetX* pt) { LineGaussSet = pt; };

    /*******************************************************************/
    // Force data manipulating.
    double& getSurfacePressure(int ID);
    double& getSurfacePressure() { return getSurfacePressure(0); };

    blVector& getBodyForceDOF(int ID);
    dbVector& getBodyForce(int ID);
    blVector& getBodyForceDOF() { return getBodyForceDOF(0); };
    dbVector& getBodyForce() { return getBodyForce(0); };

    blVector& getTractionDOF(int ID);
    dbVector& getTraction(int ID);
    blVector& getTractionDOF() { return getTractionDOF(0); };
    dbVector& getTraction() { return getTraction(0); };

    blVector& getLineForceDOF(int ID);
    dbVector& getLineForce(int ID);
    blVector& getLineForceDOF() { return getLineForceDOF(0); };
    dbVector& getLineForce() { return getLineForce(0); };

    // moment data
    blVector& getBodyMomentDOF(int ID);
    dbVector& getBodyMoment(int ID);
    blVector& getBodyMomentDOF() { return getBodyMomentDOF(0); };
    dbVector& getBodyMoment() { return getBodyMoment(0); };

    blVector& getSurfaceMomentDOF(int ID);
    dbVector& getSurfaceMoment(int ID);
    blVector& getSurfaceMomentDOF() { return getSurfaceMomentDOF(0); };
    dbVector& getSurfaceMoment() { return getSurfaceMoment(0); };

    // electric charge data
    blVector& getSurfaceElectricChargeDOF(int ID);
    dbVector& getSurfaceElectricCharge(int ID);
    blVector& getSurfaceElectricChargeDOF() { return getSurfaceElectricChargeDOF(0); };
    dbVector& getSurfaceElectricCharge() { return getSurfaceElectricCharge(0); };

    blVector& getBodyElectricChargeDOF(int ID);
    dbVector& getBodyElectricCharge(int ID);
    blVector& getBodyElectricChargeDOF() { return getBodyElectricChargeDOF(0); };
    dbVector& getBodyElectricCharge() { return getBodyElectricCharge(0); };

    // Surface normal.
    dbVector& getSurfaceNormal(int ID);
    dbVector& getSurfaceNormal() { return getSurfaceNormal(0); };

    // deformation boundary conditions
    blVector& getLineDeformationBoundDOF(int ID);
    dbVector& getLineDeformationBoundConds(int ID);
    blVector& getLineDeformationBoundDOF() { return getLineDeformationBoundDOF(0); };
    dbVector& getLineDeformationBoundConds() { return getLineDeformationBoundConds(0); };

    blVector& getSurfaceDeformationBoundDOF(int ID);
    dbVector& getSurfaceDeformationBoundConds(int ID);
    blVector& getSurfaceDeformationBoundDOF() { return getSurfaceDeformationBoundDOF(0); };
    dbVector& getSurfaceDeformationBoundConds() { return getSurfaceDeformationBoundConds(0); };

    // electric boundary condition (electric potential)
    blVector& getLineElectricBoundDOF(int ID);
    dbVector& getLineElectricBoundConds(int ID);
    blVector& getLineElectricBoundDOF() { return getLineElectricBoundDOF(0); };
    dbVector& getLineElectricBoundConds() { return getLineElectricBoundConds(0); };

    blVector& getSurfaceElectricBoundDOF(int ID);
    dbVector& getSurfaceElectricBoundConds(int ID);
    blVector& getSurfaceElectricBoundDOF() { return getSurfaceElectricBoundDOF(0); };
    dbVector& getSurfaceElectricBoundConds() { return getSurfaceElectricBoundConds(0); };

    // depolarisation boundary condition (depolarisation time)
    blVector& getLineDepolarisationBoundDOF(int ID);
    dbVector& getLineDepolarisationBoundConds(int ID);
    blVector& getLineDepolarisationBoundDOF() { return getLineDepolarisationBoundDOF(0); };
    dbVector& getLineDepolarisationBoundConds() { return getLineDepolarisationBoundConds(0); };

    blVector& getSurfaceDepolarisationBoundDOF(int ID);
    dbVector& getSurfaceDepolarisationBoundConds(int ID);
    blVector& getSurfaceDepolarisationBoundDOF() { return getSurfaceDepolarisationBoundDOF(0); };
    dbVector& getSurfaceDepolarisationBoundConds() { return getSurfaceDepolarisationBoundConds(0); };

    // micro boundary condition
    blVector& getLineMicroBoundDOF(int ID);
    dbVector& getLineMicroBoundConds(int ID);
    blVector& getLineMicroBoundDOF() { return getLineMicroBoundDOF(0); };
    dbVector& getLineMicroBoundConds() { return getLineMicroBoundConds(0); };

    blVector& getSurfaceMicroBoundDOF(int ID);
    dbVector& getSurfaceMicroBoundConds(int ID);
    blVector& getSurfaceMicroBoundDOF() { return getSurfaceMicroBoundDOF(0); };
    dbVector& getSurfaceMicroBoundConds() { return getSurfaceMicroBoundConds(0); };

  private:
    int elementType;
    int elementOrder;
    int globalID;
    int motherElementID;
    int materialID;
    int numOfNodeDOF;
    intVector nodes;
    intVector surfaceElems;
    dbMatrix surfaceNormals;

    // integration points
    intVector volumeIntegrationPts;
    intVector surfaceIntegrationPts;
    intVector lineIntegrationPts;

    // loading conditions
    dbVector surfacePressureLoads;
    blMatrix bodyForceDOF;
    dbMatrix bodyForceLoads;
    blMatrix tractionDOF;
    dbMatrix tractionLoads;
    blMatrix lineForceDOF;
    dbMatrix lineForceLoads;

    dbMatrix bodyMomentLoads;
    blMatrix bodyMomentDOF;
    dbMatrix surfaceMomentLoads;
    blMatrix surfaceMomentDOF;

    dbMatrix surfaceElectricChargeLoads;
    blMatrix surfaceElectricChargeDOF;
    dbMatrix bodyElectricChargeLoads;
    blMatrix bodyElectricChargeDOF;

    // deformation boundary conditions
    blMatrix lineDefBoundDOF;
    dbMatrix lineDefBoundConds;
    blMatrix surfaceDefBoundDOF;
    dbMatrix surfaceDefBoundConds;

    // electric boundary conditions (electric potential)
    blMatrix lineElectricBoundDOF;
    dbMatrix lineElectricBoundConds;
    blMatrix surfaceElectricBoundDOF;
    dbMatrix surfaceElectricBoundConds;

    // depolarisation boundary conditions (depolarisation time)
    blMatrix lineDepolarisationBoundDOF;
    dbMatrix lineDepolarisationBoundConds;
    blMatrix surfaceDepolarisationBoundDOF;
    dbMatrix surfaceDepolarisationBoundConds;

    // micro boundary conditions
    blMatrix lineMicroBoundDOF;
    dbMatrix lineMicroBoundConds;
    blMatrix surfaceMicroBoundDOF;
    dbMatrix surfaceMicroBoundConds;

    // approximation tools
    dbMatrix* volumeShapeFuncOrds;
    dbMatrix3* volumeShapeFuncDerivOrds;
    dbMatrix* surfaceShapeFuncOrds;
    dbMatrix* lineShapeFuncOrds;

	dbMatrix3 jacobian; // for all element integration points

    ElementTemplateX* FEVolumeSet;
    ElementTemplateX* FESurfaceSet;
    ElementTemplateX* FELineSet;

    GaussPointSetX* VolumeGaussSet;
    GaussPointSetX* SurfaceGaussSet;
    GaussPointSetX* LineGaussSet;
};

#endif /* FEMELEMENTX_H_ */
