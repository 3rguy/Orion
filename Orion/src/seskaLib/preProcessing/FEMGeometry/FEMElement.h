/// Buffers all properties of a single FEM-Element.

#ifndef FEMElement_h_
#define FEMElement_h_

#include <vector>

#include "commonTypedefs.h"
#include "ElementTemplate.h"
#include "GaussPointSets.h"

class FEMElement : public ElementTemplate {

  public:
    FEMElement(int usedDOF);
    ~FEMElement() {};  

    /// FEM element type.
    int& getElemType() { return elementType; };
    int& getElemOrder() { return elementOrder; };

    /// Global identifier manipulating.
    void setGlobalID(int globalID);
    int& getGlobalID() { return globalID; };

    /// Return the global ID of the corresponding mother element
    int& getMotherElementID() { return motherElementID; };

    /// Material 'membership' manipulating.
    int& getMaterialID() { return materialID; };

    /// Manipulating of its connected nodes.
    intVector& getNodes() { return nodes; };
    int& getNode(int idx) { return nodes[idx]; };

    /// shape functions
    ElementTemplate* getVolumeElementTemplate() { return FEVolumeSet; };
    ElementTemplate* getSurfaceElementTemplate() { return FESurfaceSet; };
    ElementTemplate* getLineElementTemplate() { return FELineSet; };

    void setVolumeElementTemplate(ElementTemplate* pt) 
      { FEVolumeSet = pt; };

    void setSurfaceElementTemplate(ElementTemplate* pt) 
      { FESurfaceSet = pt; };

    void setLineElementTemplate(ElementTemplate* pt) 
      { FELineSet = pt; };

    dbMatrix& getVolumeShapeFuncOrds() { return *volumeShapeFuncOrds; };
    dbMatrix& getSurfaceShapeFuncOrds() { return *surfaceShapeFuncOrds; };
    dbMatrix& getLineShapeFuncOrds() { return *lineShapeFuncOrds; };
    void setVolumeShapeFuncOrds(dbMatrix& ords) { volumeShapeFuncOrds = &ords; };
    void setSurfaceShapeFuncOrds(dbMatrix& ords) { surfaceShapeFuncOrds = &ords; };
    void setLineShapeFuncOrds(dbMatrix& ords) { lineShapeFuncOrds = &ords; };

    dbMatrix3& getVolumeShapeFuncDerivOrds() { return *volumeShapeFuncDerivOrds; };
    void setVolumeShapeFuncDerivOrds(dbMatrix3& ords) { volumeShapeFuncDerivOrds = &ords; };

    dbMatrix3& getJacobian(std::vector<Particle>& particles,std::ofstream& logFile);

    /// integration points
    intVector& getVolumeIntegrationPts() { return volumeIntegrationPts; };
    intVector& getSurfaceIntegrationPts() { return surfaceIntegrationPts; };
    intVector& getLineIntegrationPts() { return lineIntegrationPts; };

    GaussPointSet* getVolumeGaussSet() { return VolumeGaussSet; };
    GaussPointSet* getSurfaceGaussSet() { return SurfaceGaussSet; };
    GaussPointSet* getLineGaussSet() { return LineGaussSet; };
    void setVolumeGaussSet(GaussPointSet* pt) { VolumeGaussSet = pt; };
    void setSurfaceGaussSet(GaussPointSet* pt) { SurfaceGaussSet = pt; };
    void setLineGaussSet(GaussPointSet* pt) { LineGaussSet = pt; };

    /*******************************************************************/
    /// Force data manipulating.
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

    blVector& getFluidVolumeFluxDOF(int ID);
    dbVector& getFluidVolumeFlux(int ID);
    blVector& getFluidVolumeFluxDOF() { return getFluidVolumeFluxDOF(0); };
    dbVector& getFluidVolumeFlux() { return getFluidVolumeFlux(0); };

    blVector& getLineForceDOF(int ID);
    dbVector& getLineForce(int ID);
    blVector& getLineForceDOF() { return getLineForceDOF(0); };
    dbVector& getLineForce() { return getLineForce(0); };

    /// moment data
    blVector& getBodyMomentDOF(int ID);
    dbVector& getBodyMoment(int ID);
    blVector& getBodyMomentDOF() { return getBodyMomentDOF(0); };
    dbVector& getBodyMoment() { return getBodyMoment(0); };

    blVector& getSurfaceMomentDOF(int ID);
    dbVector& getSurfaceMoment(int ID);
    blVector& getSurfaceMomentDOF() { return getSurfaceMomentDOF(0); };
    dbVector& getSurfaceMoment() { return getSurfaceMoment(0); };

    /// electric charge data
    blVector& getSurfaceElectricChargeDOF(int ID);
    dbVector& getSurfaceElectricCharge(int ID);
    blVector& getSurfaceElectricChargeDOF() { return getSurfaceElectricChargeDOF(0); };
    dbVector& getSurfaceElectricCharge() { return getSurfaceElectricCharge(0); };

    blVector& getBodyElectricChargeDOF(int ID);
    dbVector& getBodyElectricCharge(int ID);
    blVector& getBodyElectricChargeDOF() { return getBodyElectricChargeDOF(0); };
    dbVector& getBodyElectricCharge() { return getBodyElectricCharge(0); };

    /// Surface normal.
    dbVector& getSurfaceNormal(int ID); 
    dbVector& getSurfaceNormal() { return getSurfaceNormal(0); }; 

    /// deformation boundary conditions
    blVector& getLineDeformationBoundDOF(int ID);
    dbVector& getLineDeformationBoundConds(int ID);
    blVector& getLineDeformationBoundDOF() { return getLineDeformationBoundDOF(0); };
    dbVector& getLineDeformationBoundConds() { return getLineDeformationBoundConds(0); };

    blVector& getSurfaceDeformationBoundDOF(int ID);
    dbVector& getSurfaceDeformationBoundConds(int ID);
    blVector& getSurfaceDeformationBoundDOF() { return getSurfaceDeformationBoundDOF(0); };
    dbVector& getSurfaceDeformationBoundConds() { return getSurfaceDeformationBoundConds(0); };

    /// electric boundary condition (electric potential)
    blVector& getLineElectricBoundDOF(int ID);
    dbVector& getLineElectricBoundConds(int ID);
    blVector& getLineElectricBoundDOF() { return getLineElectricBoundDOF(0); };
    dbVector& getLineElectricBoundConds() { return getLineElectricBoundConds(0); };

    blVector& getSurfaceElectricBoundDOF(int ID);
    dbVector& getSurfaceElectricBoundConds(int ID);
    blVector& getSurfaceElectricBoundDOF() { return getSurfaceElectricBoundDOF(0); };
    dbVector& getSurfaceElectricBoundConds() { return getSurfaceElectricBoundConds(0); };

    /// TPM boundary condition

    blVector& getPorePressureBoundDOF(int ID);
    dbVector& getPorePressureBoundConds(int ID);
    blVector& getPorePressureBoundDOF() { return getPorePressureBoundDOF(0); };
    dbVector& getPorePressureBoundConds() { return getPorePressureBoundConds(0); };

    /// depolarisation boundary condition (depolarisation time)
    blVector& getLineDepolarisationBoundDOF(int ID);
    dbVector& getLineDepolarisationBoundConds(int ID);
    blVector& getLineDepolarisationBoundDOF() { return getLineDepolarisationBoundDOF(0); };
    dbVector& getLineDepolarisationBoundConds() { return getLineDepolarisationBoundConds(0); };

    blVector& getSurfaceDepolarisationBoundDOF(int ID);
    dbVector& getSurfaceDepolarisationBoundConds(int ID);
    blVector& getSurfaceDepolarisationBoundDOF() { return getSurfaceDepolarisationBoundDOF(0); };
    dbVector& getSurfaceDepolarisationBoundConds() { return getSurfaceDepolarisationBoundConds(0); };

    /// micro boundary condition
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
    dbMatrix surfaceNormals;

    /// integration points
    intVector volumeIntegrationPts;
    intVector surfaceIntegrationPts;
    intVector lineIntegrationPts;

    /// loading conditions
    dbVector surfacePressureLoads;
    blMatrix bodyForceDOF;
    dbMatrix bodyForceLoads;
    blMatrix tractionDOF;
    dbMatrix tractionLoads;
    blMatrix lineForceDOF;
    dbMatrix lineForceLoads;
    blMatrix fluidVolumeFluxDOF;
    dbMatrix fluidVolumeFluxLoads;

    dbMatrix bodyMomentLoads;
    blMatrix bodyMomentDOF;
    dbMatrix surfaceMomentLoads;
    blMatrix surfaceMomentDOF;

    dbMatrix surfaceElectricChargeLoads;
    blMatrix surfaceElectricChargeDOF;
    dbMatrix bodyElectricChargeLoads;
    blMatrix bodyElectricChargeDOF;

    /// deformation boundary conditions
    blMatrix lineDefBoundDOF;
    dbMatrix lineDefBoundConds;
    blMatrix surfaceDefBoundDOF;
    dbMatrix surfaceDefBoundConds;

    /// electric boundary conditions (electric potential)
    blMatrix lineElectricBoundDOF;
    dbMatrix lineElectricBoundConds;
    blMatrix surfaceElectricBoundDOF;
    dbMatrix surfaceElectricBoundConds;

    /// TPM boundary conditions
    blMatrix porePressureBoundDOF;
    dbMatrix porePressureBoundConds;

    /// depolarisation boundary conditions (depolarisation time)
    blMatrix lineDepolarisationBoundDOF;
    dbMatrix lineDepolarisationBoundConds;
    blMatrix surfaceDepolarisationBoundDOF;
    dbMatrix surfaceDepolarisationBoundConds;

    /// micro boundary conditions
    blMatrix lineMicroBoundDOF;
    dbMatrix lineMicroBoundConds;
    blMatrix surfaceMicroBoundDOF;
    dbMatrix surfaceMicroBoundConds;

    /// approximation tools
    dbMatrix* volumeShapeFuncOrds;
    dbMatrix3* volumeShapeFuncDerivOrds;
    dbMatrix* surfaceShapeFuncOrds;
    dbMatrix* lineShapeFuncOrds;

    dbMatrix3 jacobian; // for all element integration points

    ElementTemplate* FEVolumeSet;
    ElementTemplate* FESurfaceSet;
    ElementTemplate* FELineSet;

    GaussPointSet* VolumeGaussSet;
    GaussPointSet* SurfaceGaussSet;
    GaussPointSet* LineGaussSet;
};

#endif
