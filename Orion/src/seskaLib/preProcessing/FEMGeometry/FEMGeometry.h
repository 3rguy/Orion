#ifndef FEMGeometry_h_
#define FEMGeometry_h_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "commonFunctions.h"
#include "debugFunctions.h"
#include "defs.h"
#include "ElementTemplate.h"
#include "FEMElement.h"
#include "GaussPoint.h"
#include "GaussPointSets.h"
#include "InputFileData.h"
#include "LineElementTemplates.h"
#include "Particle.h"
#include "SurfaceElementTemplates.h"
#include "VolumeElementTemplates.h"

#include "mpi.h"
#include "petscvec.h"
#include "petscsys.h"

class FEMGeometry {

  public:

    FEMGeometry(InputFileData* InputData,
		std::map<std::string,double>& modelData,
		std::ifstream& meshFile,
		std::ofstream& logFile,PetscViewer& viewerMPI,
		PetscViewer& viewerSEQ);
    ~FEMGeometry();

    /*******************************************************************/
    // Return particles' stuff.
    int getNodesNum() { return particles.size(); };
    std::vector<Particle>& getNodesVec() { return particles; };
    std::vector<int> getOldIdx() { return oldIdx; };
    std::vector<int> getNewIdx() { return newIdx; };

    // Return FEM volume elements' staff.
    //int getGlobalElemNum() { return globalElemNum; };
    int getLocalElemNum() { return nodesElements.size(); };
    std::vector<FEMElement>& getNodesElemsVec() { return nodesElements; };
    intMatrix& getBodyForceElemIdx() { return bodyForceElemIdx; };
    intVector& getGlobalLocalElemIdx() { return globalLocalElemIdx; };

    intVector& getElementRootList() { return elementRootList; };
    intVector& getNewLocalElemIdx() { return newLocalElemIdx; };
    intMatrix& getElemGaussIdx() { return elemGaussIdx; };

    // Return FEM boundary elements' staff.
    std::vector<FEMElement>& getSurfaceNodesElemsVec() 
      { return surfaceNodesElems; };
    std::vector<FEMElement>& getLineNodesElemsVec()
      { return lineNodesElems; };

    // Return volume gauss points' stuff.
    int getGlobalGaussPtsNum() { return globalGaussPtsNum; };
    int getLocalGaussPtsNum() { return gaussPoints.size(); };

    // Return vector containing all local gauss point data.
    std::vector<GaussPoint>& getGaussPointsVec() { return gaussPoints; };

    // Return boundary gauss points' stuff.
    const int& getGlobalBGaussPtsNum() { return globalBGaussPtsNum; };

    // Return gauss points indices to that forces are applied on.
    intMatrix& getBodyForceGaussPtsIdx() 
      { return bodyForceGaussPtsIdx; };

    intMatrix& getTractionBoundGaussPtsIdx() {
      return tractionBoundGaussPtsIdx; };

    intMatrix& getSurfacePressureBoundGaussPtsIdx() {
      return surfacePressureBoundGaussPtsIdx; };

    intMatrix& getLineForceBoundGaussPtsIdx() {
      return lineForceBoundGaussPtsIdx; };

    // Return gauss points indices to that moments are applied on.
    intMatrix& getBodyMomentGaussPtsIdx() 
      { return bodyMomentGaussPtsIdx; };

    intMatrix& getSurfaceMomentBoundGaussPtsIdx() 
      { return surfaceMomentBoundGaussPtsIdx; };

    // Return gauss points indices which have electric charge loads 
    // are applied.
    intMatrix& getSurfaceElectricChargeBoundGaussPtsIdx() 
      { return surfaceElectricChargeBoundGaussPtsIdx; };

    // Return gauss points indices which have electric charge loads 
    // are applied.
    intMatrix& getBodyElectricChargeGaussPtsIdx() 
      { return bodyElectricChargeGaussPtsIdx; };

    // Return indices of boundary Gauss points with deformation boundary
    // conditions applied
    intMatrix& getSurfaceDispBoundGaussPtsIdx() {
      return surfaceDispBoundGaussPtsIdx; };

    intMatrix& getLineDispBoundGaussPtsIdx() {
      return lineDispBoundGaussPtsIdx; };

    intMatrix& getSurfaceRotBoundGaussPtsIdx() {
      return surfaceRotBoundGaussPtsIdx; };

    intMatrix& getLineRotBoundGaussPtsIdx() {
      return lineRotBoundGaussPtsIdx; };

    // Return indices of boundary Gauss points with electric boundary
    // conditions applied
    intMatrix& getSurfaceElectricBoundGaussPtsIdx() {
      return surfaceElectricBoundGaussPtsIdx; };

    intMatrix& getLineElectricBoundGaussPtsIdx() {
      return lineElectricBoundGaussPtsIdx; };

    // Return indices of boundary Gauss points with depolarisation boundary
    // conditions applied
    intMatrix& getSurfaceDepolarisationBoundGaussPtsIdx() {
      return surfaceDepolarisationBoundGaussPtsIdx; };

    intMatrix& getLineDepolarisationBoundGaussPtsIdx() {
      return lineDepolarisationBoundGaussPtsIdx; };

    // Return indices of boundary Gauss points with micro boundary
    // conditions applied
    intMatrix& getSurfaceMicroBoundGaussPtsIdx() {
      return surfaceMicroBoundGaussPtsIdx; };

    intMatrix& getLineMicroBoundGaussPtsIdx() {
      return lineMicroBoundGaussPtsIdx; };

    // Return all surface integration Gauss point indices.
    intMatrix& getAllSurfaceGaussPtsIdx(std::ofstream& logFile);

    // Return all line integration Gauss point indices.
    intMatrix& getAllLineGaussPtsIdx(std::ofstream& logFile);

    //void getAllLoadBoundGaussPtsIdx(intVector& mat,
    //std::ofstream& logFile);

    // void getAllDefBoundGaussPtsIdx(intVector& mat,
    //				  std::ofstream& logFile);

    // Return vector containing all boundary gauss point data.
    std::vector<GaussPoint>& getBoundGaussPtsVec() 
      { return boundGaussPoints; };

    /*******************************************************************/
    // Particle integration stuff

    // Calculate for all particles their particle weight 
    // (Nyström Integration) and surface and volume loads.
    void setPtcleWeightsLoads(InputFileData* InputData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile,
			      PetscViewer& viewerMPI);

    // Calculate the weights for boundary particles and store loads if
    // applied.
    void setBoundPtcleWeightsLoads(InputFileData* InputData,
				   std::map<std::string,double>& modelData,
				   std::ofstream& logFile,
				   PetscViewer& viewerMPI);

    // Store the boundary point integration points essential boundary
    // conditions and any point forces if applied.
    void setBoundPointIntegration(InputFileData* InputData,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile,
				  PetscViewer& viewerMPI);

    // Determine the integration weights for these particles at which 
    // a line boundary has to be integrated and store the line loads if any 
    // are applied.
    void setLinePtcleWeightsLoads(InputFileData* InputData,
				  std::map<std::string,double>& modelData,
				  std::ofstream& logFile,
				  PetscViewer& viewerMPI);

    // Determine weights and surface normals for these particles at which 
    // the boundary has to be integrated and store the surface loads if any 
    // is applied.
    void setSurfacePtcleWeightsLoads(InputFileData* InputData,
				     std::map<std::string,double>& modelData,
				     std::ofstream& logFile);

    // Return the indices of all particles any kind of forces are applied 
    // on.
    intMatrix& getBodyForcePtcleIdx() { return bodyForcePtcleIdx; };

    intMatrix& getTractionBoundPtcleIdx() { 
      return tractionBoundPtcleIdx; };

    intMatrix& getLineForceBoundPtcleIdx() { 
      return lineForceBoundPtcleIdx; };

    intMatrix& getSurfacePressureBoundPtcleIdx() { 
      return surfacePressureBoundPtcleIdx; };

    intMatrix& getPointForceBoundPtcleIdx() { 
      return pointForceBoundPtcleIdx; };

    // Indices of particles deformation boundary conditions are 
    // applied.
    intMatrix& getPointDispBoundPtcleIdx() { 
      return pointDispBoundPtcleIdx; };

    intMatrix& getLineDispBoundPtcleIdx() { 
      return lineDispBoundPtcleIdx; };

    intMatrix& getSurfaceDispBoundPtcleIdx() { 
      return surfaceDispBoundPtcleIdx; };

    intMatrix& getPointRotBoundPtcleIdx() { 
      return pointRotBoundPtcleIdx; };

    intMatrix& getPointElectricBoundPtcleIdx() { 
      return pointElectricBoundPtcleIdx; };

    intMatrix& getPointDepolarisationBoundPtcleIdx() { 
      return pointDepolarisationBoundPtcleIdx; };

    /*******************************************************************/
    /*******************************************************************/

  private:

    int globalBGaussPtsNum,globalGaussPtsNum;

    std::vector<Particle> particles;
    std::vector<FEMElement> nodesElements;
    std::vector<FEMElement> surfaceNodesElems;
    std::vector<FEMElement> lineNodesElems;
    std::vector<GaussPoint> gaussPoints;
    std::vector<GaussPoint> boundGaussPoints;

    /*******************************************************************/
    // Gauss quadrature

    intMatrix allLineBoundGaussPtsIdx;
    intMatrix allSurfaceBoundGaussPtsIdx;

    // Indices of local gauss points any kind of forces are applied on.
    intMatrix bodyForceGaussPtsIdx;
    intMatrix lineForceBoundGaussPtsIdx;
    intMatrix tractionBoundGaussPtsIdx;
    intMatrix surfacePressureBoundGaussPtsIdx;

    // Indices of local gauss points any kind of moments are applied.
    intMatrix lineMomentBoundGaussPtsIdx;
    intMatrix surfaceMomentBoundGaussPtsIdx;
    intMatrix bodyMomentGaussPtsIdx;

    // Indices of local gauss points any kind of electric charge loads 
    // are applied.
    intMatrix surfaceElectricChargeBoundGaussPtsIdx;    
    intMatrix bodyElectricChargeGaussPtsIdx;

    // Indices of all boundary Gauss points any deformation boundary 
    // conditions are applied
    intMatrix surfaceDispBoundGaussPtsIdx;
    intMatrix lineDispBoundGaussPtsIdx;

    intMatrix surfaceRotBoundGaussPtsIdx;
    intMatrix lineRotBoundGaussPtsIdx;

    // Indices of all boundary Gauss points any electric boundary 
    // conditions are applied
    intMatrix surfaceElectricBoundGaussPtsIdx;
    intMatrix lineElectricBoundGaussPtsIdx;

    // Indices of all boundary Gauss points any depolarisation boundary 
    // conditions are applied
    intMatrix surfaceDepolarisationBoundGaussPtsIdx;
    intMatrix lineDepolarisationBoundGaussPtsIdx;

    // Indices of all boundary Gauss points any micro boundary 
    // conditions are applied
    intMatrix surfaceMicroBoundGaussPtsIdx;
    intMatrix lineMicroBoundGaussPtsIdx;

    /*******************************************************************/
    // particle integration

    // Indices of particles any kind of forces are applied on.
    intMatrix bodyForcePtcleIdx;
    intMatrix tractionBoundPtcleIdx;
    intMatrix lineForceBoundPtcleIdx;
    intMatrix surfacePressureBoundPtcleIdx;
    intMatrix pointForceBoundPtcleIdx;

    // Indices of particles any kind of moments are applied on.


    // Indices of particles deformation boundary conditions are 
    // applied.
    intMatrix pointDispBoundPtcleIdx;
    intMatrix lineDispBoundPtcleIdx;
    intMatrix surfaceDispBoundPtcleIdx;

    intMatrix pointRotBoundPtcleIdx;

    intMatrix pointElectricBoundPtcleIdx;

    intMatrix pointDepolarisationBoundPtcleIdx;

    /*******************************************************************/
    // New node number i has original read number oldIdx[i].
    intVector oldIdx;

    // Original read node number i has new number newIdx[i].
    intVector newIdx;

    // connectivity of global to local element ID
    intVector globalLocalElemIdx;

    intVector newLocalElemIdx; // old global idx -> new local idx
    intVector elementRootList;
    intMatrix elemGaussIdx; // new local elem idx -> all its local GPts

    // element indices loads are applied.
    intMatrix bodyForceElemIdx;
    intMatrix tractionElemIdx;
    intMatrix surfacePressureElemIdx;
    intMatrix lineForceElemIdx;

    intMatrix bodyMomentElemIdx;
    intMatrix surfaceMomentElemIdx;
    intMatrix lineMomentElemIdx;
    intMatrix pointMomentElemIdx;

    intMatrix surfaceElectricChargeElemIdx;
    intMatrix bodyElectricChargeElemIdx;

    // element indices deformation boundary conditions are applied.
    intMatrix surfaceDispBoundElemIdx;
    intMatrix lineDispBoundElemIdx;

    intMatrix surfaceRotBoundElemIdx;
    intMatrix lineRotBoundElemIdx;

    intMatrix surfaceElectricBoundElemIdx;
    intMatrix lineElectricBoundElemIdx;

    intMatrix surfaceDepolarisationBoundElemIdx;
    intMatrix lineDepolarisationBoundElemIdx;

    intMatrix surfaceMicroBoundElemIdx;
    intMatrix lineMicroBoundElemIdx;

    /*******************************************************************/
    // FEM element shape functions ordinates of all nodes at all
    // elements Gauss points.
    std::vector<std::map<int,int> > volumeFEShapeFuncSetID;
    dbMatrix3 volumeFEShapeFuncSets;
    std::vector<std::map<int,int> > volumeFEShapeFuncDerivSetID;
    dbMatrix4 volumeFEShapeFuncDerivSets;

    std::vector<std::map<int,int> > surfaceFEShapeFuncSetID;
    dbMatrix3 surfaceFEShapeFuncSets;

    std::vector<std::map<int,int> > lineFEShapeFuncSetID;
    dbMatrix3 lineFEShapeFuncSets;

    // Calculate a set of volume shape functions used to determine the Gauss 
    // points of a volume element.
    dbMatrix& getVolumeFEShapeFunctions(InputFileData* InputData,
					FEMElement& elem,
					std::ofstream& logFile);

    void setVolumeFEShapeFunctions(InputFileData* InputData,
				   FEMElement& elem,
				   dbMatrix& sFuncs,
				   std::ofstream& logFile);

    // Calculate a set of volume shape function derivatives of a volume 
    // element with respect to the local element coordinates.
    dbMatrix3& getVolumeFEShapeFuncDerivs(InputFileData* InputData,
					  FEMElement& elem,
					  std::ofstream& logFile);

    void setVolumeFEShapeFuncDerivs(InputFileData* InputData,
				    FEMElement& elem,
				    dbMatrix3& sFuncs,
				    std::ofstream& logFile);

    // Calculate a set of surface shape functions used to determine the Gauss 
    // points for a surface element.
    dbMatrix& getSurfaceFEShapeFunctions(InputFileData* InputData,
					 FEMElement& elem,
					 std::ofstream& logFile);

    void setSurfaceFEShapeFunctions(InputFileData* InputData,
				    FEMElement& elem,
				    dbMatrix& sFuncs,
				    std::ofstream& logFile);

    // Calculate a set of line shape functions used to determine the Gauss 
    // points for a line element.
    dbMatrix& getLineFEShapeFunctions(InputFileData* InputData,
				      FEMElement& elem,
				      std::ofstream& logFile);

    void setLineFEShapeFunctions(InputFileData* InputData,
				 FEMElement& elem,
				 dbMatrix& sFuncs,
				 std::ofstream& logFile);

    // Compute the shape function derivatives at all Gauss points with 
    // respect to global coordinates.
    void setGlobalFEShapeFuncDerivs(InputFileData* InputData,
				    std::map<std::string,double>& modelData,
				    std::ofstream& logFile,
				    PetscViewer& viewerMPI);
      
    /*******************************************************************/
    // standard element and Gauss point templates
    std::vector<std::map<int,int> > lineElemTemplateID;
    std::vector<std::map<int,int> > surfaceElemTemplateID;
    std::vector<std::map<int,int> > volumeElemTemplateID;
    std::vector<ElementTemplate*> lineElemTemplates;
    std::vector<ElementTemplate*> surfaceElemTemplates;
    std::vector<ElementTemplate*> volumeElemTemplates;

    // Return a pointer to tools and properties of a standard element
    ElementTemplate* getVolumeElementTemplate(InputFileData* InputData,
					      FEMElement& elem,
					      std::ofstream& logFile);
    
    ElementTemplate* getSurfaceElementTemplate(InputFileData* InputData,
					       FEMElement& elem,
					       std::ofstream& logFile);
    
    ElementTemplate* getLineElementTemplate(InputFileData* InputData,
					    FEMElement& elem,
					    std::ofstream& logFile);

    //-------------------------------------------------------------------
    std::vector<std::map<int,int> > volumeGaussPtTemplateID;
    std::vector<std::map<int,int> > surfaceGaussPtTemplateID;
    std::vector<std::map<int,int> > lineGaussPtTemplateID;
    std::vector<GaussPointSet*> volumeGaussPtTemplates;
    std::vector<GaussPointSet*> surfaceGaussPtTemplates;
    std::vector<GaussPointSet*> lineGaussPtTemplates;

    //std::vector<GaussPointSet> VolumeGaussSets;
    //std::vector<GaussPointSet> SurfaceGaussSets;
    //std::vector<GaussPointSet> LineGaussSets;

    GaussPointSet* getVolumeGaussSet(InputFileData* InputData,
				     FEMElement& elem,
				     std::ofstream& logFile);

    GaussPointSet* getSurfaceGaussSet(InputFileData* InputData,
				      FEMElement& elem,
				      std::ofstream& logFile);

    GaussPointSet* getLineGaussSet(InputFileData* InputData,
				   FEMElement& elem,
				   std::ofstream& logFile);


    /*******************************************************************/

    void readMeshFile(InputFileData* InputData,
		      std::map<std::string,double>& modelData,
		      std::ifstream& meshFile,
		      std::ofstream& logFile);

    // Sort the nodes according to the chosen coordinate and then store
    // a portion of elements to each processor.
    void storeElements(InputFileData* InputData,
		       std::vector<FEMElement>& allNodesElements,
		       std::map<std::string,double>& modelData,
		       std::ofstream& logFile);

    // Sort the elements to have a similar geometrical ordering as the 
    // particles.
    void matchElemNodeOrder(InputFileData* InputData,
			    std::vector<FEMElement>& allNodesElements,
			    std::vector<int>& oldGlobElemIdx,
			    std::vector<int>& newGlobElemIdx,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile);
    
    // arrange the elements in the new ordering and store only the
    // local elements
    void rearrangeElements(InputFileData* InputData,
			   std::vector<FEMElement>& allNodesElements,
			   intVector& oldGlobElemIdx,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile);

    // Transfer element info from volume to surface or line element.
    bool swapVolToFaceElemInfo(InputFileData* InputData,
			       FEMElement& sElem,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile);

    // Rearrange the element-node configuration --> deactivated.
    void rearrangeMeshConfiguration(InputFileData* InputData,
				    std::map<std::string,double>& modelData,
				    std::ofstream& logFile,
				    PetscViewer& viewerSEQ);

    // Create boundary elements.
    void createBoundElems(InputFileData* InputData,
			  std::map<std::string,double>& modelData,
			  std::ofstream& logFile); 

    // Create surface elements and store their surface load.
    void createSurfaceElems(InputFileData* InputData,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile); 

    // Create line elements and store their line load.
    void createLineElems(InputFileData* InputData,
			 std::map<std::string,double>& modelData,
			 std::ofstream& logFile); 

    // Ensure that all elements have the necessary approximation and
    // integration tools assigned.
    void setApproxAndIntTools(InputFileData* InputData,
			      std::map<std::string,double>& modelData,
			      std::ofstream& logFile); 
    
    // Calculate for all particles their particle weight 
    // ((needed for Nyström Integration and window function norming).
    void setAllPtclsWeights(InputFileData* InputData,
			    std::ofstream& logFile);

    // Determine for all local volume Gauss points their coordinates.
    void setVolumeGaussPoints(InputFileData* InputData,
			   std::map<std::string,double>& modelData,
			   std::ofstream& logFile,
			   PetscViewer& viewerMPI);

    // Create all boundary Gauss points.
    void setBoundGaussPoints(InputFileData* InputData,
			     std::map<std::string,double>& modelData,
			     std::ofstream& logFile,
			     PetscViewer& viewerMPI);

    // Calculate the surface gauss points' coordinates and store the 
    // surface load if applied.
    void setSurfaceGaussPoints(InputFileData* InputData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerMPI);

    // Calculate the line Gauss points' coordinates and store the 
    // line loads if applied.
    void setLineGaussPoints(InputFileData* InputData,
			    std::map<std::string,double>& modelData,
			    std::ofstream& logFile,
			    PetscViewer& viewerMPI);

    // Calculate from line elements surface gauss points coordinate's and 
    // integration weight. Store the line force as surface load if applied.
    void setLineSurfaceGPoints(InputFileData* InputData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile,
			       PetscViewer& viewerMPI);

    // Calucate for all local volume Gauss points their weight.
    void setVolumeGaussWeights(InputFileData* InputData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile);

    // Calculate for all surface Gauss points their weight and their 
    // surface normal vector.
    void setSGaussWeightsNormals(InputFileData* InputData,
				 std::map<std::string,double>& modelData,
				 std::ofstream& logFile);

    // Calculate for all line gauss points their weight.
    void setLineGaussPtsWeights(InputFileData* InputData,
				std::map<std::string,double>& modelData,
				std::ofstream& logFile);

};

#endif
