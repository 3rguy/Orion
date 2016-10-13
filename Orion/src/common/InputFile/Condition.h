/// Stores all properties of a loading or Dirichlet condition applied to a
/// set of elements or nodes.
/// Note: The current magnitude of the condition is computed as
/// p(t) = lambda(t)*p0 where p0 denotes the reference value and lambda
/// the current condition factor. This factor is incremented in time.

/// Displacement-Point-Constraint,Displacement-Line-Constraint,Displacement-Surface-Constraint
/// Rotation-Point-Constraint,Rotation-Line-Constraint,Rotation-Surface-Constraint
/// Stress-Point-Constraint,Stress-Line-Constraint,Stress-Surface-Constraint
/// Micro-Point-Constraints,Micro-Line-Constraints,Micro-Surface-Constraints
/// Electric-Point-Constraint,Electric-Line-Constraint,Electric-Surface-Constraint,Pore-Pressure-Constraint
/// Depolarisation-Time-Point-Constraint,Depolarisation-Time-Line-Constraint,Depolarisation-Time-Surface-Constraint

/// Point-Force-Loading,Line-Force-Loading,Traction-Loading,Surface-Pressure-Loading,Body-Force-Loading
/// Point-Moment-Loading,Line-Moment-Loading,Surface-Moment-Loading,Body-Moment-Loading
/// Electric-Surface-Charge-Loading,Electric-Body-Charge-Loading,Fluid-Volume-Flux
/// Pre-Stress-Loading

/// Elastic-Point-Force, Elastic-Line-Force, Elastic-Surface-Force

/// Particle-Displacement-Control-Constraint,Cavity-Volume-Control-Constraint

/// condition-value(t) = factor(t)*ref-condition-value
/// with e.g. factor(t) = a1*dt + a0

#ifndef Condition_h_
#define Condition_h_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "ConditionElement.h"
#include "ConditionParticle.h"

class Condition {

  public:

    Condition();
    ~Condition() {}

    bool setParam(std::string name,std::string value) {
      using namespace std;
      bool info = false;
      if(name == "type") {
        type = value;
        info = true;
      }

      return info;
    };

    bool setParam(std::string name,double value) {
      using namespace std;
      bool info = false;
      if(name == "ID") {
        ID = value;
        info = true;
      }
      else if(name == "controlMode") {
        controlMode = value;
        info = true;
      }
      else if(name == "surfaceApplicationType") {
        surfaceApplicationType = value;
        info = true;
      }
      else if(name == "function") {
        function = value;
        info = true;
      }
      else if(name == "functionParam") {
        functionParams[0] = value;
        info = true;
      }
      else if(name == "conditionTime") {
        conditionTime = value;
        info = true;
      }
      else if(name == "previousConditionTime") {
        previousConditionTime = value;
        info = true;
      }
      else if(name == "conditionTimeIncrement") {
        conditionTimeIncrement = value;
        info = true;
      }
      else if(name == "minConditionTimeIncrement") {
        minConditionTimeIncrement = value;
        info = true;
      }
      else if(name == "maxConditionTimeIncrement") {
        maxConditionTimeIncrement = value;
        info = true;
      }
      else if(name == "conditionApplicationTimePeriod") {
        conditionApplicationTimePeriod = value;
        info = true;
      }
      else if(name == "factor") {
        factor = value;
        info = true;
      }
      else if(name == "previousFactor") {
        previousFactor = value;
        info = true;
      }
      else if(name == "maxFactor") {
        maxFactor = value;
        info = true;
      }
      else if(name == "minFactor") {
        minFactor = value;
        info = true;
      }
      else if(name == "maxInverseFactor") {
        maxInverseFactor = value;
        info = true;
      }
      else if(name == "referenceFactor") {
        referenceFactor = value;
        info = true;
      }
      else if(name == "increment") {
        increment = value;
        info = true;
      }
      else if(name == "dirichletControlConditionID") {
        dirichletControlConditionID = value-1;
        info = true;
      }
      else if(name == "cardiacBeatNumber") {
        cardiacBeatNumber = value;
        info = true;
      }
      else if(name == "cardiacCyclePhase") {
        cardiacCyclePhase = value;
        info = true;
      }
      else if(name == "endDiastolicCavityVolume") {
        endDiastolicCavityVolume = value;
        info = true;
      }
      else if(name == "endSystolicCavityVolume") {
        endSystolicCavityVolume = value;
        info = true;
      }
      else if(name == "endRelaxationCavityVolume") {
        endRelaxationCavityVolume = value;
        info = true;
      }
      else if(name == "currentStrokeVolume") {
        currentStrokeVolume = value;
        info = true;
      }
      else if(name == "diastolicFillingFunction") {
        diastolicFillingFunction = value;
        info = true;
      }
      else if(name == "endDiastolicTime") {
        endDiastolicTime = value;
        info = true;
      }
      else if(name == "endIVRTime") {
        endIVRTime = value;
        info = true;
      }
      else if(name == "cardiacResidualPressure") {
        cardiacResidualPressure = value;
        info = true;
      }
      else if(name == "endDiastolicPressure") {
        endDiastolicPressure = value;
        info = true;
      }
      else if(name == "endIVCPressure") {
        endIVCPressure = value;
        info = true;
      }
      else if(name == "endSystolicPressure") {
        endSystolicPressure = value;
        info = true;
      }
      else if(name == "endIVRPressure") {
        endIVRPressure = value;
        info = true;
      }
      else if(name == "rateOfChangeOfPressure") {
        rateOfChangeOfPressure = value;
        info = true;
      }
      else if(name == "previousRateOfChangeOfPressure") {
        previousRateOfChangeOfPressure = value;
        info = true;
      }
      else if(name == "flowRate") {
        flowRate = value;
        info = true;
      }
      else if(name == "previousFlowRate") {
        previousFlowRate = value;
        info = true;
      }
      else if(name == "flowRateGradient") {
        flowRateGradient = value;
        info = true;
      }
      else if(name == "flowResistance") {
        flowResistance = value;
        info = true;
      }
      else if(name == "arterialCompliance") {
        arterialCompliance = value;
        info = true;
      }
      else if(name == "peripheralResistance") {
        peripheralResistance = value;
        info = true;
      }
      else {
        std::cerr << "In Condition::setParam parameter "<< name << " does not exist." << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      return info;
    };



    int& getID() { return ID; };

    std::string& getType() { return type; };

    int& getFunction() { return function; };
    dbVector& getFunctionParams() { return functionParams; };
    blVector& getConditionDOFs() { return DOFs; };
    dbVector& getConditionValues() { return conditionValues; };
    double& getConditionTime() { return conditionTime; };
    double& getPreviousConditionTime() { return previousConditionTime; };
    double& getConditionTimeIncrement() { return conditionTimeIncrement; };
    double& getMinConditionTimeIncrement() { return minConditionTimeIncrement; };
    double& getMaxConditionTimeIncrement() { return maxConditionTimeIncrement; };

    double& getFactor() { return factor; };
    double& getPreviousFactor() { return previousFactor; };
    double& getMaxFactor();
    double& getMinFactor();
    double& getMaxInverseFactor() { return maxInverseFactor; };
    double& getMinInverseFactor() { return minInverseFactor; };
    double& getIncrement() { return increment; }; // increment of factor

    // set the maxFactor such that the condition value matches the given
    // maxValue
    void setMaxFactor(double maxValue);
    double getCurrentCondition();


    double& getConditionApplicationTimePeriod() { return conditionApplicationTimePeriod; };

    std::vector<ConditionElement>& getElements() { return elements; };
    std::vector<ConditionParticle>& getParticles() { return particles; };
    intVector& getNodes(std::ofstream& logFile);
    intVector& getGaussPoints() { return gaussPoints; };
    intVector& getLocalBoundPtcls() { return localBoundPtcls; };
    dbVector& getSurfaceNormal() { return surfaceNormal; };

    // surface orientation-dependent condition application
    int& getSurfaceApplicationType() { return surfaceApplicationType; };

    /// displacement and cavity controlled simulation
    int& getControlMode() { return controlMode; };
    int& getDirichletControlConditionID() { return dirichletControlConditionID; }; // loads only!
    int& getControlNode();
    int& getControlDOF();
    double& getReferenceFactor() { return referenceFactor; };
    double& getNewReferenceFactor() { return newReferenceFactor; };

    dbVector& getExternalDirichletControlDeformations() { return externalDirichletControlDeformations; };
    double& getStepDirichletControlDeformation() { return stepDirichletControlDeformation; };
    double& getDeltaDirichletControlDeformation() { return deltaDirichletControlDeformation; };

    double& getReferenceCavityVolume() { return referenceCavityVolume; };
    double& getCurrentReferenceCavityVolume() { return currentReferenceCavityVolume; };
    double& getTargetCavityVolume() { return targetCavityVolume; };
    double& getCurrentCavityVolume() { return currentCavityVolume; };
    double& getPreviousCavityVolume() { return previousCavityVolume; };
    double& getEndDiastolicCavityVolume() { return endDiastolicCavityVolume; };
    double& getEndSystolicCavityVolume() { return endSystolicCavityVolume; };
    double& getEndRelaxationCavityVolume() { return endRelaxationCavityVolume; };

    /// rotating conditions
    dbVector& getRotationAxis() { return rotationAxis; };
    dbVector& getRotationCenter() { return rotationCenter; }

    /// cardiac mechanics
    int& getCardiacBeatNumber() { return cardiacBeatNumber; };
    int& getCardiacCyclePhase() { return cardiacCyclePhase; };
    int& getDiastolicFillingFunction() { return diastolicFillingFunction; };
    double& getEndDiastolicTime() { return endDiastolicTime; };
    double& getEndIVRTime() { return endIVRTime; };
    double& getEndDiastolicPressure() { return endDiastolicPressure; };
    double& getEndSystolicPressure() { return endSystolicPressure; };
    double& getEndIVCPressure() { return  endIVCPressure; };
    double& getEndIVRPressure() { return  endIVRPressure; };
    // residual filling pressure
    double& getCardiacResidualPressure() { return cardiacResidualPressure; };

    double& getRateOfChangeOfPressure() { return rateOfChangeOfPressure; };
    double& getPreviousRateOfChangeOfPressure() { return previousRateOfChangeOfPressure; };
    double& getFlowRate() { return flowRate; };
    double& getPreviousFlowRate() { return previousFlowRate; };
    double& getFlowRateGradient() { return flowRateGradient; };
    double& getRequiredStrokeVolume() { return requiredStrokeVolume; };
    double& getCurrentStrokeVolume() { return currentStrokeVolume; };

    /// get the three-element windkessel parameters
    double& getArterialCompliance() { return arterialCompliance; };
    double& getPeripheralResistance() { return peripheralResistance; };
    double& getFlowResistance() { return flowResistance; };
    double& getWindkesselODEResidual() { return windkesselODEResidual; };
    double& getWindkesselODEResidualTolerence() { return windkesselODEResidualTolerence; };
    double& getPressureIncrement() { return cavityPressureIncrement; };
    int& getLinkedControlConditionID() { return linkedControlConditionID; };
    
    /// return the condition time increment of the next simulation step
    /// using the specified condition function
    double getNewConditionTimeIncrement(
        double incrementFactor,
        std::map<std::string,double>& calcData,
        std::map<std::string,double>& modelData,std::ofstream& logFile);

    // Check whether the condition time increment is within the
    // user-specified limits
    void checkConditionTimeIncrementLimits(double& newConditionTimeIncrement,
                                           std::ofstream& logFile);

    /// increase incrementally the condition factor via the condition time
    /// increment
    void setFactor(double incrementFactor,
                   std::map<std::string,double>& calcData,
                   std::map<std::string,double>& modelData,
                   std::ofstream& logFile);

    /// update the condition factor and the corresponding condition time
    /// increment
    void updateFactor(double& factorIncrement,
                      std::map<std::string,double>& calcData,
                      std::map<std::string,double>& modelData,
                      std::ofstream& logFile);

    /// determine the condition factor for given condition time
    double computeFactor(double t,std::map<std::string,double>& calcData,
                         std::map<std::string,double>& modelData,
                         std::ofstream& logFile);

    /// check whether the inverse simulation is active
    bool checkInverseSimulationActive(std::map<std::string,double>& calcData,
                                      std::map<std::string,double>& modelData,
                                      std::ofstream& logFile);

    /// check whether specified condition type is admissible
    bool isDirichletType(std::string& t);
    bool isPointDirichletCondition();
    bool isLineDirichletCondition();
    bool isSurfaceDirichletCondition();
    bool isLoadingType(std::string& t);
    bool isPointLoadingCondition();
    bool isLineLoadingCondition();
    bool isSurfaceLoadingCondition();

    /// check whether the condition is a Dirichlet-control condition.
    bool isDirichletControlType(std::string& t);
    bool isSurfaceDirichletControlCondition();
    bool isDirichletControlCondition();

  protected:
    
    int ID;
    std::string type;

    std::vector<std::string> admissibleDirichletTypes;
    std::vector<std::string> admissiblePointDirichletTypes;
    std::vector<std::string> admissibleLineDirichletTypes;
    std::vector<std::string> admissibleSurfaceDirichletTypes;
    std::vector<std::string> admissibleLoadingTypes;
    std::vector<std::string> admissiblePointLoadingTypes;
    std::vector<std::string> admissibleLineLoadingTypes;
    std::vector<std::string> admissibleSurfaceLoadingTypes;

    std::vector<std::string> admissibleDirichletControlTypes;
    std::vector<std::string> admissibleSurfaceDirichletControlTypes;

    int function;
    double conditionTimeIncrement; // usually equal to problem time increment
    double minConditionTimeIncrement,maxConditionTimeIncrement;
    double conditionTime,previousConditionTime; // usually equal to problem time
    dbVector functionParams;
    double conditionApplicationTimePeriod; // needed for dynamics

    double factor;
    double previousFactor;
    double maxFactor,minFactor;
    double maxInverseFactor,minInverseFactor;

    double increment;

    double deltaFactor,stepDeltaFactor;

    /// needed for displacement/volume control mode
    double newReferenceFactor;
    double referenceFactor;

    /// condition
    blVector DOFs;
    dbVector conditionValues;

    /// applied to elements, nodes, Gauss points or boundary particles
    std::vector<ConditionElement> elements;
    std::vector<ConditionParticle> particles;
    intVector nodes;
    intVector gaussPoints;
    intVector localBoundPtcls;

    dbVector surfaceNormal;

    // surface orientation-dependent condition application
    int surfaceApplicationType;

    /// rotating conditions
    dbVector rotationAxis;
    dbVector rotationCenter;

    /// ------------------------------------------------------------------
    /// loading application direct (1), displacement (2), cavity-volume (3)
    /// or windkessel-controlled (4)
    int controlMode;
    int dirichletControlConditionID; // ID of assembly of loading condition(s) controlled together
    double stepDirichletControlDeformation,deltaDirichletControlDeformation;
    dbVector externalDirichletControlDeformations;

    /// displacement control
    int controlNode;
    int controlDOF;

    /// cavity volume control
    double referenceCavityVolume;  // unloaded
    double currentReferenceCavityVolume; // cavity volume at outset of current simulation step
    double currentCavityVolume; // cavity volume at current simulation and iteration step
    double targetCavityVolume; // target cavity volume of current simulation step
    double previousCavityVolume; // cavity volume at previous simulation and iteration step
    double endDiastolicCavityVolume,endSystolicCavityVolume,endRelaxationCavityVolume;

    /// 3-element Windkessel
    int linkedControlConditionID;
    double cavityPressureIncrement;
    double windkesselODEResidual,windkesselODEResidualTolerence;
    double arterialCompliance,peripheralResistance,flowResistance;
    double requiredStrokeVolume,currentStrokeVolume;

    /// cardiac mechanics
    ///
    /// needs cavity Dirichlet conditions applied, in case of pressure-controlled
    /// diastole additionally set input variable
    /// 'diastolePressureControlled' = 1

    bool inverseSimulationActive;
    int cardiacCyclePhase,cardiacBeatNumber;
    int diastolicFillingFunction; // slow/fast filling
    double endDiastolicTime,endIVRTime;
    double endDiastolicPressure,endSystolicPressure,
        endIVCPressure,endIVRPressure;
    // residual filling pressure
    double cardiacResidualPressure;
    double rateOfChangeOfPressure,previousRateOfChangeOfPressure;
    double flowRate,previousFlowRate,flowRateGradient;


};

#endif
