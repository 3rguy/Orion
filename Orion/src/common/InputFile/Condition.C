/// Stores all properties of a loading or Dirichlet condition applied to a
/// set of elements or nodes.
///
/// Note: The current magnitude of the condition is computed as
/// p(t) = lambda(t)*p0 where p0 denotes the reference value and lambda
/// the current condition factor. This factor is incremented in time.

#include "Condition.h"

Condition::Condition() :
    ID( -1), type("none"), function(0), conditionTime(0),
    previousConditionTime(0), conditionTimeIncrement(0),
    minConditionTimeIncrement(0), maxConditionTimeIncrement(0),
    conditionApplicationTimePeriod(0), factor(0), previousFactor(0),
    maxFactor(0), minFactor(0), maxInverseFactor(0), minInverseFactor(0),
    deltaFactor(0), stepDeltaFactor(0), increment(0), controlMode(1),
    controlNode( -1), controlDOF( -1), dirichletControlConditionID(0),
    rateOfChangeOfPressure(0), previousRateOfChangeOfPressure(0), flowRate(0),
    previousFlowRate(0), flowRateGradient(0), newReferenceFactor(0),
    referenceFactor(0), inverseSimulationActive(0), referenceCavityVolume(0),
    currentReferenceCavityVolume(0), currentCavityVolume(0),
    cardiacBeatNumber(1), cardiacCyclePhase(1), diastolicFillingFunction(0),
    endDiastolicTime(0), endIVRTime(0), endDiastolicPressure(0),
    endSystolicPressure(0), endIVCPressure(0), endIVRPressure(0),
    cardiacResidualPressure(0), windkesselODEResidual(0),
    windkesselODEResidualTolerence(0), arterialCompliance(0),
    peripheralResistance(0), flowResistance(0), previousCavityVolume(0),
    cavityPressureIncrement(0), linkedControlConditionID(0),
    targetCavityVolume(0), requiredStrokeVolume(0), endDiastolicCavityVolume(0),
    endSystolicCavityVolume(0), endRelaxationCavityVolume(0),
    currentStrokeVolume(0), deltaDirichletControlDeformation(0),
    stepDirichletControlDeformation(0), surfaceApplicationType(0) {

  using namespace std;

  // point Dirichlet
  pushBackVector(admissiblePointDirichletTypes,
                 (string) "Displacement-Point-Constraint");
  pushBackVector(admissiblePointDirichletTypes,
                 (string) "Rotation-Point-Constraint");
  pushBackVector(admissiblePointDirichletTypes,
                 (string) "Micro-Point-Constraint");
  pushBackVector(admissiblePointDirichletTypes,
                 (string) "Stress-Point-Constraint");
  pushBackVector(admissiblePointDirichletTypes,
                 (string) "Electric-Point-Constraint");
  pushBackVector(admissiblePointDirichletTypes,
                 (string) "Depolarisation-Time-Point-Constraint");
  pushBackVector(admissibleDirichletTypes,admissiblePointDirichletTypes);

  // line Dirichlet
  pushBackVector(admissibleLineDirichletTypes,
                 (string) "Displacement-Line-Constraint");
  pushBackVector(admissibleLineDirichletTypes,
                 (string) "Rotation-Line-Constraint");
  pushBackVector(admissibleLineDirichletTypes,(string) "Micro-Line-Constraint");
  pushBackVector(admissibleLineDirichletTypes,
                 (string) "Electric-Line-Constraint");
  pushBackVector(admissibleLineDirichletTypes,
                 (string) "Stress-Line-Constraint");
  pushBackVector(admissibleLineDirichletTypes,
                 (string) "Depolarisation-Time-Line-Constraint");
  pushBackVector(admissibleDirichletTypes,admissibleLineDirichletTypes);

  // surface Dirichlet
  pushBackVector(admissibleSurfaceDirichletTypes,
                 (string) "Displacement-Surface-Constraint");
  pushBackVector(admissibleSurfaceDirichletTypes,
                 (string) "Rotation-Surface-Constraint");
  pushBackVector(admissibleSurfaceDirichletTypes,
                 (string) "Electric-Surface-Constraint");
  pushBackVector(admissibleSurfaceDirichletTypes,
                  (string) "Pore-Pressure-Constraint");
  pushBackVector(admissibleSurfaceDirichletTypes,
                 (string) "Micro-Surface-Constraint");
  pushBackVector(admissibleSurfaceDirichletTypes,
                 (string) "Stress-Surface-Constraint");
  pushBackVector(admissibleSurfaceDirichletTypes,
                 (string) "Depolarisation-Time-Surface-Constraint");
  pushBackVector(admissibleDirichletTypes,
                 (string) "Cavity-SurfaceVolume-Control-Constraint");
  pushBackVector(admissibleDirichletTypes,admissibleSurfaceDirichletTypes);

  // point loading
  pushBackVector(admissiblePointLoadingTypes,(string) "Point-Force-Loading");
  pushBackVector(admissiblePointLoadingTypes,(string) "Point-Moment-Loading");
  pushBackVector(admissiblePointLoadingTypes,(string) "Elastic-Point-Force");
  pushBackVector(admissibleLoadingTypes,admissiblePointLoadingTypes);

  // line loading
  pushBackVector(admissibleLineLoadingTypes,(string) "Line-Force-Loading");
  pushBackVector(admissibleLineLoadingTypes,(string) "Line-Moment-Loading");
  pushBackVector(admissibleLineLoadingTypes,(string) "Elastic-Line-Force");
  pushBackVector(admissibleLoadingTypes,admissibleLineLoadingTypes);

  // surface loading
  pushBackVector(admissibleSurfaceLoadingTypes,(string) "Traction-Loading");
  pushBackVector(admissibleSurfaceLoadingTypes,
                 (string) "Surface-Pressure-Loading");
  pushBackVector(admissibleSurfaceLoadingTypes,
                 (string) "Surface-Moment-Loading");
  pushBackVector(admissibleSurfaceLoadingTypes,
                 (string) "Elastic-Surface-Force");
  pushBackVector(admissibleSurfaceLoadingTypes,
                 (string) "Electric-Surface-Charge-Loading");
  pushBackVector(admissibleSurfaceLoadingTypes,(string) "Fluid-Volume-Flux");
  pushBackVector(admissibleLoadingTypes,admissibleSurfaceLoadingTypes);

  // body loading
  pushBackVector(admissibleLoadingTypes,(string) "Body-Force-Loading");
  pushBackVector(admissibleLoadingTypes,(string) "Body-Moment-Loading");
  pushBackVector(admissibleLoadingTypes,
                 (string) "Electric-Body-Charge-Loading");
  pushBackVector(admissibleLoadingTypes, (string) "Fluid-Volume-Flux");
  pushBackVector(admissibleLoadingTypes, (string) "Pre-Stress-Loading");

  // Dirichlet control conditions
  pushBackVector(admissibleDirichletControlTypes,
                 (string) "Particle-Displacement-Control-Constraint");
  pushBackVector(admissibleSurfaceDirichletControlTypes,
                 (string) "Cavity-Volume-Control-Constraint");
  pushBackVector(admissibleDirichletControlTypes,
                 admissibleSurfaceDirichletControlTypes);

  /*********************************************************************/
  // set default function if not specified by user
  if(functionParams.size() == 0) {
    resizeArray(functionParams,2);
    functionParams[0] = 0.0;
    functionParams[1] = 1.0;
  }

}

/***********************************************************************/
/***********************************************************************/
/// check whether specified Dirichlet condition type is admissible
bool Condition::isDirichletType(std::string& t) {

  using namespace std;

  bool flag = false;

  if(find(admissibleDirichletTypes.begin(),admissibleDirichletTypes.end(),t)
    != admissibleDirichletTypes.end()) flag = true;

  return flag;
}

/// check whether condition type is a point Dirichlet condition
bool Condition::isPointDirichletCondition() {

  using namespace std;

  return (contains(admissiblePointDirichletTypes,type));
}

/// check whether condition type is a line Dirichlet condition
bool Condition::isLineDirichletCondition() {

  using namespace std;

  return (contains(admissibleLineDirichletTypes,type));
}

/// check whether condition type is a surface Dirichlet condition
bool Condition::isSurfaceDirichletCondition() {

  using namespace std;

  return (contains(admissibleSurfaceDirichletTypes,type));
}

/***********************************************************************/
/***********************************************************************/
/// check whether specified loading condition type is admissible
bool Condition::isLoadingType(std::string& t) {

  using namespace std;

  bool flag = false;

  if(find(admissibleLoadingTypes.begin(),admissibleLoadingTypes.end(),t)
    != admissibleLoadingTypes.end()) flag = true;

  return flag;
}

/// check whether condition type is a point loading
bool Condition::isPointLoadingCondition() {

  using namespace std;

  return (contains(admissiblePointLoadingTypes,type));
}

/// check whether condition type is a line loading
bool Condition::isLineLoadingCondition() {

  using namespace std;

  return (contains(admissibleLineLoadingTypes,type));
}

/// check whether condition type is a surface loading
bool Condition::isSurfaceLoadingCondition() {

  using namespace std;

  return (contains(admissibleSurfaceLoadingTypes,type));
}

/***********************************************************************/
/***********************************************************************/
/// check whether the condition is a Dirichlet-control condition.
bool Condition::isDirichletControlType(std::string& t) {

  using namespace std;

  bool flag = false;

  if(find(admissibleDirichletControlTypes.begin(),
          admissibleDirichletControlTypes.end(),t)
    != admissibleDirichletControlTypes.end()) flag = true;

  return flag;
}

bool Condition::isDirichletControlCondition() {

  using namespace std;

  return (contains(admissibleDirichletControlTypes,type));
}

/// check whether condition type is a surface Dirichlet-control condition
bool Condition::isSurfaceDirichletControlCondition() {

  using namespace std;

  return (contains(admissibleSurfaceDirichletControlTypes,type));
}

/***********************************************************************/
/***********************************************************************/
intVector& Condition::getNodes(std::ofstream& logFile) {

  using namespace std;

  if(nodes.size() == 0) {

    for(int i = 0;i < elements.size();i++) {

      ConditionElement& elem = elements[i];
      intVector& elemNodes = elem.getNodes();

      nodes.insert(nodes.end(),elemNodes.begin(),elemNodes.end());

    }

    int nSize = nodes.size();
    removeRedundantEntries(nodes,0,nSize,logFile);
  }

  return (nodes);
}

/***********************************************************************/
/***********************************************************************/
/// displacement control
int& Condition::getControlNode() {

  using namespace std;

  if(controlNode == -1) {

    if(nodes.size() != 1) {
      cerr
          << "In Condition::getControlNode each displacement-control condition\n"
          << "can only be associated with one node." << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    else controlNode = nodes[0];
  }

  return controlNode;
}
int& Condition::getControlDOF() {

  using namespace std;

  if(controlDOF == -1) {

    int n = 0;
    for(int i = 0;i < DOFs.size();i++) {
      if(DOFs[i]) {
        controlDOF = i;
        n++;
      }
    }
    if(n != 1) {
      cerr
          << "In Condition::getControlNode each displacement-control condition\n"
          << "can only be associated with one node." << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

  return controlDOF;
}

/***********************************************************************/
/***********************************************************************/
/// return the condition time increment of the next simulation step
/// using the specified condition function
double Condition::getNewConditionTimeIncrement(
    double incrementFactor,std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int dynamicsType = (int) calcData["dynamicSimulation"];

  double newFactor,newIncrement;

  // t_n+2 = t_n+1 + c*dt
  double currentTime = conditionTime + conditionTimeIncrement * incrementFactor;

  // set the admissible condition time increment dt = t_n+2 - t_n+1
  switch(function) {

  // f = Df
  case 0:
    newIncrement = 0;
    break;

    /********************************************************************/
    // f(s) = a1*s + a0 <= fmax
    // with smax = (fmax-a0)/a1
  case 1: {

    newFactor = functionParams[1] * currentTime + functionParams[0];

    // f < f0
    if(fabs(newFactor)
      < fabs(maxFactor) || fabs(newFactor-maxFactor) <= DBL_EPSILON) newIncrement =
      currentTime - conditionTime;

    // f = f0
    else if(fabs(newFactor) > fabs(maxFactor)) newIncrement = (maxFactor
      - functionParams[0]) / functionParams[1] - conditionTime;

    else newIncrement = 0;

    break;
  }
    /********************************************************************/
    // f(s) = -a1*s + a0 >= fmin
    // with smax = -(fmin-a0)/a1
  case 2: {

    newFactor = -functionParams[1] * currentTime + functionParams[0];

    // f < f0
    if(newFactor > minFactor || fabs(newFactor - minFactor) <= DBL_EPSILON) newIncrement =
      currentTime - conditionTime;

    // f = f0
    else if(newFactor < minFactor) newIncrement = -(minFactor
      - functionParams[0]) / functionParams[1] - conditionTime;

    else newIncrement = 0;

    break;
  }
    /********************************************************************/
    // f(s) = a1*s + a0 for s <= s0
    // f(s) = f0        for s > s0
  case 5: {

    // s < s0
    if(conditionApplicationTimePeriod - currentTime > DBL_EPSILON * 1.0e+03) newIncrement =
      currentTime - conditionTime;

    // s = s0
    else if(conditionApplicationTimePeriod - currentTime
      < -DBL_EPSILON * 1.0e+03
      && conditionApplicationTimePeriod - conditionTime > DBL_EPSILON * 1.0e+03) newIncrement =
      conditionApplicationTimePeriod - conditionTime;

    // s > s0
    else newIncrement = 0;

    break;
  }
    /********************************************************************/
    // f(s) = a1*s + a0                 for s <= s0  and
    // f(s) = -a1*(s-s0) + (a1*s0 + a0) for s0 < s < 2*s0
    // f(s) = 0                         for  >= 2*s0
  case 6: {

    // s <= s0
    if(conditionApplicationTimePeriod - currentTime >= -DBL_EPSILON * 1.0e+03) newIncrement =
      currentTime - conditionTime;

    // s = s0
    else if(conditionApplicationTimePeriod - currentTime
      < -DBL_EPSILON * 1.0e+03
      && conditionApplicationTimePeriod - conditionTime > DBL_EPSILON * 1.0e+03) newIncrement =
      conditionApplicationTimePeriod - conditionTime;

    // s0 < s <= 2*s0
    else if(conditionApplicationTimePeriod - currentTime
      < -DBL_EPSILON * 1.0e+03
      && 2.0 * conditionApplicationTimePeriod - currentTime
        > DBL_EPSILON * 1.0e+03) newIncrement = currentTime - conditionTime;

    // s > 2*s0
    else newIncrement = 0;

    break;
  }
  default:
    logFile << "In Condition::getIncrement condition function " << function
        << "\n" << "is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return (newIncrement);
}

/************************************************************************/
/************************************************************************/
// Check whether the condition time increment is within the
// user-specified limits
void Condition::checkConditionTimeIncrementLimits(
    double& newConditionTimeIncrement,std::ofstream& logFile) {

  using namespace std;

  // keep the simulation increment within a given range
  // adjust to user-defined minimum time increment
  if(fabs(newConditionTimeIncrement) < fabs(minConditionTimeIncrement)) {

    newConditionTimeIncrement = minConditionTimeIncrement;

    cerr << "WARNING: minimum condition increment=" << minConditionTimeIncrement
        << " enforced" << endl;
    logFile << "WARNING: minimum condition increment="
        << minConditionTimeIncrement << " enforced" << endl;

  }

  // adjust to user-defined maximum time increment
  else if(fabs(newConditionTimeIncrement) > fabs(maxConditionTimeIncrement)) {

    newConditionTimeIncrement = maxConditionTimeIncrement;

    cerr << "WARNING: maximum condition increment=" << maxConditionTimeIncrement
        << " enforced" << endl;
    logFile << "WARNING: maximum condition increment="
        << maxConditionTimeIncrement << " enforced" << endl;

  }

}

/***********************************************************************/
/***********************************************************************/
// increase incrementally the condition factor via the condition time
// increment
void Condition::setFactor(double incrementFactor,
                          std::map<std::string,double>& calcData,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile) {

  using namespace std;

  bool extend = (bool) calcData["restartFileID"];
  int dynamicsType = (int) calcData["dynamicSimulation"];
  double tMax = calcData["maxSimulationTime"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(calcData["iterationStep"] == 0
    && (calcData["simulationStep"] > 1 || extend)) {

    previousFactor = factor;
    previousConditionTime = conditionTime;
  }

  double fullFactor,fullTime,fullConditionTime;

  // (normal) forward simulation
  if( !inverseSimulationActive) {

    if(conditionTimeIncrement > conditionApplicationTimePeriod) {
      logFile << "In Condition::setFactor " << type
          << " 'condition time increment'=" << conditionTimeIncrement << "\n"
          << "must be smaller or equal to than 'condition application time period'="
          << conditionApplicationTimePeriod << "." << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Choose a loading function and set loading factor.
    switch(function) {

    // f = Df
    case 0:

      break;

      /********************************************************************/
      // f(s) = a1*s + a0 <= fmax
      // with smax = (fmax-a0)/a1
    case 1: {

      if((int) calcData["simulationStep"] == 1
        && fabs(factor) > fabs(maxFactor)) {
        logFile << "In Condition::setFactor " << type
            << " 'max. condition factor'=" << maxFactor << "\n"
            << "must be larger than current condition factor=" << factor << "."
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      if(maxFactor <= 0) {
        logFile << "In Condition::setFactor max. condition factor\n" << "of "
            << type << " ID=" << ID << " must be larger than zero." << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      // check whether simulation has finished - statics
      if(dynamicsType == 0) {

        // s_n+1
        conditionTimeIncrement *= incrementFactor;
        conditionTime = previousConditionTime + conditionTimeIncrement;

        // lambda at s_n+1
        factor = functionParams[1] * conditionTime + functionParams[0];

        // check whether simulation has finished
        if(fabs(factor - maxFactor) < DBL_EPSILON * 1.0e+03
          || fabs(factor) > fabs(maxFactor)) {
          factor = maxFactor;

          // s_n+1
          conditionTime = (factor - functionParams[0]) / functionParams[1];
          conditionTimeIncrement = conditionTime - previousConditionTime;

          calcData["simulationCompleted"] = 1;
        }

      }

      // check whether simulation has finished - dynamics
      else if(dynamicsType > 0) {

        // s_n
        conditionTime += 0.5 * conditionTimeIncrement;

        // s_n+1/2
        conditionTimeIncrement *= incrementFactor;
        conditionTime += 0.5 * conditionTimeIncrement;

        // lambda at s_n+1/2
        factor = functionParams[1] * conditionTime + functionParams[0];

        // check whether simulation has finished
        fullFactor = functionParams[1]
          * (conditionTime + 0.5 * conditionTimeIncrement) + functionParams[0];

        if(fabs(fullFactor - maxFactor) < DBL_EPSILON * 1.0e+03
          || fabs(fullFactor) > fabs(maxFactor)) {

          // s_n+1
          fullConditionTime = (maxFactor - functionParams[0])
            / functionParams[1];
          conditionTimeIncrement = fullConditionTime - previousConditionTime;

          // s_n+1/2 (in case of dynamics mid-point rule)
          conditionTime = previousConditionTime + 0.5 * conditionTimeIncrement;

          // lambda at s_n+1/2
          factor = functionParams[1] * conditionTime + functionParams[0];

          calcData["simulationCompleted"] = 1;
        }

      }

      break;
    }
      /********************************************************************/
      // f(s) = -a1*s + a0 >= fmin
      // with smax = -(fmin-a0)/a1
    case 2: {

      if(functionParams[0] <= minFactor) {
        logFile << "In Condition::setFactor min. condition factor\n" << "of "
            << type << " ID=" << ID
            << " must be larger than function parameter 'a0'." << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }
      if((int) calcData["simulationStep"] == 1 && factor > functionParams[0]) {
        logFile << "In Condition::setFactor condition factor\n" << "of " << type
            << " ID=" << ID << " must be larger than function parameter 'a0'."
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      // check whether simulation has finished - statics
      if(dynamicsType == 0) {

        // s_n+1
        conditionTimeIncrement *= incrementFactor;
        conditionTime = previousConditionTime + conditionTimeIncrement;

        // lambda at s_n+1
        factor = -functionParams[1] * conditionTime + functionParams[0];

        // check whether simulation has finished
        if(fabs(factor - minFactor) < DBL_EPSILON * 1.0e+03
          || factor < minFactor) {
          factor = minFactor;

          // s_n+1
          conditionTime = -(factor - functionParams[0]) / functionParams[1];
          conditionTimeIncrement = conditionTime - previousConditionTime;

          calcData["simulationCompleted"] = 1;
        }

      }

      // check whether simulation has finished - dynamics
      else if(dynamicsType > 0) {
        logFile << "In Condition::setFactor loading function " << function
            << "\n" << "of " << type << " ID=" << ID
            << "is not supported for dynamics!" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      break;
    }
      /********************************************************************/
      // f(s) = a1*s + a0 for s <= s0
      // f(s) = f0        for s > s0
    case 5: {

      // s_n
      conditionTime += 0.5 * conditionTimeIncrement;

      // s_n+1/2
      conditionTimeIncrement *= incrementFactor;
      conditionTime += 0.5 * conditionTimeIncrement;

      fullConditionTime = conditionTime + 0.5 * conditionTimeIncrement;

      // lambda at s_n+1/2
      if(conditionApplicationTimePeriod - fullTime >= -DBL_EPSILON * 1.0e+03) factor =
        functionParams[1] * conditionTime + functionParams[0];

      else {

        // s_n+1
        fullConditionTime = (maxFactor - functionParams[0]) / functionParams[1];
        conditionTimeIncrement = fullConditionTime - previousConditionTime;

        // s_n+1/2 (in case of dynamics mid-point rule)
        conditionTime = previousConditionTime + 0.5 * conditionTimeIncrement;

        // lambda at s_n+1/2
        factor = functionParams[1] * conditionTime + functionParams[0];
      }

      // check whether simulation has finished
      if(tMax - fullTime < -DBL_EPSILON * 1.0e+03) calcData["simulationCompleted"] =
        1;

      break;
    }
      /********************************************************************/
      // f(s) = a1*s + a0                 for s <= s0  and
      // f(s) = -a1*(s-s0) + (a1*s0 + a0) for s0 < s < 2*s0
      // f(s) = 0
    case 6: {

      // s_n
      conditionTime += 0.5 * conditionTimeIncrement;

      // s_n+1/2

      conditionTimeIncrement *= incrementFactor;
      conditionTime += 0.5 * conditionTimeIncrement;

      fullTime = conditionTime + 0.5 * conditionTimeIncrement;

      // --------------------------------------
      // lambda at s_n+1/2

      // s <= s0
      if(conditionApplicationTimePeriod - fullTime >= -DBL_EPSILON * 1.0e+03) factor =
        functionParams[1] * conditionTime + functionParams[0];

      // s0 < s <= 2*s0
      else if(conditionApplicationTimePeriod - fullTime < -DBL_EPSILON * 1.0e+03
        && 2.0 * conditionApplicationTimePeriod - fullTime
          > DBL_EPSILON * 1.0e+03) {
        factor = -functionParams[1] * conditionTime
          + 2 * functionParams[1] * conditionApplicationTimePeriod
          + functionParams[0];
      }

      // s > 2*s0
      else {

        conditionTimeIncrement = 0;

        // s_n+1/2
        conditionTime = 0;

        // lambda at s_n+1/2

        factor = 0;
      }

      // ----------------------------------------
      // check whether simulation has finished
      if(tMax - fullTime < -DBL_EPSILON * 1.0e+03) calcData["simulationCompleted"] =
        1;

      break;
    }
    default:
      logFile << "In Condition::updateFactor loading function " << function
          << "\n" << "of " << type << " ID=" << ID << "is not supported!"
          << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }
  /*********************************************************************/
  // inverse simulation (to determine unknown unloaded configuration)
  else {

    if(fabs(factor) > fabs(maxInverseFactor)
      && (bool) calcData["inverseLoadingControlledSimulation"]) {
      logFile << "In Condition::setFactor " << type << " ID=" << ID
          << " max. inverse condition factor=" << maxInverseFactor << "\n"
          << "must be larger than current condition factor=" << factor << "."
          << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    // Choose a loading function and set loading factor.
    switch(function) {

    // f = Df
    case 0:

      break;

      /********************************************************************/
      // f(s) = a1*s + a0 <= fmax
      // with smax = (fmax-a0)/a1
    case 1: {

      if(maxInverseFactor <= 0) {
        logFile << "In Condition::updateFactor max. inverse condition factor\n"
            << "of " << type << " ID=" << ID << " must be larger than zero."
            << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      // check whether simulation has finished - statics
      if(dynamicsType == 0) {

        // s_n+1
        conditionTimeIncrement *= incrementFactor;
        conditionTime = previousConditionTime + conditionTimeIncrement;

        // lambda at s_n+1
        factor = functionParams[1] * conditionTime + functionParams[0];

        // check whether simulation has finished
        if(fabs(factor - maxInverseFactor) < DBL_EPSILON * 1.0e+03
          || fabs(factor) > fabs(maxInverseFactor)) {

          factor = maxInverseFactor;

          // s_n+1
          conditionTime = (factor - functionParams[0]) / functionParams[1];
          conditionTimeIncrement = conditionTime - previousConditionTime;

          calcData["simulationCompleted"] = 1;
        }

      }

      // check whether simulation has finished - dynamics
      else if(dynamicsType > 0) {

        // s_n
        conditionTime += 0.5 * conditionTimeIncrement;

        // s_n+1/2
        conditionTimeIncrement *= incrementFactor;
        conditionTime += 0.5 * conditionTimeIncrement;

        // lambda at s_n+1/2
        factor = functionParams[1] * conditionTime + functionParams[0];

        // check whether simulation has finished
        fullFactor = functionParams[1]
          * (conditionTime + 0.5 * conditionTimeIncrement) + functionParams[0];

        if(fabs(fullFactor - maxInverseFactor) < DBL_EPSILON * 1.0e+03
          || fabs(fullFactor) > fabs(maxInverseFactor)) {

          // s_n+1
          fullConditionTime = (maxInverseFactor - functionParams[0])
            / functionParams[1];
          conditionTimeIncrement = fullConditionTime - previousConditionTime;

          // s_n+1/2 (in case of dynamics mid-point rule)
          conditionTime = previousConditionTime + 0.5 * conditionTimeIncrement;

          // lambda at s_n+1/2
          factor = functionParams[1] * conditionTime + functionParams[0];

          calcData["simulationCompleted"] = 1;
        }

      }

      break;
    }
    default:
      logFile << "In Condition::updateFactor loading function " << function
          << "\n" << "of " << type << " ID=" << ID << " is not supported!"
          << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  // update increment
  increment = factor - previousFactor;

  logFile << type << " ID=" << ID << " factor=" << factor << " increment="
      << increment << endl;
  if(rank == 0) cout << type << " ID=" << ID << " factor=" << factor
      << " increment=" << increment << endl;

}

/***********************************************************************/
/***********************************************************************/
/// update the condition factor and the corresponding condition time
/// increment
void Condition::updateFactor(double& factorIncrement,
                             std::map<std::string,double>& calcData,
                             std::map<std::string,double>& modelData,
                             std::ofstream& logFile) {

  using namespace std;

  bool extend = (bool) calcData["restartFileID"];
  int dynamicsType = (int) calcData["dynamicSimulation"];
  double tMax = calcData["maxSimulationTime"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(calcData["iterationStep"] == 0
    && (calcData["simulationStep"] > 1 || extend)) {

    previousFactor = factor;
    previousConditionTime = conditionTime;
  }

  // update the factor
  increment = factorIncrement;
  factor += increment;

  /*********************************************************************/
  // update the condition time
  //
  // Choose a loading function and set loading factor.
  switch(function) {

  // f = Df
  case 0:

    break;

    /********************************************************************/
    // f(s) = a1*s + a0 <= fmax,
    //
    // i.e.: s = (f(s)-a0)/a1
  case 1: {

    // s_n+1 or s_n+1/2 (in case of dynamics mid-point rule)
    conditionTime = (factor - functionParams[0]) / functionParams[1];

    break;
  }
    /********************************************************************/
    // f(s) = -a1*s + a0 >= fmin,
    //
    // i.e.: s = -(f(s)-a0)/a1
  case 2: {

    // s_n+1 or s_n+1/2 (in case of dynamics mid-point rule)
    conditionTime = -(factor - functionParams[0]) / functionParams[1];

    break;
  }
    /********************************************************************/
    // f(s) = a1*s + a0 for s <= s0
    // f(s) = f0        for s > s0
    //
    // i.e.: s = (f(s)-a0)/a1
  case 5: {

    // s_n+1/2
    conditionTime = (factor - functionParams[0]) / functionParams[1];

    break;
  }
    /********************************************************************/
    // f(s) = a1*s + a0                 for s <= s0  and
    // f(s) = -a1*(s-s0) + (a1*s0 + a0) for s0 < s < 2*s0
    // f(s) = 0
    //
    // i.e.: s = (f(s)-a0)/a1 and
    //       s = (2*a1*s0 + a0 - f(s))/a1
  case 6: {

    // s < s0
    if(fabs(factor) > fabs(previousFactor))

    // s_n+1/2
    conditionTime = (factor - functionParams[0]) / functionParams[1];

    else

    // s_n+1/2
    conditionTime = (2.0 * functionParams[1] * conditionApplicationTimePeriod
      + functionParams[0] - factor) / functionParams[1];

    break;
  }
  default:
    logFile << "In Condition::updateFactor loading function " << function
        << "\n" << "is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  conditionTimeIncrement = conditionTime - previousConditionTime;

  logFile << type << " ID=" << ID << " factor=" << factor << " increment="
      << increment << endl;
  if(rank == 0) cout << type << " ID=" << ID << " factor=" << factor
      << " increment=" << increment << endl;

}

/***********************************************************************/
/***********************************************************************/
/// determine the condition factor for given condition time
double Condition::computeFactor(double t,std::map<std::string,double>& calcData,
                                std::map<std::string,double>& modelData,
                                std::ofstream& logFile) {

  using namespace std;

  bool extend = (bool) calcData["restartFileID"];
  int dynamicsType = (int) calcData["dynamicSimulation"];
  double tMax = calcData["maxSimulationTime"];

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double lambda = 0;

  // Choose a loading function and set loading factor.
  switch(function) {

  // f = Df
  case 0:

    break;

    /********************************************************************/
    // f(s) = a1*s + a0 <= fmax
    // with smax = (fmax-a0)/a1
  case 1:

    lambda = functionParams[1] * t + functionParams[0];
    break;

    /********************************************************************/
    // f(s) = -a1*s + a0 >= fmin
    // with smax = -(fmin-a0)/a1
  case 2:

    lambda = -functionParams[1] * t + functionParams[0];
    break;

    /********************************************************************/
    // f(s) = a1*s + a0 for s <= s0
    // f(s) = f0        for s > s0
  case 5: {

    if(conditionApplicationTimePeriod - t >= -DBL_EPSILON * 1.0e+03) lambda =
      functionParams[1] * conditionTime + functionParams[0];

    else lambda = maxFactor;

    break;
  }
    /********************************************************************/
    // f(s) = a1*s + a0                 for s <= s0  and
    // f(s) = -a1*(s-s0) + (a1*s0 + a0) for s0 < s < 2*s0
    // f(s) = 0
  case 6: {

    // s <= s0
    if(conditionApplicationTimePeriod - t >= -DBL_EPSILON * 1.0e+03) lambda =
      functionParams[1] * conditionTime + functionParams[0];

    // s0 < s <= 2*s0
    else if(conditionApplicationTimePeriod - t < -DBL_EPSILON * 1.0e+03
      && 2.0 * conditionApplicationTimePeriod - t > DBL_EPSILON * 1.0e+03) lambda =
      -functionParams[1] * conditionTime
        + 2 * functionParams[1] * conditionApplicationTimePeriod
        + functionParams[0];

    // s > 2*s0
    else lambda = 0;

    break;
  }
  default:
    logFile << "In Condition::computeFactor loading function " << function
        << "\n" << "is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return lambda;

}

/***********************************************************************/
/***********************************************************************/
double& Condition::getMaxFactor() {

  using namespace std;

  // (normal) forward simulation
  if( !inverseSimulationActive) return maxFactor;

  // inverse simulation (to determine unknown unloaded configuration)
  else return maxInverseFactor;

}
double& Condition::getMinFactor() {

  using namespace std;

  // (normal) forward simulation
  if( !inverseSimulationActive) return minFactor;

  // inverse simulation (to determine unknown unloaded configuration)
  else return minInverseFactor;

}

/***********************************************************************/
/***********************************************************************/
// set the maxFactor such that the condition value matches the given
// maxValue
void Condition::setMaxFactor(double maxValue) {

  using namespace std;

  double value = 0;

  for(int i = 0;i < conditionValues.size();i++)

    if(DOFs[i]) value += pow(conditionValues[i],2);

  value = sqrt(value);

  // (normal) forward simulation
  if( !inverseSimulationActive) maxFactor = maxValue / value;

  // inverse simulation (to determine unknown unloaded configuration)
  else maxInverseFactor = maxValue / value;

}

/***********************************************************************/
/***********************************************************************/
/// compute the current values of the condition,
/// i.e. condition-value(t) = factor(t)*ref-condition-value
double Condition::getCurrentCondition() {

  using namespace std;

  double value = 0;

  for(int i = 0;i < conditionValues.size();i++)

    if(DOFs[i]) value += pow(conditionValues[i],2);

  value = factor * sqrt(value);
  return value;
}

/***********************************************************************/
/***********************************************************************/
/// check whether the inverse simulation is active
bool Condition::checkInverseSimulationActive(
    std::map<std::string,double>& calcData,
    std::map<std::string,double>& modelData,std::ofstream& logFile) {

  using namespace std;

  int dynamicsType = (int) calcData["dynamicSimulation"];

  double currentFactor;

  inverseSimulationActive = false;

  // Choose a loading function and set loading factor.
  switch(function) {

  // f = Df
  case 0:
    //inverseSimulationActive = false;
    break;

    /********************************************************************/
    // f(s) = a1*s + a0 <= fmax
    // with tmax = (fmax-a0)/a1
  case 1: {

    // check whether simulation has finished - statics
    if(dynamicsType == 0) {

      // lambda at s_n
      currentFactor = functionParams[1] * conditionTime + functionParams[0];

      // check whether simulation has finished
      if(fabs(currentFactor) < fabs(maxInverseFactor)) inverseSimulationActive =
        true;

    }

    // check whether simulation has finished - dynamics
    else if(dynamicsType > 0) {

      // s_n
      conditionTime += 0.5 * conditionTimeIncrement;

      // lambda at s_n
      currentFactor = functionParams[1] * conditionTime + functionParams[0];

      if(fabs(currentFactor) < fabs(maxInverseFactor)) inverseSimulationActive =
        true;

    }

    break;
  }

  default:
    logFile << "In Condition::checkInverseSimulationActive loading function "
        << function << "\n" << "is not supported!" << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return inverseSimulationActive;
}
