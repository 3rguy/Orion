#include "ParticleDistribution.h"
#include <functional>


/**********************************************************************/
/**********************************************************************/
// Determine for each process a portion of particles and set a root list,
// which process each particle belongs to.
void ParticleDistribution::setPtcleProcDistrib(InputFileData* InputData,
                                               std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int startIdx,endIdx;
  int ptclePortion = (int)ceil((double)particlesNum/size);

  if(ptclePortion*(rank+1) <= particlesNum) {
    startIdx = ptclePortion*rank;
    endIdx = ptclePortion*(rank+1);
  }
  else if(ptclePortion*rank <= particlesNum) {
    startIdx = ptclePortion*rank;
    endIdx = particlesNum;
  }
  else {
    startIdx = particlesNum;
    endIdx = particlesNum;
  }
    
  ptclePortion = endIdx-startIdx;

  globalLocalPtcleIdx.resize(particlesNum);
  globalLocalPtcleIdx.assign(globalLocalPtcleIdx.size(),-1);

  intVector localRootList(ptclePortion);
  exclusiveLocalPtcls.resize(ptclePortion);
  int m=0;

  for(int i=startIdx;i<endIdx;i++) {
    localRootList[m] = rank;
    globalLocalPtcleIdx[i] = m;
    exclusiveLocalPtcls[m] = i;
    m++;
  }

  /*********************************************************************/
  // Assemble the global list.

  ptcleRootList.resize(particlesNum);

  int* recvCounts = new int[size];
  int* displs = new int[size];

  MPI_Allgather(&ptclePortion,1,MPI_INT,recvCounts,1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  MPI_Allgatherv(&localRootList[0],recvCounts[rank],MPI_INT,
		 &ptcleRootList[0],recvCounts,displs,MPI_INT,
		 MPI_COMM_WORLD);

  delete[] displs,recvCounts;

//#ifdef _geometryDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"################ set particle distribution ##############"<<endl;
  logFile<<"particle portion: "<<ptclePortion<<endl;
  logFile<<"startIdx: "<<startIdx<<"; endIdx: "<<endIdx<<"; ptclePortion: "
	 <<ptclePortion<<endl;
  logFile<<"******************* local root list *********************"<<endl;
  for(int i=0;i<ptclePortion;i++)
    logFile<<"PARTICLE "<<i<<" "<<localRootList[i]<<endl;
  logFile<<"********************* root list *************************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"PARTICLE "<<i<<" "<<ptcleRootList[i]<<endl;
  }
//#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine the necessary influence radii for all particles.
void ParticleDistribution::setInfluenceRadii(InputFileData* InputData,
                                             std::map<std::string,double>& modelData,
                                             std::ofstream& logFile,
                                             PetscViewer& viewerMPI,
                                             PetscViewer& viewerSEQ) {

  using namespace std;

  double multiplier = InputData->getValue("influenceRadiusMultiplier");
  int mode = (int)InputData->getValue("radiusDeterminationAlgorithm");
  int minPtcleSupport = (int)InputData->getValue("minParticleSupport");

  if(minPtcleSupport > particlesNum) {
    logFile<<"In ParticleDistribution::setInfluenceRadii too less "
	   <<"particles in the system ("<<particlesNum<<") to "
	   <<"satisfy the minimum particle support ("
	   <<minPtcleSupport<<")!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  dbVector globalInfluenceRadii(particlesNum*3);

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef _influenceRadiusDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"************ Calculation of influence radii **********"<<endl;
  logFile<<"rank: "<<rank<<endl;
  logFile<<"******************** rootlist ************************"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"PARTICLE "<<i<<": "<<ptcleRootList[i]<<endl;
  logFile<<"******************************************************"<<endl;
#endif

  dbVector localInfRadii;
  dbVector localInfluenceRadii;

  switch(mode) {

    // determine particle influence radii correspondingly to the distances
    // to their neighbouring particles
  case 1:

    setPtcleRadsDistDepend(InputData,modelData,logFile);

    break;

    //--------------------------------------------------------------------
    // Calculated the particles' influence radius support dependent.
    //case 2:

    //setPtcleRadsSupportDepend(InputData,modelData,logFile);

    //for(int i=0;i<particlesNum;i++) {

    //localInfRadii = particles[i].getRadii();
    //dbVector& globalInfRadii = particles[i].getRadii();

    //MPI_Allreduce(&localInfRadii[0],&globalInfRadii[0],3,
    //	    MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    //}

    //break;

    //-------------------------------------------------------------------
    // Calculated the particles' influence radius weight dependent.
  case 3:

    for(int i=0;i<particlesNum;i++) {

      dbVector& radii = particles[i].getRadii();

      for(int j=0;j<radii.size();j++)

	radii[j] = particles[i].getWeight();

    }

    break;

    //-------------------------------------------------------------------
    // determine particle influence radii correspondingly to the distances
    // to their neighbouring particles, but with different values for
    // negative and positive coordinate direction
  case 5:

    setPtcleRadsDistDepend2(InputData,modelData,logFile);

    break;

  default:

    logFile<<"In ParticleDistribution::setInfluenceRadii Particle "
	   <<"influence radii have not been determined.\n"
	   <<"Set input datum 'radiusDeterminationAlgorithm'<1,2,3>"
	   <<"is not set."<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


  /**********************************************************************/
  // Post-process all particles' influence radii.

  postProcInfRads(InputData,modelData,logFile);

}

/************************************************************************/
/************************************************************************/
// determine particle influence radii correspondingly to the distances
// to their neighbouring particles
//
// controlling input parameter 
// -> minDirectionalPtcleSupport: minimum number of neighbour particles
//    for each coordinate directions (plus/minus); the number should 
//    not be larger than maximal possible support everywhere
// -> minDirectPtcleSuppReduction: support deficit in plus or minus 
//    direction is usually added to the opposite direction (plus/minus)
//    but can be reduced by this parameter
void ParticleDistribution::setPtcleRadsDistDepend(InputFileData* InputData,
						  std::vector<Particle>& ptcls,
                                                  std::map<std::string,double>& modelData,
                                                  std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int numOfPtcls = ptcls.size();

  for(int i=0;i<numOfPtcls;i++) {

    dbVector& radii = ptcls[i].getRadii();
    radii.resize(usedDims);

  }


#ifdef _influenceRadiusDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"****** default influence zone determination ********"<<endl;
#endif

  dbMatrix delta = getKroneckerSymbol(usedDims);
  dbMatrix g(usedDims,dbVector(usedDims));

  double searchAngle =
    InputData->getValue("influenceRadiusDeterminationAngle");

  // the search directions are the cartesian coordinate axes
  if(searchAngle != 45) {

    searchAngle = 0;
    g = delta;

  }

  // diagonal search: rotate each of the cartesian coordinate axes by 45
  // e_1 = [  1/sqrt(3), -1/sqrt(3), 1/sqrt(3) ]
  // e_2 = [ -1/sqrt(3),  1/sqrt(3), 1/sqrt(3) ]
  // e_3 = [ -1/sqrt(3), -1/sqrt(3), 1/sqrt(3) ]

  else {

    g[0][0] =  1.0/sqrt(3.0);
    g[0][1] = -1.0/sqrt(3.0);
    g[0][2] =  1.0/sqrt(3.0);

    g[1][0] = -1.0/sqrt(3.0);
    g[1][1] =  1.0/sqrt(3.0);
    g[1][2] =  1.0/sqrt(3.0);

    g[2][0] = -1.0/sqrt(3.0);
    g[2][1] = -1.0/sqrt(3.0);
    g[2][2] =  1.0/sqrt(3.0);

  }


#ifdef _influenceRadiusDebugMode_
  logFile<<"----------------------------------------------------"<<endl;
  logFile<<"search directions:"<<endl;
  for(int i=0;i<g.size();i++) {
    logFile<<"e["<<i<<"]: ";
    for(int j=0;j<g[i].size();j++)
      logFile<<g[i][j]<<" ";
    logFile<<endl;
  }
#endif

  // set the search sector respectively angles

  double pi = 3.141592654;
  searchAngle = searchAngle*pi/180.0;

  double minAngle = searchAngle; // alpha_0
  double maxAngle = pi/6.0+searchAngle; // 30 degree + alpha_0

#ifdef _influenceRadiusDebugMode_
  logFile<<"----------------------------------------------------"<<endl;
  logFile<<"searchAngle: "<<minAngle<<" - "<<maxAngle
	 <<endl;
#endif

  /*********************************************************************/
  // choice between two modes:
  // 1) particle 'i' gets influence radius to support enough neighbour 
  //    particles (active support)
  // 2) closest neighbour particles 'j' get suitable influence radius
  //    such that particle 'i' has enough supporting neighbour 
  //    particles (passive support)
  
  bool activeSupport = false;


  int minPtcleSupport = (int)InputData->getValue("minDirectionalPtcleSupport");
  int supportReduction = (int)InputData->getValue("minDirectPtcleSuppReduction");

  double distance,absDist,maxAbsDist;
  dbVector absDists(numOfPtcls);
  dbMatrix maxDists(numOfPtcls,dbVector(usedDims));
  dbMatrix dists(numOfPtcls,dbVector(usedDims));

  int round = 0;
  double tol = 1.0e-08;

  double angle,dirCos,weight;
  int minusCounter,plusCounter;
  intVector idxPlus(numOfPtcls);
  intVector idxMinus(numOfPtcls);

  dbVector weightedRadsPlus(numOfPtcls);
  dbVector weightedRadsMinus(numOfPtcls);

  intVector neededPtcls(usedDims);
  dbVector neededRads(usedDims);

  // ensure that each particle 'i' has got the minimum number of
  // supporting particles
  for(int i=0;i<numOfPtcls;i++) {

#ifdef _influenceRadiusDebugMode_
    logFile<<"**************************************************"<<endl;
    logFile<<"PARTICLE "<<i<<" coords: "<<ptcls[i].getCoord(0)
	   <<" "<<ptcls[i].getCoord(1)<<" "<<ptcls[i].getCoord(2)<<endl;
#endif

    clearArray(absDists);

    // Loop over all particles and determine the distance to the
    // neighbouring particles.
    for(int j=0;j<numOfPtcls;j++) {

      // do not include current particle 'i' itself as a neighbour
      if(i == j)
	continue;

#ifdef _influenceRadiusDebugMode_
      logFile<<"********"<<endl;
      logFile<<"neighbour particle: "<<j<<endl;
#endif

      for(int k=0;k<usedDims;k++) {

	dists[j][k] =
	  ptcls[i].getCoord(k)-ptcls[j].getCoord(k);

	absDists[j] += pow(dists[j][k],2);


#ifdef _influenceRadiusDebugMode_
	logFile<<"distance["<<k<<"]: "<<dists[j][k]<<" maxDist: "<<maxDists[i][k]<<endl;
#endif

	if(maxDists[i][k] < fabs(dists[j][k]))
	  maxDists[i][k] = fabs(dists[j][k]);

      }

      absDists[j] = sqrt(absDists[j]);


#ifdef _influenceRadiusDebugMode_
      logFile<<"absDist: "<<absDists[j]<<endl;
#endif

    }


    round = 0;

    clearArray(neededPtcls);
    clearArray(neededRads);

    // determine separately for each direction the closest particles
    for(int dof=0;dof<usedDims;dof++) {

      // 10 degree*round + 30 degree + alpha_0; for each fail the search
      // angle is increased by 10 degree
      maxAngle = pi/18*round + pi/6.0 + searchAngle;

#ifdef _influenceRadiusDebugMode_
      logFile<<"--------------------------------------------------"<<endl;
      logFile<<"DOF "<<dof<<" round "<<round<<": "<<endl;
#endif

      plusCounter = 0;
      minusCounter = 0;

      for(int j=0;j<numOfPtcls;j++) {

	// do not include current particle 'i' itself as a neighbour
	if(i == j)
	  continue;

	dirCos = 0;

	// determine the directional cosine foreach search direction
	// separately
	switch(dof) {

	  // search direction 1
	case 0:

	  for(int k=0;k<usedDims;k++)

	    dirCos += g[0][k]*dists[j][k];

	  break;

	  // search direction 2
	case 1:

	  for(int k=0;k<usedDims;k++)

	    dirCos += g[1][k]*dists[j][k];

	  break;

	  // search direction 3
	case 2:

	  for(int k=0;k<usedDims;k++)

	    dirCos += g[2][k]*dists[j][k];

	  break;
	}


	if(absDists[j]> 0)
	  dirCos = dirCos/absDists[j];

	else
	  dirCos = 0;

	// determine the angle included by the search direction and
	// neighbour particles direction with respect to particle 'i'
	angle = acos(fabs(dirCos));

	if(angle < minAngle || angle > maxAngle) {

#ifdef _influenceRadiusDebugMode_
	  logFile<<"neigh ptcle "<<j
		 <<" dirCos "<<dirCos<<" angle: "<<angle<<" dists ";
	  for(int k=0;k<usedDims;k++)
	    logFile<<dists[j][k]<<" ";
	  logFile<<"dof "<<dof<<" SKIPPED (not in search sector)"<<endl;
#endif


	  continue;

	}

#ifdef _influenceRadiusDebugMode_
	else {

	  logFile<<plusCounter+minusCounter<<".) neigh ptcle "<<j
		 <<" dirCos "<<dirCos<<" angle: "<<angle<<" dists ";
	  for(int k=0;k<usedDims;k++)
	    logFile<<dists[j][k]<<" ";
	  logFile<<"dof "<<dof<<endl;

	}
#endif

	// plus direction
	if(dirCos >= 0 && fabs(dists[j][dof]) > tol) {

	  idxPlus[plusCounter] = j;
	  weightedRadsPlus[plusCounter] = fabs(dists[j][dof]);

	  plusCounter++;
	}

	// minus direction
	else if(dirCos < 0 && fabs(dists[j][dof]) > tol) {

	  idxMinus[minusCounter] = j;
	  weightedRadsMinus[minusCounter] = fabs(dists[j][dof]);

	  minusCounter++;
	}

      }

      sortValuesIdx(weightedRadsPlus,idxPlus,0,plusCounter-1);
      sortValuesIdx(weightedRadsMinus,idxMinus,0,minusCounter-1);

      // remove redundant entries
      int m=0;

      for(int j=1;j<plusCounter;j++) {

	if(fabs(weightedRadsPlus[m]-weightedRadsPlus[j]) > tol) {
	  ++m;
	  weightedRadsPlus[m] = weightedRadsPlus[j];
	  idxPlus[m] = idxPlus[j];
	}

      }

      if(plusCounter > 0)
	plusCounter = m+1;

      m=0;

      for(int j=1;j<minusCounter;j++) {

	if(fabs(weightedRadsMinus[m]-weightedRadsMinus[j]) > tol) {
	  ++m;
	  weightedRadsMinus[m] = weightedRadsMinus[j];
	  idxMinus[m] = idxMinus[j];
	}

      }

      if(minusCounter > 0)
	minusCounter = m+1;

      // ----------------------------------------------------------------
      // set the array index of 'idxPlus' and 'idxMinus'

      int plusID,minusID;

      // available support in positive coordinate direction more than
      // minimum particle support
      if(plusCounter >= minPtcleSupport)
	plusID = minPtcleSupport-1;

      // available support in positive coordinate direction less than
      // minimum particle support
      else if(plusCounter < minPtcleSupport && plusCounter > 0)
	plusID = plusCounter-1;

      // no available support in positive coordinate direction
      else
	plusID = -1;

      // available support in negative coordinate direction more than
      // minimum particle support
      if(minusCounter >= minPtcleSupport)
	minusID = minPtcleSupport-1;

      // available support in negative coordinate direction less than
      // minimum particle support
      else if(minusCounter < minPtcleSupport && minusCounter > 0)
	minusID = minusCounter-1;

      // no available support in negative coordinate direction
      else
	minusID = -1;

      // ----------------------------------------------------------------
      // calculate the support deficit in positive coordinate direction
      // and increase the support in negative direction correspondingly

      // mininum support in plus direction insufficient
      if(plusCounter < minPtcleSupport) {

	int deficit = minPtcleSupport - plusCounter;
	deficit -= supportReduction;

	// add the reduced deficit to the minus direction and
	// neglect plus direction
	if(minusID+deficit < minusCounter) {

	  if(deficit > 0)
	    minusID += deficit;

	  plusID = -1;
	}

	// also insufficient support in minus direction -> plus
	// direction must be respected, even if bigger than minus one
	else if(plusID >= 0 && minusID >= 0) {

	  if(weightedRadsMinus[minusID] < weightedRadsPlus[plusID] &&
	     fabs(weightedRadsMinus[minusID]-weightedRadsPlus[plusID]) > tol)

	    minusID = -1;

	  else

	    plusID = -1;

	}

	// add no deficit and neglect plus direction
	else {
	  minusID = minusCounter-1;
	  plusID = -1;
	}

      }

      // ----------------------------------------------------------------
      // alternatively calculate the support deficit in positive
      // coordinate direction and increase the support in positive
      // direction correspondingly

      // minimum support in minus direction insufficient
      else if(minusCounter < minPtcleSupport) {

	int deficit = minPtcleSupport - minusCounter;
	deficit -= supportReduction;

	// add the reduced deficit to the plus direction and
	// neglect minus direction
	if(plusID+deficit < plusCounter) {

	  if(deficit > 0)
	    plusID += deficit;

	  minusID = -1;
	}

	// also insufficient support in plus direction -> minus
	// direction must be respected even if bigger than plus one
	else if(plusID >= 0 && minusID >= 0) {

	  if(weightedRadsPlus[plusID] < weightedRadsMinus[minusID] &&
	     fabs(weightedRadsPlus[plusID]-weightedRadsMinus[minusID]) > tol)

	    plusID = -1;

	  else

	    minusID = -1;

	}

	// add no deficit and neglect minus direction
	else {
	  plusID = plusCounter-1;
	  minusID = -1;
	}

      }


#ifdef _influenceRadiusDebugMode_
      logFile<<"------------------------------------------------"<<endl;
      logFile<<"plusCounter = "<<plusCounter<<" minusCounter = "
	     <<minusCounter<<endl;
      logFile<<"plusID = "<<plusID<<" minusID = "<<minusID<<endl;
      for(int j=0;j<plusCounter;j++) {
	logFile<<j<<".) neigh ptcle: "<<idxPlus[j]
	       <<" weightedRadPlus["<<j<<"] = "
	       <<weightedRadsPlus[j]<<" dists: ";
	for(int k=0;k<usedDims;k++)
	  logFile<<dists[idxPlus[j]][k]<<" ";
	logFile<<endl;
      }
      for(int j=0;j<minusCounter;j++) {
	logFile<<j<<".) neigh ptcle: "<<idxMinus[j]
	       <<" weightedRadMinus["<<j<<"] = "
	       <<weightedRadsMinus[j]<<" dists: ";
	for(int k=0;k<usedDims;k++)
	  logFile<<dists[idxMinus[j]][k]<<" ";
	logFile<<endl;
      }
#endif

      if(activeSupport) {
	
	dbVector& radii = ptcls[i].getRadii();
	
	// enough neighbours in plus and minus direction available
	if(plusID >= 0 && minusID >= 0) {
	  
	  if(weightedRadsPlus[plusID] <
	     weightedRadsMinus[minusID])
	    
	    radii[dof] = weightedRadsPlus[plusID];
	  
	  else
	    
	    radii[dof] = weightedRadsMinus[minusID];
	  
	  round = 0;
	}
	
	// enough neighbours only plus direction available
	else if(plusID >= 0 && minusID < 0) {
	  
	  radii[dof] = weightedRadsPlus[plusID];
	  round = 0;
	}
	
	// enough neighbours only minus direction available
	else if(plusID < 0 && minusID >= 0) {
	  
	  radii[dof] = weightedRadsMinus[minusID];
	  round = 0;
	}
	
	// no neighbour found -> search angle is increased by 10 degree
	// and the search algorithm for current DOF is executed again
	else if(maxAngle < minAngle + pi*15.0/18.0) {
	  
	  round += 1;
	  dof -= 1;
	  
	}
	else {
	  logFile<<"In ParticleDistribution::setPtcleRadsDistDepend particle\n "
		 <<i<<" DOF "<<dof<<" has not enough neighbour particles!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
	
      }
      
      else {
	
	// enough neighbours in plus and minus direction available
	if(plusID >= 0 && minusID >= 0) {
	  
	  if(weightedRadsPlus[plusID] <
	     weightedRadsMinus[minusID]) {
	    
	    neededPtcls[dof] = idxPlus[plusID];
	    neededRads[dof] = weightedRadsPlus[plusID];
	  }
	  else {
	    
	    neededPtcls[dof] = idxMinus[minusID];
	    neededRads[dof] = weightedRadsMinus[minusID];
	  }
	  
	  round = 0;
	}
	
	// enough neighbours only plus direction available
	else if(plusID >= 0 && minusID < 0) {
	  
	  neededPtcls[dof] = idxPlus[plusID];
	  neededRads[dof] = weightedRadsPlus[plusID];
	  round = 0;
	}
	
	// enough neighbours only minus direction available
	else if(plusID < 0 && minusID >= 0) {
	  
	  neededPtcls[dof] = idxMinus[minusID];
	  neededRads[dof] = weightedRadsMinus[minusID];
	  round = 0;
	}
	
	// no neighbour found -> search angle is increased by 10 degree
	// and the search algorithm for current DOF is executed again
	else if(maxAngle < minAngle + pi*15.0/18.0) {
	  
	  round += 1;
	  dof -= 1;
	  
	}
	else {
	  logFile<<"In ParticleDistribution::setPtcleRadsDistDepend particle\n "
		 <<i<<" DOF "<<dof<<" has not enough neighbour particles!"
		 <<endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
	
      }
      
    }
    
    if(!activeSupport) {
      
#ifdef _influenceRadiusDebugMode_
      logFile<<"*************************************************"<<endl;
      logFile<<"*********** needed influence radii **************"<<endl;
      for(int dof=0;dof<neededPtcls.size();dof++) {
	logFile<<"NEIGHBOUR "<<neededPtcls[dof]<<": dist["<<dof<<"]="
	       <<neededRads[dof]<<endl;
      }
#endif

      // determine the distance of those particles found to be required to 
      // support particle 'i' and ensure that all other particles 'j' with 
      // equal or less distance also support particle 'i'
      maxAbsDist = 0;
      
      for(int dof=0;dof<neededPtcls.size();dof++) {
	
	if(absDists[neededPtcls[dof]] > maxAbsDist)
	  maxAbsDist = absDists[neededPtcls[dof]];
	
      }
      
      for(int j=0;j<numOfPtcls;j++) {
	
	if(absDists[j] <= maxAbsDist) {
	  
	  dbVector& radii = ptcls[j].getRadii();
	  
	  for(int k=0;k<dists[j].size();k++) {
	    
	    if(radii[k] < fabs(dists[j][k]))
	      radii[k] = fabs(dists[j][k]);
	    
	  }
	  
	}

      }
	
      
#ifdef _influenceRadiusDebugMode_
      intVector idx(numOfPtcls);
      for(int j=0;j<numOfPtcls;j++)
	idx[j] = j;
      logFile<<"maxAbsDist="<<maxAbsDist<<endl;
      logFile<<"*********** distance sorted neighbours ***********"<<endl;
      sortValuesIdx(absDists,idx,0,absDists.size()-1);
      for(int j=0;j<numOfPtcls;j++) {
	dbVector& nRads = ptcls[idx[j]].getRadii();
	logFile<<"NEIGHBOUR "<<idx[j]<<": dist="<<absDists[j]<<" radii: ";
	for(int k=0;k<nRads.size();k++)
	  logFile<<nRads[k]<<" ";
	logFile<<"(dists: ";
	for(int k=0;k<dists[j].size();k++)
	  logFile<<dists[j][k]<<" " ;
	logFile<<")"<<endl;
      }
#endif

    }
    
  }
  
#ifdef _influenceRadiusDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"****************** temp influence radii **************"<<endl;
  for(int i=0;i<numOfPtcls;i++) {
    dbVector& coords = ptcls[i].getCoords();
    logFile<<"PARTICLE "<<i<<"(";
    for(int j=0;j<coords.size();j++)
      logFile<<coords[j]<<" ";
    logFile<<") radii: ";
    dbVector& radii = ptcls[i].getRadii();
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<endl;
  }
#endif
  
}

/************************************************************************/
/************************************************************************/
// determine particle influence radii correspondingly to the distances
// to their neighbouring particles - different values in negative and
// positive coordinate direction
//
// controlling input parameter 
// -> minDirectionalPtcleSupport: minimum number of neighbour particles
//    for each coordinate directions (plus/minus); the number should 
//    not be larger than maximal possible support everywhere
// -> minDirectPtcleSuppReduction: support deficit in plus or minus 
//    direction is usually added to the opposite direction (plus/minus)
//    but can be reduced by this parameter
void ParticleDistribution::setPtcleRadsDistDepend2(InputFileData* InputData,
                                                   std::map<std::string,double>& modelData,
                                                   std::ofstream& logFile) {
  
  using namespace std;
  
  InputData->setValue("shapefunctionType",5.0);
  //InputData->setValue("windowFunctionType",1.0);
  InputData->setValue("plusMinusDirectionDependentRadius",1.0);

  int usedDims = (int)modelData["usedDimensions"];
  int minPtcleSupport = (int)InputData->getValue("minDirectionalPtcleSupport");
  int supportReduction = (int)InputData->getValue("minDirectPtcleSuppReduction");

  int numOfPtcls = particles.size();
  vector<Particle>& ptcls = particles;

  for(int i=0;i<numOfPtcls;i++) {

    dbVector& radii = ptcls[i].getRadii();
    radii.resize(usedDims*2);

  }

#ifdef _influenceRadiusDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"****** default influence zone determination ********"<<endl;
#endif

  dbMatrix delta = getKroneckerSymbol(usedDims);
  dbMatrix g(usedDims,dbVector(usedDims));

  g = delta;
  double searchAngle = 0;

  // set the search sector respectively angles

  double pi = 3.141592654;
  searchAngle = searchAngle*pi/180.0;

  double minAngle = searchAngle; // alpha_0
  double maxAngle = pi/6.0+searchAngle; // 30 degree + alpha_0

#ifdef _influenceRadiusDebugMode_
  logFile<<"----------------------------------------------------"<<endl;
  logFile<<"searchAngle: "<<minAngle<<" - "<<maxAngle
	 <<endl;
#endif

  /********************************************************************/
  // loop over all directions and particles to determine foreach
  // particle a default influence zone

  double distance;

  dbVector absDists(numOfPtcls);
  dbMatrix maxDists(numOfPtcls,dbVector(usedDims));
  dbMatrix dists(numOfPtcls,dbVector(usedDims));

  int round = 0;
  double tol = 1.0e-08;
  double absDist,angle,dirCos,weight;
  int minusCounter,plusCounter;
  intVector idxPlus(numOfPtcls);
  intVector idxMinus(numOfPtcls);

  dbVector weightedRadsPlus(numOfPtcls);
  dbVector weightedRadsMinus(numOfPtcls);


  for(int i=0;i<numOfPtcls;i++) {

#ifdef _influenceRadiusDebugMode_
    logFile<<"**************************************************"<<endl;
    logFile<<"PARTICLE "<<i<<" coords: "<<ptcls[i].getCoord(0)
	   <<" "<<ptcls[i].getCoord(1)<<" "<<ptcls[i].getCoord(2)<<endl;
#endif

    clearArray(absDists);

    // Loop over all particles and determine the distance to the
    // neighbouring particles.
    for(int j=0;j<numOfPtcls;j++) {

      // do not include current particle 'i' itself as a neighbour
      if(i == j)
	continue;

#ifdef _influenceRadiusDebugMode_
      logFile<<"********"<<endl;
      logFile<<"neighbour particle: "<<j<<" ";
#endif

      for(int k=0;k<usedDims;k++) {

	dists[j][k] =
	  ptcls[j].getCoord(k)-ptcls[i].getCoord(k);

	absDists[j] += pow(dists[j][k],2);


#ifdef _influenceRadiusDebugMode_
	logFile<<"distance: "<<dists[j][k]<<" maxDist: "<<maxDists[i][k]<<endl;
#endif

	if(maxDists[i][k] < fabs(dists[j][k]))
	  maxDists[i][k] = fabs(dists[j][k]);

      }

      absDists[j] = sqrt(absDists[j]);

#ifdef _influenceRadiusDebugMode_
      logFile<<" absDist: "<<absDists[j]<<endl;
#endif

    }


    round = 0;

    // determine separately for each direction the closest particles
    for(int dof=0;dof<usedDims;dof++) {

      // 10 degree*round + 30 degree + alpha_0; for each fail the search
      // angle is increased by 10 degree
      maxAngle = pi/18*round + pi/6.0 + searchAngle;

#ifdef _influenceRadiusDebugMode_
      logFile<<"--------------------------------------------------"<<endl;
      logFile<<"DOF "<<dof<<" round "<<round<<": "<<endl;
#endif

      plusCounter = 0;
      minusCounter = 0;

      for(int j=0;j<numOfPtcls;j++) {

	// do not include current particle 'i' itself as a neighbour
	if(i == j)
	  continue;

	dirCos = 0;

	// determine the directional cosine foreach search direction
	// separately
	switch(dof) {

	  // search direction 1
	case 0:

	  for(int k=0;k<usedDims;k++)

	    dirCos += g[0][k]*dists[j][k];

	  break;

	  // search direction 2
	case 1:

	  for(int k=0;k<usedDims;k++)

	    dirCos += g[1][k]*dists[j][k];

	  break;

	  // search direction 3
	case 2:

	  for(int k=0;k<usedDims;k++)

	    dirCos += g[2][k]*dists[j][k];

	  break;
	}


	if(absDists[j]> 0)
	  dirCos = dirCos/absDists[j];

	else
	  dirCos = 0;

	// determine the angle included by the search direction and
	// neighbour particles direction with respect to particle 'i'
	angle = acos(fabs(dirCos));

	if(angle < minAngle || angle > maxAngle) {

#ifdef _influenceRadiusDebugMode_
	  logFile<<"neigh ptcle "<<j
		 <<" dirCos "<<dirCos<<" angle: "<<angle<<" dists ";
	  for(int k=0;k<usedDims;k++)
	    logFile<<dists[j][k]<<" ";
	  logFile<<"dof "<<dof<<" SKIPPED"<<endl;
#endif


	  continue;

	}

#ifdef _influenceRadiusDebugMode_
	else {

	  logFile<<plusCounter+minusCounter<<".) neigh ptcle "<<j
		 <<" dirCos "<<dirCos<<" angle: "<<angle<<" dists ";
	  for(int k=0;k<usedDims;k++)
	    logFile<<dists[j][k]<<" ";
	  logFile<<"dof "<<dof<<endl;

	}
#endif

	// plus direction
	if(dirCos >= 0 && fabs(dists[j][dof]) > tol) {

	  idxPlus[plusCounter] = j;
	  weightedRadsPlus[plusCounter] = fabs(dists[j][dof]);

	  plusCounter++;
	}

	// minus direction
	else if(dirCos < 0 && fabs(dists[j][dof]) > tol) {

	  idxMinus[minusCounter] = j;
	  weightedRadsMinus[minusCounter] = fabs(dists[j][dof]);

	  minusCounter++;
	}

      }

      sortValuesIdx(weightedRadsPlus,idxPlus,0,plusCounter-1);
      sortValuesIdx(weightedRadsMinus,idxMinus,0,minusCounter-1);

      // remove redundant entries
      int m=0;

      for(int j=1;j<plusCounter;j++) {

	if(fabs(weightedRadsPlus[m]-weightedRadsPlus[j]) > tol) {
	  ++m;
	  weightedRadsPlus[m] = weightedRadsPlus[j];
	  idxPlus[m] = idxPlus[j];
	}

      }

      if(plusCounter > 0)
	plusCounter = m+1;

      m=0;

      for(int j=1;j<minusCounter;j++) {

	if(fabs(weightedRadsMinus[m]-weightedRadsMinus[j]) > tol) {
	  ++m;
	  weightedRadsMinus[m] = weightedRadsMinus[j];
	  idxMinus[m] = idxMinus[j];
	}

      }

      if(minusCounter > 0)
	minusCounter = m+1;

      // ----------------------------------------------------------------
      // set the array index of 'idxPlus' and 'idxMinus'

      int plusID,minusID;

      // available support in positive and negative coordinate direction
      // more or equal than minimum particle support
      if(plusCounter >= minPtcleSupport
	 && minusCounter >= minPtcleSupport) {

	plusID = minPtcleSupport-1;
	minusID = minPtcleSupport-1;

      }

      // available support in positive coordinate direction less
      // and in negative coordinate direction more or equal than
      // minimum particle support
      else if(plusCounter < minPtcleSupport) {

	int deficit = minPtcleSupport - plusCounter;
	deficit -= supportReduction;

	if(deficit < 0)
	  deficit = 0;

	minusCounter < minPtcleSupport+deficit ? minusID = minusCounter-1 :
	  minusID = minPtcleSupport+deficit-1;

	plusID = plusCounter-1;

      }

      // available support in negative coordinate direction more or equal
      // and in positive coordinate direction less than minimum particle
      // support
      else if(minusCounter < minPtcleSupport) {

	int deficit = minPtcleSupport - minusCounter;
	deficit -= supportReduction;

	if(deficit < 0)
	  deficit = 0;

	plusCounter < minPtcleSupport+deficit ? plusID = plusCounter-1 :
	  plusID = minPtcleSupport+deficit-1;

	minusID = minusCounter-1;

      }


#ifdef _influenceRadiusDebugMode_
      logFile<<"------------------------------------------------"<<endl;
      logFile<<"plusCounter = "<<plusCounter<<" minusCounter = "
	     <<minusCounter<<endl;
      logFile<<"plusID = "<<plusID<<" minusID = "<<minusID<<endl;
      for(int j=0;j<plusCounter;j++) {
	logFile<<j<<".) neigh ptcle: "<<idxPlus[j]
	       <<" weightedRadPlus["<<j<<"] = "
	       <<weightedRadsPlus[j]<<" dists: ";
	for(int k=0;k<usedDims;k++)
	  logFile<<dists[idxPlus[j]][k]<<" ";
	logFile<<endl;
      }
      for(int j=0;j<minusCounter;j++) {
	logFile<<j<<".) neigh ptcle: "<<idxMinus[j]
	       <<" weightedRadMinus["<<j<<"] = "
	       <<weightedRadsMinus[j]<<" dists: ";
	for(int k=0;k<usedDims;k++)
	  logFile<<dists[idxMinus[j]][k]<<" ";
	logFile<<endl;
      }
#endif


      // enough neighbours in plus and minus direction available
      if(plusID >= 0 && minusID >= 0) {

	dbVector& radii = ptcls[i].getRadii();

	// set positive and negative radius of particle 'i'
	if(radii[dof] < weightedRadsPlus[plusID])
	  radii[dof] = weightedRadsPlus[plusID];

	if(radii[usedDims+dof] < weightedRadsMinus[minusID])
	  radii[usedDims+dof] = weightedRadsMinus[minusID];

	// set radius of the closest positive and negative neighbours
	for(int k=0;k<=plusID;k++) {

	  dbVector& radii2 = ptcls[idxPlus[k]].getRadii();

	  if(radii2[usedDims+dof] < weightedRadsPlus[k])
	    radii2[usedDims+dof] = weightedRadsPlus[k];

	  dbVector& radii3 = ptcls[idxMinus[k]].getRadii();

	  if(radii3[dof] < weightedRadsMinus[k])
	    radii3[dof] = weightedRadsMinus[k];

	}


	round = 0;

#ifdef _influenceRadiusDebugMode_
	logFile<<"ptcle "<<i<<" set rads: ";
	for(int r=0;r<radii.size();r++)
	  logFile<<radii[r]<<" ";
	logFile<<endl;
	for(int k=0;k<=plusID;k++) {
	  dbVector& radii2 = ptcls[idxPlus[k]].getRadii();
	  logFile<<"ptcle "<<idxPlus[k]<<" set rads: ";
	  for(int r=0;r<radii2.size();r++)
	    logFile<<radii2[r]<<" ";
	  logFile<<endl;
	}
	for(int k=0;k<=minusID;k++) {
	  dbVector& radii3 = ptcls[idxMinus[k]].getRadii();
	  logFile<<"ptcle "<<idxMinus[k]<<" set rads: ";
	  for(int r=0;r<radii3.size();r++)
	    logFile<<radii3[r]<<" ";
	  logFile<<endl;
	}
#endif

      }

      // enough neighbours only in plus direction available
      else if(plusID >= 0 && minusID < 0) {

	// set radius of particle 'i'
	dbVector& radii = ptcls[i].getRadii();

	if(radii[dof] < weightedRadsPlus[plusID])
	  radii[dof] = weightedRadsPlus[plusID];

	// set radius of the closest positive neighbours
	for(int k=0;k<=plusID;k++) {

	  dbVector& radii2 = ptcls[idxPlus[k]].getRadii();

	  if(radii2[usedDims+dof] < weightedRadsPlus[k])
	    radii2[usedDims+dof] = weightedRadsPlus[k];

	}

	round = 0;

#ifdef _influenceRadiusDebugMode_
	logFile<<"ptcle "<<i<<" set rads: ";
	for(int r=0;r<radii.size();r++)
	  logFile<<radii[r]<<" ";
	logFile<<endl;
	for(int k=0;k<=plusID;k++) {
	  dbVector& radii2 = ptcls[idxPlus[k]].getRadii();
	  logFile<<"ptcle "<<idxPlus[k]<<" set rads: ";
	  for(int r=0;r<radii2.size();r++)
	    logFile<<radii2[r]<<" ";
	  logFile<<endl;
	}
#endif

      }

      // enough neighbours only in minus direction available
      else if(plusID < 0 && minusID >= 0) {

	// set radius of particle 'i'
	dbVector& radii = ptcls[i].getRadii();

	if(radii[usedDims+dof] < weightedRadsMinus[minusID])
	  radii[usedDims+dof] = weightedRadsMinus[minusID];

	// set radius of the closest positive and negative neighbours
	for(int k=0;k<=minusID;k++) {

	  dbVector& radii2 = ptcls[idxMinus[k]].getRadii();

	  if(radii2[dof] < weightedRadsMinus[k])
	    radii2[dof] = weightedRadsMinus[k];

	}

	round = 0;

#ifdef _influenceRadiusDebugMode_
	logFile<<"ptcle "<<i<<" set rads: ";
	for(int r=0;r<radii.size();r++)
	  logFile<<radii[r]<<" ";
	logFile<<endl;
	for(int k=0;k<=minusID;k++) {
	  dbVector& radii2 = ptcls[idxMinus[k]].getRadii();
	  logFile<<"ptcle "<<idxMinus[k]<<" set rads: ";
	  for(int r=0;r<radii2.size();r++)
	    logFile<<radii2[r]<<" ";
	  logFile<<endl;
	}
#endif

      }

      // no neighbour found -> search angle is increased by 10 degree
      // and the search algorithm for current DOF is executed again
      else if(maxAngle < minAngle + pi*15.0/18.0) {

	round += 1;
	dof -= 1;

      }
      else {
	logFile<<"In ParticleDistribution::setPtcleRadsDistDepend particle\n "
	       <<i<<" DOF "<<dof<<" has not enough neighbour particles!"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

    }

  }


#ifdef _influenceRadiusDebugMode_
  logFile<<"#################################################"<<endl;
  logFile<<"************* temp influence radii **************"<<endl;
  for(int i=0;i<numOfPtcls;i++) {
    logFile<<"PARTICLE "<<i<<": radii ";
    dbVector& radii = ptcls[i].getRadii();
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<"coords: ";
    for(int j=0;j<usedDims;j++)
      logFile<<ptcls[i].getCoord(j)<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine the influence radius of all particles support dependent.
void ParticleDistribution::setPtcleRadsSupportDepend(InputFileData* InputData,
                                                     std::vector<GaussPoint>& gPoints,
                                                     std::vector<GaussPoint>& bGPoints,
                                                     std::map<std::string,double>& modelData,
                                                     std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];

  int minPtcleSupport = (int)InputData->getValue("minParticleSupport");

  int numOfPtcls = particles.size();
  vector<Particle>& ptcls = particles;


  int counter;
  double oldRadius;

  intVector radiiIdx(numOfPtcls);
  dbVector radii(numOfPtcls);
  dbMatrix dists(numOfPtcls,dbVector(usedDims));

  // Determine the distances of a portion of particles to its neighbours
  // and sort them
  for(int i=0;i<exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];

    clearArray(radii);

#ifdef _influenceRadiusDebugMode_
    logFile<<"PARTICLE: "<<ptcle<<":"<<endl;
#endif

    // Loop over all particles.
    for(int j=0;j<numOfPtcls;j++) {

#ifdef _influenceRadiusDebugMode_
      logFile<<"********"<<endl;
      logFile<<"neighbour particle: "<<j<<" ";
#endif

      radiiIdx[j] = j;

      for(int dof=0;dof<usedDims;dof++) {

	dists[j][dof] =
	  ptcls[i].getCoord(dof)-ptcls[j].getCoord(dof);

	radii[j] += pow(dists[j][dof],2);

#ifdef _influenceRadiusDebugMode_
	logFile<<"DOF "<<dof<<" distance = "<<dists[j][dof]<<"| ";
#endif

      }

      radii[j] = sqrt(radii[j]);

#ifdef _influenceRadiusDebugMode_
      logFile<<" radius "<<radii[j]<<endl;
#endif

    }

    // sort the distances
    sortValuesIdx(radii,radiiIdx,0,numOfPtcls-1);

#ifdef _influenceRadiusDebugMode_
    logFile<<"*************** all particle distances *************"<<endl;
    for(int j=0;j<numOfPtcls;j++) {
      logFile<<j<<".) neigh ptcle "<<radiiIdx[j]<<" distances: ";
      for(int dof=0;dof<usedDims;dof++)
	logFile<<dists[radiiIdx[j]][dof]<<" ";
      logFile<<radii[j]<<endl;
    }
#endif

    counter = 0;
    oldRadius = 1.0e+16;

    // set the necessary influence radius of the closest particles so
    // that minimum support is given
    for(int j=0;j<numOfPtcls;j++) {

      if(j >= minPtcleSupport && radii[j] > oldRadius)

	break;

      dbVector& infRadii = ptcls[radiiIdx[j]].getRadii();

      for(int k=0;k<infRadii.size();k++)

	if(infRadii[k] < fabs(dists[radiiIdx[j]][k]))
	  infRadii[k] = fabs(dists[radiiIdx[j]][k]);

      oldRadius = radii[j];

    }

  }

#ifdef _influenceRadiusDebugMode_
  logFile<<"******* tmp influenceRadii(with ptcle support) ********"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"Particle "<<i<<"(old: "<<oldIdx[i]<<"): "
	   <<particles[i].getRadius(0)<<" "
	   <<particles[i].getRadius(1)<<" "
	   <<particles[i].getRadius(2)<<endl;
#endif

  /*********************************************************************/
  // Determine the distances of the local Gauss points to their neighbour
  // particles and sort them.

  for(int i=0;i<gPoints.size();i++) {

    clearArray(radii);

#ifdef _influenceRadiusDebugMode_
    logFile<<"GAUSS POINT: "<<i<<":"<<endl;
#endif

    // Loop over all particles.
    for(int j=0;j<numOfPtcls;j++) {

#ifdef _influenceRadiusDebugMode_
      logFile<<"********"<<endl;
      logFile<<"neighbour particle: "<<j<<" ";
#endif

      radiiIdx[j] = j;

      for(int dof=0;dof<usedDims;dof++) {

	dists[j][dof] =
	  gPoints[i].getCoord(dof)-ptcls[j].getCoord(dof);

	radii[j] += pow(dists[j][dof],2);

#ifdef _influenceRadiusDebugMode_
	logFile<<"DOF "<<dof<<" distance = "<<dists[j][dof]<<"| ";
#endif

      }

      radii[j] = sqrt(radii[j]);

#ifdef _influenceRadiusDebugMode_
      logFile<<" radius "<<radii[j]<<endl;
#endif

    }

    // sort the distances
    sortValuesIdx(radii,radiiIdx,0,numOfPtcls-1);

#ifdef _influenceRadiusDebugMode_
    logFile<<"******** all Gauss pts - particle distances ********"<<endl;
    for(int j=0;j<particlesNum;j++) {
      logFile<<j<<".) neigh ptcle "<<radiiIdx[j]<<" distances: ";
      for(int dof=0;dof<usedDims;dof++)
	logFile<<dists[radiiIdx[j]][dof]<<" ";
      logFile<<radii[j]<<endl;
    }
#endif

#ifdef _influenceRadiusDebugMode_
    logFile<<"supporting ptcls : ";
#endif

    oldRadius = 1.0e+16;

    // set the necessary influence radius of the closest particles so
    // that minimum support is given
    for(int j=0;j<numOfPtcls;j++) {

      if(j >= minPtcleSupport && radii[j] > oldRadius)

	break;

      dbVector& infRadii = ptcls[radiiIdx[j]].getRadii();

      for(int k=0;k<infRadii.size();k++)

	if(infRadii[k] < fabs(dists[radiiIdx[j]][k]))
	  infRadii[k] = fabs(dists[radiiIdx[j]][k]);

      oldRadius = radii[j];

#ifdef _influenceRadiusDebugMode_
      logFile<<radiiIdx[j]<<" ";
#endif

    }

#ifdef _influenceRadiusDebugMode_
    logFile<<endl;
#endif

  }

#ifdef _influenceRadiusDebugMode_
  logFile<<"******* tmp influenceRadii(with gPts support) ********"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"Particle "<<i<<"(old: "<<oldIdx[i]<<"): "
	   <<particles[i].getRadius(0)<<" "
	   <<particles[i].getRadius(1)<<" "
	   <<particles[i].getRadius(2)<<endl;
#endif

  /*********************************************************************/
  // Determine the distances of all boundary Gauss points to their
  // neighbour particles and sort them.

  for(int i=0;i<bGPoints.size();i++) {

    clearArray(radii);

#ifdef _influenceRadiusDebugMode_
    logFile<<"BOUND GAUSS POINT: "<<i<<":"<<endl;
#endif

    // Loop over all particles.
    for(int j=0;j<numOfPtcls;j++) {

#ifdef _influenceRadiusDebugMode_
      logFile<<"********"<<endl;
      logFile<<"neighbour particle: "<<j<<" ";
#endif

      radiiIdx[j] = j;

      for(int dof=0;dof<usedDims;dof++) {

	dists[j][dof] =
	  bGPoints[i].getCoord(dof)-ptcls[j].getCoord(dof);

	radii[j] += pow(dists[j][dof],2);

#ifdef _influenceRadiusDebugMode_
	logFile<<"DOF "<<dof<<" distance = "<<dists[j][dof]<<"| ";
#endif

      }

      radii[j] = sqrt(radii[j]);

#ifdef _influenceRadiusDebugMode_
      logFile<<" radius "<<radii[j]<<endl;
#endif

    }

    // sort the distances
    sortValuesIdx(radii,radiiIdx,0,numOfPtcls-1);

#ifdef _influenceRadiusDebugMode_
    logFile<<"******** all Gauss pt - particle distances *********"<<endl;
    for(int j=0;j<particlesNum;j++) {
      logFile<<j<<".) neigh ptcle "<<radiiIdx[j]<<" distances: ";
      for(int dof=0;dof<usedDims;dof++)
	logFile<<dists[radiiIdx[j]][dof]<<" ";
      logFile<<radii[j]<<endl;
    }
#endif

#ifdef _influenceRadiusDebugMode_
    logFile<<"supporting ptcls : ";
#endif

    oldRadius = 1.0e+16;

    // set the necessary influence radius of the closest particles so
    // that minimum support is given
    for(int j=0;j<numOfPtcls;j++) {

      if(j >= minPtcleSupport && radii[j] > oldRadius)

	break;

      dbVector& infRadii = ptcls[radiiIdx[j]].getRadii();

      for(int k=0;k<infRadii.size();k++)

	if(infRadii[k] < fabs(dists[radiiIdx[j]][k]))
	  infRadii[k] = fabs(dists[radiiIdx[j]][k]);

      oldRadius = radii[j];

#ifdef _influenceRadiusDebugMode_
      logFile<<radiiIdx[j]<<" ";
#endif
    }

#ifdef _influenceRadiusDebugMode_
    logFile<<endl;
#endif

  }

#ifdef _influenceRadiusDebugMode_
  logFile<<"******* tmp influenceRadii(with bGPts support) ********"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"Particle "<<i<<"(old: "<<oldIdx[i]<<"): "
	   <<particles[i].getRadius(0)<<" "
	   <<particles[i].getRadius(1)<<" "
	   <<particles[i].getRadius(2)<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Post-process all particles' influence radii.
void ParticleDistribution::postProcInfRads(InputFileData* InputData,
					   std::vector<Particle>& ptcls,
                                           std::map<std::string,double>& modelData,
                                           std::ofstream& logFile) {


  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int minPtcleSupport = (int)InputData->getValue("minDirectionalPtcleSupport");
  bool plusMinusDependent =
    (bool)InputData->getValue("plusMinusDirectionDependentRadius");
  double multiplier = InputData->getValue("influenceRadiusMultiplier");

  int numOfPtcls = ptcls.size();

  // multiply all determined influence radii by a constant

  for(int i=0;i<numOfPtcls;i++) {

    dbVector& radii = ptcls[i].getRadii();

    for(int k=0;k<radii.size();k++)

      radii[k] *= multiplier;

  }

#ifdef _influenceRadiusDebugMode_
  logFile<<"multiplier="<<multiplier<<endl;
  logFile<<"*********** factored global influenceRadii ***********"<<endl;
  for(int i=0;i<numOfPtcls;i++) {
    dbVector& radii = ptcls[i].getRadii();
    logFile<<"Particle "<<i<<": ";
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<endl;
  }
#endif

  /**********************************************************************/
  // Check if there 'zero-radii'
  double minRadius;
  dbVector& dummy = ptcls[0].getRadii();
  vector<bool> zeros(dummy.size());

  double radTol = 1.0e-07;

  int mode = (int)InputData->getValue("radiusDeterminationAlgorithm");
  double radiusRatio = InputData->getValue("positiveNegativeRadiusRatio");
  double minRadiusTol = InputData->getValue("mininumInfluenceRadiusTolerance");

#ifdef _influenceRadiusDebugMode_
  logFile<<"********* global influenceRadii zeros check **********"<<endl;
#endif

  // Loop over all particles influence radii.
  for(int i=0;i<numOfPtcls;i++) {

    dbVector& radii = ptcls[i].getRadii();

    minRadius = 0;
    zeros.assign(zeros.size(),false);

    // Loop over all three coordinate directions to determine the zeros.
    for(int j=0;j<radii.size();j++) {

      if(radii[j] < minRadiusTol)
	zeros[j] = true;

      else if(radii[j] < minRadius || minRadius == 0)
	minRadius = radii[j];

    }

    if(minRadius == 0) {
      logFile<<"Particle "<<i<<" has got a zero influence zone in \n"
	     <<"(ParticleDistribution::postProcInfRads)!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

#ifdef _influenceRadiusDebugMode_
    logFile<<"Particle "<<i<<": minRadius = "
	   <<minRadius<<" radii: ";
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" -> "<<zeros[j]<<" | ";
    logFile<<endl;
#endif

    if(mode == 5 || mode == 6) {

      for(int j=0;j<usedDims;j++) {

	// radius of positive side zero but negative not
	if(zeros[j] && !zeros[j+usedDims])

	  radii[j] = radii[j+usedDims];

	// radius of negative side zero but positive not
	else if(zeros[j+usedDims] && !zeros[j])

	  radii[j+usedDims] = radii[j];

      }

    }

    else if(mode != 5 && mode != 6) {

      for(int j=0;j<radii.size();j++)

	if(zeros[j])
	  radii[j] = minRadius;

    }

  }



#ifdef _influenceRadiusDebugMode_
  logFile<<"*********** global influenceRadii ******************"<<endl;
  for(int i=0;i<numOfPtcls;i++) {
    dbVector& radii = ptcls[i].getRadii();
    logFile<<"Particle "<<i<<": ";
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<endl;
  }
#endif


  /**********************************************************************/
  // Set all particles the same influence radius for a certain direction.
  double x1Radius = InputData->getValue("x1InfluenceRadius");
  double x2Radius = InputData->getValue("x2InfluenceRadius");
  double x3Radius = InputData->getValue("x3InfluenceRadius");

  if(x1Radius!= 0 && !plusMinusDependent) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();
      radii[0] = x1Radius;

    }

  }
  else if(x1Radius != 0 && plusMinusDependent) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();
      radii[0] = x1Radius;
      radii[usedDims] = x1Radius;
    }

  }

  if(x2Radius != 0 && !plusMinusDependent) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();
      radii[1] = x2Radius;

    }

  }
  else if(x2Radius != 0 && plusMinusDependent) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();
      radii[1] = x2Radius;
      radii[usedDims+1] = x2Radius;
    }

  }

  if(x3Radius != 0 && !plusMinusDependent) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();
      radii[2] = x3Radius;

    }

  }
  else if(x3Radius != 0 && plusMinusDependent) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();
      radii[2] = x3Radius;
      radii[usedDims+2] = x3Radius;
    }

  }


  /**********************************************************************/
  // Determine minimum and maximum influence radius.
  minRadius = 1.0e+16;
  double maxRadius = 0;

  for(int i=0;i<numOfPtcls;i++) {

    dbVector& radii = ptcls[i].getRadii();

    for(int j=0;j<radii.size();j++) {

      if(radii[j] > maxRadius)
	maxRadius = radii[j];

      if(radii[j] < minRadius && radii[j] > DBL_EPSILON)
	minRadius = radii[j];

    }

  }

  /**********************************************************************/
  // Set for all particles the same influence radius
  int uniformRadii = (int)InputData->getValue("uniformInfluenceRadii");

  if(uniformRadii == 1) {

    for(int i=0;i<numOfPtcls;i++) {

      dbVector& radii = ptcls[i].getRadii();

      for(int j=0;j<radii.size();j++)

	radii[j] = maxRadius;

    }

    minRadius = maxRadius;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(rank == 0) {
    logFile<<"influence radii "<<minRadius<<" - "<<maxRadius<<endl;
    cout<<"influence radii "<<minRadius<<" - "<<maxRadius<<endl;
  }

#ifdef _geometryDebugMode_
  logFile<<"######################################################"<<endl;
  logFile<<"*********** stored influenceRadii ******************"<<endl;
  for(int i=0;i<numOfPtcls;i++) {
    logFile<<"Particle "<<i<<"(old: "<<oldIdx[i]<<"): ";
    dbVector& radii = ptcls[i].getRadii();
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine for all particles their neighbour particles.
void ParticleDistribution::setPtclePtclsConn(InputFileData* InputData,
                                             std::map<std::string,double>& modelData,
                                             std::ofstream& logFile) {

  using namespace std;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //   int startIdx = exclusiveLocalPtcls.front();
  //   int endIdx = exclusiveLocalPtcls.back()+1;

  bool supported;
  int supportSize,minPtcleSupport;
  int localMinSupport;

  dbVector radius(3);

  int localMaxSupport = 0;

#ifdef _geometryDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"********** calculation of ptcle-ptcle connlist **********"<<endl;
#elif defined _forceVecModificationDebugMode_
  logFile<<"#########################################################"<<endl;
  logFile<<"********** calculation of ptcle-ptcle connlist **********"<<endl;
#endif

  // Loop over a local portion of particles and create for each a support list.
  for(int i=0;i<exclusiveLocalPtcls.size();i++) {
    int& ptcle = exclusiveLocalPtcls[i];

    intVector& suppPtcls = particles[ptcle].getSupportPtcls();
    suppPtcls.resize(localMaxSupport);

    supportSize=0;

#ifdef _geometryDebugMode_
    logFile<<"PARTICLE "<<i<<"**********"<<endl;
#endif

    // Loop over all neighbour particles and determine which support particle
    // 'i'.
    for(int j=0;j<particlesNum;j++) {

      // Check if current particle 'i' is supported by particle 'j'.
      if(particles[j].querySupported(InputData,particles[ptcle].getCoords(),
				     modelData,logFile)) {

	if(supportSize < suppPtcls.size())
	  suppPtcls[supportSize] = j;

	else
	  suppPtcls.push_back(j);

	supportSize++;
      }

#ifdef _geometryDebugMode_
      logFile<<"neighbour "<<j<<": coords: "<<particles[j].getCoord(0)
	     <<" "<<particles[j].getCoord(1)
	     <<" "<<particles[j].getCoord(2)<<endl;
      logFile<<"radii: "<<particles[j].getRadius(0)
	     <<" "<<particles[j].getRadius(1)
	     <<" "<<particles[j].getRadius(2)
	     <<"; dists: "
	     <<fabs(particles[ptcle].getCoord(0)-particles[j].getCoord(0))
	     <<" "<<fabs(particles[ptcle].getCoord(1)-particles[j].getCoord(1))
             <<" "<<fabs(particles[ptcle].getCoord(2)-particles[j].getCoord(2))
	     <<"; supportSize = "<<supportSize<<endl;
#endif

    }

    resizeArray(suppPtcls,supportSize);


    // -------------------------------------------------------------------

    if(supportSize > localMaxSupport)
      localMaxSupport = supportSize;

    if(i == 0)
      localMinSupport = supportSize;

    else if(supportSize < localMinSupport)
      localMinSupport = supportSize;

#ifdef _geometryDebugMode_
    logFile<<"localMinSupport = "<<localMinSupport<<endl;
    logFile<<"localMaxSupport = "<<localMaxSupport<<endl;
#endif

  }

  /**********************************************************************/
  // broadcast the locally determined particle support lists to the other
  // processes

  int globalSize;

  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();

    globalSize = suppPtcls.size();

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&globalSize,1,MPI_INT,ptcleRootList[i],
	      MPI_COMM_WORLD);

    if(rank != ptcleRootList[i])
      suppPtcls.resize(globalSize);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&suppPtcls[0],globalSize,MPI_INT,ptcleRootList[i],
	      MPI_COMM_WORLD);
  }

  // Determine global maximum number of supporting particles for a
  // particle.
  MPI_Allreduce(&localMaxSupport,&maxPtcleSupport,1,
		MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  // Determine global minimum number of supporting particles for a
  // particle.
  MPI_Allreduce(&localMinSupport,&minPtcleSupport,1,
		MPI_INT,MPI_MIN,MPI_COMM_WORLD);


  logFile<<"minPtcleSupport = "<<minPtcleSupport<<endl;
  logFile<<"maxPtcleSupport = "<<maxPtcleSupport<<endl;

  if(rank == 0)
    cout<<"ptcleSupport: "<<minPtcleSupport<<" - "<<maxPtcleSupport<<endl;

#ifdef _geometryDebugMode_
  logFile<<"********** particle - particle support list **********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    logFile<<"Ptcle "<<i<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
#elif defined _forceVecModificationDebugMode_
  logFile<<"********** particle - particle support list **********"<<endl;
  for(int i=0;i<particlesNum;i++) {
    intVector& suppPtcls = particles[i].getSupportPtcls();
    logFile<<"Ptcle "<<i<<": ";
    for(int j=0;j<suppPtcls.size();j++)
      logFile<<suppPtcls[j]<<" ";
    logFile<<endl;
  }
#endif

}

/**********************************************************************/
/**********************************************************************/
// Determine for given point their neighbour particles.
void ParticleDistribution::setPointPtclsConn(InputFileData* InputData,
                                             dbVector& coords,
                                             int& supportSize,
                                             intVector& supportingPtcls,
                                             std::map<std::string,double>& modelData,
                                             std::ofstream& logFile) {

  using namespace std;

  supportSize=0;
  supportingPtcls.resize(maxPtcleSupport);

  // Loop over all particles.
  for(int i=0;i<particlesNum;i++) {

#ifdef _supportDebugMode_
    logFile<<"ptcle "<<i<<endl;
#endif

    // Check if current point is supported.
    if(particles[i].querySupported(InputData,coords,modelData,logFile)) {

      if(supportSize < supportingPtcls.size())
	supportingPtcls[supportSize] = i;

      else
	supportingPtcls.push_back(i);

      supportSize++;

#ifdef _supportDebugMode_
      logFile<<"support"<<endl;
#endif
    }


#ifdef _supportDebugMode_
    logFile<<endl;
#endif

  }

  resizeArray(supportingPtcls,supportSize);
}

/**********************************************************************/
/**********************************************************************/
// Determine for given point their neighbour particles and the minimum
// distance with the supporting particle.
void ParticleDistribution::setPointPtclsConn(InputFileData* InputData,
                                             dbVector& coords,
                                             int& supportSize,
                                             intVector& supportingPtcls,
                                             double& localMinPtcleDist,
                                             std::map<std::string,double>& modelData,
                                             std::ofstream& logFile) {

  using namespace std;

  dbVector radius(3), ptcleDistVec;

  double sum=0,ptcleDist, minPtcleDist = 99999;
  bool rad_set = false;

  bool duplexPtcle;

  supportSize=0;
  supportingPtcls.resize(maxPtcleSupport);

  // Loop over all particles.
  for(int i=0;i<particlesNum;i++) {

    duplexPtcle = false;

    // Check if current point is supported.
    if(particles[i].querySupported(InputData,coords,ptcleDist,modelData,logFile)) {

      if(supportSize < supportingPtcls.size())
	supportingPtcls[supportSize] = i;

      else
	supportingPtcls.push_back(i);


      if(ptcleDist != 0)
	ptcleDistVec.push_back(ptcleDist);

      // ptcle 'i' is not the particle itself which needs support
      if(fabs(particles[i].getCoord(0) - coords[0]) > DBL_EPSILON &&
	 fabs(particles[i].getCoord(1) - coords[1]) > DBL_EPSILON &&
	 fabs(particles[i].getCoord(2) - coords[2]) > DBL_EPSILON) {

	// particle 'i' closer than the ones already checked
	if(minPtcleDist > ptcleDist && ptcleDist > DBL_EPSILON)
	  
	  minPtcleDist = ptcleDist;
	
      }
      else if(!duplexPtcle)

	duplexPtcle = true;

      else {
	logFile<<"In ParticleDistribution::setPointPtclsConn particle "
	       <<i<<"(GiD node: "<<oldIdx[i]+1<<") with coords "
	       <<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<"\n"
	       <<"appears with different nodal ID more than once in\n" 
	       <<"'mesh.dat'. Check GiD model e.g. for duplex surfaces"
	       <<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    

      supportSize++;
    }

  }

  for (int s=0;s<ptcleDistVec.size();s++ ){
    sum += ptcleDistVec[s];
  }

  //localMinPtcleDist = minPtcleDist;
  localMinPtcleDist = sum/ptcleDistVec.size();

  resizeArray(supportingPtcls,supportSize);
}

/************************************************************************/
/************************************************************************/
// Update the particle coordinates.
void ParticleDistribution::updateParticleCoords(InputFileData* InputData,
                                                dbVector& allDOF,
                                                std::map<std::string,double>& modelData,
                                                std::ofstream& logFile) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int usedDims = (int)modelData["usedDimensions"];
  int numOfDOF = particlesNum*usedDOF;

  if(allDOF.size() < (unsigned int)numOfDOF) {
    logFile <<"To less degrees of freedom increments to update "
	    <<"particle coordinates in "
	    <<"RKPMDiscretising::updateParticleCoords!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifdef _calculationDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"************* old particle coordinates *************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Particle "<<i<<": ";
    for(int j=0;j<usedDims;j++)
      logFile<<particles[i].getCoord(j)<<", ";
    logFile<<endl;
  }
#endif

  // Loop over all particles to update their coordinates.
  for(int i=0;i<particlesNum;i++) {

    // Loop over all particle dimensions.
    for(int j=0;j<usedDims;j++) {
      double& coord = particles[i].getCoord(j);
      coord += allDOF[i*usedDOF+j];
    }

  }

#ifdef _calculationDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"********** updated particle coordinates ************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    logFile<<"Particle "<<i<<": ";
    for(int j=0;j<usedDims;j++)
      logFile<<particles[i].getCoord(j)<<", ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Set all particle degree of freedom to null.
void ParticleDistribution::clearAllPtcleDOF(InputFileData* InputData,
                                            std::map<std::string,double>& calcData,
                                            std::map<std::string,double>& modelData,
                                            std::ofstream& logFile) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];
  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  double value = 0;

  // Loop over all particles
  for(int i=0;i<particlesNum;i++)

    clearArray(particles[i].getDOF());

#ifdef _calculationDebugMode_
  logFile<<"####################################################"<<endl;
  logFile<<"************** cleared particle DOF ****************"<<endl;
  for(int i=0;i<particlesNum;i++) {
    dbVector& DOF = particles[i].getDOF();
    logFile<<"Particle "<<i<<": ";
    for(int j=0;j<DOF.size();j++)
      logFile<<DOF[j]<<", ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Create a list of the exclusively local particles from all processors 
// in consecutive ordering and establish the mapping of global ptcle ID to
// to this sequential ordering.
//
// allExclLocalPtcls: consecutive ordering of all exclusively local 
//                particles
// allGlobalExclLocalPtcleIdx: mapping of global particle IDs to the all 
//                             processor's consecutive exclusively local 
//                             particle ordering
// localPtclsStartIdx: start position of processor's exclusively local 
//                     particle portion within allExclLocalPtcls
void ParticleDistribution::setAllExclLocalPtcleOrder(InputFileData* InputData,
                                                     std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  allExclLocalPtclsStartIdx.resize(size);
  allExclLocalPtcls.resize(particlesNum);
  allGlobalExclLocalPtcleIdx.resize(particlesNum);


  int sendCount = exclusiveLocalPtcls.size();
  intVector recvCounts(size);

  MPI_Allgather(&sendCount,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		MPI_COMM_WORLD);

  allExclLocalPtclsStartIdx[0] = 0;

  for(int i=1;i<size;i++)

    allExclLocalPtclsStartIdx[i] =
      allExclLocalPtclsStartIdx[i-1]+recvCounts[i-1];


  MPI_Allgatherv(&exclusiveLocalPtcls[0],recvCounts[rank],MPI_INT,
		 &allExclLocalPtcls[0],&recvCounts[0],&
		 allExclLocalPtclsStartIdx[0],MPI_INT,MPI_COMM_WORLD);

  // ---------------------------------------------------------------------
  // idx linkage of all local consecutive particles to global particle ID

  for(int i=0;i<allExclLocalPtcls.size();i++)

    allGlobalExclLocalPtcleIdx[allExclLocalPtcls[i]] = i;


}

/************************************************************************/
/************************************************************************/
// Merge the local portions of entries stored at particles
void ParticleDistribution::mergeLocalPtcleVectors(InputFileData* InputData,
                                                  dbVector& localVec,int ptcleEntries,
                                                  dbVector& globalVec,
                                                  std::map<std::string,double>& calcData,
                                                  std::map<std::string,double>& modelData,
                                                  std::ofstream& logFile) {

  using namespace std;

  intVector& allGlobalLocalIdx = allGlobalExclLocalPtcleIdx;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  intVector recvCounts(size);
  intVector displs(size);
    
  int sendCount = localVec.size();

  MPI_Allgather(&sendCount,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		MPI_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];


  dbVector tmpVec(particlesNum*ptcleEntries);

  MPI_Allgatherv(&localVec[0],recvCounts[rank],MPI_DOUBLE,
		 &tmpVec[0],&recvCounts[0],&displs[0],
		 MPI_DOUBLE,MPI_COMM_WORLD);

#ifdef _geometryDebugMode_
  logFile<<"******************** local vector ********************"<<endl;
  for(int i=0,k=0;i<localVec.size();i+=ptcleEntries,k++) {
    logFile<<k<<".) ";
    for(int j=0;j<ptcleEntries;j++)
      logFile<<localVec[i+j]<<" ";
    logFile<<endl;
  }
  logFile<<"**************** tmp global vector *******************"<<endl;
  for(int i=0,k=0;i<tmpVec.size();i+=ptcleEntries,k++) {
    logFile<<k<<".) ";
    for(int j=0;j<ptcleEntries;j++)
      logFile<<tmpVec[i+j]<<" ";
    logFile<<endl;
  }
#endif

  resizeArray(localVec,0);


  // ---------------------------------------------------------------------

  if(globalVec.size() < particlesNum*ptcleEntries)

    globalVec.resize(particlesNum*ptcleEntries);

  clearArray(globalVec);


  // Loop over all global particles to put it in the right order.
  for(int i=0;i<allGlobalLocalIdx.size();i++) {

    // Loop over all entries of a particle.
    for(int k=0;k<ptcleEntries;k++)

      globalVec[i*ptcleEntries+k] =
	tmpVec[allGlobalLocalIdx[i]*ptcleEntries+k];

  }

#ifdef _geometryDebugMode_
  logFile<<"******************* global vector ********************"<<endl;
  for(int i=0,k=0;i<globalVec.size();i+=ptcleEntries,k++) {
    logFile<<k<<".) ";
    for(int j=0;j<ptcleEntries;j++)
      logFile<<globalVec[i+j]<<" ";
    logFile<<endl;
  }
#endif

}

/************************************************************************/
/************************************************************************/
// Deallocate all class arrays of particles which are not locally needed.
void ParticleDistribution::deleteNonlocalPtcls(InputFileData* InputData,
                                               std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  intVector tmpVec(globalLocalPtcleIdx);

  // Loop over the local portion of particles and determine those
  // local particles which are supported by non-local particles
  for(int i=0;i<ptcleRootList.size();i++) {

    if(ptcleRootList[i] == rank) {

      intVector& sPtcls = particles[i].getSupportPtcls();

      // Loop over all supporting particles of current particle point.
      for(int j=0;j<sPtcls.size();j++)

	if(globalLocalPtcleIdx[sPtcls[j]] == -1)
	  tmpVec[sPtcls[j]] = -2;

    }

  }

#ifdef _geometryDebugMode_
  logFile<<"******  global->local ptcle indices to delete ********"<<endl;
  for(int i=0;i<tmpVec.size();i++)
    logFile<<"global ptcle number: "<<i<<" local one: "
	   <<tmpVec[i]<<endl;
#endif

  allocateArray(extendedLocalGlobalPtcls,particlesNum);
  int m=0;

  for(int i=0;i<particlesNum;i++) {

    // locally not needed particles
    if(tmpVec[i] == -1) {

      dbVector& coords = particles[i].getCoords();
      resizeArray(coords,0);

      dbVector& influenceRadii = particles[i].getRadii();
      resizeArray(influenceRadii,0);

      dbVector& degreesOfFreedom = particles[i].getDOF();
      resizeArray(degreesOfFreedom,0);

      intVector& supportingPtcls = particles[i].getSupportPtcls();
      resizeArray(supportingPtcls,0);

      intVector& influencingSpheres = particles[i].getInflSpheres();
      resizeArray(influencingSpheres,0);

      dbVector& shapeFuncs = particles[i].getShapeFuncs();
      dbMatrix& firstDerivShapeFuncs = particles[i].getFirstDerivShapes();
      dbMatrix& secondDerivShapeFuncs = particles[i].getSecondDerivShapes();
      resizeArray(shapeFuncs,0);
      resizeArray(firstDerivShapeFuncs,0);
      resizeArray(secondDerivShapeFuncs,0);

      intVector& elems = particles[i].getElems();
      resizeArray(elems,0);
    }
    // indirectly locally needed particles
    else if(tmpVec[i] == -2) {

      extendedLocalGlobalPtcls[m] = i;
      m++;

      dbVector& coords = particles[i].getCoords();
      resizeArray(coords,0);

      dbVector& influenceRadii = particles[i].getRadii();
      resizeArray(influenceRadii,0);

      intVector& supportingPtcls = particles[i].getSupportPtcls();
      resizeArray(supportingPtcls,0);

      intVector& influencingSpheres = particles[i].getInflSpheres();
      resizeArray(influencingSpheres,0);

      dbVector& shapeFuncs = particles[i].getShapeFuncs();
      dbMatrix& firstDerivShapeFuncs = particles[i].getFirstDerivShapes();
      dbMatrix& secondDerivShapeFuncs = particles[i].getSecondDerivShapes();
      resizeArray(shapeFuncs,0);
      resizeArray(firstDerivShapeFuncs,0);
      resizeArray(secondDerivShapeFuncs,0);

      intVector& elems = particles[i].getElems();
      resizeArray(elems,0);
    }

    //particles[i].clearArrays();
  }

  resizeArray(extendedLocalGlobalPtcls,m);

#ifdef _geometryDebugMode_
  logFile<<"************** extendedLocalGlobalPtcls **************"<<endl;
  for(int i=0;i<extendedLocalGlobalPtcls.size();i++)
    logFile<<i<<".) "<<extendedLocalGlobalPtcls[i]<<endl;
#endif

}

/**********************************************************************/
/**********************************************************************/
// Create a vector containing all particle coordinates.
void ParticleDistribution::getAllPtcleCoords(InputFileData* InputData,
					     dbVector& allCoords,
					     std::map<std::string,double>& calcData,
					     std::map<std::string,double>& modelData,
					     std::ofstream& logFile) {
  
  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];

  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  dbVector localCoords;
  localCoords.resize(exclusiveLocalPtcls.size()*usedDims);

  int m=0;

  // Loop over the local portion of particles.
  for(int i=0;i<ptcleRootList.size();i++) {

    if(ptcleRootList[i] == rank) {
      dbVector& coords = particles[i].getCoords();
      
      for(int k=0;k<coords.size();k++) {
	localCoords[m] = coords[k];
	m++;
      }

    }

  }

  mergeLocalPtcleVectors(InputData,localCoords,usedDims,allCoords,
			 calcData,modelData,logFile);

}
