#include "VolumeElementTemplates.h"




/***********************************************************************/
/***********************************************************************/

Tetra4ElementTemplate::Tetra4ElementTemplate() {

  allocateArray(nodalCoords,4,4);

  nodalCoords[0][0] = 1;
  nodalCoords[0][1] = 0;
  nodalCoords[0][2] = 0;
  nodalCoords[0][3] = 0;
 
  nodalCoords[1][0] = 0;
  nodalCoords[1][1] = 1;
  nodalCoords[1][2] = 0;
  nodalCoords[1][3] = 0;

  nodalCoords[2][0] = 0;
  nodalCoords[2][1] = 0;
  nodalCoords[2][2] = 1;
  nodalCoords[2][3] = 0;
 
  nodalCoords[3][0] = 0; 
  nodalCoords[3][1] = 0;
  nodalCoords[3][2] = 0;
  nodalCoords[3][3] = 1;

  nodalWeights.resize(4);

  nodalWeights[0] = 1.0;
  nodalWeights[1] = 1.0;
  nodalWeights[2] = 1.0;
  nodalWeights[3] = 1.0;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    intVector idx(4);

    idx[0] = 1;
    idx[1] = 2;
    idx[2] = 3;
    idx[3] = 0;

    reorderVector(nodalCoords,idx); 
    reorderVector(nodalWeights,idx);

  }

}

/***********************************************************************/
double Tetra4ElementTemplate::N(int func,dbVector& coords) {

  using namespace std;

  double result;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    switch(func) {
    case 1:
      result = N1(coords);
      break;
    case 2:
      result = N2(coords);
      break;
    case 3:
      result = N3(coords);
      break;
    case 0:
      result = N4(coords);
      break;
    default:
      cerr <<"Tetra4ElementTemplate has no function N"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  else {

    switch(func) {
    case 0:
      result = N1(coords);
      break;
    case 1:
      result = N2(coords);
      break;
    case 2:
      result = N3(coords);
      break;
    case 3:
      result = N4(coords);
      break;
    default:
      cerr <<"Tetra4ElementTemplate has no function N"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }


  }

  return result;
}

double Tetra4ElementTemplate::N1(dbVector& coords) {
  return (coords[0]);
}
double Tetra4ElementTemplate::N2(dbVector& coords) {
  return (coords[1]);
}
double Tetra4ElementTemplate::N3(dbVector& coords) {
  return (coords[2]);
}
double Tetra4ElementTemplate::N4(dbVector& coords) {
  return (coords[3]);
}

/***********************************************************************/
double Tetra4ElementTemplate::dN(int func,int deriv,dbVector& coords) {

  using namespace std;

  double result;

  switch(func) {
  case 0:
    result = dN1(deriv,coords);
    break;
  case 1:
    result = dN2(deriv,coords);
    break;
  case 2:
    result = dN3(deriv,coords);
    break;
  case 3:
    result = dN4(deriv,coords);
    break;
  default:
    cerr <<"Tetra4ElementTemplate has no function dN"<<func+1<<"!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  return result;
}

double Tetra4ElementTemplate::dN1(int deriv,dbVector& coords) {

  if(deriv == 0)
    return (1.0);
  else
    return (0.0);

}
double Tetra4ElementTemplate::dN2(int deriv,dbVector& coords) {

  if(deriv == 1)
    return (1.0);
  else
    return (0.0);
}
double Tetra4ElementTemplate::dN3(int deriv,dbVector& coords) {

  if(deriv == 2)
    return (1.0);
  else
    return (0.0);
}
double Tetra4ElementTemplate::dN4(int deriv,dbVector& coords) {

  if(deriv == 3)
    return (1.0);
  else
    return (0.0);
}

// dL_j/dx_i needed for  dN_j/dx_k = dN_j/dL_r dL_r/dx_k

/***********************************************************************/
// Calculate the length of a edge of volume element to integrate
// over this edge.
void Tetra4ElementTemplate::lineMetric(std::vector<Particle>& particles,
				       intVector& nodes,
				       dbVector& coords,double& length,
				       std::ofstream& logFile) {
  
  // length = sqrt( |dx|^2 + |dy|^2 + |dz|^2 )
  double deltaX = particles[nodes[0]-1].getCoord(0) -
    particles[nodes[1]-1].getCoord(0);

  double deltaY = particles[nodes[0]-1].getCoord(1) -
    particles[nodes[1]-1].getCoord(1);

  double deltaZ = particles[nodes[0]-1].getCoord(2) -
    particles[nodes[1]-1].getCoord(2);

  length = sqrt(pow(deltaX,2)+pow(deltaY,2)+pow(deltaZ,2));

}

/***********************************************************************/
// Calculate surface normal vector on volume element's side and its 
// length, which is the absolute value of the surface, in order to 
// integrate over this volume element's side.
void Tetra4ElementTemplate::surfaceMetric(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  dbVector& surfaceNormal,
					  double& surface,
					  std::ofstream& logFile) {

  using namespace std;

  // Calculate the surface normal vector
  surfaceNormal[0] = 
    (particles[nodes[0]-1].getCoord(1)*
     particles[nodes[1]-1].getCoord(2) - 
     particles[nodes[0]-1].getCoord(2)*
     particles[nodes[1]-1].getCoord(1)) +
     
    (particles[nodes[1]-1].getCoord(1)*
     particles[nodes[2]-1].getCoord(2) - 
     particles[nodes[1]-1].getCoord(2)*
     particles[nodes[2]-1].getCoord(1)) +
     
    (particles[nodes[2]-1].getCoord(1)*
     particles[nodes[0]-1].getCoord(2) - 
     particles[nodes[2]-1].getCoord(2)*
     particles[nodes[0]-1].getCoord(1));
  
  surfaceNormal[1] = 
    (particles[nodes[0]-1].getCoord(2)*
     particles[nodes[1]-1].getCoord(0) - 
     particles[nodes[0]-1].getCoord(0)*
     particles[nodes[1]-1].getCoord(2)) +
     
    (particles[nodes[1]-1].getCoord(2)*
     particles[nodes[2]-1].getCoord(0) - 
     particles[nodes[1]-1].getCoord(0)*
     particles[nodes[2]-1].getCoord(2)) +
     
    (particles[nodes[2]-1].getCoord(2)*
     particles[nodes[0]-1].getCoord(0) - 
     particles[nodes[2]-1].getCoord(0)*
     particles[nodes[0]-1].getCoord(2));
  
  surfaceNormal[2] = 
    (particles[nodes[0]-1].getCoord(0)*
     particles[nodes[1]-1].getCoord(1) - 
     particles[nodes[0]-1].getCoord(1)*
     particles[nodes[1]-1].getCoord(0)) +
     
    (particles[nodes[1]-1].getCoord(0)*
     particles[nodes[2]-1].getCoord(1) - 
     particles[nodes[1]-1].getCoord(1)*
     particles[nodes[2]-1].getCoord(0)) +
     
    (particles[nodes[2]-1].getCoord(0)*
     particles[nodes[0]-1].getCoord(1) - 
     particles[nodes[2]-1].getCoord(1)*
     particles[nodes[0]-1].getCoord(0));
  
  // Calculate the surface
  surface = 0.5*sqrt(pow(surfaceNormal[0],2)+pow(surfaceNormal[1],2)
		     +pow(surfaceNormal[2],2));

  surfaceNormal[0] /= (2*surface);
  surfaceNormal[1] /= (2*surface);
  surfaceNormal[2] /= (2*surface);

}

/***********************************************************************/
// Calculate the transformation and metric factor for the integration
// of a volume element's volume.
double Tetra4ElementTemplate::metricFactorVolume(std::vector<Particle>& particles,
						 intVector& nodes,
						 dbVector& coords,
						 std::ofstream& logFile) {

  using namespace std;

  double det;
  dbMatrix mat(4,dbVector(4));

  mat[0][0] = 1;
  mat[1][0] = 1;
  mat[2][0] = 1;
  mat[3][0] = 1;

  mat[0][1] = particles[nodes[0]-1].getCoord(0);
  mat[1][1] = particles[nodes[1]-1].getCoord(0); 
  mat[2][1] = particles[nodes[2]-1].getCoord(0);
  mat[3][1] = particles[nodes[3]-1].getCoord(0);

  mat[0][2] = particles[nodes[0]-1].getCoord(1);
  mat[1][2] = particles[nodes[1]-1].getCoord(1); 
  mat[2][2] = particles[nodes[2]-1].getCoord(1);
  mat[3][2] = particles[nodes[3]-1].getCoord(1);

  mat[0][3] = particles[nodes[0]-1].getCoord(2);
  mat[1][3] = particles[nodes[1]-1].getCoord(2); 
  mat[2][3] = particles[nodes[2]-1].getCoord(2);
  mat[3][3] = particles[nodes[3]-1].getCoord(2);


#ifdef _commonDebugMode_
  logFile<<"********** determinant calculation ******************"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"nodal coords: "<<endl;
  for(int i=0;i<nodes.size();i++) {
    dbVector& coordinates = particles[nodes[i]-1].getCoords();
    for(int j=0;j<coordinates.size();j++)
      logFile<< coordinates[j]<<" ";
    logFile<<endl;
  }
#endif

  calcDetDoubleDense(mat,det,logFile);

#ifdef _FEdebugMode_
  logFile<<"determinante = "<<det<<endl;
#endif

  det = fabs(det/6.0);

  if (det < DBL_EPSILON) {
    logFile<<"Tetra4ElementTemplate::metricFactorVolume integration "
	   <<"weight to small!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  return(det);

}


/************************************************************************/
/************************************************************************/
Cube8ElementTemplate::Cube8ElementTemplate() {

  allocateArray(nodalCoords,8,3);


  // default Zienkiewicz nodal numbering
  nodalCoords[0][0] = 1;
  nodalCoords[0][1] = 1;
  nodalCoords[0][2] = 1;

  nodalCoords[1][0] = -1;
  nodalCoords[1][1] = 1;
  nodalCoords[1][2] = 1;

  nodalCoords[2][0] = -1;
  nodalCoords[2][1] = -1;
  nodalCoords[2][2] = 1;

  nodalCoords[3][0] = 1; 
  nodalCoords[3][1] = -1; 
  nodalCoords[3][2] = 1;

  nodalCoords[4][0] = 1;
  nodalCoords[4][1] = 1; 
  nodalCoords[4][2] = -1;

  nodalCoords[5][0] = -1;
  nodalCoords[5][1] = 1;
  nodalCoords[5][2] = -1;

  nodalCoords[6][0] = -1;
  nodalCoords[6][1] = -1;
  nodalCoords[6][2] = -1;

  nodalCoords[7][0] = 1;
  nodalCoords[7][1] = -1;
  nodalCoords[7][2] = -1;

  nodalWeights.resize(8);

  nodalWeights[0] = 1.0;
  nodalWeights[1] = 1.0;
  nodalWeights[2] = 1.0;
  nodalWeights[3] = 1.0;
  nodalWeights[4] = 1.0;
  nodalWeights[5] = 1.0;
  nodalWeights[6] = 1.0;
  nodalWeights[7] = 1.0;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    intVector idx(8);

    idx[0] = 6;
    idx[1] = 7;
    idx[2] = 4;
    idx[3] = 5;
    idx[4] = 2;
    idx[5] = 3;
    idx[6] = 0;
    idx[7] = 1;

    reorderVector(nodalCoords,idx); 
    reorderVector(nodalWeights,idx);

  }

}

/************************************************************************/
double Cube8ElementTemplate::N(int func,dbVector& coords) {

  using namespace std;

  double result;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    switch(func) {
    case 6:
      result = N1(coords);
      break;
    case 7:
      result = N2(coords);
      break;
    case 4:
      result = N3(coords);
      break;
    case 5:
      result = N4(coords);
      break;
    case 2:
      result = N5(coords);
      break;
    case 3:
      result = N6(coords);
      break;
    case 0:
      result = N7(coords);
      break;
    case 1:
      result = N8(coords);
      break;
    default:
      cerr <<"Cube8ElementTemplate has no function N"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }
  // default Zienkiewicz nodal numbering
  else {

    switch(func) {
    case 0:
      result = N1(coords);
      break;
    case 1:
      result = N2(coords);
      break;
    case 2:
      result = N3(coords);
      break;
    case 3:
      result = N4(coords);
      break;
    case 4:
      result = N5(coords);
      break;
    case 5:
      result = N6(coords);
      break;
    case 6:
      result = N7(coords);
      break;
    case 7:
      result = N8(coords);
      break;
    default:
      cerr <<"Cube8ElementTemplate has no function N"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  return result;
}

double Cube8ElementTemplate::N1(dbVector& coords) {
  return (1.0/8.0*(1+coords[0])*(1+coords[1])*(1+coords[2]));
}
double Cube8ElementTemplate::N2(dbVector& coords) {
  return (1.0/8.0*(1-coords[0])*(1+coords[1])*(1+coords[2]));
}
double Cube8ElementTemplate::N3(dbVector& coords) {
  return (1.0/8.0*(1-coords[0])*(1-coords[1])*(1+coords[2]));
}
double Cube8ElementTemplate::N4(dbVector& coords) {
  return (1.0/8.0*(1+coords[0])*(1-coords[1])*(1+coords[2]));
}
double Cube8ElementTemplate::N5(dbVector& coords) {
  return (1.0/8.0*(1+coords[0])*(1+coords[1])*(1-coords[2]));
}
double Cube8ElementTemplate::N6(dbVector& coords) {
  return (1.0/8.0*(1-coords[0])*(1+coords[1])*(1-coords[2]));
}
double Cube8ElementTemplate::N7(dbVector& coords) {
  return (1.0/8.0*(1-coords[0])*(1-coords[1])*(1-coords[2]));
}
double Cube8ElementTemplate::N8(dbVector& coords) {
  return (1.0/8.0*(1+coords[0])*(1-coords[1])*(1-coords[2]));
}

/************************************************************************/
double Cube8ElementTemplate::dN(int func,int deriv,dbVector& coords) {

  using namespace std;

  double result;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    switch(func) {
    case 6:
      result = dN1(deriv,coords);
      break;
    case 7:
      result = dN2(deriv,coords);
      break;
    case 4:
      result = dN3(deriv,coords);
      break;
    case 5:
      result = dN4(deriv,coords);
      break;
    case 2:
      result = dN5(deriv,coords);
      break;
    case 3:
      result = dN6(deriv,coords);
      break;
    case 0:
      result = dN7(deriv,coords);
      break;
    case 1:
      result = dN8(deriv,coords);
      break;

    default:
      cerr <<"Cube8ElementTemplate has no function dN"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }
  // default Zienkiewicz nodal numbering
  else {

    switch(func) {
    case 0:
      result = dN1(deriv,coords);
      break;
    case 1:
      result = dN2(deriv,coords);
      break;
    case 2:
      result = dN3(deriv,coords);
      break;
    case 3:
      result = dN4(deriv,coords);
      break;
    case 4:
      result = dN5(deriv,coords);
      break;
    case 5:
      result = dN6(deriv,coords);
      break;
    case 6:
      result = dN7(deriv,coords);
      break;
    case 7:
      result = dN8(deriv,coords);
      break;
    default:
      cerr <<"Cube8ElementTemplate has no function dN"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  return result;
}

double Cube8ElementTemplate::dN1(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (1.0/8.0*(1+coords[1])*(1+coords[2]));
    break;

  case 1:
    return (1.0/8.0*(1+coords[0])*(1+coords[2]));
    break;

  case 2:
    return (1.0/8.0*(1+coords[0])*(1+coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN1/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube8ElementTemplate::dN2(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (-1.0/8.0*(1+coords[1])*(1+coords[2]));
    break;

  case 1:
    return (1.0/8.0*(1-coords[0])*(1+coords[2]));
    break;

  case 2:
    return (1.0/8.0*(1-coords[0])*(1+coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN2/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube8ElementTemplate::dN3(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (-1.0/8.0*(1-coords[1])*(1+coords[2]));
    break;

  case 1:
    return (-1.0/8.0*(1-coords[0])*(1+coords[2]));
    break;

  case 2:
    return (1.0/8.0*(1-coords[0])*(1-coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN3/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube8ElementTemplate::dN4(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (1.0/8.0*(1-coords[1])*(1+coords[2]));
    break;

  case 1:
    return (-1.0/8.0*(1+coords[0])*(1+coords[2]));
    break;

  case 2:
    return (1.0/8.0*(1+coords[0])*(1-coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN4/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
  
}
double Cube8ElementTemplate::dN5(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (1.0/8.0*(1+coords[1])*(1-coords[2]));
    break;

  case 1:
    return (1.0/8.0*(1+coords[0])*(1-coords[2]));
    break;

  case 2:
    return (-1.0/8.0*(1+coords[0])*(1+coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN5/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
  
}
double Cube8ElementTemplate::dN6(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (-1.0/8.0*(1+coords[1])*(1-coords[2]));
    break;

  case 1:
    return (1.0/8.0*(1-coords[0])*(1-coords[2]));
    break;

  case 2:
    return (-1.0/8.0*(1-coords[0])*(1+coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN6/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
  
}
double Cube8ElementTemplate::dN7(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (-1.0/8.0*(1-coords[1])*(1-coords[2]));
    break;

  case 1:
    return (-1.0/8.0*(1-coords[0])*(1-coords[2]));
    break;

  case 2:
    return (-1.0/8.0*(1-coords[0])*(1-coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN7/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
  
}
double Cube8ElementTemplate::dN8(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (1.0/8.0*(1-coords[1])*(1-coords[2]));
    break;

  case 1:
    return (-1.0/8.0*(1+coords[0])*(1-coords[2]));
    break;

  case 2:
    return (-1.0/8.0*(1+coords[0])*(1-coords[1]));
    break;

  default:
    cerr <<"In Cube8ElementTemplate::dN8/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }
  
}

/***********************************************************************/
// Calculate the length of a edge of volume element to integrate
// over this edge (eta = zeta = 1).
void Cube8ElementTemplate::lineMetric(std::vector<Particle>& particles,
				      intVector& nodes,
				      dbVector& coords,double& length,
				      std::ofstream& logFile) {
  
  // length = sqrt( |dx|^2 + |dy|^2 + |dz|^2 )
  double deltaX = 
    0.5*particles[nodes[0]-1].getCoord(0)-
    0.5*particles[nodes[1]-1].getCoord(0);

  double deltaY = 
    0.5*particles[nodes[0]-1].getCoord(1)-
    0.5*particles[nodes[1]-1].getCoord(1);

  double deltaZ = 
    0.5*particles[nodes[0]-1].getCoord(2)-
    0.5*particles[nodes[1]-1].getCoord(2);

  length = sqrt(pow(deltaX,2)+pow(deltaY,2)+pow(deltaZ,2));

}

/************************************************************************/
// Compute Jacobian dx_i/dL_j
void Cube8ElementTemplate::jacobian(std::vector<Particle>& particles,
				    intVector& nodes,
				    dbVector& coords,
				    dbMatrix& Jmat,
				    std::ofstream& logFile) {

  using namespace std;

  if(Jmat.size() == 0)
    allocateArray(Jmat,3,3);

  clearArray(Jmat);

  // dx_i/dL_j = dN^I/dL_j x^I_i


  for(int i=0;i<3;i++) {

    for(int j=0;j<3;j++) {
      
      for(int k=0;k<nodes.size();k++) { 

	Jmat[i][j] += dN(k,j,coords)*particles[nodes[k]-1].getCoord(i);
      }

    }

  }

#ifdef _FEdebugMode_
  dbMatrix Jmat2;
  allocateArray(Jmat2,3,3);
  Jmat2[0][0] = 
    1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(0)
    -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(0)
    -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(0)
    +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(0)
    +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(0)
    -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(0)
    -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(0)
    +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(0);
  
  Jmat2[0][1] = 
    1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(0)
    +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(0)
    -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(0)
    -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(0)
    +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(0)
    +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(0)
    -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(0)
    -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(0);

  Jmat2[0][2] = 
    1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(0)
    +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(0)
    +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(0)
    +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(0)
    -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(0)
    -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(0)
    -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(0)
    -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(0);

  Jmat2[1][0] = 
    1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(1)
    -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(1)
    -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(1)
    +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(1)
    +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(1)
    -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(1)
    -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(1)
    +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(1);

  Jmat2[1][1] = 
    1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(1)
    +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(1)
    -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(1)
    -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(1)
    +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(1)
    +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(1)
    -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(1)
    -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(1);

  Jmat2[1][2] = 
    1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(1)
    +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(1)
    +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(1)
    +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(1)
    -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(1)
    -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(1)
    -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(1)
    -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(1);

  Jmat2[2][0] = 
    1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(2)
    -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(2)
    -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(2)
    +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(2)
    +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(2)
    -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(2)
    -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(2)
    +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(2);

  Jmat2[2][1] = 
    1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(2)
    +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(2)
    -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(2)
    -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(2)
    +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(2)
    +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(2)
    -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(2)
    -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(2);
  
  Jmat2[2][2] = 
    1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(2)
    +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(2)
    +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(2)
    +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(2)
    -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(2)
    -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(2)
    -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(2)
    -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(2);

  logFile<<"------------------------------------------------------"<<endl;
  logFile<<"------------------ jacobi matrix ---------------------"<<endl;
  for(int i=0;i<Jmat.size();i++)
    for(int j=0;j<Jmat[i].size();j++)
      logFile<<"J["<<i<<"]["<<j<<"]="<<Jmat[i][j]<<" =? "
	     <<Jmat2[i][j]<<endl;
#endif

}


/************************************************************************/
// Calculate surface normal vector on volume element's side and its 
// length, which is the absolute value of the surface, in order to 
// integrate over this volume element's side.
void Cube8ElementTemplate::surfaceMetric(std::vector<Particle>& particles,
					 intVector& nodes,
					 dbVector& coords,
					 dbVector& surfaceNormal,
					 double& surface,
					 std::ofstream& logFile) {

  using namespace std;

  // Calculate the surface normal vector
  surfaceNormal[0] = 
    ((1.0/4.0*(1-coords[1])*particles[nodes[0]-1].getCoord(1)
     -1.0/4.0*(1-coords[1])*particles[nodes[1]-1].getCoord(1)
     -1.0/4.0*(1+coords[1])*particles[nodes[2]-1].getCoord(1)
     +1.0/4.0*(1+coords[1])*particles[nodes[3]-1].getCoord(1))*

    (-1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(2)
     -1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(2)
     +1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(2)
     +1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(2))) -

    ((1.0/4.0*(1-coords[1])*particles[nodes[0]-1].getCoord(2)
     -1.0/4.0*(1-coords[1])*particles[nodes[1]-1].getCoord(2)
     -1.0/4.0*(1+coords[1])*particles[nodes[2]-1].getCoord(2)
     +1.0/4.0*(1+coords[1])*particles[nodes[3]-1].getCoord(2))*

    (-1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(1)
     -1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(1)
     +1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(1)
     +1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(1)));

  surfaceNormal[1] = 
    ((1.0/4.0*(1-coords[1])*particles[nodes[0]-1].getCoord(2)
     -1.0/4.0*(1-coords[1])*particles[nodes[1]-1].getCoord(2)
     -1.0/4.0*(1+coords[1])*particles[nodes[2]-1].getCoord(2)
     +1.0/4.0*(1+coords[1])*particles[nodes[3]-1].getCoord(2))*

    (-1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(0)
     -1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(0)
     +1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(0)
     +1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(0))) -

    ((1.0/4.0*(1-coords[1])*particles[nodes[0]-1].getCoord(0)
     -1.0/4.0*(1-coords[1])*particles[nodes[1]-1].getCoord(0)
     -1.0/4.0*(1+coords[1])*particles[nodes[2]-1].getCoord(0)
     +1.0/4.0*(1+coords[1])*particles[nodes[3]-1].getCoord(0))*

    (-1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(2)
     -1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(2)
     +1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(2)
     +1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(2)));

  surfaceNormal[2] = 
    ((1.0/4.0*(1-coords[1])*particles[nodes[0]-1].getCoord(0)
     -1.0/4.0*(1-coords[1])*particles[nodes[1]-1].getCoord(0)
     -1.0/4.0*(1+coords[1])*particles[nodes[2]-1].getCoord(0)
     +1.0/4.0*(1+coords[1])*particles[nodes[3]-1].getCoord(0))*

    (-1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(1)
     -1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(1)
     +1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(1)
     +1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(1))) -

    ((1.0/4.0*(1-coords[1])*particles[nodes[0]-1].getCoord(1)
     -1.0/4.0*(1-coords[1])*particles[nodes[1]-1].getCoord(1)
     -1.0/4.0*(1+coords[1])*particles[nodes[2]-1].getCoord(1)
     +1.0/4.0*(1+coords[1])*particles[nodes[3]-1].getCoord(1))*

    (-1.0/4.0*(1+coords[0])*particles[nodes[0]-1].getCoord(0)
     -1.0/4.0*(1-coords[0])*particles[nodes[1]-1].getCoord(0)
     +1.0/4.0*(1-coords[0])*particles[nodes[2]-1].getCoord(0)
     +1.0/4.0*(1+coords[0])*particles[nodes[3]-1].getCoord(0)));

  // Calculate the surface
  surface = sqrt(pow(surfaceNormal[0],2)+pow(surfaceNormal[1],2)
		 +pow(surfaceNormal[2],2));

  surfaceNormal[0] = surfaceNormal[0]/surface;
  surfaceNormal[1] = surfaceNormal[1]/surface;
  surfaceNormal[2] = surfaceNormal[2]/surface;

}

/***********************************************************************/
// Calculate the transformation and metric factor for the integration
// of a volume element's volume.
double Cube8ElementTemplate::metricFactorVolume(std::vector<Particle>& particles,
						intVector& nodes,
						dbVector& coords,
						std::ofstream& logFile) {


  using namespace std;

  double det;
  dbMatrix Jmat;

  jacobian(particles,nodes,coords,Jmat,logFile);
  int info =  calcDetDoubleDense(Jmat,det,logFile);

  det = fabs(det);

  if (det < DBL_EPSILON || info != 0) {
    logFile<<"Cube8ElementTemplate::metricFactorVolume integration "
	   <<"weight to small!"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

#ifdef _FEdebugMode_
  double det2 =
   ((1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(0)
     -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(0)
     -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(0)
     +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(0)
     +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(0)
     -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(0)
     -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(0)
     +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(0))
    *(1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(1)
      +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(1)
      -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(1)
      -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(1)
      +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(1)
      +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(1)
      -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(1)
      -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(1))
    *(1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(2)
      +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(2)
      +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(2)
      +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(2)
      -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(2)
      -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(2)
      -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(2)
      -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(2)))

   +((1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(1)
      -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(1)
      -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(1)
      +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(1)
      +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(1)
      -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(1)
      -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(1)
      +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(1))
     *(1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(2)
       +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(2)
       -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(2)
       -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(2)
       +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(2)
       +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(2)
       -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(2)
       -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(2))
     *(1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(0)
       +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(0)))
   
   +((1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(2)
      -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(2)
      -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(2)
      +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(2)
      +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(2)
      -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(2)
      -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(2)
      +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(2))
     *(1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(0)
       +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(0))
     *(1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(1)
       +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(1)
       +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(1)
       +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(1)
       -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(1)
       -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(1)
       -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(1)
       -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(1)))
   
   -((1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(2)
      -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(2)
      -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(2)
      +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(2)
      +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(2)
      -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(2)
      -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(2)
      +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(2))
     *(1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(1)
       +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(1)
       -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(1)
       -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(1)
       +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(1)
       +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(1)
       -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(1)
       -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(1))
     *(1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(0)
       +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(0)))
   
   -((1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(0)
      -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(0)
      -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(0)
      +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(0)
      +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(0)
      -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(0)
      -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(0)
      +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(0))
     *(1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(2)
       +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(2)
       -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(2)
       -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(2)
       +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(2)
       +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(2)
       -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(2)
       -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(2))
     *(1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(1)
       +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(1)
       +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(1)
       +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(1)
       -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(1)
       -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(1)
       -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(1)
       -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(1)))
   
   -((1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[0]-1].getCoord(1)
      -1.0/8.0*(1+coords[1])*(1+coords[2])*particles[nodes[1]-1].getCoord(1)
      -1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[2]-1].getCoord(1)
      +1.0/8.0*(1-coords[1])*(1+coords[2])*particles[nodes[3]-1].getCoord(1)
      +1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[4]-1].getCoord(1)
      -1.0/8.0*(1+coords[1])*(1-coords[2])*particles[nodes[5]-1].getCoord(1)
      -1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[6]-1].getCoord(1)
      +1.0/8.0*(1-coords[1])*(1-coords[2])*particles[nodes[7]-1].getCoord(1))
     *(1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[0]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[1]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1+coords[2])*particles[nodes[2]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1+coords[2])*particles[nodes[3]-1].getCoord(0)
       +1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[4]-1].getCoord(0)
       +1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[5]-1].getCoord(0)
       -1.0/8.0*(1-coords[0])*(1-coords[2])*particles[nodes[6]-1].getCoord(0)
       -1.0/8.0*(1+coords[0])*(1-coords[2])*particles[nodes[7]-1].getCoord(0))
     *(1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[0]-1].getCoord(2)
       +1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[1]-1].getCoord(2)
       +1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[2]-1].getCoord(2)
       +1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[3]-1].getCoord(2)
       -1.0/8.0*(1+coords[0])*(1+coords[1])*particles[nodes[4]-1].getCoord(2)
       -1.0/8.0*(1-coords[0])*(1+coords[1])*particles[nodes[5]-1].getCoord(2)
       -1.0/8.0*(1-coords[0])*(1-coords[1])*particles[nodes[6]-1].getCoord(2)
       -1.0/8.0*(1+coords[0])*(1-coords[1])*particles[nodes[7]-1].getCoord(2)));

  logFile<<"det="<<det<<" ?= "<<det2<<endl;
#endif

  return(det);
}

/************************************************************************/
/************************************************************************/
Cube27ElementTemplate::Cube27ElementTemplate() {

  allocateArray(nodalCoords,27,3);

  nodalCoords[0][0] = 1; // 6
  nodalCoords[0][1] = 1;
  nodalCoords[0][2] = 1;

  nodalCoords[1][0] = -1; // 7
  nodalCoords[1][1] = 1;
  nodalCoords[1][2] = 1;

  nodalCoords[2][0] = -1; // 4
  nodalCoords[2][1] = -1;
  nodalCoords[2][2] = 1;

  nodalCoords[3][0] = 1;  // 5
  nodalCoords[3][1] = -1;
  nodalCoords[3][2] = 1;

  nodalCoords[4][0] = 1; // 2
  nodalCoords[4][1] = 1;
  nodalCoords[4][2] = -1;

  nodalCoords[5][0] = -1; // 3
  nodalCoords[5][1] = 1;
  nodalCoords[5][2] = -1;

  nodalCoords[6][0] = -1; // 0
  nodalCoords[6][1] = -1;
  nodalCoords[6][2] = -1;

  nodalCoords[7][0] = 1; // 1
  nodalCoords[7][1] = -1;
  nodalCoords[7][2] = -1;

  nodalCoords[8][0] = 0; // 18
  nodalCoords[8][1] = 1;
  nodalCoords[8][2] = 1;

  nodalCoords[9][0] = -1; // 19
  nodalCoords[9][1] = 0;
  nodalCoords[9][2] = 1;

  nodalCoords[10][0] = 0;
  nodalCoords[10][1] = -1; // 16
  nodalCoords[10][2] = 1;

  nodalCoords[11][0] = 1;
  nodalCoords[11][1] = 0; // 17
  nodalCoords[11][2] = 1;

  nodalCoords[12][0] = 1; // 14
  nodalCoords[12][1] = 1;
  nodalCoords[12][2] = 0;

  nodalCoords[13][0] = -1; // 15
  nodalCoords[13][1] = 1;
  nodalCoords[13][2] = 0;

  nodalCoords[14][0] = -1;
  nodalCoords[14][1] = -1; // 12
  nodalCoords[14][2] = 0;

  nodalCoords[15][0] = 1;
  nodalCoords[15][1] = -1; // 13
  nodalCoords[15][2] = 0;

  nodalCoords[16][0] = 0; // 10
  nodalCoords[16][1] = 1;
  nodalCoords[16][2] = -1;

  nodalCoords[17][0] = -1; // 11
  nodalCoords[17][1] = 0;
  nodalCoords[17][2] = -1;

  nodalCoords[18][0] = 0; // 8
  nodalCoords[18][1] = -1;
  nodalCoords[18][2] = -1;

  nodalCoords[19][0] = 1; // 9
  nodalCoords[19][1] = 0;
  nodalCoords[19][2] = -1;

  nodalCoords[20][0] = 0; // 25
  nodalCoords[20][1] = 0;
  nodalCoords[20][2] = 1;

  nodalCoords[21][0] = 0; // 23
  nodalCoords[21][1] = 1;
  nodalCoords[21][2] = 0;

  nodalCoords[22][0] = -1; // 24
  nodalCoords[22][1] = 0;
  nodalCoords[22][2] = 0;

  nodalCoords[23][0] = 0; // 21
  nodalCoords[23][1] = -1;
  nodalCoords[23][2] = 0;

  nodalCoords[24][0] = 1; // 22
  nodalCoords[24][1] = 0;
  nodalCoords[24][2] = 0;

  nodalCoords[25][0] = 0; // 20
  nodalCoords[25][1] = 0;
  nodalCoords[25][2] = -1;

  nodalCoords[26][0] = 0; // 26
  nodalCoords[26][1] = 0;
  nodalCoords[26][2] = 0;

  nodalWeights.resize(27);

  nodalWeights[0] = 0.125;
  nodalWeights[1] =  0.125;
  nodalWeights[2] =  0.125;
  nodalWeights[3] =  0.125;
  nodalWeights[4] =  0.125;
  nodalWeights[5] =  0.125;
  nodalWeights[6] =  0.125;
  nodalWeights[7] =  0.125;

  nodalWeights[8] = 0.25;
  nodalWeights[9] =  0.25;
  nodalWeights[10] =  0.25;
  nodalWeights[11] =  0.25;
  nodalWeights[12] =  0.25;
  nodalWeights[13] =  0.25;
  nodalWeights[14] =  0.25;
  nodalWeights[15] =  0.25;
  nodalWeights[16] =  0.25;
  nodalWeights[17] =  0.25;
  nodalWeights[18] =  0.25;
  nodalWeights[19] =  0.25;

  nodalWeights[20] = 0.5;
  nodalWeights[21] = 0.5;
  nodalWeights[22] = 0.5;
  nodalWeights[23] = 0.5;
  nodalWeights[24] = 0.5;
  nodalWeights[25] = 0.5;

  nodalWeights[26] = 1.0;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    intVector idx(27);

    // idx[new] = old(Zienkiewicz)
    idx[0] = 6;
    idx[1] = 7;
    idx[2] = 4;
    idx[3] = 5;
    idx[4] = 2;
    idx[5] = 3;
    idx[6] = 0;
    idx[7] = 1;  
    idx[8] = 18;
    idx[9] = 19;
    idx[10] = 16;
    idx[11] = 17;
    idx[12] = 14;
    idx[13] = 15;
    idx[14] = 12;
    idx[15] = 13;
    idx[16] = 10;
    idx[17] = 11;
    idx[18] = 8;
    idx[19] = 9;
    idx[20] = 25;
    idx[21] = 23;
    idx[22] = 24;
    idx[23] = 21;
    idx[24] = 22;
    idx[25] = 20;
    idx[26] = 26;

    reorderVector(nodalCoords,idx); 
    reorderVector(nodalWeights,idx);
  }

}

/************************************************************************/
double Cube27ElementTemplate::N(int func,dbVector& coords) {

  using namespace std;

  double result;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    switch(func) {
    case 6:
      result = N1(coords);
      break;
    case 7:
      result = N2(coords);
      break;
    case 4:
      result = N3(coords);
      break;
    case 5:
      result = N4(coords);
      break;
    case 2:
      result = N5(coords);
      break;
    case 3:
      result = N6(coords);
      break;
    case 0:
      result = N7(coords);
      break;
    case 1:
      result = N8(coords);
      break;
    case 18:
      result = N9(coords);
      break;
    case 19:
      result = N10(coords);
      break;
    case 16:
      result = N11(coords);
      break;
    case 17:
      result = N12(coords);
      break;
    case 14:
      result = N13(coords);
      break;
    case 15:
      result = N14(coords);
      break;
    case 12:
      result = N15(coords);
      break;
    case 13:
      result = N16(coords);
      break;
    case 10:
      result = N17(coords);
      break;
    case 11:
      result = N18(coords);
      break;
    case 8:
      result = N19(coords);
      break;
    case 9:
      result = N20(coords);
      break;
    case 25:
      result = N21(coords);
      break;
    case 23:
      result = N22(coords);
      break;
    case 24:
      result = N23(coords);
      break;
    case 21:
      result = N24(coords);
      break;
    case 22:
      result = N25(coords);
      break;
    case 20:
      result = N26(coords);
      break;
    case 26:
      result = N27(coords);
      break;
    default:
      cerr <<"Cube27ElementTemplate has no function N"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  // default Zienkiewicz nodal numbering
  else {

    switch(func) {
    case 0:
      result = N1(coords);
      break;
    case 1:
      result = N2(coords);
      break;
    case 2:
      result = N3(coords);
      break;
    case 3:
      result = N4(coords);
      break;
    case 4:
      result = N5(coords);
      break;
    case 5:
      result = N6(coords);
      break;
    case 6:
      result = N7(coords);
      break;
    case 7:
      result = N8(coords);
      break;
    case 8:
      result = N9(coords);
      break;
    case 9:
      result = N10(coords);
      break;
    case 10:
      result = N11(coords);
      break;
    case 11:
      result = N12(coords);
      break;
    case 12:
      result = N13(coords);
      break;
    case 13:
      result = N14(coords);
      break;
    case 14:
      result = N15(coords);
      break;
    case 15:
      result = N16(coords);
      break;
    case 16:
      result = N17(coords);
      break;
    case 17:
      result = N18(coords);
      break;
    case 18:
      result = N19(coords);
      break;
    case 19:
      result = N20(coords);
      break;
    case 20:
      result = N21(coords);
      break;
    case 21:
      result = N22(coords);
      break;
    case 22:
      result = N23(coords);
      break;
    case 23:
      result = N24(coords);
      break;
    case 24:
      result = N25(coords);
      break;
    case 25:
      result = N26(coords);
      break;
    case 26:
      result = N27(coords);
      break;
    default:
      cerr <<"Cube27ElementTemplate has no function N"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  return result;
}

double Cube27ElementTemplate::N1(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N2(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N3(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N4(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N5(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N6(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)-coords[2])
);
}
double Cube27ElementTemplate::N7(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N8(dbVector& coords) {
  return (
	  1.0/8.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N9(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N10(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N11(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N12(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N13(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N14(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N15(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N16(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N17(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N18(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N19(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N20(dbVector& coords) {
  return (
	  -1.0/4.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N21(dbVector& coords) {
  return (
	  1.0/2.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)+coords[2])
	  );
}
double Cube27ElementTemplate::N22(dbVector& coords) {
  return (
	  1.0/2.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)+coords[1])
	  *(pow(coords[2],2)-1.0)
);
}
double Cube27ElementTemplate::N23(dbVector& coords) {
  return (
	  1.0/2.0*(pow(coords[0],2)-coords[0])
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N24(dbVector& coords) {
  return (
	  1.0/2.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)-coords[1])
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N25(dbVector& coords) {
  return (
	  1.0/2.0*(pow(coords[0],2)+coords[0])
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)-1.0)
	  );
}
double Cube27ElementTemplate::N26(dbVector& coords) {
  return (
	  1.0/2.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)-coords[2])
	  );
}
double Cube27ElementTemplate::N27(dbVector& coords) {
  return (
	  -1.0*(pow(coords[0],2)-1.0)
	  *(pow(coords[1],2)-1.0)
	  *(pow(coords[2],2)-1.0)
	  );
}

/************************************************************************/
double Cube27ElementTemplate::dN(int func,int deriv,dbVector& coords) {

  using namespace std;

  double result;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    switch(func) {
    case 6:
      result = dN1(deriv,coords);
      break;
    case 7:
      result = dN2(deriv,coords);
      break;
    case 4:
      result = dN3(deriv,coords);
      break;
    case 5:
      result = dN4(deriv,coords);
      break;
    case 2:
      result = dN5(deriv,coords);
      break;
    case 3:
      result = dN6(deriv,coords);
      break;
    case 0:
      result = dN7(deriv,coords);
      break;
    case 1:
      result = dN8(deriv,coords);
      break;
    case 18:
      result = dN9(deriv,coords);
      break;
    case 19:
      result = dN10(deriv,coords);
      break;
    case 16:
      result = dN11(deriv,coords);
      break;
    case 17:
      result = dN12(deriv,coords);
      break;
    case 14:
      result = dN13(deriv,coords);
      break;
    case 15:
      result = dN14(deriv,coords);
      break;
    case 12:
      result = dN15(deriv,coords);
      break;
    case 13:
      result = dN16(deriv,coords);
      break;
    case 10:
      result = dN17(deriv,coords);
      break;
    case 11:
      result = dN18(deriv,coords);
      break;
    case 8:
      result = dN19(deriv,coords);
      break;
    case 9:
      result = dN20(deriv,coords);
      break;
    case 25:
      result = dN21(deriv,coords);
      break;
    case 23:
      result = dN22(deriv,coords);
      break;
    case 24:
      result = dN23(deriv,coords);
      break;
    case 21:
      result = dN24(deriv,coords);
      break;
    case 22:
      result = dN25(deriv,coords);
      break;
    case 20:
      result = dN26(deriv,coords);
      break;
    case 26:
      result = dN27(deriv,coords);
      break;
    default:
      cerr <<"Cube27ElementTemplate has no function dN"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  // default Zienkiewicz nodal numbering
  else {

    switch(func) {
    case 0:
      result = dN1(deriv,coords);
      break;
    case 1:
      result = dN2(deriv,coords);
      break;
    case 2:
      result = dN3(deriv,coords);
      break;
    case 3:
      result = dN4(deriv,coords);
      break;
    case 4:
      result = dN5(deriv,coords);
      break;
    case 5:
      result = dN6(deriv,coords);
      break;
    case 6:
      result = dN7(deriv,coords);
      break;
    case 7:
      result = dN8(deriv,coords);
      break;
    case 8:
      result = dN9(deriv,coords);
      break;
    case 9:
      result = dN10(deriv,coords);
      break;
    case 10:
      result = dN11(deriv,coords);
      break;
    case 11:
      result = dN12(deriv,coords);
      break;
    case 12:
      result = dN13(deriv,coords);
      break;
    case 13:
      result = dN14(deriv,coords);
      break;
    case 14:
      result = dN15(deriv,coords);
      break;
    case 15:
      result = dN16(deriv,coords);
      break;
    case 16:
      result = dN17(deriv,coords);
      break;
    case 17:
      result = dN18(deriv,coords);
      break;
    case 18:
      result = dN19(deriv,coords);
      break;
    case 19:
      result = dN20(deriv,coords);
      break;
    case 20:
      result = dN21(deriv,coords);
      break;
    case 21:
      result = dN22(deriv,coords);
      break;
    case 22:
      result = dN23(deriv,coords);
      break;
    case 23:
      result = dN24(deriv,coords);
      break;
    case 24:
      result = dN25(deriv,coords);
      break;
    case 25:
      result = dN26(deriv,coords);
      break;
    case 26:
      result = dN27(deriv,coords);
      break;
    default:
      cerr <<"Cube8ElementTemplate has no function dN"<<func+1<<"!"<< endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

  }

  return result;
}

double Cube27ElementTemplate::dN1(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN1/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN2(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN2/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN3(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN3/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN4(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN4/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN5(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN5/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN6(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN6/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN7(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN7/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN8(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/8.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/8.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN8/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN9(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/2.0*(coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN9/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN10(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/2.0*(pow(coords[0],2)-coords[0])
	    *(coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN10/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN11(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/2.0*(coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN11/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN12(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/2.0*(pow(coords[0],2)+coords[0])
	    *(coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN12/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN13(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)+coords[0])
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    -1.0/2.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN13/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN14(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)-coords[0])
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    -1.0/2.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN14/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN15(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)-coords[0])
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    -1.0/2.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN15/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN16(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)+coords[0])
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    -1.0/2.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN16/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN17(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/2.0*(coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN17/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN18(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/2.0*(pow(coords[0],2)-coords[0])
	    *(coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN18/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN19(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/2.0*(coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN19/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN20(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -1.0/4.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    -1.0/2.0*(pow(coords[0],2)+coords[0])
	    *(coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    -1.0/4.0*(pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN20/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN21(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    (coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 1:
    return (
	    (pow(coords[0],2)-1.0)
	    *(coords[1])
	    *(pow(coords[2],2)+coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/2.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)-1.0)
	    *(2.0*coords[2]+1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN21/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN22(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    (coords[0])
	    *(pow(coords[1],2)+coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    1.0/2.0*(pow(coords[0],2)-1.0)
	    *(2.0*coords[1]+1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    (pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)+coords[1])
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN22/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN23(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/2.0*(2.0*coords[0]-1.0)
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    (pow(coords[0],2)-coords[0])
	    *(coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    (pow(coords[0],2)-coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN23/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN24(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    (coords[0])
	    *(pow(coords[1],2)-coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    1.0/2.0*(pow(coords[0],2)-1.0)
	    *(2.0*coords[1]-1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    (pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)-coords[1])
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN24/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN25(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    1.0/2.0*(2.0*coords[0]+1.0)
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    (pow(coords[0],2)+coords[0])
	    *(coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    (pow(coords[0],2)+coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN25/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN26(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    (coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 1:
    return (
	    (pow(coords[0],2)-1.0)
	    *(coords[1])
	    *(pow(coords[2],2)-coords[2])
	    );
    break;

  case 2:
    return (
	    1.0/2.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)-1.0)
	    *(2.0*coords[2]-1.0)
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN26/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}
double Cube27ElementTemplate::dN27(int deriv,dbVector& coords) {

  using namespace std;

  switch(deriv) {

  case 0:
    return (
	    -2.0*(coords[0])
	    *(pow(coords[1],2)-1.0)
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 1:
    return (
	    -2.0*(pow(coords[0],2)-1.0)
	    *(coords[1])
	    *(pow(coords[2],2)-1.0)
	    );
    break;

  case 2:
    return (
	    -2.0*(pow(coords[0],2)-1.0)
	    *(pow(coords[1],2)-1.0)
	    *(coords[2])
	    );
    break;

  default:
    cerr <<"In Cube27ElementTemplate::dN27/dx"<<deriv<<" not available!"
	 <<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

}

/************************************************************************/
// Compute Jacobi matrix dx_i/dL_j
void Cube27ElementTemplate::jacobian(std::vector<Particle>& particles,
				     intVector& nodes,
				     dbVector& coords,
				     dbMatrix& Jmat,
				     std::ofstream& logFile) {
  
  using namespace std;

  if(Jmat.size() == 0)
    allocateArray(Jmat,3,3);

  clearArray(Jmat);

  // dx_i/dL_j = dN^I/dL_j x^I_i

  for(int i=0;i<3;i++) 

    for(int j=0;j<3;j++) 

      // loop over all nodes
      for(int k=0;k<nodalCoords.size();k++)

	Jmat[i][j] += dN(k,j,coords)
	  *particles[nodes[k]-1].getCoord(i);

#ifdef _FEdebugMode_
  logFile<<"********************** Jacobian *********************"<<endl;
  logFile<<"*****************************************************"<<endl;
  for(int i=0;i<Jmat.size();i++)
    for(int j=0;j<Jmat[i].size();j++)
      logFile<<"J["<<i<<"]["<<j<<"] = "<<Jmat[i][j]<<endl;
#endif
  
}

/***********************************************************************/
// Calculate the length of a edge of volume element to integrate
// over this edge (eta = zeta = 1).
void Cube27ElementTemplate::lineMetric(std::vector<Particle>& particles,
				      intVector& nodes,
				      dbVector& coords,double& length,
				      std::ofstream& logFile) {
  
  // length = sqrt( |dx|^2 + |dy|^2 + |dz|^2 )

  double xi = coords[0];

  double x1 = particles[nodes[0]-1].getCoord(0);
  double x2 = particles[nodes[1]-1].getCoord(0);
  double x9 = particles[nodes[2]-1].getCoord(0);

  double y1 = particles[nodes[0]-1].getCoord(1);
  double y2 = particles[nodes[1]-1].getCoord(1);
  double y9 = particles[nodes[2]-1].getCoord(1);
  
  double z1 = particles[nodes[0]-1].getCoord(2);
  double z2 = particles[nodes[1]-1].getCoord(2);
  double z9 = particles[nodes[2]-1].getCoord(2);

  length = sqrt(-16.0*z2*xi*xi*z9-16.0*z1*xi*xi*z9+8.0*z1*xi*xi*z2+8.0*
		y2*xi*y9-8.0*y1*xi*y9-16.0*y2*xi*xi*y9-16.0*y1*xi*xi*y9+
		8.0*y1*xi*xi*y2+8.0*x2*xi*x9-8.0*x1*xi*x9-16.0*x2*xi*xi*
		x9-16.0*x1*xi*xi*x9+8.0*x1*xi*xi*x2+8.0*z2*xi*z9-8.0*z1*
		xi*z9+16.0*xi*xi*z9*z9+z2*z2-4.0*z2*z2*xi+4.0*z2*z2*xi*
		xi-2.0*z1*z2+z1*z1+4.0*z1*z1*xi+4.0*z1*z1*xi*xi+16.0*xi*
		xi*y9*y9+y2*y2-4.0*y2*y2*xi+4.0*y2*y2*xi*xi-2.0*y1*y2+y1*
		y1+4.0*y1*y1*xi+4.0*y1*y1*xi*xi+16.0*xi*xi*x9*x9+x2*x2-
		4.0*x2*x2*xi+4.0*x2*x2*xi*xi-2.0*x1*x2+x1*x1+4.0*x1*x1*
		xi+4.0*x1*x1*xi*xi)/2.0;

}

/************************************************************************/
// Calculate surface normal vector on volume element's side and its 
// length, which is the absolute value of the surface, in order to 
// integrate over this volume element's side.
void Cube27ElementTemplate::surfaceMetric(std::vector<Particle>& particles,
					  intVector& nodes,
					  dbVector& coords,
					  dbVector& surfaceNormal,
					  double& surface,
					  std::ofstream& logFile) {

  using namespace std;

  // Calculate the surface normal vector

  double xi = coords[0];
  double eta = coords[1];

  double x1 = particles[nodes[0]-1].getCoord(0);
  double x2 = particles[nodes[1]-1].getCoord(0);
  double x3 = particles[nodes[2]-1].getCoord(0);
  double x4 = particles[nodes[3]-1].getCoord(0);
  double x9 = particles[nodes[4]-1].getCoord(0);
  double x10 = particles[nodes[5]-1].getCoord(0);
  double x11 = particles[nodes[6]-1].getCoord(0);
  double x12 = particles[nodes[7]-1].getCoord(0);
  double x21 = particles[nodes[8]-1].getCoord(0);

  double y1 = particles[nodes[0]-1].getCoord(1);
  double y2 = particles[nodes[1]-1].getCoord(1);
  double y3 = particles[nodes[2]-1].getCoord(1);
  double y4 = particles[nodes[3]-1].getCoord(1);
  double y9 = particles[nodes[4]-1].getCoord(1);
  double y10 = particles[nodes[5]-1].getCoord(1);
  double y11 = particles[nodes[6]-1].getCoord(1);
  double y12 = particles[nodes[7]-1].getCoord(1);
  double y21 = particles[nodes[8]-1].getCoord(1);
  
  double z1 = particles[nodes[0]-1].getCoord(2);
  double z2 = particles[nodes[1]-1].getCoord(2);
  double z3 = particles[nodes[2]-1].getCoord(2);
  double z4 = particles[nodes[3]-1].getCoord(2);
  double z9 = particles[nodes[4]-1].getCoord(2);
  double z10 = particles[nodes[5]-1].getCoord(2);
  double z11 = particles[nodes[6]-1].getCoord(2);
  double z12 = particles[nodes[7]-1].getCoord(2);
  double z21 = particles[nodes[8]-1].getCoord(2);


  surfaceNormal[0] = 
    (2.0*(xi/4.0+1.0/8.0)*(eta*eta+eta)*y1+2.0*(xi/4.0-1.0/8.0)*
     (eta*eta+eta)*y2+2.0*(xi/4.0-1.0/8.0)*(eta*eta-eta)*y3+2.0*
     (xi/4.0+1.0/8.0)*(eta*eta-eta)*y4-xi*(eta*eta+eta)*y9+2.0*
     (-xi/2.0+1.0/4.0)*(eta*eta-1.0)*y10-xi*(eta*eta-eta)*y11+2.0*
     (-xi/2.0-1.0/4.0)*(eta*eta-1.0)*y12+2.0*xi*(eta*eta-1.0)*y21)*
    (2.0*(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*z1+2.0*(xi*xi/8.0-xi/8.0)*
     (2.0*eta+1.0)*z2+2.0*(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*z3+2.0*
     (xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*z4+2.0*(-xi*xi/4.0+1.0/4.0)*
     (2.0*eta+1.0)*z9+4.0*(-xi*xi/4.0+xi/4.0)*eta*z10+2.0*
     (-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*z11+4.0*(-xi*xi/4.0-xi/4.0)*eta*
     z12+4.0*(xi*xi/2.0-1.0/2.0)*eta*z21)-
    (2.0*(xi/4.0+1.0/8.0)*(eta*eta+eta)*z1+2.0*(xi/4.0-1.0/8.0)*
     (eta*eta+eta)*z2+2.0*(xi/4.0-1.0/8.0)*(eta*eta-eta)*z3+2.0*
     (xi/4.0+1.0/8.0)*(eta*eta-eta)*z4-xi*(eta*eta+eta)*z9+2.0*
     (-xi/2.0+1.0/4.0)*(eta*eta-1.0)*z10-xi*(eta*eta-eta)*z11+2.0*
     (-xi/2.0-1.0/4.0)*(eta*eta-1.0)*z12+2.0*xi*(eta*eta-1.0)*z21)*
    (2.0*(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*y1+2.0*(xi*xi/8.0-xi/8.0)*
     (2.0*eta+1.0)*y2+2.0*(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*y3+2.0*
     (xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*y4+2.0*(-xi*xi/4.0+1.0/4.0)*
     (2.0*eta+1.0)*y9+4.0*(-xi*xi/4.0+xi/4.0)*eta*y10+2.0*
     (-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*y11+4.0*(-xi*xi/4.0-xi/4.0)*
     eta*y12+4.0*(xi*xi/2.0-1.0/2.0)*eta*y21); 

  surfaceNormal[1] = 
    (2.0*(xi/4.0+1.0/8.0)*(eta*eta+eta)*z1+2.0*(xi/4.0-1.0/8.0)*
     (eta*eta+eta)*z2+2.0*(xi/4.0-1.0/8.0)*(eta*eta-eta)*z3+2.0*
     (xi/4.0+1.0/8.0)*(eta*eta-eta)*z4-xi*(eta*eta+eta)*z9+2.0*
     (-xi/2.0+1.0/4.0)*(eta*eta-1.0)*z10-xi*(eta*eta-eta)*z11+2.0*
     (-xi/2.0-1.0/4.0)*(eta*eta-1.0)*z12+2.0*xi*(eta*eta-1.0)*z21)*
    (2.0*(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*x1+2.0*(xi*xi/8.0-xi/8.0)*
     (2.0*eta+1.0)*x2+2.0*(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*x3+2.0*
     (xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*x4+2.0*(-xi*xi/4.0+1.0/4.0)*
     (2.0*eta+1.0)*x9+4.0*(-xi*xi/4.0+xi/4.0)*eta*x10+2.0*
     (-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*x11+4.0*(-xi*xi/4.0-xi/4.0)*
     eta*x12+4.0*(xi*xi/2.0-1.0/2.0)*eta*x21)-
    (2.0*(xi/4.0+1.0/8.0)*(eta*eta+eta)*x1+2.0*(xi/4.0-1.0/8.0)*
     (eta*eta+eta)*x2+2.0*(xi/4.0-1.0/8.0)*(eta*eta-eta)*x3+2.0*
     (xi/4.0+1.0/8.0)*(eta*eta-eta)*x4-xi*(eta*eta+eta)*x9+2.0*
     (-xi/2.0+1.0/4.0)*(eta*eta-1.0)*x10-xi*(eta*eta-eta)*x11+2.0*
     (-xi/2.0-1.0/4.0)*(eta*eta-1.0)*x12+2.0*xi*(eta*eta-1.0)*x21)*
    (2.0*(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*z1+2.0*(xi*xi/8.0-xi/8.0)*
     (2.0*eta+1.0)*z2+2.0*(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*z3+2.0*
     (xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*z4+2.0*(-xi*xi/4.0+1.0/4.0)*
     (2.0*eta+1.0)*z9+4.0*(-xi*xi/4.0+xi/4.0)*eta*z10+2.0*
     (-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*z11+4.0*(-xi*xi/4.0-xi/4.0)*
     eta*z12+4.0*(xi*xi/2.0-1.0/2.0)*eta*z21);

  surfaceNormal[2] = 
    (2.0*(xi/4.0+1.0/8.0)*(eta*eta+eta)*x1+2.0*(xi/4.0-1.0/8.0)*
     (eta*eta+eta)*x2+2.0*(xi/4.0-1.0/8.0)*(eta*eta-eta)*x3+2.0*
     (xi/4.0+1.0/8.0)*(eta*eta-eta)*x4-xi*(eta*eta+eta)*x9+2.0*
     (-xi/2.0+1.0/4.0)*(eta*eta-1.0)*x10-xi*(eta*eta-eta)*x11+2.0*
     (-xi/2.0-1.0/4.0)*(eta*eta-1.0)*x12+2.0*xi*(eta*eta-1.0)*x21)*
    (2.0*(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*y1+2.0*(xi*xi/8.0-xi/8.0)*
     (2.0*eta+1.0)*y2+2.0*(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*y3+2.0*
     (xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*y4+2.0*(-xi*xi/4.0+1.0/4.0)*
     (2.0*eta+1.0)*y9+4.0*(-xi*xi/4.0+xi/4.0)*eta*y10+2.0*
     (-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*y11+4.0*(-xi*xi/4.0-xi/4.0)*
     eta*y12+4.0*(xi*xi/2.0-1.0/2.0)*eta*y21)-
    (2.0*(xi/4.0+1.0/8.0)*(eta*eta+eta)*y1+2.0*(xi/4.0-1.0/8.0)*
     (eta*eta+eta)*y2+2.0*(xi/4.0-1.0/8.0)*(eta*eta-eta)*y3+2.0*
     (xi/4.0+1.0/8.0)*(eta*eta-eta)*y4-xi*(eta*eta+eta)*y9+2.0*
     (-xi/2.0+1.0/4.0)*(eta*eta-1.0)*y10-xi*(eta*eta-eta)*y11+2.0*
     (-xi/2.0-1.0/4.0)*(eta*eta-1.0)*y12+2.0*xi*(eta*eta-1.0)*y21)*
    (2.0*(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*x1+2.0*(xi*xi/8.0-xi/8.0)*
     (2.0*eta+1.0)*x2+2.0*(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*x3+2.0*
     (xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*x4+2.0*(-xi*xi/4.0+1.0/4.0)*
     (2.0*eta+1.0)*x9+4.0*(-xi*xi/4.0+xi/4.0)*eta*x10+2.0*
     (-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*x11+4.0*(-xi*xi/4.0-xi/4.0)*
     eta*x12+4.0*(xi*xi/2.0-1.0/2.0)*eta*x21); 

  // Calculate the surface
  surface = sqrt(pow(surfaceNormal[0],2)+pow(surfaceNormal[1],2)
		 +pow(surfaceNormal[2],2));

  surfaceNormal[0] = (-1)*surfaceNormal[0]/surface;
  surfaceNormal[1] = (-1)*surfaceNormal[1]/surface;
  surfaceNormal[2] = (-1)*surfaceNormal[2]/surface;

}

/***********************************************************************/
// Calculate the transformation and metric factor for the integration
// of a volume element's volume.
double Cube27ElementTemplate::metricFactorVolume(std::vector<Particle>& particles,
						 intVector& nodes,
						 dbVector& coords,
						 std::ofstream& logFile) {

  using namespace std;

  double det;
  dbMatrix Jmat;

  jacobian(particles,nodes,coords,Jmat,logFile);
  int info =  calcDetDoubleDense(Jmat,det,logFile);

  det = fabs(det);

  if (det < DBL_EPSILON || info != 0) {
    logFile<<"Cube27ElementTemplate::metricFactorVolume integration "
	   <<"weight to small!"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


#ifdef _FEdebugMode_
  logFile<<"********************** Jacobian *********************"<<endl;
  logFile<<"*****************************************************"<<endl;
  double det2;
  dbMatrix Jmat2;
  allocateArray(Jmat2,3,3);

  double xi = coords[0];
  double eta = coords[1];
  double zeta = coords[2];

  double x1 = particles[nodes[0]-1].getCoord(0);
  double x2 = particles[nodes[1]-1].getCoord(0);
  double x3 = particles[nodes[2]-1].getCoord(0);
  double x4 = particles[nodes[3]-1].getCoord(0);
  double x5 = particles[nodes[4]-1].getCoord(0);
  double x6 = particles[nodes[5]-1].getCoord(0);
  double x7 = particles[nodes[6]-1].getCoord(0);
  double x8 = particles[nodes[7]-1].getCoord(0);
  double x9 = particles[nodes[8]-1].getCoord(0);
  double x10 = particles[nodes[9]-1].getCoord(0);
  double x11 = particles[nodes[10]-1].getCoord(0);
  double x12 = particles[nodes[11]-1].getCoord(0);
  double x13 = particles[nodes[12]-1].getCoord(0);
  double x14 = particles[nodes[13]-1].getCoord(0);
  double x15 = particles[nodes[14]-1].getCoord(0);
  double x16 = particles[nodes[15]-1].getCoord(0);
  double x17 = particles[nodes[16]-1].getCoord(0);
  double x18 = particles[nodes[17]-1].getCoord(0);
  double x19 = particles[nodes[18]-1].getCoord(0);
  double x20 = particles[nodes[19]-1].getCoord(0);
  double x21 = particles[nodes[20]-1].getCoord(0);
  double x22 = particles[nodes[21]-1].getCoord(0);
  double x23 = particles[nodes[22]-1].getCoord(0);
  double x24 = particles[nodes[23]-1].getCoord(0);
  double x25 = particles[nodes[24]-1].getCoord(0);
  double x26 = particles[nodes[25]-1].getCoord(0);
  double x27 = particles[nodes[26]-1].getCoord(0);

  double y1 = particles[nodes[0]-1].getCoord(1);
  double y2 = particles[nodes[1]-1].getCoord(1);
  double y3 = particles[nodes[2]-1].getCoord(1);
  double y4 = particles[nodes[3]-1].getCoord(1);
  double y5 = particles[nodes[4]-1].getCoord(1);
  double y6 = particles[nodes[5]-1].getCoord(1);
  double y7 = particles[nodes[6]-1].getCoord(1);
  double y8 = particles[nodes[7]-1].getCoord(1);
  double y9 = particles[nodes[8]-1].getCoord(1);
  double y10 = particles[nodes[9]-1].getCoord(1);
  double y11 = particles[nodes[10]-1].getCoord(1);
  double y12 = particles[nodes[11]-1].getCoord(1);
  double y13 = particles[nodes[12]-1].getCoord(1);
  double y14 = particles[nodes[13]-1].getCoord(1);
  double y15 = particles[nodes[14]-1].getCoord(1);
  double y16 = particles[nodes[15]-1].getCoord(1);
  double y17 = particles[nodes[16]-1].getCoord(1);
  double y18 = particles[nodes[17]-1].getCoord(1);
  double y19 = particles[nodes[18]-1].getCoord(1);
  double y20 = particles[nodes[19]-1].getCoord(1);
  double y21 = particles[nodes[20]-1].getCoord(1);
  double y22 = particles[nodes[21]-1].getCoord(1);
  double y23 = particles[nodes[22]-1].getCoord(1);
  double y24 = particles[nodes[23]-1].getCoord(1);
  double y25 = particles[nodes[24]-1].getCoord(1);
  double y26 = particles[nodes[25]-1].getCoord(1);
  double y27 = particles[nodes[26]-1].getCoord(1);

  double z1 = particles[nodes[0]-1].getCoord(2);
  double z2 = particles[nodes[1]-1].getCoord(2);
  double z3 = particles[nodes[2]-1].getCoord(2);
  double z4 = particles[nodes[3]-1].getCoord(2);
  double z5 = particles[nodes[4]-1].getCoord(2);
  double z6 = particles[nodes[5]-1].getCoord(2);
  double z7 = particles[nodes[6]-1].getCoord(2);
  double z8 = particles[nodes[7]-1].getCoord(2);
  double z9 = particles[nodes[8]-1].getCoord(2);
  double z10 = particles[nodes[9]-1].getCoord(2);
  double z11 = particles[nodes[10]-1].getCoord(2);
  double z12 = particles[nodes[11]-1].getCoord(2);
  double z13 = particles[nodes[12]-1].getCoord(2);
  double z14 = particles[nodes[13]-1].getCoord(2);
  double z15 = particles[nodes[14]-1].getCoord(2);
  double z16 = particles[nodes[15]-1].getCoord(2);
  double z17 = particles[nodes[16]-1].getCoord(2);
  double z18 = particles[nodes[17]-1].getCoord(2);
  double z19 = particles[nodes[18]-1].getCoord(2);
  double z20 = particles[nodes[19]-1].getCoord(2);
  double z21 = particles[nodes[20]-1].getCoord(2);
  double z22 = particles[nodes[21]-1].getCoord(2);
  double z23 = particles[nodes[22]-1].getCoord(2);
  double z24 = particles[nodes[23]-1].getCoord(2);
  double z25 = particles[nodes[24]-1].getCoord(2);
  double z26 = particles[nodes[25]-1].getCoord(2);
  double z27 = particles[nodes[26]-1].getCoord(2);


  // define the Jacobian (computed by Matlab)

Jmat2[0][0] = -2.0*xi*(eta*eta-1.0)*(zeta*zeta-1.0)*x27+xi*(eta*eta-1.0)*(zeta*zeta-zeta)*x26+(xi/4.0+1.0/8.0)*(eta*eta+eta)*(zeta*zeta-zeta)*x5+(xi/4.0-1.0/8.0)*(eta*eta+eta)*(zeta*zeta-zeta)*x6+(-xi/2.0+1.0/4.0)*(eta*eta-eta)*(zeta*zeta-1.0)*x15+xi*(eta*eta-1.0)*(zeta*zeta+zeta)*x21+(-xi/2.0-1.0/4.0)*(eta*eta-1.0)*(zeta*zeta+zeta)*x12+(xi/4.0-1.0/8.0)*(eta*eta-eta)*(zeta*zeta-zeta)*x7-xi*(eta*eta+eta)*(zeta*zeta-zeta)*x17/2.0-xi*(eta*eta+eta)*(zeta*zeta+zeta)*x9/2.0+(-xi/2.0-1.0/4.0)*(eta*eta-eta)*(zeta*zeta-1.0)*x16+(xi-1.0/2.0)*(eta*eta-1.0)*(zeta*zeta-1.0)*x23+(-xi/2.0-1.0/4.0)*(eta*eta+eta)*(zeta*zeta-1.0)*x13+(xi/4.0+1.0/8.0)*(eta*eta+eta)*(zeta*zeta+zeta)*x1+(-xi/2.0+1.0/4.0)*(eta*eta-1.0)*(zeta*zeta-zeta)*x18+(xi/4.0-1.0/8.0)*(eta*eta-eta)*(zeta*zeta+zeta)*x3+(-xi/2.0+1.0/4.0)*(eta*eta+eta)*(zeta*zeta-1.0)*x14+(xi+1.0/2.0)*(eta*eta-1.0)*(zeta*zeta-1.0)*x25+(-xi/2.0-1.0/4.0)*(eta*eta-1.0)*(zeta*zeta-zeta)*x20+(xi/4.0-1.0/8.0)*(eta*eta+eta)*(zeta*zeta+zeta)*x2+xi*(eta*eta+eta)*(zeta*zeta-1.0)*x22-xi*(eta*eta-eta)*(zeta*zeta-zeta)*x19/2.0+xi*(eta*eta-eta)*(zeta*zeta-1.0)*x24-xi*(eta*eta-eta)*(zeta*zeta+zeta)*x11/2.0+(-xi/2.0+1.0/4.0)*(eta*eta-1.0)*(zeta*zeta+zeta)*x10+(xi/4.0+1.0/8.0)*(eta*eta-eta)*(zeta*zeta-zeta)*x8+(xi/4.0+1.0/8.0)*(eta*eta-eta)*(zeta*zeta+zeta)*x4;      

Jmat2[0][1] = (-xi*xi/4.0+xi/4.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*x15+2.0*(xi*xi/2.0-1.0/2.0)*eta*(zeta*zeta+zeta)*x21+2.0*(-xi*xi/4.0-xi/4.0)*eta*(zeta*zeta+zeta)*x12+(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*x7+(-xi*xi/4.0+1.0/4.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*x17+(-xi*xi/4.0+1.0/4.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*x9+(-xi*xi/4.0-xi/4.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*x16+2.0*(xi*xi/2.0-xi/2.0)*eta*(zeta*zeta-1.0)*x23+2.0*(xi*xi/2.0+xi/2.0)*eta*(zeta*zeta-1.0)*x25+2.0*(-xi*xi/4.0-xi/4.0)*eta*(zeta*zeta-zeta)*x20+(xi*xi/8.0-xi/8.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*x2+(xi*xi/2.0-1.0/2.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*x22+(-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*x19+(xi*xi/2.0-1.0/2.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*x24+(-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*x11+2.0*(-xi*xi/4.0+xi/4.0)*eta*(zeta*zeta+zeta)*x10+(xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*x8+(-xi*xi/4.0-xi/4.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*x13+(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*x1+2.0*(-xi*xi/4.0+xi/4.0)*eta*(zeta*zeta-zeta)*x18+(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*x3+(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*x5+(xi*xi/8.0-xi/8.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*x6+2.0*(-xi*xi+1.0)*eta*(zeta*zeta-1.0)*x27+2.0*(xi*xi/2.0-1.0/2.0)*eta*(zeta*zeta-zeta)*x26+(xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*x4+(-xi*xi/4.0+xi/4.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*x14;      

Jmat2[0][2] = 2.0*(-xi*xi/4.0-xi/4.0)*(eta*eta-eta)*zeta*x16+2.0*(xi*xi/2.0-xi/2.0)*(eta*eta-1.0)*zeta*x23+(xi*xi/8.0+xi/8.0)*(eta*eta+eta)*(2.0*zeta-1.0)*x5+(xi*xi/8.0-xi/8.0)*(eta*eta+eta)*(2.0*zeta-1.0)*x6+2.0*(-xi*xi/4.0+xi/4.0)*(eta*eta-eta)*zeta*x15+(xi*xi/2.0-1.0/2.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*x21+(-xi*xi/4.0-xi/4.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*x12+(xi*xi/8.0-xi/8.0)*(eta*eta-eta)*(2.0*zeta-1.0)*x7+(-xi*xi/4.0+1.0/4.0)*(eta*eta+eta)*(2.0*zeta-1.0)*x17+(-xi*xi/4.0+1.0/4.0)*(eta*eta+eta)*(2.0*zeta+1.0)*x9+2.0*(xi*xi/2.0-1.0/2.0)*(eta*eta+eta)*zeta*x22+(-xi*xi/4.0+1.0/4.0)*(eta*eta-eta)*(2.0*zeta-1.0)*x19+2.0*(xi*xi/2.0-1.0/2.0)*(eta*eta-eta)*zeta*x24+(-xi*xi/4.0+1.0/4.0)*(eta*eta-eta)*(2.0*zeta+1.0)*x11+2.0*(-xi*xi/4.0+xi/4.0)*(eta*eta+eta)*zeta*x14+2.0*(xi*xi/2.0+xi/2.0)*(eta*eta-1.0)*zeta*x25+(-xi*xi/4.0-xi/4.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*x20+(xi*xi/8.0-xi/8.0)*(eta*eta+eta)*(2.0*zeta+1.0)*x2+(-xi*xi/4.0+xi/4.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*x10+(xi*xi/8.0+xi/8.0)*(eta*eta-eta)*(2.0*zeta-1.0)*x8+2.0*(-xi*xi/4.0-xi/4.0)*(eta*eta+eta)*zeta*x13+(xi*xi/8.0+xi/8.0)*(eta*eta+eta)*(2.0*zeta+1.0)*x1+(-xi*xi/4.0+xi/4.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*x18+(xi*xi/8.0-xi/8.0)*(eta*eta-eta)*(2.0*zeta+1.0)*x3+2.0*(-xi*xi+1.0)*(eta*eta-1.0)*zeta*x27+(xi*xi/2.0-1.0/2.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*x26+(xi*xi/8.0+xi/8.0)*(eta*eta-eta)*(2.0*zeta+1.0)*x4;      

Jmat2[1][0] = -xi*(eta*eta+eta)*(zeta*zeta+zeta)*y9/2.0+(-xi/2.0+1.0/4.0)*(eta*eta+eta)*(zeta*zeta-1.0)*y14+(xi/4.0-1.0/8.0)*(eta*eta-eta)*(zeta*zeta-zeta)*y7-xi*(eta*eta+eta)*(zeta*zeta-zeta)*y17/2.0-2.0*xi*(eta*eta-1.0)*(zeta*zeta-1.0)*y27-xi*(eta*eta-eta)*(zeta*zeta-zeta)*y19/2.0+(-xi/2.0-1.0/4.0)*(eta*eta-1.0)*(zeta*zeta-zeta)*y20+(-xi/2.0+1.0/4.0)*(eta*eta-1.0)*(zeta*zeta-zeta)*y18+xi*(eta*eta-1.0)*(zeta*zeta+zeta)*y21+(-xi/2.0+1.0/4.0)*(eta*eta-eta)*(zeta*zeta-1.0)*y15+(-xi/2.0-1.0/4.0)*(eta*eta+eta)*(zeta*zeta-1.0)*y13+xi*(eta*eta+eta)*(zeta*zeta-1.0)*y22+xi*(eta*eta-1.0)*(zeta*zeta-zeta)*y26-xi*(eta*eta-eta)*(zeta*zeta+zeta)*y11/2.0+xi*(eta*eta-eta)*(zeta*zeta-1.0)*y24+(xi/4.0+1.0/8.0)*(eta*eta+eta)*(zeta*zeta+zeta)*y1+(xi/4.0-1.0/8.0)*(eta*eta+eta)*(zeta*zeta+zeta)*y2+(xi/4.0-1.0/8.0)*(eta*eta-eta)*(zeta*zeta+zeta)*y3+(xi/4.0+1.0/8.0)*(eta*eta-eta)*(zeta*zeta+zeta)*y4+(xi-1.0/2.0)*(eta*eta-1.0)*(zeta*zeta-1.0)*y23+(xi+1.0/2.0)*(eta*eta-1.0)*(zeta*zeta-1.0)*y25+(xi/4.0+1.0/8.0)*(eta*eta+eta)*(zeta*zeta-zeta)*y5+(-xi/2.0-1.0/4.0)*(eta*eta-1.0)*(zeta*zeta+zeta)*y12+(xi/4.0-1.0/8.0)*(eta*eta+eta)*(zeta*zeta-zeta)*y6+(-xi/2.0-1.0/4.0)*(eta*eta-eta)*(zeta*zeta-1.0)*y16+(-xi/2.0+1.0/4.0)*(eta*eta-1.0)*(zeta*zeta+zeta)*y10+(xi/4.0+1.0/8.0)*(eta*eta-eta)*(zeta*zeta-zeta)*y8;      

Jmat2[1][1] = (-xi*xi/4.0-xi/4.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*y16+2.0*(-xi*xi/4.0-xi/4.0)*eta*(zeta*zeta-zeta)*y20+2.0*(-xi*xi/4.0+xi/4.0)*eta*(zeta*zeta-zeta)*y18+2.0*(-xi*xi/4.0+xi/4.0)*eta*(zeta*zeta+zeta)*y10+(xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*y8+(-xi*xi/4.0+1.0/4.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*y9+(-xi*xi/4.0+xi/4.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*y14+(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*y7+(-xi*xi/4.0+1.0/4.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*y17+2.0*(-xi*xi+1.0)*eta*(zeta*zeta-1.0)*y27+(-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*y19+2.0*(xi*xi/2.0+xi/2.0)*eta*(zeta*zeta-1.0)*y25+(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*y5+2.0*(-xi*xi/4.0-xi/4.0)*eta*(zeta*zeta+zeta)*y12+(xi*xi/8.0-xi/8.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*y6+2.0*(xi*xi/2.0-1.0/2.0)*eta*(zeta*zeta-zeta)*y26+(-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*y11+(xi*xi/2.0-1.0/2.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*y24+2.0*(xi*xi/2.0-xi/2.0)*eta*(zeta*zeta-1.0)*y23+(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*y3+(xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*y4+2.0*(xi*xi/2.0-1.0/2.0)*eta*(zeta*zeta+zeta)*y21+(-xi*xi/4.0+xi/4.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*y15+(-xi*xi/4.0-xi/4.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*y13+(xi*xi/2.0-1.0/2.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*y22+(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*y1+(xi*xi/8.0-xi/8.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*y2;      

Jmat2[1][2] = 2.0*(-xi*xi+1.0)*(eta*eta-1.0)*zeta*y27+(-xi*xi/4.0+1.0/4.0)*(eta*eta-eta)*(2.0*zeta-1.0)*y19+(-xi*xi/4.0-xi/4.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*y20+(-xi*xi/4.0+xi/4.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*y18+2.0*(-xi*xi/4.0+xi/4.0)*(eta*eta+eta)*zeta*y14+(xi*xi/8.0-xi/8.0)*(eta*eta-eta)*(2.0*zeta-1.0)*y7+(-xi*xi/4.0+1.0/4.0)*(eta*eta+eta)*(2.0*zeta-1.0)*y17+2.0*(xi*xi/2.0+xi/2.0)*(eta*eta-1.0)*zeta*y25+(xi*xi/8.0+xi/8.0)*(eta*eta+eta)*(2.0*zeta-1.0)*y5+(-xi*xi/4.0-xi/4.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*y12+(xi*xi/8.0-xi/8.0)*(eta*eta+eta)*(2.0*zeta-1.0)*y6+2.0*(-xi*xi/4.0-xi/4.0)*(eta*eta-eta)*zeta*y16+(-xi*xi/4.0+xi/4.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*y10+(xi*xi/8.0+xi/8.0)*(eta*eta-eta)*(2.0*zeta-1.0)*y8+(-xi*xi/4.0+1.0/4.0)*(eta*eta+eta)*(2.0*zeta+1.0)*y9+2.0*(xi*xi/2.0-1.0/2.0)*(eta*eta-eta)*zeta*y24+2.0*(xi*xi/2.0-xi/2.0)*(eta*eta-1.0)*zeta*y23+2.0*(xi*xi/2.0-1.0/2.0)*(eta*eta+eta)*zeta*y22+(xi*xi/2.0-1.0/2.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*y26+(-xi*xi/4.0+1.0/4.0)*(eta*eta-eta)*(2.0*zeta+1.0)*y11+(xi*xi/2.0-1.0/2.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*y21+2.0*(-xi*xi/4.0+xi/4.0)*(eta*eta-eta)*zeta*y15+2.0*(-xi*xi/4.0-xi/4.0)*(eta*eta+eta)*zeta*y13+(xi*xi/8.0+xi/8.0)*(eta*eta+eta)*(2.0*zeta+1.0)*y1+(xi*xi/8.0-xi/8.0)*(eta*eta+eta)*(2.0*zeta+1.0)*y2+(xi*xi/8.0-xi/8.0)*(eta*eta-eta)*(2.0*zeta+1.0)*y3+(xi*xi/8.0+xi/8.0)*(eta*eta-eta)*(2.0*zeta+1.0)*y4;      

Jmat2[2][0] = (xi/4.0-1.0/8.0)*(eta*eta-eta)*(zeta*zeta-zeta)*z7+(xi/4.0+1.0/8.0)*(eta*eta+eta)*(zeta*zeta+zeta)*z1+(xi/4.0-1.0/8.0)*(eta*eta+eta)*(zeta*zeta+zeta)*z2+(xi/4.0-1.0/8.0)*(eta*eta-eta)*(zeta*zeta+zeta)*z3+(xi/4.0+1.0/8.0)*(eta*eta-eta)*(zeta*zeta+zeta)*z4+(xi/4.0+1.0/8.0)*(eta*eta+eta)*(zeta*zeta-zeta)*z5+(xi/4.0-1.0/8.0)*(eta*eta+eta)*(zeta*zeta-zeta)*z6+xi*(eta*eta+eta)*(zeta*zeta-1.0)*z22-xi*(eta*eta-eta)*(zeta*zeta-zeta)*z19/2.0+(xi/4.0+1.0/8.0)*(eta*eta-eta)*(zeta*zeta-zeta)*z8+xi*(eta*eta-1.0)*(zeta*zeta-zeta)*z26+(-xi/2.0-1.0/4.0)*(eta*eta-1.0)*(zeta*zeta-zeta)*z20+xi*(eta*eta-1.0)*(zeta*zeta+zeta)*z21+(-xi/2.0+1.0/4.0)*(eta*eta-1.0)*(zeta*zeta-zeta)*z18+(-xi/2.0-1.0/4.0)*(eta*eta-eta)*(zeta*zeta-1.0)*z16-2.0*xi*(eta*eta-1.0)*(zeta*zeta-1.0)*z27+xi*(eta*eta-eta)*(zeta*zeta-1.0)*z24+(-xi/2.0+1.0/4.0)*(eta*eta-eta)*(zeta*zeta-1.0)*z15+(-xi/2.0+1.0/4.0)*(eta*eta+eta)*(zeta*zeta-1.0)*z14-xi*(eta*eta+eta)*(zeta*zeta+zeta)*z9/2.0+(-xi/2.0+1.0/4.0)*(eta*eta-1.0)*(zeta*zeta+zeta)*z10-xi*(eta*eta-eta)*(zeta*zeta+zeta)*z11/2.0+(-xi/2.0-1.0/4.0)*(eta*eta-1.0)*(zeta*zeta+zeta)*z12+(-xi/2.0-1.0/4.0)*(eta*eta+eta)*(zeta*zeta-1.0)*z13-xi*(eta*eta+eta)*(zeta*zeta-zeta)*z17/2.0+(xi-1.0/2.0)*(eta*eta-1.0)*(zeta*zeta-1.0)*z23+(xi+1.0/2.0)*(eta*eta-1.0)*(zeta*zeta-1.0)*z25;      

Jmat2[2][1] = (xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*z7+(xi*xi/8.0-xi/8.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*z2+(xi*xi/8.0-xi/8.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*z3+(xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*z4+(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*z5+(xi*xi/8.0-xi/8.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*z6+2.0*(-xi*xi/4.0-xi/4.0)*eta*(zeta*zeta+zeta)*z12+(-xi*xi/4.0-xi/4.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*z13+(-xi*xi/4.0+1.0/4.0)*(2.0*eta+1.0)*(zeta*zeta-zeta)*z17+2.0*(xi*xi/2.0-xi/2.0)*eta*(zeta*zeta-1.0)*z23+2.0*(xi*xi/2.0+xi/2.0)*eta*(zeta*zeta-1.0)*z25+(xi*xi/8.0+xi/8.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*z1+2.0*(-xi*xi/4.0-xi/4.0)*eta*(zeta*zeta-zeta)*z20+2.0*(xi*xi/2.0-1.0/2.0)*eta*(zeta*zeta+zeta)*z21+2.0*(-xi*xi/4.0+xi/4.0)*eta*(zeta*zeta-zeta)*z18+(-xi*xi/4.0+xi/4.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*z15+(-xi*xi/4.0+xi/4.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*z14+2.0*(xi*xi/2.0-1.0/2.0)*eta*(zeta*zeta-zeta)*z26+(-xi*xi/4.0-xi/4.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*z16+2.0*(-xi*xi+1.0)*eta*(zeta*zeta-1.0)*z27+(xi*xi/2.0-1.0/2.0)*(2.0*eta-1.0)*(zeta*zeta-1.0)*z24+(xi*xi/2.0-1.0/2.0)*(2.0*eta+1.0)*(zeta*zeta-1.0)*z22+(-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*z19+(xi*xi/8.0+xi/8.0)*(2.0*eta-1.0)*(zeta*zeta-zeta)*z8+(-xi*xi/4.0+1.0/4.0)*(2.0*eta+1.0)*(zeta*zeta+zeta)*z9+2.0*(-xi*xi/4.0+xi/4.0)*eta*(zeta*zeta+zeta)*z10+(-xi*xi/4.0+1.0/4.0)*(2.0*eta-1.0)*(zeta*zeta+zeta)*z11;      

Jmat2[2][2] = (xi*xi/8.0-xi/8.0)*(eta*eta+eta)*(2.0*zeta-1.0)*z6+(xi*xi/8.0-xi/8.0)*(eta*eta-eta)*(2.0*zeta-1.0)*z7+2.0*(xi*xi/2.0+xi/2.0)*(eta*eta-1.0)*zeta*z25+(xi*xi/8.0+xi/8.0)*(eta*eta+eta)*(2.0*zeta+1.0)*z1+(xi*xi/8.0-xi/8.0)*(eta*eta+eta)*(2.0*zeta+1.0)*z2+(xi*xi/8.0-xi/8.0)*(eta*eta-eta)*(2.0*zeta+1.0)*z3+(xi*xi/8.0+xi/8.0)*(eta*eta-eta)*(2.0*zeta+1.0)*z4+(xi*xi/8.0+xi/8.0)*(eta*eta+eta)*(2.0*zeta-1.0)*z5+(-xi*xi/4.0-xi/4.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*z12+2.0*(-xi*xi/4.0-xi/4.0)*(eta*eta+eta)*zeta*z13+(-xi*xi/4.0+1.0/4.0)*(eta*eta+eta)*(2.0*zeta-1.0)*z17+2.0*(xi*xi/2.0-xi/2.0)*(eta*eta-1.0)*zeta*z23+2.0*(xi*xi/2.0-1.0/2.0)*(eta*eta+eta)*zeta*z22+(-xi*xi/4.0+1.0/4.0)*(eta*eta-eta)*(2.0*zeta-1.0)*z19+(xi*xi/8.0+xi/8.0)*(eta*eta-eta)*(2.0*zeta-1.0)*z8+(-xi*xi/4.0+1.0/4.0)*(eta*eta+eta)*(2.0*zeta+1.0)*z9+(-xi*xi/4.0+xi/4.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*z10+(-xi*xi/4.0+1.0/4.0)*(eta*eta-eta)*(2.0*zeta+1.0)*z11+2.0*(-xi*xi/4.0+xi/4.0)*(eta*eta-eta)*zeta*z15+2.0*(-xi*xi/4.0+xi/4.0)*(eta*eta+eta)*zeta*z14+(xi*xi/2.0-1.0/2.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*z26+(-xi*xi/4.0-xi/4.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*z20+(xi*xi/2.0-1.0/2.0)*(eta*eta-1.0)*(2.0*zeta+1.0)*z21+(-xi*xi/4.0+xi/4.0)*(eta*eta-1.0)*(2.0*zeta-1.0)*z18+2.0*(-xi*xi/4.0-xi/4.0)*(eta*eta-eta)*zeta*z16+2.0*(-xi*xi+1.0)*(eta*eta-1.0)*zeta*z27+2.0*(xi*xi/2.0-1.0/2.0)*(eta*eta-eta)*zeta*z24;

  for(int i=0;i<Jmat.size();i++)
    for(int j=0;j<Jmat[i].size();j++)
      logFile<<"J["<<i<<"]["<<j<<"] = "<<Jmat[i][j]<<"=?"
	     <<Jmat2[i][j]<<endl;
  info =  calcDetDoubleDense(Jmat2,det2,logFile);
  logFile<<"det="<<det<<"=?"<<det2<<endl;
#endif


  return(det);
}
  
