// Calculate the window function and its derivations at a requested 
// point.

#include "WindowFunctionSets.h"


/***********************************************************************/
// closed Gauss spline
//
// W(r,h) = e^{-1/(1-r^2)}*(1-r^2)^8

/***********************************************************************/
// Calculate the window function in 3-D.
void CubicSplineWindowFunc::calcWinFunc(double& x,double& y,double& z,
					double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 0.5) {
    if(yabs <= 0.5) {
      if(zabs <= 0.5) {
	wFunc = W1(x)*W1(y)*W1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	wFunc = W1(x)*W1(y)*W2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else if(yabs > 0.5 && yabs <= 1.0) {
      if(zabs <= 0.5) {
	wFunc = W1(x)*W2(y)*W1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	wFunc = W1(x)*W2(y)*W2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else {
      cerr<<"While calculating window functions chosen point is out of " 
	  <<"particle range!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else if(xabs > 0.5 && xabs <= 1.0) {
    if(yabs <= 0.5) {
      if(zabs <= 0.5) {
	wFunc = W2(x)*W1(y)*W1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	wFunc = W2(x)*W1(y)*W2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else if(yabs > 0.5 && yabs <= 1.0) {
      if(zabs <= 0.5) {
	wFunc = W2(x)*W2(y)*W1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	wFunc = W2(x)*W2(y)*W2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else {
      cerr<<"While calculating window functions chosen point is out of " 
	  <<"particle range!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void CubicSplineWindowFunc::calcWinFunc1stDerivs(double& x,double& y,
						 double& z,double& xDWFunc,
						 double& yDWFunc,
						 double& zDWFunc) {

  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 0.5) {
    if(yabs <= 0.5) {
      if(zabs <= 0.5) {
	xDWFunc = dW1(x)*W1(y)*W1(z);
	yDWFunc = W1(x)*dW1(y)*W1(z);
	zDWFunc = W1(x)*W1(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xDWFunc = dW1(x)*W1(y)*W2(z);	
	yDWFunc = W1(x)*dW1(y)*W2(z);	
	zDWFunc = W1(x)*W1(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else if(yabs > 0.5 && yabs <= 1.0) {
      if(zabs<= 0.5) {
	xDWFunc = dW1(x)*W2(y)*W1(z);
	yDWFunc = W1(x)*dW2(y)*W1(z);
	zDWFunc = W1(x)*W2(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xDWFunc = dW1(x)*W2(y)*W2(z);
	yDWFunc = W1(x)*dW2(y)*W2(z);
	zDWFunc = W1(x)*W2(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else {
      cerr<<"While calculating window functions chosen point is out of " 
	  <<"particle range!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else if(xabs > 0.5 && xabs <= 1.0) {
    if(yabs <= 0.5) {
      if(zabs <= 0.5) {
	xDWFunc = dW2(x)*W1(y)*W1(z);
	yDWFunc = W2(x)*dW1(y)*W1(z);
	zDWFunc = W2(x)*W1(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xDWFunc = dW2(x)*W1(y)*W2(z);
	yDWFunc = W2(x)*dW1(y)*W2(z);
	zDWFunc = W2(x)*W1(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else if(yabs > 0.5 && yabs <= 1.0) {
      if(zabs <= 0.5) {
	xDWFunc = dW2(x)*W2(y)*W1(z);
	yDWFunc = W2(x)*dW2(y)*W1(z);
	zDWFunc = W2(x)*W2(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xDWFunc = dW2(x)*W2(y)*W2(z);
	yDWFunc = W2(x)*dW2(y)*W2(z);
	zDWFunc = W2(x)*W2(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else {
      cerr<<"While calculating window functions chosen point is out of " 
	  <<"particle range!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void CubicSplineWindowFunc::calcWinFunc2ndDerivs(double& x,double& y,
						 double& z,
						 double& xxDWFunc,
						 double& yyDWFunc,
						 double& zzDWFunc,
						 double& xyDWFunc,
						 double& yzDWFunc,
						 double& zxDWFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 0.5) {
    if(yabs <= 0.5) {
      if(zabs <= 0.5) {
	xxDWFunc = d2W1(x)*W1(y)*W1(z);
	yyDWFunc = W1(x)*d2W1(y)*W1(z);
	zzDWFunc = W1(x)*W1(y)*d2W1(z);

	xyDWFunc = dW1(x)*dW1(y)*W1(z);
	yzDWFunc = W1(x)*dW1(y)*dW1(z);
	zxDWFunc = dW1(x)*W1(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xxDWFunc = d2W1(x)*W1(y)*W2(z);	
	yyDWFunc = W1(x)*d2W1(y)*W2(z);	
	zzDWFunc = W1(x)*W1(y)*d2W2(z);

	xyDWFunc = dW1(x)*dW1(y)*W2(z);	
	yzDWFunc = W1(x)*dW1(y)*dW2(z);	
	zxDWFunc = dW1(x)*W1(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else if(yabs > 0.5 && yabs <= 1.0) {
      if(zabs<= 0.5) {
	xxDWFunc = d2W1(x)*W2(y)*W1(z);
	yyDWFunc = W1(x)*d2W2(y)*W1(z);
	zzDWFunc = W1(x)*W2(y)*d2W1(z);

	xyDWFunc = dW1(x)*dW2(y)*W1(z);
	yzDWFunc = W1(x)*dW2(y)*dW1(z);
	zxDWFunc = dW1(x)*W2(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xxDWFunc = d2W1(x)*W2(y)*W2(z);
	yyDWFunc = W1(x)*d2W2(y)*W2(z);
	zzDWFunc = W1(x)*W2(y)*d2W2(z);

	xyDWFunc = dW1(x)*dW2(y)*W2(z);
	yzDWFunc = W1(x)*dW2(y)*dW2(z);
	zxDWFunc = dW1(x)*W2(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else {
      cerr<<"While calculating window functions chosen point is out of " 
	  <<"particle range!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else if(xabs > 0.5 && xabs <= 1.0) {
    if(yabs <= 0.5) {
      if(zabs <= 0.5) {
	xxDWFunc = d2W2(x)*W1(y)*W1(z);
	yyDWFunc = W2(x)*d2W1(y)*W1(z);
	zzDWFunc = W2(x)*W1(y)*d2W1(z);

	xyDWFunc = dW2(x)*dW1(y)*W1(z);
	yzDWFunc = W2(x)*dW1(y)*dW1(z);
	zxDWFunc = dW2(x)*W1(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xxDWFunc = d2W2(x)*W1(y)*W2(z);
	yyDWFunc = W2(x)*d2W1(y)*W2(z);
	zzDWFunc = W2(x)*W1(y)*d2W2(z);

	xyDWFunc = dW2(x)*dW1(y)*W2(z);
	yzDWFunc = W2(x)*dW1(y)*dW2(z);
	zxDWFunc = dW2(x)*W1(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else if(yabs > 0.5 && yabs <= 1.0) {
      if(zabs <= 0.5) {
	xxDWFunc = d2W2(x)*W2(y)*W1(z);
	yyDWFunc = W2(x)*d2W2(y)*W1(z);
	zzDWFunc = W2(x)*W2(y)*d2W1(z);

	xyDWFunc = dW2(x)*dW2(y)*W1(z);
	yzDWFunc = W2(x)*dW2(y)*dW1(z);
	zxDWFunc = dW2(x)*W2(y)*dW1(z);
      }
      else if(zabs > 0.5 && zabs <= 1.0) {
	xxDWFunc = d2W2(x)*W2(y)*W2(z);
	yyDWFunc = W2(x)*d2W2(y)*W2(z);
	zzDWFunc = W2(x)*W2(y)*d2W2(z);

	xyDWFunc = dW2(x)*dW2(y)*W2(z);
	yzDWFunc = W2(x)*dW2(y)*dW2(z);
	zxDWFunc = dW2(x)*W2(y)*dW2(z);
      }
      else {
	cerr<<"While calculating window functions chosen point is out of " 
	    <<"particle range!"<<endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }
    }
    else {
      cerr<<"While calculating window functions chosen point is out of " 
	  <<"particle range!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}


/***********************************************************************/
// The cubic spline.
double CubicSplineWindowFunc::W1(double& rad) {

  double rabs = fabs(rad);

  return (2.0/3.0-4.0*pow(rabs,2)+4.0*pow(rabs,3));
}

double CubicSplineWindowFunc::W2(double& rad) {

  double rabs = fabs(rad);

  return (4.0/3.0-4.0*rabs+4.0*pow(rabs,2)-4.0/3.0*pow(rabs,3));
}


// First order derivation of the cubic spline.
double CubicSplineWindowFunc::dW1(double& rad) {

  double rabs = fabs(rad);

  return(-8.0*rad+12.0*rabs*rad);
  
}

double CubicSplineWindowFunc::dW2(double& rad) {

  double rabs = fabs(rad);

  double sign;

  if(rabs > 0)
    sign = rad/rabs;

  else
    sign = 1.0;

  return (-4.0*sign+8.0*rad-4.0*rabs*rad);

}

// Second order derivation of the cubic spline.
double CubicSplineWindowFunc::d2W1(double& rad) {

  double rabs = fabs(rad);

  return (-8.0+24.0*rabs);
}
double CubicSplineWindowFunc::d2W2(double& rad) {

  double rabs = fabs(rad);

  return (8.0-8.0*rabs);

}

/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void QuarticSplineWindowFunc::calcWinFunc(double& x,double& y,double& z,
					  double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1)

    wFunc = W(x)*W(y)*W(z);
 
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void QuarticSplineWindowFunc::calcWinFunc1stDerivs(double& x,double& y,
						   double& z,double& xDWFunc,
						   double& yDWFunc,
						   double& zDWFunc) {

  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xDWFunc = dW(x)*W(y)*W(z);
    yDWFunc = W(x)*dW(y)*W(z);
    zDWFunc = W(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void QuarticSplineWindowFunc::calcWinFunc2ndDerivs(double& x,double& y,
						   double& z,
						   double& xxDWFunc,
						   double& yyDWFunc,
						   double& zzDWFunc,
						   double& xyDWFunc,
						   double& yzDWFunc,
						   double& zxDWFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xxDWFunc = d2W(x)*W(y)*W(z);
    yyDWFunc = W(x)*d2W(y)*W(z);
    zzDWFunc = W(x)*W(y)*d2W(z);
    
    xyDWFunc = dW(x)*dW(y)*W(z);
    yzDWFunc = W(x)*dW(y)*dW(z);
    zxDWFunc = dW(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}


/***********************************************************************/
// The quartic spline.
double QuarticSplineWindowFunc::W(double& rad) {
  
  double rabs = fabs(rad);

  return(1.0-6.0*pow(rabs,2)+8.0*pow(rabs,3)-3.0*pow(rabs,4));
}

// First order derivation of the quartic spline.
double QuarticSplineWindowFunc::dW(double& rad) {

  double rabs = fabs(rad);

  return(-12.0*rad+24.0*rabs*rad-12.0*pow(rad,2)*rad);
}

// Second order derivation of the quartic spline.
double QuarticSplineWindowFunc::d2W(double& rad) {

  double rabs = fabs(rad);

  return (-12.0+48.0*rabs-36.0*pow(rabs,2));
}

/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void ExponentialWindowFunc::calcWinFunc(double& x,double& y,double& z,
					  double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs < 1 && yabs < 1 && zabs < 1)

    wFunc = W(x)*W(y)*W(z);

  else if(xabs == 1 || yabs == 1 || zabs == 1)

    wFunc = 0;
 
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void ExponentialWindowFunc::calcWinFunc1stDerivs(double& x,double& y,
						 double& z,double& xDWFunc,
						 double& yDWFunc,
						 double& zDWFunc) {

  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs < 1 && yabs < 1 && zabs < 1) {
    xDWFunc = dW(x)*W(y)*W(z);
    yDWFunc = W(x)*dW(y)*W(z);
    zDWFunc = W(x)*W(y)*dW(z);
  }

  else if(xabs == 1 || yabs == 1 || zabs == 1) {
    xDWFunc = 0;
    yDWFunc = 0;
    zDWFunc = 0;
  }

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void ExponentialWindowFunc::calcWinFunc2ndDerivs(double& x,double& y,
						   double& z,
						   double& xxDWFunc,
						   double& yyDWFunc,
						   double& zzDWFunc,
						   double& xyDWFunc,
						   double& yzDWFunc,
						   double& zxDWFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs < 1 && yabs < 1 && zabs < 1) {
    xxDWFunc = d2W(x)*W(y)*W(z);
    yyDWFunc = W(x)*d2W(y)*W(z);
    zzDWFunc = W(x)*W(y)*d2W(z);
    
    xyDWFunc = dW(x)*dW(y)*W(z);
    yzDWFunc = W(x)*dW(y)*dW(z);
    zxDWFunc = dW(x)*W(y)*dW(z);
  }
  else if(xabs == 1 || yabs == 1 || zabs == 1) {
    xxDWFunc = 0;
    yyDWFunc = 0;
    zzDWFunc = 0;
    
    xyDWFunc = 0;
    yzDWFunc = 0;
    zxDWFunc = 0;
  }

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}


/***********************************************************************/
// The exponential spline.
double ExponentialWindowFunc::W(double& rad) {
  
  double factor = 5.0;
  double eTerm = exp((-1.0)*pow(factor*rad,2));

  return(eTerm);
}

// First order derivation of the exponential spline.
double ExponentialWindowFunc::dW(double& rad) {

  double factor = 5.0;
  double eTerm = exp((-1.0)*pow(factor*rad,2));

  return((-2.0)*factor*rad*eTerm);
}

// Second order derivation of the exponential spline.
double ExponentialWindowFunc::d2W(double& rad) {

  double factor = 5.0;
  double eTerm = exp((-1.0)*pow(factor*rad,2));

  return(2.0*factor*eTerm*(2.0*factor*pow(rad,2)-1.0));
}


/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void GaussWindowFunc::calcWinFunc(double& x,double& y,double& z,
				  double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs < 1 && yabs < 1 && zabs < 1)

    wFunc = W(x)*W(y)*W(z);

  else if(xabs == 1 || yabs == 1 || zabs == 1)

    wFunc = 0;
 
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void GaussWindowFunc::calcWinFunc1stDerivs(double& x,double& y,
					   double& z,double& xDWFunc,
					   double& yDWFunc,
					   double& zDWFunc) {

  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs < 1 && yabs < 1 && zabs < 1) {
    xDWFunc = dW(x)*W(y)*W(z);
    yDWFunc = W(x)*dW(y)*W(z);
    zDWFunc = W(x)*W(y)*dW(z);
  }

  else if(xabs == 1 || yabs == 1 || zabs == 1) {

    xDWFunc = 0;
    yDWFunc = 0;
    zDWFunc = 0;
  }

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void GaussWindowFunc::calcWinFunc2ndDerivs(double& x,double& y,
					   double& z,
					   double& xxDWFunc,
					   double& yyDWFunc,
					   double& zzDWFunc,
					   double& xyDWFunc,
					   double& yzDWFunc,
					   double& zxDWFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs < 1 && yabs < 1 && zabs < 1) {
    xxDWFunc = d2W(x)*W(y)*W(z);
    yyDWFunc = W(x)*d2W(y)*W(z);
    zzDWFunc = W(x)*W(y)*d2W(z);
    
    xyDWFunc = dW(x)*dW(y)*W(z);
    yzDWFunc = W(x)*dW(y)*dW(z);
    zxDWFunc = dW(x)*W(y)*dW(z);
  }

  else if(xabs == 1 || yabs == 1 || zabs == 1) {

    xxDWFunc = 0;
    yyDWFunc = 0;
    zzDWFunc = 0;

    xyDWFunc = 0;
    yzDWFunc = 0;
    zxDWFunc = 0;
  }

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
// The gauss spline.
double GaussWindowFunc::W(double& rad) {
  
  double c = 3.2;
  double eTerm = pow(c,2)*exp(c/(pow(rad,2)-1));
  
  return(eTerm);

}

// First order derivation of the gauss spline.
double GaussWindowFunc::dW(double& rad) {

  double c = 3.2;
  double eTerm = pow(c,2)*exp(c/(pow(rad,2)-1));
  double pTerm = pow(rad,2)-1;

  return(-2.0*pow(c,3)/pow(pTerm,2)*rad*eTerm);
}

// Second order derivation of the gauss spline.
double GaussWindowFunc::d2W(double& rad) {

  double c = 3.2;  
  double eTerm = pow(c,2)*exp(c/(pow(rad,2)-1));
  double pTerm = pow(rad,2)-1;

  return(2.0*pow(c,3)*eTerm/pow(pTerm,2)*
	 (4.0/pTerm*pow(rad,2)-1.0+2.0*c/pow(pTerm,2)*pow(rad,2)));
}

/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void GaussWindowFunc2::calcWinFunc(double& x,double& y,double& z,
				   double& c,Particle& ptcle,
				   double& wFunc) {

  using namespace std;

  double k = 1.0;

  dbVector& radius = ptcle.getRadii();

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= radius[0] && yabs <= radius[1] && zabs <= radius[2])

    wFunc = W(x,radius[0],c,k)*W(y,radius[1],c,k)*W(z,radius[2],c,k);

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void GaussWindowFunc2::calcWinFunc1stDerivs(double& x,double& y,
					    double& z,double& c,
					    Particle& ptcle,
					    double& xDWFunc,
					    double& yDWFunc,
					    double& zDWFunc) {

  using namespace std;

  double result;

  double k = 1.0;
  dbVector& radius = ptcle.getRadii();

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= radius[0] && yabs <= radius[1] && zabs <= radius[2]) {

    xDWFunc = dW(x,radius[0],c,k)*W(y,radius[1],c,k)*W(z,radius[2],c,k);
    yDWFunc = W(x,radius[0],c,k)*dW(y,radius[1],c,k)*W(z,radius[2],c,k);
    zDWFunc = W(x,radius[0],c,k)*W(y,radius[1],c,k)*dW(z,radius[2],c,k);
  }

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void GaussWindowFunc2::calcWinFunc2ndDerivs(double& x,double& y,
					    double& z,double& c,
					    Particle& ptcle,
					    double& xxDWFunc,
					    double& yyDWFunc,
					    double& zzDWFunc,
					    double& xyDWFunc,
					    double& yzDWFunc,
					    double& zxDWFunc) {
  
  using namespace std;

  double k = 1.0;
  dbVector& radius = ptcle.getRadii();

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs <= radius[0] && yabs <= radius[1] && zabs <= radius[2]) {

    xxDWFunc = d2W(x,radius[0],c,k)*W(y,radius[1],c,k)
      *W(z,radius[2],c,k);
    yyDWFunc = W(x,radius[0],c,k)*d2W(y,radius[1],c,k)
      *W(z,radius[2],c,k);
    zzDWFunc = W(x,radius[0],c,k)*W(y,radius[1],c,k)
      *d2W(z,radius[2],c,k);
    
    xyDWFunc = dW(x,radius[0],c,k)*dW(y,radius[1],c,k)
      *W(z,radius[2],c,k);
    yzDWFunc = W(x,radius[0],c,k)*dW(y,radius[1],c,k)
      *dW(z,radius[2],c,k);
    zxDWFunc = dW(x,radius[0],c,k)*W(y,radius[1],c,k)
      *dW(z,radius[2],c,k);
  }

  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
// The gauss spline.
double GaussWindowFunc2::W(double& r,double& radius,double& c,double& k) {

  // [ e-((x-xI)/c)^2k ] [ 1 - e-(r/c)^2k ]^(-1) -
  // [ e-(r/c)^2k ] [ 1 - e-(r/c)^2k ]^(-1)

  double pTerm1 = pow(r/c,2.0*k);
  double pTerm2 = pow(radius/c,2.0*k);

  double eTerm1 = exp(-pTerm1);
  double eTerm2 = exp(-pTerm2);

  return((eTerm1-eTerm2)/(1.0-eTerm2));

}

// First order derivation of the gauss spline.
double GaussWindowFunc2::dW(double& r,double& radius,double& c,double& k) {

  // [ e-((x-xI)/c)^2k (-1)((x-xI)/c)^(2k-1) 2k/c ] [ 1 - e-(r/c)^2k ]^(-1)

  double pTerm1 = pow(r/c,2.0*k);
  double pTerm2 = pow(radius/c,2.0*k);

  double pTerm3 = -2.0*k/c*pow(r/c,2.0*k-1.0);

  double eTerm1 = exp(-pTerm1);
  double eTerm2 = exp(-pTerm2);

  return(eTerm1*pTerm3/(1.0-eTerm2));
}

// Second order derivation of the gauss spline.
double GaussWindowFunc2::d2W(double& r,double& radius,double& c,double& k) {
  
  // [ e-((x-xI)/c)^2k ((x-xI)/c)^(4k-2) (2k/c)^2 + 
  //   e-((x-xI)/c)^2k (-1)((x-xI)/c)^(2k-2) 2k/c^2 (2k-1) ] 
  // [ 1 - e-(r/c)^2k ]^(-1)

  double pTerm1 = pow(r/c,2.0*k);
  double pTerm2 = pow(radius/c,2.0*k);

  double pTerm3 = pow(2.0*k/c,2.0)*pow(r/c,4.0*k-2.0);
  double pTerm4 = (-1.0)*pow(r/c,2.0*k-2.0)*2.0*k/pow(c,2.0)*(2.0*k-1.0);

  
  double eTerm1 = exp(-pTerm1);
  double eTerm2 = exp(-pTerm2);
  
  
  return((eTerm1*pTerm3 + eTerm1*pTerm4)/(1.0-eTerm2));
}

/************************************************************************/
/************************************************************************/
// Calculate the window function in 3-D.
void TenthOrderSplineWinFunc::calcWinFunc(double& x,double& y,double& z,
					  double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1)

    wFunc = W(x)*W(y)*W(z);
 
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void TenthOrderSplineWinFunc::calcWinFunc1stDerivs(double& x,double& y,
						   double& z,double& xDWFunc,
						   double& yDWFunc,
						   double& zDWFunc) {

  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xDWFunc = dW(x)*W(y)*W(z);
    yDWFunc = W(x)*dW(y)*W(z);
    zDWFunc = W(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

// Calculate the second order derivations of the window function in 3-D.
void TenthOrderSplineWinFunc::calcWinFunc2ndDerivs(double& x,double& y,
						   double& z,
						   double& xxDWFunc,
						   double& yyDWFunc,
						   double& zzDWFunc,
						   double& xyDWFunc,
						   double& yzDWFunc,
						   double& zxDWFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xxDWFunc = d2W(x)*W(y)*W(z);
    yyDWFunc = W(x)*d2W(y)*W(z);
    zzDWFunc = W(x)*W(y)*d2W(z);
    
    xyDWFunc = dW(x)*dW(y)*W(z);
    yzDWFunc = W(x)*dW(y)*dW(z);
    zxDWFunc = dW(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}


/***********************************************************************/
// Tenth order spline.
double TenthOrderSplineWinFunc::W(double& rad) {
  
  double a = 0.98;
  double term = (a - a*pow(rad,2));

  return(pow(term,5));
}

// First order derivation of the tenth order spline.
double TenthOrderSplineWinFunc::dW(double& rad) {

  double a = 0.98;
  double term = (a - a*pow(rad,2));

  return(-10.0*a*rad*pow(term,4));
}

// Second order derivation of the tenth order spline.
double TenthOrderSplineWinFunc::d2W(double& rad) {

  double a = 0.98;
  double term = (a - a*pow(rad,2));

  return (80.0*pow(a*rad,2)*pow(term,3) - 10.0*a*pow(term,4));
}

/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void QuarticSplineWindowFunc2::calcWinFunc(double& x,double& y,double& z,
					  double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1)

    wFunc = W(x)*W(y)*W(z);
 
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void QuarticSplineWindowFunc2::calcWinFunc1stDerivs(double& x,double& y,
						   double& z,double& xDWFunc,
						   double& yDWFunc,
						   double& zDWFunc) {

  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xDWFunc = dW(x)*W(y)*W(z);
    yDWFunc = W(x)*dW(y)*W(z);
    zDWFunc = W(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void QuarticSplineWindowFunc2::calcWinFunc2ndDerivs(double& x,double& y,
						   double& z,
						   double& xxDWFunc,
						   double& yyDWFunc,
						   double& zzDWFunc,
						   double& xyDWFunc,
						   double& yzDWFunc,
						   double& zxDWFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xxDWFunc = d2W(x)*W(y)*W(z);
    yyDWFunc = W(x)*d2W(y)*W(z);
    zzDWFunc = W(x)*W(y)*d2W(z);
    
    xyDWFunc = dW(x)*dW(y)*W(z);
    yzDWFunc = W(x)*dW(y)*dW(z);
    zxDWFunc = dW(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
// The quartic spline.
double QuarticSplineWindowFunc2::W(double& rad) {

  using namespace std;

  dbVector a(5);
  dbVector b(5);
  dbVector c(5);
  dbVector d(5);
  dbVector e(5);
  dbVector f(5);

  a[0] = 27.0/32.0;
  a[1] = 27.0/8.0;
  a[2] = 81.0/16.0;
  a[3] = 27.0/8.0;
  a[4] = 27.0/32.0;

  b[0] = 17.0/96.0;
  b[1] = -5.0/8.0;
  b[2] = -63.0/16.0;
  b[3] = -45.0/8.0;
  b[4] = -81.0/32.0;

  c[0] = 11.0/48.0;
  c[1] = 0.0;
  c[2] = -9.0/8.0;
  c[3] = 0.0;
  c[4] = 27.0/16.0;

  d[0] = 11.0/48.0;
  d[1] = 0.0;
  d[2] = -9.0/8.0;
  d[3] = 0.0;
  d[4] = 27.0/16.0;

  e[0] = 17.0/96.0;
  e[1] = 5.0/8.0;
  e[2] = -63.0/16.0;
  e[3] = 45.0/8.0;
  e[4] = -81.0/32.0;

  f[0] = 27.0/32.0;
  f[1] = -27.0/8.0;
  f[2] = 81.0/16.0;
  f[3] = -27.0/8.0;
  f[4] = 27.0/32.0;


  double value;
  
  if(rad >= -1.0 && rad < -2.0/3.0)

    value =  a[0] + a[1]*rad + a[2]*pow(rad,2) + a[3]*pow(rad,3) + a[4]*pow(rad,4);

  else if(rad >= -2.0/3.0 && rad < -1.0/3.0)

    value =  b[0] + b[1]*rad + b[2]*pow(rad,2) + b[3]*pow(rad,3) + b[4]*pow(rad,4);

  else if(rad >= -1.0/3.0 && rad < 0.0)

    value =  c[0] + c[1]*rad + c[2]*pow(rad,2) + c[3]*pow(rad,3) + c[4]*pow(rad,4);

  else if(rad >= 0.0 && rad < 1.0/3.0)

    value =  d[0] + d[1]*rad + d[2]*pow(rad,2) + d[3]*pow(rad,3) + d[4]*pow(rad,4);

  else if(rad >= 1.0/3.0 && rad < 2.0/3.0)

    value =  e[0] + e[1]*rad + e[2]*pow(rad,2) + e[3]*pow(rad,3) + e[4]*pow(rad,4);

  else if(rad >= 2.0/3.0 && rad <= 1.0)

    value =  f[0] + f[1]*rad + f[2]*pow(rad,2) + f[3]*pow(rad,3) + f[4]*pow(rad,4);


  return(value);
}

// First order derivation of the quartic spline.
double QuarticSplineWindowFunc2::dW(double& rad) {

  using namespace std;

  dbVector a(5);
  dbVector b(5);
  dbVector c(5);
  dbVector d(5);
  dbVector e(5);
  dbVector f(5);

  a[0] = 27.0/32.0;
  a[1] = 27.0/8.0;
  a[2] = 81.0/16.0;
  a[3] = 27.0/8.0;
  a[4] = 27.0/32.0;

  b[0] = 17.0/96.0;
  b[1] = -5.0/8.0;
  b[2] = -63.0/16.0;
  b[3] = -45.0/8.0;
  b[4] = -81.0/32.0;

  c[0] = 11.0/48.0;
  c[1] = 0.0;
  c[2] = -9.0/8.0;
  c[3] = 0.0;
  c[4] = 27.0/16.0;

  d[0] = 11.0/48.0;
  d[1] = 0.0;
  d[2] = -9.0/8.0;
  d[3] = 0.0;
  d[4] = 27.0/16.0;

  e[0] = 17.0/96.0;
  e[1] = 5.0/8.0;
  e[2] = -63.0/16.0;
  e[3] = 45.0/8.0;
  e[4] = -81.0/32.0;

  f[0] = 27.0/32.0;
  f[1] = -27.0/8.0;
  f[2] = 81.0/16.0;
  f[3] = -27.0/8.0;
  f[4] = 27.0/32.0;

  double value;
  
  if(rad >= -1.0 && rad < -2.0/3.0)

    value =  a[1] + 2.0*a[2]*rad + 3.0*a[3]*pow(rad,2) + 4.0*a[4]*pow(rad,3);

  else if(rad >= -2.0/3.0 && rad < -1.0/3.0)

    value =  b[1] + 2.0*b[2]*rad + 3.0*b[3]*pow(rad,2) + 4.0*b[4]*pow(rad,3);

  else if(rad >= -1.0/3.0 && rad < 0.0)

    value =  c[1] + 2.0*c[2]*rad + 3.0*c[3]*pow(rad,2) + 4.0*c[4]*pow(rad,3);

  else if(rad >= 0.0 && rad < 1.0/3.0)

    value =  d[1] + 2.0*d[2]*rad + 3.0*d[3]*pow(rad,2) + 4.0*d[4]*pow(rad,3);

  else if(rad >= 1.0/3.0 && rad < 2.0/3.0)

    value =  e[1] + 2.0*e[2]*rad + 3.0*e[3]*pow(rad,2) + 4.0*e[4]*pow(rad,3);

  else if(rad >= 2.0/3.0 && rad <= 1.0)

    value =  f[1] + 2.0*f[2]*rad + 3.0*f[3]*pow(rad,2) + 4.0*f[4]*pow(rad,3);


  return(value);
}

// Second order derivation of the quartic spline.
double QuarticSplineWindowFunc2::d2W(double& rad) {

  using namespace std;

  dbVector a(5);
  dbVector b(5);
  dbVector c(5);
  dbVector d(5);
  dbVector e(5);
  dbVector f(5);

  a[0] = 27.0/32.0;
  a[1] = 27.0/8.0;
  a[2] = 81.0/16.0;
  a[3] = 27.0/8.0;
  a[4] = 27.0/32.0;

  b[0] = 17.0/96.0;
  b[1] = -5.0/8.0;
  b[2] = -63.0/16.0;
  b[3] = -45.0/8.0;
  b[4] = -81.0/32.0;

  c[0] = 11.0/48.0;
  c[1] = 0.0;
  c[2] = -9.0/8.0;
  c[3] = 0.0;
  c[4] = 27.0/16.0;

  d[0] = 11.0/48.0;
  d[1] = 0.0;
  d[2] = -9.0/8.0;
  d[3] = 0.0;
  d[4] = 27.0/16.0;

  e[0] = 17.0/96.0;
  e[1] = 5.0/8.0;
  e[2] = -63.0/16.0;
  e[3] = 45.0/8.0;
  e[4] = -81.0/32.0;

  f[0] = 27.0/32.0;
  f[1] = -27.0/8.0;
  f[2] = 81.0/16.0;
  f[3] = -27.0/8.0;
  f[4] = 27.0/32.0;


  double value;
  
  if(rad >= -1.0 && rad < -2.0/3.0)

    value =  2.0*a[2] + 6.0*a[3]*rad + 12.0*a[4]*pow(rad,2);

  else if(rad >= -2.0/3.0 && rad < -1.0/3.0)

    value =  2.0*b[2] + 6.0*b[3]*rad + 12.0*b[4]*pow(rad,2);

  else if(rad >= -1.0/3.0 && rad < 0.0)

    value =  2.0*c[2] + 6.0*c[3]*rad + 12.0*c[4]*pow(rad,2);

  else if(rad >= 0.0 && rad < 1.0/3.0)

    value =  2.0*d[2] + 6.0*d[3]*rad + 12.0*d[4]*pow(rad,2);

  else if(rad >= 1.0/3.0 && rad < 2.0/3.0)

    value =  2.0*e[2] + 6.0*e[3]*rad + 12.0*e[4]*pow(rad,2);

  else if(rad >= 2.0/3.0 && rad <= 1.0)

    value =  2.0*f[2] + 6.0*f[3]*rad + 12.0*f[4]*pow(rad,2);


  return(value);
}

/***********************************************************************/
/***********************************************************************/
// Calculate the window function in 3-D.
void QuinticSplineWindowFunc::calcWinFunc(double& x,double& y,double& z,
					  double& wFunc) {

  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1)

    wFunc = W(x)*W(y)*W(z);
 
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/***********************************************************************/
// Calculate the first order derivations the window function in 3-D.
void QuinticSplineWindowFunc::calcWinFunc1stDerivs(double& x,double& y,
						   double& z,double& xDWFunc,
						   double& yDWFunc,
						   double& zDWFunc) {
  
  using namespace std;

  double result;

  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);

  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xDWFunc = dW(x)*W(y)*W(z);
    yDWFunc = W(x)*dW(y)*W(z);
    zDWFunc = W(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
}

/************************************************************************/
// Calculate the second order derivations of the window function in 3-D.
void QuinticSplineWindowFunc::calcWinFunc2ndDerivs(double& x,double& y,
						   double& z,
						   double& xxDWFunc,
						   double& yyDWFunc,
						   double& zzDWFunc,
						   double& xyDWFunc,
						   double& yzDWFunc,
						   double& zxDWFunc) {
  
  using namespace std;

  double result;
  double xabs = fabs(x);
  double yabs = fabs(y);
  double zabs = fabs(z);


  if(xabs <= 1 && yabs <= 1 && zabs <= 1) {
    xxDWFunc = d2W(x)*W(y)*W(z);
    yyDWFunc = W(x)*d2W(y)*W(z);
    zzDWFunc = W(x)*W(y)*d2W(z);
    
    xyDWFunc = dW(x)*dW(y)*W(z);
    yzDWFunc = W(x)*dW(y)*dW(z);
    zxDWFunc = dW(x)*W(y)*dW(z);
  }
  else {
    cerr<<"While calculating window functions chosen point is out of " 
	<<"particle range!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
// The quartic spline.
double QuinticSplineWindowFunc::W(double& rad) {

  using namespace std;

  dbVector a(6);
  dbVector b(6);
  dbVector c(6);
  dbVector d(6);
  dbVector e(6);

  a[0] = -0.25431315104168916719;
  a[1] = -7.3750813802085168902;
  a[2] = -26.957194010417296681;
  a[3] = -39.164225260417538266;
  a[4] = -25.68562825520888282;
  a[5] = -6.3578287760417984131;

  b[0] = 0.63557942708333836634;
  b[1] = 0.040690104166682949938;
  b[2] = -2.2379557291666523078;
  b[3] = 2.0345052083335293247;
  b[4] = 8.6466471354170071351;
  b[5] = 5.0862630208334982385;

  c[0] = 0.63395182291666940522;
  c[1] = 1.043609643147647148e-14;
  c[2] = -2.6448567708333459159;
  c[3] = 6.1580370432542011902e-14;
  c[4] = 3.5603841145833388104;
  c[5] = -1.7000290064572679233e-13;

  d[0] = 0.63557942708333503568;
  d[1] = -0.04069010416665364005;
  d[2] = -2.2379557291667055985;
  d[3] = -2.0345052083331389703;
  d[4] = 8.6466471354163392249;
  d[5] = -5.0862630208331713888;

  e[0] = -0.25431315104164509133;
  e[1] = 7.3750813802081296444;
  e[2] = -26.957194010416003493;
  e[3] = 39.164225260415726382;
  e[4] = -25.685628255207724635;
  e[5] = 6.3578287760415177488;

  double value;
  
  if(rad >= -1.0 && rad < -0.6)

    value =  a[0] + a[1]*rad + a[2]*pow(rad,2) + a[3]*pow(rad,3) 
      + a[4]*pow(rad,4) + a[5]*pow(rad,5);

  else if(rad >= -0.6 && rad < -0.2)

    value =  b[0] + b[1]*rad + b[2]*pow(rad,2) + b[3]*pow(rad,3) 
      + b[4]*pow(rad,4) + b[5]*pow(rad,5);

  else if(rad >= -0.2 && rad < 0.2)

    value =  c[0] + c[1]*rad + c[2]*pow(rad,2) + c[3]*pow(rad,3) 
      + c[4]*pow(rad,4) + c[5]*pow(rad,5);

  else if(rad >= 0.2 && rad < 0.6)

    value =  d[0] + d[1]*rad + d[2]*pow(rad,2) + d[3]*pow(rad,3) 
      + d[4]*pow(rad,4) + d[5]*pow(rad,5);

  else

    value =  e[0] + e[1]*rad + e[2]*pow(rad,2) + e[3]*pow(rad,3) 
      + e[4]*pow(rad,4) + e[5]*pow(rad,5);


  return(value);
}

// First order derivation of the quartic spline.
double QuinticSplineWindowFunc::dW(double& rad) {

  using namespace std;

  dbVector a(6);
  dbVector b(6);
  dbVector c(6);
  dbVector d(6);
  dbVector e(6);


  a[0] = -0.25431315104168916719;
  a[1] = -7.3750813802085168902;
  a[2] = -26.957194010417296681;
  a[3] = -39.164225260417538266;
  a[4] = -25.68562825520888282;
  a[5] = -6.3578287760417984131;

  b[0] = 0.63557942708333836634;
  b[1] = 0.040690104166682949938;
  b[2] = -2.2379557291666523078;
  b[3] = 2.0345052083335293247;
  b[4] = 8.6466471354170071351;
  b[5] = 5.0862630208334982385;

  c[0] = 0.63395182291666940522;
  c[1] = 1.043609643147647148e-14;
  c[2] = -2.6448567708333459159;
  c[3] = 6.1580370432542011902e-14;
  c[4] = 3.5603841145833388104;
  c[5] = -1.7000290064572679233e-13;

  d[0] = 0.63557942708333503568;
  d[1] = -0.04069010416665364005;
  d[2] = -2.2379557291667055985;
  d[3] = -2.0345052083331389703;
  d[4] = 8.6466471354163392249;
  d[5] = -5.0862630208331713888;

  e[0] = -0.25431315104164509133;
  e[1] = 7.3750813802081296444;
  e[2] = -26.957194010416003493;
  e[3] = 39.164225260415726382;
  e[4] = -25.685628255207724635;
  e[5] = 6.3578287760415177488;

  double value;
  
  if(rad >= -1.0 && rad < -0.6)

    value =  a[1] + 2.0*a[2]*rad + 3.0*a[3]*pow(rad,2) + 4.0*a[4]*pow(rad,3) + 5.0*a[5]*pow(rad,4);

  else if(rad >= -0.6 && rad < -0.2)

    value =  b[1] + 2.0*b[2]*rad + 3.0*b[3]*pow(rad,2) + 4.0*b[4]*pow(rad,3) + 5.0*b[5]*pow(rad,4);

  else if(rad >= -0.2 && rad < 0.2)

    value =  c[1] + 2.0*c[2]*rad + 3.0*c[3]*pow(rad,2) + 4.0*c[4]*pow(rad,3) + 5.0*c[5]*pow(rad,4);

  else if(rad >= 0.2 && rad < 0.6)

    value =  d[1] + 2.0*d[2]*rad + 3.0*d[3]*pow(rad,2) + 4.0*d[4]*pow(rad,3) + 5.0*d[5]*pow(rad,4);

  else

    value =  e[1] + 2.0*e[2]*rad + 3.0*e[3]*pow(rad,2) + 4.0*e[4]*pow(rad,3) + 5.0*e[5]*pow(rad,4);


  return(value);

}

// Second order derivation of the quartic spline.
double QuinticSplineWindowFunc::d2W(double& rad) {

  using namespace std;

  dbVector a(6);
  dbVector b(6);
  dbVector c(6);
  dbVector d(6);
  dbVector e(6);

  a[0] = -0.25431315104168916719;
  a[1] = -7.3750813802085168902;
  a[2] = -26.957194010417296681;
  a[3] = -39.164225260417538266;
  a[4] = -25.68562825520888282;
  a[5] = -6.3578287760417984131;

  b[0] = 0.63557942708333836634;
  b[1] = 0.040690104166682949938;
  b[2] = -2.2379557291666523078;
  b[3] = 2.0345052083335293247;
  b[4] = 8.6466471354170071351;
  b[5] = 5.0862630208334982385;

  c[0] = 0.63395182291666940522;
  c[1] = 1.043609643147647148e-14;
  c[2] = -2.6448567708333459159;
  c[3] = 6.1580370432542011902e-14;
  c[4] = 3.5603841145833388104;
  c[5] = -1.7000290064572679233e-13;

  d[0] = 0.63557942708333503568;
  d[1] = -0.04069010416665364005;
  d[2] = -2.2379557291667055985;
  d[3] = -2.0345052083331389703;
  d[4] = 8.6466471354163392249;
  d[5] = -5.0862630208331713888;

  e[0] = -0.25431315104164509133;
  e[1] = 7.3750813802081296444;
  e[2] = -26.957194010416003493;
  e[3] = 39.164225260415726382;
  e[4] = -25.685628255207724635;
  e[5] = 6.3578287760415177488;

  double value;
  
  if(rad >= -1.0 && rad < -0.6)

    value =  2.0*a[2] + 6.0*a[3]*rad + 12.0*a[4]*pow(rad,2) + 20.0*a[5]*pow(rad,3);

  else if(rad >= -0.6 && rad < -0.2)

    value =  2.0*b[2] + 6.0*b[3]*rad + 12.0*b[4]*pow(rad,2) + 20.0*b[5]*pow(rad,3);

  else if(rad >= -0.2 && rad < 0.2)

    value =  2.0*c[2] + 6.0*c[3]*rad + 12.0*c[4]*pow(rad,2) + 20.0*c[5]*pow(rad,3);

  else if(rad >= 0.2 && rad < 0.6)

    value =  2.0*d[2] + 6.0*d[3]*rad + 12.0*d[4]*pow(rad,2) + 20.0*d[5]*pow(rad,3);

  else

    value =  2.0*e[2] + 6.0*e[3]*rad + 12.0*e[4]*pow(rad,2) + 20.0*e[5]*pow(rad,3);


  return(value);
}
