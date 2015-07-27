#include "debugFunctions.h"

/***********************************************************************/
/************************************************************************/
// Test routine for kontinuum mechanics' stuff.
void testContinuumMechanics(std::ofstream& logFile) {

  using namespace std;

  logFile<<"#####################################################"<<endl;
  logFile<<"######### test kontinuum mechanics's stuff ##########"<<endl;

  dbVector u(3);
  dbVector v(3);

  for(int i=0;i<3;i++) {
    u[i] = i+1.2;
    v[i] = i+1.5;
  }

  double scalar = 0;

  // vector scalar product.
  for(int i=0;i<3;i++)
    scalar += u[i]*v[i];

  logFile<<"vector scalar product: "<<scalar<<" bzw. ";

  scalarProduct(u,v,scalar,logFile);
  logFile<<scalar<<endl;

  // permutation
  intMatrix permutations = getPermutations(3);

  logFile<<"permutations"<<endl;
  for(int i=0;i<3;i++)
    logFile<<" +1 fr epsilon["<<permutations[i][0]<<"]["
	   <<permutations[i][1]<<"]["<<permutations[i][2]<<"]"<<endl;

  for(int i=3;i<6;i++)
    logFile<<" -1 fr epsilon["<<permutations[i][0]<<"]["
	   <<permutations[i][1]<<"]["<<permutations[i][2]<<"]"<<endl;

  // vector cross product.
  dbVector w(3);

  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];

  logFile<<"vector cross product:"<<endl;
  for(int i=0;i<3;i++)
    logFile<<"w["<<i<<"] = "<<w[i]<<endl;

  scalar = 1;
  dbVector x;
  crossProduct(u,v,x,scalar,logFile);
  for(int i=0;i<3;i++)
    logFile<<"w["<<i<<"] = "<<x[i]<<endl;

  // inner Tensor product
  dbMatrix T(3,dbVector(3));

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      T[i][j] = i/0.1 + j/0.2;

  clearArray(w);

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      w[i] += T[i][j]*u[j];

  logFile<<"inner Tensor product w = Tu:"<<endl;
  for(int i=0;i<3;i++)
    logFile<<"w["<<i<<"] = "<<w[i]<<endl;

  dbVector a;
  innerTensorProduct(T,u,a,false,logFile);
  for(int i=0;i<3;i++)
    logFile<<"w["<<i<<"] = "<<a[i]<<endl;

  clearArray(w);

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      w[i] += T[j][i]*u[j];

  logFile<<"inner Tensor product w = T_t u:"<<endl;
  for(int i=0;i<3;i++)
    logFile<<"w["<<i<<"] = "<<w[i]<<endl;

  innerTensorProduct(T,u,a,true,logFile);
  for(int i=0;i<3;i++)
    logFile<<"w["<<i<<"] = "<<a[i]<<endl;

  dbMatrix S(3,dbVector(3));

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      S[i][j] = i/0.1 + j/0.3;

  dbMatrix R(3,dbVector(3));

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
	R[i][j] += T[i][k]*S[k][j];

  logFile<<"inner Tensor product R = TS:"<<endl;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<R[i][j]<<endl;

  dbMatrix W;
  innerTensorProduct(T,S,W,false,false,logFile);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<W[i][j]<<endl;

  clearArray(R);

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
	R[i][j] += T[i][k]*S[j][k];

  logFile<<"inner Tensor product R = TS_t:"<<endl;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<R[i][j]<<endl;

  innerTensorProduct(T,S,W,false,true,logFile);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<W[i][j]<<endl;

  clearArray(R);

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
	R[i][j] += T[k][i]*S[j][k];

  logFile<<"inner Tensor product R = T_t S_t:"<<endl;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<R[i][j]<<endl;

  innerTensorProduct(T,S,W,true,true,logFile);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<W[i][j]<<endl;

  clearArray(R);

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
	R[i][j] += T[k][i]*S[k][j];

  logFile<<"inner Tensor product R = T_t S_t:"<<endl;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<R[i][j]<<endl;

  innerTensorProduct(T,S,W,true,false,logFile);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      logFile<<"R["<<i<<"]["<<j<<"] = "<<W[i][j]<<endl;

}

/***********************************************************************/
/***********************************************************************/
// Plot a cubic spline and its first and second order derivations.
void testSpline(std::ofstream& logFile) {

  using namespace std;

  int splineType;
  double xPlotCoord,yPlotCoord,zPlotCoord;

  ifstream inputFile("spline-test-input.dat");

  // Read input file
  if(inputFile) {

    string dummyArray;

    inputFile >> dummyArray >> splineType;
    inputFile >> dummyArray >> xPlotCoord;
    inputFile >> dummyArray >> yPlotCoord;
    inputFile >> dummyArray >> zPlotCoord;

  }
  else {
    cout<<"enter desired window function type (1 = cubic spline, \n"
	<<"2 = quartic spline, 3 = exponential, 4 = gauss spline,\n"
	<<"5 = 10th order spline): "<<endl;
    cin >> splineType;
    cout<<"enter x-plot coordinate [-1.0;1.0] (default = -100): ";
    cin >> xPlotCoord;
    cout<<endl;
    cout<<"enter y-plot coordinate [-1.0;1.0] (default = -100): ";
    cin >> yPlotCoord;
    cout<<endl;
    cout<<"enter z-plot coordinate [-1.0;1.0] (default = -100): ";
    cin >> zPlotCoord;
    cout<<endl;
  }

  /*********************************************************************/
  dbVector x(1001);
  dbVector y(1001);
  dbVector z(1001);

  dbVector windowFuncs(1001);
  dbVector xDerivWinFuncs(1001);
  dbVector yDerivWinFuncs(1001);
  dbVector zDerivWinFuncs(1001);

  dbVector xxDerivWinFuncs(1001);
  dbVector yyDerivWinFuncs(1001);
  dbVector zzDerivWinFuncs(1001);

  dbVector xyDerivWinFuncs(1001);
  dbVector yzDerivWinFuncs(1001);
  dbVector zxDerivWinFuncs(1001);

  dbVector plotOrd(1001);

  // Loop over all plot points.
  for(int i=0;i<=1000;i++) {

    if(xPlotCoord == -100)
      x[i] = (-500.0 + i)/500.0;

    else
      x[i] = xPlotCoord;

    if(yPlotCoord == -100)
      y[i] = (-500.0 + i)/500.0;

    else
      y[i] = yPlotCoord;

    if(zPlotCoord == -100)
      z[i] = (-500.0 + i)/500.0;

    else
      z[i] = zPlotCoord;

    if(i <= 500)
      plotOrd[i] = -sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
    else
      plotOrd[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);

    if(splineType == 1) {

      using namespace CubicSplineWindowFunc;

      calcWinFunc(x[i],y[i],z[i],windowFuncs[i]);
      //windowFuncs[i] /= pow(500.0,3);

      calcWinFunc1stDerivs(x[i],y[i],z[i],xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);
      //xDerivWinFuncs[i] /= pow(500.0,4);
      //yDerivWinFuncs[i] /= pow(500.0,4);
      //zDerivWinFuncs[i] /= pow(500.0,4);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);

      //xxDerivWinFuncs[i] /= pow(500.0,5);
      //yyDerivWinFuncs[i] /= pow(500.0,5);
      //zzDerivWinFuncs[i] /= pow(500.0,5);
      //xyDerivWinFuncs[i] /= pow(500.0,5);
      //yzDerivWinFuncs[i] /= pow(500.0,5);
      //zxDerivWinFuncs[i] /= pow(500.0,5);
    }
    else if (splineType == 2) {

      using namespace QuarticSplineWindowFunc;

      calcWinFunc(x[i],y[i],z[i],windowFuncs[i]);

      calcWinFunc1stDerivs(x[i],y[i],z[i],xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);

    }
    else if(splineType == 3) {

      using namespace ExponentialWindowFunc;

      calcWinFunc(x[i],y[i],z[i],windowFuncs[i]);

      calcWinFunc1stDerivs(x[i],y[i],z[i],xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);
    }
    else if(splineType == 4) {

      using namespace GaussWindowFunc;

      calcWinFunc(x[i],y[i],z[i],windowFuncs[i]);

      calcWinFunc1stDerivs(x[i],y[i],z[i],xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);
    }

    else if(splineType == 5) {

      using namespace TenthOrderSplineWinFunc;

      calcWinFunc(x[i],y[i],z[i],windowFuncs[i]);

      calcWinFunc1stDerivs(x[i],y[i],z[i],xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);
    }

    else if(splineType == 6) {

      using namespace QuarticSplineWindowFunc2;

      calcWinFunc(x[i],y[i],z[i],windowFuncs[i]);

      calcWinFunc1stDerivs(x[i],y[i],z[i],xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);
    }

    else if(splineType == 7) {

      using namespace GaussWindowFunc2;

      double value = 1.0;
      double maxDist = 1.0;

      Particle ptcle(3);
      ptcle.setID(0);
      double& weight = ptcle.getWeight();
      weight = 1.0;
      ptcle.setCoords(0,0,0);
      ptcle.setRadii(value,value,value);

      calcWinFunc(x[i],y[i],z[i],maxDist,ptcle,windowFuncs[i]);

      calcWinFunc1stDerivs(x[i],y[i],z[i],maxDist,ptcle,
			   xDerivWinFuncs[i],yDerivWinFuncs[i],
			   zDerivWinFuncs[i]);

      calcWinFunc2ndDerivs(x[i],y[i],z[i],maxDist,ptcle,
			   xxDerivWinFuncs[i],
			   yyDerivWinFuncs[i],zzDerivWinFuncs[i],
			   xyDerivWinFuncs[i],yzDerivWinFuncs[i],
			   zxDerivWinFuncs[i]);
    }


    else
      MPI_Abort(MPI_COMM_WORLD,1);

  }


  ofstream shp,d0shp,d1shp,d2shp,d00shp,
    d11shp,d22shp,d01shp,d12shp,d20shp;
  shp.open("spline.gra");
  d0shp.open("d0spline.gra");
  d1shp.open("d1spline.gra");
  d2shp.open("d2spline.gra");

  d00shp.open("d00spline.gra");
  d11shp.open("d11spline.gra");
  d22shp.open("d22spline.gra");

  d01shp.open("d01spline.gra");
  d12shp.open("d12spline.gra");
  d20shp.open("d20spline.gra");
    


  /* Write the header.*/
  shp<<"# Graf: \"Spline function\"\n"
     <<"#\n"
     <<"# X: \"0-direction-coordinate\" Y: \"splinefunction-value\""
     <<endl;
  d0shp<<"# Graf: \"X derivation of the spline function\"\n"
       <<"#\n"
       <<"# X: \"XYZ-direction-coordinate\" Y: "
       <<"\"X-derivation-splinefunction-value\""
       <<endl;
  d1shp<<"# Graf: \"Y derivation of the spline function\"\n"
       <<"#\n"
       <<"# X: \"XYZ-direction-coordinate\" Y: "
       <<"\"Y-derivation-splinefunction-value\""
       <<endl;
  d2shp<<"# Graf: \"Z derivation of the spline function\"\n"
       <<"#\n"
       <<"# X: \"XYZ-direction-coordinate\" Y: "
       <<"\"Z-derivation-splinefunction-value\""
       <<endl;

  d00shp<<"# Graf: \"XX derivation of the spline function\"\n"
	<<"#\n"
        <<"# X: \"XYZ-direction-coordinate\" Y: "
	<<"\"XX-derivation-splinefunction-value\""
	<<endl;
  d11shp<<"# Graf: \"YY derivation of the spline function\"\n"
	<<"#\n"
        <<"# X: \"XYZ-direction-coordinate\" Y: "
	<<"\"YY-derivation-splinefunction-value\""
	<<endl;
  d22shp<<"# Graf: \"ZZ derivation of the spline function\"\n"
	<<"#\n"
        <<"# X: \"XYZ-direction-coordinate\" Y: "
	<<"\"ZZ-derivationsplinefunction-value\""
	<<endl;

  d01shp<<"# Graf: \"XY derivation of the spline function\"\n"
	<<"#\n"
        <<"# X: \"XYZ-direction-coordinate\" Y: "
	<<"\"XY-derivatedSplinefunction-value\""
	<<endl;
  d12shp<<"# Graf: \"YZ derivation of the spline function\"\n"
	<<"#\n"
        <<"# X: \"XYZ-direction-coordinate\" Y: "
	<<"\"YZ-derivatedSplinefunction-value\""
	<<endl;
  d20shp<<"# Graf: \"ZX derivation of the spline function\"\n"
	<<"#\n"
        <<"# X: \"XYZ-direction-coordinate\" Y: "
	<<"\"ZX-derivatedSplinefunction-value\""
	<<endl;

  // Loop over all graph data and plot their ordinates.
  for(int i=0;i<=1000;i++) {
    shp<<plotOrd[i]<<" "<<windowFuncs[i]<<endl;

    d0shp<<plotOrd[i]<<" "<<xDerivWinFuncs[i]<<endl;
    d1shp<<plotOrd[i]<<" "<<yDerivWinFuncs[i]<<endl;
    d2shp<<plotOrd[i]<<" "<<zDerivWinFuncs[i]<<endl;

    d00shp<<plotOrd[i]<<" "<<xxDerivWinFuncs[i]<<endl;
    d11shp<<plotOrd[i]<<" "<<yyDerivWinFuncs[i]<<endl;
    d22shp<<plotOrd[i]<<" "<<zzDerivWinFuncs[i]<<endl;

    d01shp<<plotOrd[i]<<" "<<xyDerivWinFuncs[i]<<endl;
    d12shp<<plotOrd[i]<<" "<<yzDerivWinFuncs[i]<<endl;
    d20shp<<plotOrd[i]<<" "<<zxDerivWinFuncs[i]<<endl;
  }

  // MPI_Abort(MPI_COMM_WORLD,1);
}

/***********************************************************************/
/***********************************************************************/
// Plot the RKPM shape functions and its first and second order 
// derivations.
void testShapes(InputFileData* InputData,std::vector<Particle>& ptcls,
                std::map<std::string,double>& modelData,
                std::ofstream& logFile,PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDims = (int)modelData["usedDimensions"];

  bool symmetricRadii;
  int numOfPltPoints,pltPtcle,startPtcle,endPtcle,pltPoint;
  double radius;
  double xPlotCoord,yPlotCoord,zPlotCoord;

  ifstream inputFile("shape-test-input.dat");


  int radiusAlgorithm =
    (int)InputData->getValue("radiusDeterminationAlgorithm");
  int windowFuncType = (int)InputData->getValue("windowFunctionType");
  int polynomialType = (int)InputData->getValue("basisPolynomType");
  int polynomialOrder = (int)InputData->getValue("basisPolynomOrder");
  int shapefunctionType = (int)InputData->getValue("shapefunctionType");

  logFile<<"######################################################"<<endl;
  logFile<<"############### test shapefunction ###################"<<endl;

  // Read input file
  if(inputFile) {

    string dummyArray;

    // particle for shapefunction plotting is the one closest to the specified
    // plot point
    inputFile >> dummyArray >> pltPoint;
    inputFile >> dummyArray >> pltPtcle;

    // plot path
    inputFile >> dummyArray >> numOfPltPoints;
    inputFile >> dummyArray >> startPtcle;
    inputFile >> dummyArray >> endPtcle;

    // check symmetrize asymmetric radii
    inputFile >> dummyArray >> symmetricRadii;


    logFile<<"pltPoint = "<<pltPoint<<endl;
    logFile<<"pltPtcle = "<<pltPtcle<<endl;
    logFile<<"numOfPltPoints = "<<numOfPltPoints<<endl;
    logFile<<"startPtcle = "<<startPtcle<<endl;
    logFile<<"endPtcle = "<<endPtcle<<endl;
    logFile<<"symmetricRadii = "<<symmetricRadii<<endl;

    logFile<<"radiusAlgorithm = "<<radiusAlgorithm<<endl;
    logFile<<"windowFunctionType = "<<windowFuncType<<endl;
    logFile<<"basisPolynomType = "<<polynomialType<<endl;
    logFile<<"basisPolynomOrder = "<<polynomialOrder<<endl;
    logFile<<"shapefunctionType = "<<shapefunctionType<<endl;
  }

  else {
    cerr <<"Can't open input file 'shape-test-input.dat'!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  /*********************************************************************/
  // set a particle vector
  int numOfPtcls = ptcls.size();

  dbVector maxCoords(3);
  dbVector minCoords(3);


  logFile<<"******************************************************"<<endl;

  for(int i=0;i<ptcls.size();i++) {
    double& weight = ptcls[i].getWeight();
    //weight = 1.0;
    dbVector& coords = ptcls[i].getCoords();
    dbVector& radii = ptcls[i].getRadii();

    logFile<<"PARTICLE "<<i<<": "<<coords[0]<<" "<<coords[1]<<" "
	   <<coords[2]<<endl;

    for(int j=0;j<3;j++) {

      if(coords[j] > maxCoords[j])

	maxCoords[j] = coords[j];

      if(i == 0)

	minCoords[j] = coords[j];

      else if(coords[j] < minCoords[j])

	minCoords[j] = coords[j];

    }

  }

  for(int j=0;j<3;j++) {

    logFile<<"maxCoords["<<j<<"] = "<<maxCoords[j]<<endl;
    logFile<<"minCoords["<<j<<"] = "<<minCoords[j]<<endl;

  }

  /*********************************************************************/
  // Set plot points.

  if(startPtcle == 0 && endPtcle == 0) {
    startPtcle = 0;
    endPtcle = numOfPtcls-1;
  }
  else {
    startPtcle -= 1;
    endPtcle -= 1;
  }


  double xStart = ptcls[startPtcle].getCoord(0);
  double yStart = ptcls[startPtcle].getCoord(1);
  double zStart = ptcls[startPtcle].getCoord(2);
  double xEnd = ptcls[endPtcle].getCoord(0);
  double yEnd = ptcls[endPtcle].getCoord(1);
  double zEnd = ptcls[endPtcle].getCoord(2);

  dbVector pltOrds(numOfPltPoints);
  dbMatrix pltCoords(numOfPltPoints,dbVector(3));

  intMatrix suppPtcls(numOfPltPoints);

  double deltaX = xEnd-xStart;
  double deltaY = yEnd-yStart;
  double deltaZ = zEnd-zStart;
  double deltaPltPath =
    sqrt(pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ,2))/(numOfPltPoints-1);

  double plotStartValue = -deltaPltPath*numOfPltPoints/2.0;


  deltaX /= (numOfPltPoints-1);
  deltaY /= (numOfPltPoints-1);
  deltaZ /= (numOfPltPoints-1);

  logFile<<"******************************************************"<<endl;
  logFile<<"determine plot points"<<endl;
  logFile<<"startPtle = "<<startPtcle<<endl;
  logFile<<"endPtcle = "<<endPtcle<<endl;
  logFile<<"xStart = "<<xStart<<" -> xEnd = "<<xEnd
	 <<"; deltaX = "<<deltaX<<endl;
  logFile<<"yStart = "<<yStart<<" -> yEnd = "<<yEnd
	 <<"; deltaY = "<<deltaY<<endl;
  logFile<<"zStart = "<<zStart<<" -> zEnd = "<<zEnd
	 <<"; deltaZ = "<<deltaZ<<endl;
  logFile<<"numOfPltPoints = "<<numOfPltPoints<<"; deltaPltPath = "
	 <<deltaPltPath<<endl;

  for(int i=0;i<numOfPltPoints;i++) {

    pltOrds[i] = plotStartValue + i*deltaPltPath;

    pltCoords[i][0] = xStart + i*deltaX;
    pltCoords[i][1] = yStart + i*deltaY;
    pltCoords[i][2] = zStart + i*deltaZ;


    logFile<<"Plot-POINT "<<i<<": xCoord = "<<pltCoords[i][0]<<" yCoord "
	   <<pltCoords[i][1]<<" zCoord = "<<pltCoords[i][2]<<" ord = "
	   <<pltOrds[i]<<endl;
  }

  /*********************************************************************/
  // set supporting particle list for all plot points

  bool determinePltPtcle = false;
  double maxDistance = 1.0e+07;


  if(pltPtcle < 0)
    determinePltPtcle =  true;

  else
    pltPtcle -= 1;

  for(int i=0;i<numOfPltPoints;i++)

    for(int j=0;j<numOfPtcls;j++) {

      // Check if current point 'i' is supported by particle 'j'.
      if(ptcls[j].querySupported(InputData,pltCoords[i],modelData,
				 logFile)) {

	suppPtcls[i].push_back(j);

	if(determinePltPtcle && pltPoint == i)

	  if(calcDistance(ptcls[j].getCoords(),pltCoords[i]) < maxDistance) {

	    maxDistance = calcDistance(ptcls[j].getCoords(),pltCoords[i]);
	    pltPtcle = j;
	  }

      }

    }


  //#ifdef _commonDebugMode_
  logFile<<"*******************************************************"<<endl;
  logFile<<"pltPtcle = "<<pltPtcle<<endl;
  logFile<<"-------------------------------------------------------"<<endl;
  for(int i=0;i<suppPtcls.size();i++) {
    logFile<<"plot point "<<i<<"("<<pltCoords[i][0]<<", "
	   <<pltCoords[i][1]<<", "<<pltCoords[i][2]<<"): ";
    for(int j=0;j<suppPtcls[i].size();j++)
      logFile<<suppPtcls[i][j]<<" ";
    logFile<<endl;
  }
  //#endif


  /*********************************************************************/
  // set custom spline and its integral

  if(radiusAlgorithm == 6 || shapefunctionType == 5) {

    //WindowFuncAsym WFuncSet;

    for(int i=0;i<ptcls.size();i++) {

#ifdef _commonDebugMode_
      logFile<<"ptcle "<<i<<":"<<endl;
#endif

      if(symmetricRadii) {
	dbVector& radii = ptcls[i].getRadii();

	for(int k=0;k<usedDims;k++) {

	  if(radii[k] > radii[k+usedDims])
	    radii[k+usedDims] = radii[k];

	  else if(radii[k] < radii[k+usedDims])
	    radii[k] = radii[k+usedDims];

	}

      }

    }

  }

  /**********************************************************************/
  // plot spline

  ofstream spline,dxSpline,dySpline,dzSpline,dxxSpline,dyySpline,
    dzzSpline,dxySpline,dyzSpline,dzxSpline;
  spline.open("spline.grf");
  dxSpline.open("dxSpline.grf");
  dySpline.open("dySpline.grf");
  dzSpline.open("dzSpline.grf");
  dxxSpline.open("dxxSpline.grf");
  dyySpline.open("dyySpline.grf");
  dzzSpline.open("dzzSpline.grf");
  dxySpline.open("dxySpline.grf");
  dyzSpline.open("dyzSpline.grf");
  dzxSpline.open("dzxSpline.grf");

  // Write the header.
  spline<<"# Graf: \"Spline function\"\n"
	<<"#\n"
        <<"# X: \"plot-coordinate\" Y: \"splinefunction-value\""
	<<endl;

  dxSpline<<"# Graf: \"Spline function x-deriv\"\n"
	  <<"#\n"
          <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	  <<endl;

  dySpline<<"# Graf: \"Spline function y-deriv\"\n"
	  <<"#\n"
          <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	  <<endl;

  dzSpline<<"# Graf: \"Spline function z-deriv\"\n"
	  <<"#\n"
          <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	  <<endl;

  dxxSpline<<"# Graf: \"Spline function xx-deriv\"\n"
	   <<"#\n"
           <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	   <<endl;

  dyySpline<<"# Graf: \"Spline function yy-deriv\"\n"
	   <<"#\n"
           <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	   <<endl;

  dzzSpline<<"# Graf: \"Spline function zz-deriv\"\n"
	   <<"#\n"
           <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	   <<endl;

  dxySpline<<"# Graf: \"Spline function xy-deriv\"\n"
	   <<"#\n"
           <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	   <<endl;

  dyzSpline<<"# Graf: \"Spline function yz-deriv\"\n"
	   <<"#\n"
           <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	   <<endl;

  dzxSpline<<"# Graf: \"Spline function zx-deriv\"\n"
	   <<"#\n"
           <<"# X: \"plot-coordinate\" Y: \"splinefunctionderiv-value\""
	   <<endl;

  int supportSize;

  dbVector W(numOfPltPoints);
  dbVector dxW(numOfPltPoints);
  dbVector dyW(numOfPltPoints);
  dbVector dzW(numOfPltPoints);

  dbVector dxxW(numOfPltPoints);
  dbVector dyyW(numOfPltPoints);
  dbVector dzzW(numOfPltPoints);

  dbVector dxyW(numOfPltPoints);
  dbVector dyzW(numOfPltPoints);
  dbVector dzxW(numOfPltPoints);

  dbVector xnorm(3);

  dbVector& radii = ptcls[pltPtcle].getRadii();


  int firstSuppPltPoint = -1;
  int lastSuppPltPoint = -1;

  logFile<<"######################################################"<<endl;
  logFile<<"**************** spline plotting *********************"<<endl;
  for(int i=0;i<radii.size();i++)
    logFile<<radii[i]<<endl;

  for(int i=0;i<numOfPltPoints;i++) {
    supportSize = suppPtcls[i].size();

    int pos = findIntVecPos(pltPtcle,0,suppPtcls[i].size(),suppPtcls[i]);

    if(shapefunctionType != 5 && pos != -1) {

      xnorm[0] = (pltCoords[i][0] - ptcls[pltPtcle].getCoord(0))/radii[0];
      xnorm[1] = (pltCoords[i][1] - ptcls[pltPtcle].getCoord(1))/radii[1];
      xnorm[2] = (pltCoords[i][2] - ptcls[pltPtcle].getCoord(2))/radii[2];

      if(windowFuncType == 1) {

	using namespace CubicSplineWindowFunc;

	calcWinFunc(xnorm[0],xnorm[1],xnorm[2],W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxW[i],dyW[i],
			     dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxxW[i],
			     dyyW[i],dzzW[i],
			     dxyW[i],dyzW[i],
			     dzxW[i]);

      }
      else if (windowFuncType == 2) {

	using namespace QuarticSplineWindowFunc;

	calcWinFunc(xnorm[0],xnorm[1],xnorm[2],W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxW[i],
			     dyW[i],dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxxW[i],
			     dyyW[i],dzzW[i],
			     dxyW[i],dyzW[i],
			     dzxW[i]);

      }
      else if(windowFuncType == 3) {

	using namespace ExponentialWindowFunc;

	calcWinFunc(xnorm[0],xnorm[1],xnorm[2],W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxW[i],
			     dyW[i],dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxxW[i],
			     dyyW[i],dzzW[i],dxyW[i],dyzW[i],
			     dzxW[i]);
      }
      else if(windowFuncType == 4) {

	using namespace GaussWindowFunc;

	calcWinFunc(xnorm[0],xnorm[1],xnorm[2],W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxW[i],
			     dyW[i],dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxxW[i],
			     dyyW[i],dzzW[i],
			     dxyW[i],dyzW[i],
			     dzxW[i]);
      }

      else if(windowFuncType == 5) {

	using namespace TenthOrderSplineWinFunc;

	calcWinFunc(xnorm[0],xnorm[1],xnorm[2],W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxW[i],
			     dyW[i],dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxxW[i],
			     dyyW[i],dzzW[i],
			     dxyW[i],dyzW[i],
			     dzxW[i]);
      }

      else if(windowFuncType == 6) {

	using namespace QuarticSplineWindowFunc2;

	calcWinFunc(xnorm[0],xnorm[1],xnorm[2],W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxW[i],
			     dyW[i],dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],dxxW[i],
			     dyyW[i],dzzW[i],
			     dxyW[i],dyzW[i],dzxW[i]);
      }

      else if(windowFuncType == 7) {

	using namespace GaussWindowFunc2;

	double value = 1.0;
	double maxDist = 1.0;

	Particle ptcle(3);
	ptcle.setID(0);
	double& weight = ptcle.getWeight();
	weight = 1.0;
	ptcle.setCoords(0,0,0);
	ptcle.setRadii(value,value,value);

	calcWinFunc(xnorm[0],xnorm[1],
		    xnorm[2],maxDist,ptcle,W[i]);

	calcWinFunc1stDerivs(xnorm[0],xnorm[1],
			     xnorm[2],maxDist,ptcle,
			     dxW[i],dyW[i],dzW[i]);

	calcWinFunc2ndDerivs(xnorm[0],xnorm[1],
			     xnorm[2],maxDist,ptcle,
			     dxxW[i],dyyW[i],dzzW[i],
			     dxyW[i],dyzW[i],dzxW[i]);
      }


      else
	MPI_Abort(MPI_COMM_WORLD,1);

    }

    else if(shapefunctionType == 5 && pos != -1) {

      xnorm[0] = pltCoords[i][0] - ptcls[pltPtcle].getCoord(0);
      xnorm[1] = pltCoords[i][1] - ptcls[pltPtcle].getCoord(1);
      xnorm[2] = pltCoords[i][2] - ptcls[pltPtcle].getCoord(2);


      WindowFuncAsym WFuncSet;

      WFuncSet.calcSplineValue(InputData,ptcls[pltPtcle],xnorm,W[i],
			       dxW[i],dyW[i],dzW[i],
			       dxxW[i],dyyW[i],dzzW[i],
			       dxyW[i],dyzW[i],dzxW[i],
			       modelData,logFile,viewerSEQ);


    }


    logFile<<"PLT-POINT "<<i<<" xnorm[0]="<<xnorm[0]
	   <<" xnorm[1]="<<xnorm[1]<<" xnorm[2]="<<xnorm[2]
	   <<" s="<<pltOrds[i]<<" W="<<W[i]<<endl;
  }


  // Loop over all graph data and plot their ordinates.
  for(int i=0;i<numOfPltPoints;i++) {
    spline<<pltOrds[i]<<" "<<W[i]<<endl;

    dxSpline<<pltOrds[i]<<" "<<dxW[i]<<endl;
    dySpline<<pltOrds[i]<<" "<<dyW[i]<<endl;
    dzSpline<<pltOrds[i]<<" "<<dzW[i]<<endl;

    dxxSpline<<pltOrds[i]<<" "<<dxxW[i]<<endl;
    dyySpline<<pltOrds[i]<<" "<<dyyW[i]<<endl;
    dzzSpline<<pltOrds[i]<<" "<<dzzW[i]<<endl;

    dxySpline<<pltOrds[i]<<" "<<dxyW[i]<<endl;
    dyzSpline<<pltOrds[i]<<" "<<dyzW[i]<<endl;
    dzxSpline<<pltOrds[i]<<" "<<dzxW[i]<<endl;
  }

  /*********************************************************************/
  // Set shape functions.
  ShapefunctionSet* ShapeSet;


  dbMatrix shapes(numOfPltPoints);

  dbMatrix xDerivShapes(numOfPltPoints);
  dbMatrix yDerivShapes(numOfPltPoints);
  dbMatrix zDerivShapes(numOfPltPoints);

  dbMatrix xxDerivShapes(numOfPltPoints);
  dbMatrix yyDerivShapes(numOfPltPoints);
  dbMatrix zzDerivShapes(numOfPltPoints);

  dbMatrix xyDerivShapes(numOfPltPoints);
  dbMatrix yzDerivShapes(numOfPltPoints);
  dbMatrix zxDerivShapes(numOfPltPoints);

  bool secondOrderDerivs = false;

  // calculation of the shapefunctions including their first order
  // derivation.
  if(!secondOrderDerivs)

    for(int i=0;i<numOfPltPoints;i++) {
      supportSize = suppPtcls[i].size();

      if(shapes[i].size() < supportSize)
	shapes[i].resize(supportSize);

      if(shapefunctionType != 1) {

	if(xDerivShapes[i].size() < supportSize) {

	  xDerivShapes[i].resize(supportSize);
	  yDerivShapes[i].resize(supportSize);
	  zDerivShapes[i].resize(supportSize);
	}

      }


      if(shapefunctionType == 1) {

	ShapeSet = new ShepardShapeFunc(InputData,supportSize,suppPtcls[i],
					ptcls,pltCoords[i][0],pltCoords[i][1],
					pltCoords[i][2],1,modelData,
					logFile,viewerSEQ);

	shapes[i] = ShapeSet->getShapefunctions();
	dbMatrix& firstDerivShapes = ShapeSet->getFirstDerivShapes();

	xDerivShapes[i] = firstDerivShapes[0];
	yDerivShapes[i] = firstDerivShapes[1];
	zDerivShapes[i] = firstDerivShapes[2];
      }
      else if(shapefunctionType == 2)

	RKPMShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				  ptcls,pltCoords[i][0],pltCoords[i][1],
				  pltCoords[i][2],shapes[i],
				  xDerivShapes[i],yDerivShapes[i],
				  zDerivShapes[i],modelData,logFile,
				  viewerSEQ);

      else if(shapefunctionType == 3)

	EFGShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				 ptcls,pltCoords[i][0],pltCoords[i][1],
				 pltCoords[i][2],shapes[i],
				 xDerivShapes[i],yDerivShapes[i],
				 zDerivShapes[i],modelData,logFile,
				 viewerSEQ);

      else if(shapefunctionType == 4)

	OrthoShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				   ptcls,pltCoords[i][0],pltCoords[i][1],
				   pltCoords[i][2],shapes[i],
				   xDerivShapes[i],yDerivShapes[i],
				   zDerivShapes[i],modelData,logFile,
				   viewerSEQ);

      else if(shapefunctionType == 5)

	AsymShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				  ptcls,pltCoords[i][0],pltCoords[i][1],
				  pltCoords[i][2],shapes[i],
				  xDerivShapes[i],yDerivShapes[i],
				  zDerivShapes[i],modelData,logFile,
				  viewerSEQ);
    }

  // calculation of the shapefunctions including their first and second
  // order derivation.
  else

    for(int i=0;i<numOfPltPoints;i++) {
      supportSize = suppPtcls[i].size();

      if(shapes[i].size() < supportSize)
	shapes[i].resize(supportSize);

      if(shapefunctionType != 1) {

	if(xDerivShapes[i].size() < supportSize) {

	  xDerivShapes[i].resize(supportSize);
	  yDerivShapes[i].resize(supportSize);
	  zDerivShapes[i].resize(supportSize);

	  xxDerivShapes[i].resize(supportSize);
	  yyDerivShapes[i].resize(supportSize);
	  zzDerivShapes[i].resize(supportSize);
	  xyDerivShapes[i].resize(supportSize);
	  yzDerivShapes[i].resize(supportSize);
	  zxDerivShapes[i].resize(supportSize);

	}

      }


      if(shapefunctionType == 1) {

	ShapeSet = new ShepardShapeFunc(InputData,supportSize,suppPtcls[i],
					ptcls,pltCoords[i][0],pltCoords[i][1],
					pltCoords[i][2],2,modelData,
					logFile,viewerSEQ);

	shapes[i] = ShapeSet->getShapefunctions();
	dbMatrix& firstDerivShapes = ShapeSet->getFirstDerivShapes();
	dbMatrix& secondDerivShapes = ShapeSet->getSecondDerivShapes();

	xDerivShapes[i] = firstDerivShapes[0];
	yDerivShapes[i] = firstDerivShapes[1];
	zDerivShapes[i] = firstDerivShapes[2];

	xxDerivShapes[i] = secondDerivShapes[0];
	yyDerivShapes[i] = secondDerivShapes[1];
	zzDerivShapes[i] = secondDerivShapes[2];
	xyDerivShapes[i] = secondDerivShapes[3];
	yzDerivShapes[i] = secondDerivShapes[4];
	zxDerivShapes[i] = secondDerivShapes[5];
      }
      else if(shapefunctionType == 2)

	RKPMShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				  ptcls,pltCoords[i][0],pltCoords[i][1],
				  pltCoords[i][2],shapes[i],
				  xDerivShapes[i],yDerivShapes[i],
				  zDerivShapes[i],xxDerivShapes[i],
				  yyDerivShapes[i],zzDerivShapes[i],
				  xyDerivShapes[i],yzDerivShapes[i],
				  zxDerivShapes[i],modelData,logFile,
				  viewerSEQ);

      else if(shapefunctionType == 3)

	EFGShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				 ptcls,pltCoords[i][0],pltCoords[i][1],
				 pltCoords[i][2],shapes[i],
				 xDerivShapes[i],yDerivShapes[i],
				 zDerivShapes[i],xxDerivShapes[i],
				 yyDerivShapes[i],zzDerivShapes[i],
				 xyDerivShapes[i],yzDerivShapes[i],
				 zxDerivShapes[i],modelData,logFile,
				 viewerSEQ);

      else if(shapefunctionType == 4)

	OrthoShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				   ptcls,pltCoords[i][0],pltCoords[i][1],
				   pltCoords[i][2],shapes[i],
				   xDerivShapes[i],yDerivShapes[i],
				   zDerivShapes[i],xxDerivShapes[i],
				   yyDerivShapes[i],zzDerivShapes[i],
				   xyDerivShapes[i],yzDerivShapes[i],
				   zxDerivShapes[i],modelData,logFile,
				   viewerSEQ);

      else if(shapefunctionType == 5)

	AsymShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				  ptcls,pltCoords[i][0],pltCoords[i][1],
				  pltCoords[i][2],shapes[i],
				  xDerivShapes[i],yDerivShapes[i],
				  zDerivShapes[i],xxDerivShapes[i],
				  yyDerivShapes[i],zzDerivShapes[i],
				  xyDerivShapes[i],yzDerivShapes[i],
				  zxDerivShapes[i],modelData,logFile,
				  viewerSEQ);
    }

  logFile<<"finished calculation shapefunctions"<<endl;


#ifdef _commonDebugMode_
  logFile<<"**************** all spline volumes ****************"<<endl;
  for(int i=0;i<ptcls.size();i++) {
    CustomSpline* spline = ptcls[i].getSpline();
    dbVector& radii = ptcls[i].getRadii();
    logFile<<"PTCLE "<<i<<" radii: ";
    for(int j=0;j<radii.size();j++)
      logFile<<radii[j]<<" ";
    logFile<<"volume="<<spline->getSplineIntegral()
	   <<" r^3="<<radii[0]*radii[1]*radii[2]<<endl;
  }
#endif

  /*********************************************************************/
  // Write graph files.


  ofstream support;
  support.open("support.gra");
    
  // Write the headers.
  support<<"# Graf: \"number of supporting particles\"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: \"support number\""
	 <<endl;


  // Loop over all plot points and plot the number of supporting
  // particles
  double ptcleSupport,totalWeight;

  double scaleSupportPlot = 1.0;

  for(int i=0;i<suppPtcls.size();i++) {

    ptcleSupport = suppPtcls[i].size()*scaleSupportPlot;

    // plot support
    support<<pltOrds[i]<<" "<<ptcleSupport<<endl;

  }

  // ====================================================================

  ofstream shp,d0shp,d1shp,d2shp,d00shp,d11shp,d22shp,d01shp,d12shp,d20shp;
  shp.open("shape.grf");
  d0shp.open("d0shape.grf");
  d1shp.open("d1shape.grf");
  d2shp.open("d2shape.grf");

  d00shp.open("d00shape.grf");
  d11shp.open("d11shape.grf");
  d22shp.open("d22shape.grf");

  d01shp.open("d01shape.grf");
  d12shp.open("d12shape.grf");
  d20shp.open("d20shape.grf");
    

  // Plot the shape function of a particular particle.
  if(pltPtcle != -1) {

    // Write the header.
    shp<<"# Graf: \"Shape function of particle "
       <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
       <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
       <<"#\n"
       <<"# X: \"plot-ordinate\" Y: \"shapefunction-value\""
       <<endl;
    d0shp<<"# Graf: \"X derivation of the shape function of particle "
	 <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	 <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"X-derivation-shapefunction-value\""
	 <<endl;
    d1shp<<"# Graf: \"Y derivation of the shape function of particle "
	 <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	 <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"Y-derivation-shapefunction-value\""
	 <<endl;
    d2shp<<"# Graf: \"Z derivation of the shape function of particle "
	 <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	 <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"Z-derivation-shapefunction-value\""
	 <<endl;

    d00shp<<"# Graf: \"XX derivation of the shape function of particle "
	  <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	  <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"XX-derivation-shapefunction-value\""
	  <<endl;
    d11shp<<"# Graf: \"YY derivation of the shape function of particle "
	  <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	  <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"YY-derivation-shapefunction-value\""
	  <<endl;
    d22shp<<"# Graf: \"ZZ derivation of the shape function of particle "
	  <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	  <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"ZZ-derivationshapefunction-value\""
	  <<endl;

    d01shp<<"# Graf: \"XY derivation of the shape function of particle "
	  <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	  <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"XY-derivatedShapefunction-value\""
	  <<endl;
    d12shp<<"# Graf: \"YZ derivation of the shape function of particle "
	  <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	  <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"YZ-derivatedShapefunction-value\""
	  <<endl;
    d20shp<<"# Graf: \"ZX derivation of the shape function of particle "
	  <<ptcls[pltPtcle].getCoord(0)<<" "<<ptcls[pltPtcle].getCoord(1)<<" "
	  <<ptcls[pltPtcle].getCoord(2)<<" \"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"ZX-derivatedShapefunction-value\""
	  <<endl;

    logFile<<"***************************************************"<<endl;
    logFile<<"**************** start plotting *******************"<<endl;

    // Loop over all plot points and plot the ordinates of a particluar shape
    // function.
    int pos;

    for(int i=0;i<numOfPltPoints;i++) {

      pos = findIntVecPos(pltPtcle,0,suppPtcls[i].size(),suppPtcls[i]);

      if(pos != -1) {

	//#ifdef _commonDebugMode_
	logFile<<"plot point "<<i<<" -> pos "<<pos<<"; suppPtcle "
	       <<suppPtcls[i][pos]<<" =? "<<pltPtcle<<endl;
	//#endif

	shp<<pltOrds[i]<<" "<<shapes[i][pos]<<endl;

	d0shp<<pltOrds[i]<<" "<<xDerivShapes[i][pos]<<endl;
	d1shp<<pltOrds[i]<<" "<<yDerivShapes[i][pos]<<endl;
	d2shp<<pltOrds[i]<<" "<<zDerivShapes[i][pos]<<endl;

	if(secondOrderDerivs) {
	  d00shp<<pltOrds[i]<<" "<<xxDerivShapes[i][pos]<<endl;
	  d11shp<<pltOrds[i]<<" "<<yyDerivShapes[i][pos]<<endl;
	  d22shp<<pltOrds[i]<<" "<<zzDerivShapes[i][pos]<<endl;

	  d01shp<<pltOrds[i]<<" "<<xyDerivShapes[i][pos]<<endl;
	  d12shp<<pltOrds[i]<<" "<<yzDerivShapes[i][pos]<<endl;
	  d20shp<<pltOrds[i]<<" "<<zxDerivShapes[i][pos]<<endl;
	}
	else {
	  d00shp<<pltOrds[i]<<" 0"<<endl;
	  d11shp<<pltOrds[i]<<" 0"<<endl;
	  d22shp<<pltOrds[i]<<" 0"<<endl;

	  d01shp<<pltOrds[i]<<" 0"<<endl;
	  d12shp<<pltOrds[i]<<" 0"<<endl;
	  d20shp<<pltOrds[i]<<" 0"<<endl;
	}

      }
      else {
	shp<<pltOrds[i]<<" 0"<<endl;

	d0shp<<pltOrds[i]<<" 0"<<endl;
	d1shp<<pltOrds[i]<<" 0"<<endl;
	d2shp<<pltOrds[i]<<" 0"<<endl;

	d00shp<<pltOrds[i]<<" 0"<<endl;
	d11shp<<pltOrds[i]<<" 0"<<endl;
	d22shp<<pltOrds[i]<<" 0"<<endl;

	d01shp<<pltOrds[i]<<" 0"<<endl;
	d12shp<<pltOrds[i]<<" 0"<<endl;
	d20shp<<pltOrds[i]<<" 0"<<endl;
      }

    }

  }
  //---------------------------------------------------------------------
  // Plot summation of all shape functions.
  else {

    // Write the header.
    shp<<"# Graf: \"summation of all shape functions\"\n"
       <<"#\n"
       <<"# X: \"plot-ordinate\" Y: \"shapefunction-value\""
       <<endl;
    d0shp<<"# Graf: \"summation of all X derivation of the shape functions\"\n"
	 <<"#\n"
	 <<"# X: \"plot-ordinate\" Y: "
	 <<"\"X-derivation-shapefunction-value\""
         <<endl;
    d1shp<<"# Graf: \"summation of all Y derivation of the shape functions\"\n"
	 <<"#\n"
	 <<"# X: \"plot-ordinate\" Y: "
	 <<"\"Y-derivation-shapefunction-value\""
         <<endl;
    d2shp<<"# Graf: \"summation of all Z derivation of the shape functions\"\n"
	 <<"#\n"
	 <<"# X: \"plot-ordinate\" Y: "
	 <<"\"Z-derivation-shapefunction-value\""
         <<endl;

    d00shp<<"# Graf: \"summation of all XX derivation of the shape functions\"\n"
	  <<"#\n"
	  <<"# X: \"plot-ordinate\" Y: "
	  <<"\"XX-derivation-shapefunction-value\""
          <<endl;
    d11shp<<"# Graf: \"summation of all YY derivation of the shape functions\"\n"
	  <<"#\n"
	  <<"# X: \"plot-ordinate\" Y: "
	  <<"\"YY-derivation-shapefunction-value\""
          <<endl;
    d22shp<<"# Graf: \"summation of all ZZ derivation of the shape functions\"\n"
	  <<"#\n"
	  <<"# X: \"plot-ordinate\" Y: "
	  <<"\"ZZ-derivationshapefunction-value\""
          <<endl;

    d01shp<<"# Graf: \"summation of all XY derivation of the shape functions\"\n"
	  <<"#\n"
	  <<"# X: \"plot-ordinate\" Y: "
	  <<"\"XY-derivatedShapefunction-value\""
          <<endl;
    d12shp<<"# Graf: \"summation of all YZ derivation of the shape functions\"\n"
	  <<"#\n"
	  <<"# X: \"plot-ordinate\" Y: "
	  <<"\"YZ-derivatedShapefunction-value\""
          <<endl;
    d20shp<<"# Graf: \"summation of all ZX derivation of the shape functions\"\n"
	  <<"#\n"
	  <<"# X: \"plot-ordinate\" Y: "
	  <<"\"ZX-derivatedShapefunction-value\""
          <<endl;

    // Loop over all plot points, sum the shape functions and plot them.
    double shape;
    double dxShape,dyShape,dzShape;
    double dxxShape,dyyShape,dzzShape,dxyShape,dyzShape,dzxShape;

    for(int i=0;i<numOfPltPoints;i++) {

      shape = dxShape = dyShape = dzShape = dxxShape = dyyShape =
	dzzShape = dxyShape = dyzShape = dzxShape = 0;

      for(int j=0;j<shapes[i].size();j++) {
	shape += shapes[i][j];
	dxShape += xDerivShapes[i][j];
	dyShape += yDerivShapes[i][j];
	dzShape += zDerivShapes[i][j];

	if(secondOrderDerivs) {
	  dxxShape += xxDerivShapes[i][j];
	  dyyShape += yyDerivShapes[i][j];
	  dzzShape += zzDerivShapes[i][j];
	  dxyShape += xyDerivShapes[i][j];
	  dyzShape += yzDerivShapes[i][j];
	  dzxShape += zxDerivShapes[i][j];
	}

      }

      shp<<pltOrds[i]<<" "<<shape<<endl;

      d0shp<<pltOrds[i]<<" "<<dxShape<<endl;
      d1shp<<pltOrds[i]<<" "<<dyShape<<endl;
      d2shp<<pltOrds[i]<<" "<<dzShape<<endl;

      d00shp<<pltOrds[i]<<" "<<dxxShape<<endl;
      d11shp<<pltOrds[i]<<" "<<dyyShape<<endl;
      d22shp<<pltOrds[i]<<" "<<dzzShape<<endl;

      d01shp<<pltOrds[i]<<" "<<dxyShape<<endl;
      d12shp<<pltOrds[i]<<" "<<dyzShape<<endl;
      d20shp<<pltOrds[i]<<" "<<dzxShape<<endl;
    }

  }

}

/************************************************************************/
/************************************************************************/
// Test the curve fitting property of a MLS shape functions and its 
// first and second order derivations within the domain and on its
// boundary.

void testMLS(InputFileData* InputData,std::vector<Particle>& ptcls,
             std::map<std::string,double>& modelData,
             std::ofstream& logFile,PetscViewer& viewerSEQ) {

  using namespace std;

  int numOfPtcls = ptcls.size();

  int usedDims = (int)modelData["usedDimensions"];
  int radiusAlgorithm =
    (int)InputData->getValue("radiusDeterminationAlgorithm");
  int windowFuncType = (int)InputData->getValue("windowFunctionType");
  int polynomialType = (int)InputData->getValue("basisPolynomType");
  int polynomialOrder = (int)InputData->getValue("basisPolynomOrder");
  int shapefunctionType = (int)InputData->getValue("shapefunctionType");

  logFile<<"******************************************************"<<endl;
  logFile<<"****************** mls test routine ******************"<<endl;
  logFile<<"DBL_EPSILON = "<<DBL_EPSILON<<endl;

  bool symmetricRadii;
  int numOfPltPoints,startPtcle,endPtcle,pltPtcle,pltPoint;
  bool absoluteResiduum;
  double minValue,maxPercentage;

  int shapeFuncPlotPosition,functionOrder,functionType;
  bool plotSupport;
  int gidPlotPoint;
  double scaleSupportPlot,scaleSupportWeightPlot,scaleSupportResiduumPlot;

  string dummyArray;
  ifstream inputFile("mls-test-input.dat");

  if(!inputFile) {
    logFile<<"Can't open input file mls-test-input.dat!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Read input file
  if(inputFile) {

    // plot path
    inputFile >> dummyArray >> numOfPltPoints;
    inputFile >> dummyArray >> startPtcle;
    inputFile >> dummyArray >> endPtcle;

    // test function
    inputFile >> dummyArray >> functionType;
    inputFile >> dummyArray >> functionOrder;

    // plot residuum

    // minimal ratio between the exact and the approximated function values
    inputFile >> dummyArray >> minValue;

    // maximal plotted percentage value for the relative residuum
    inputFile >> dummyArray >> maxPercentage;

    // plot the absolute values of the residuum
    inputFile >> dummyArray >> absoluteResiduum;

    // plot support
    inputFile >> dummyArray >> plotSupport;
    inputFile >> dummyArray >> scaleSupportPlot;
    inputFile >> dummyArray >> scaleSupportWeightPlot;

    // plot 3-D figure of the support of particle plot point
    inputFile >> dummyArray >> gidPlotPoint;

    // check symmetrize asymmetric radii
    inputFile >> dummyArray >> symmetricRadii;


    logFile<<"numOfPltPoints = "<<numOfPltPoints<<endl;
    logFile<<"startPtcle = "<<startPtcle<<endl;
    logFile<<"endPtcle = "<<endPtcle<<endl;
    logFile<<"radiusAlgorithm = "<<radiusAlgorithm<<endl;

    logFile<<"windowFunctionType = "<<windowFuncType<<endl;
    logFile<<"basisPolynomType = "<<polynomialType<<endl;
    logFile<<"basisPolynomOrder = "<<polynomialOrder<<endl;
    logFile<<"shapefunctionType = "<<shapefunctionType<<endl;

    logFile<<"functionType = "<<functionType<<endl;
    logFile<<"functionOrder = "<<functionOrder<<endl;

    logFile<<"minValue = "<<minValue<<endl;
    logFile<<"maxPercentage = "<<maxPercentage<<endl;
    logFile<<"absoluteResiduum = "<<absoluteResiduum<<endl;

    logFile<<"plotSupport = "<<plotSupport<<endl;
    logFile<<"scaleSupportPlot = "<<scaleSupportPlot<<endl;
    logFile<<"scaleSupportWeightPlot = "<<scaleSupportWeightPlot<<endl;
    logFile<<"gidPlotPoint = "<<gidPlotPoint<<endl; ;
    logFile<<"symmetricRadii = "<<symmetricRadii<<endl;
  }

  else {
    cerr <<"Can't open input file 'mls-test-input.dat'!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


  /*********************************************************************/
  // set a particle vector
  bool secondOrderDerivs = true;

  dbVector maxCoords(3);
  dbVector minCoords(3);


  logFile<<"******************************************************"<<endl;

  for(int i=0;i<numOfPtcls;i++) {
    double& weight = ptcls[i].getWeight();
    weight = 1.0;
    dbVector& coords = ptcls[i].getCoords();
    dbVector& radii = ptcls[i].getRadii();

    logFile<<"PARTICLE "<<i<<": "<<coords[0]<<" "<<coords[1]<<" "
	   <<coords[2]<<endl;

    for(int j=0;j<3;j++) {

      if(coords[j] > maxCoords[j])

	maxCoords[j] = coords[j];

      if(i == 0)

	minCoords[j] = coords[j];

      else if(coords[j] < minCoords[j])

	minCoords[j] = coords[j];

    }

  }

  for(int j=0;j<3;j++) {

    logFile<<"maxCoords["<<j<<"] = "<<maxCoords[j]<<endl;
    logFile<<"minCoords["<<j<<"] = "<<minCoords[j]<<endl;

  }

  /*********************************************************************/
  // Set plot points.

  if(startPtcle > ptcls.size() || endPtcle > ptcls.size()) {
    cerr <<"In testMLS startPtcle or endPtcle nonexisting!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  if(startPtcle == 0 && endPtcle == 0) {
    startPtcle = 0;
    endPtcle = numOfPtcls-1;
  }
  else {
    startPtcle -= 1;
    endPtcle -= 1;
  }
  pltPtcle -= 1;

  double xStart = ptcls[startPtcle].getCoord(0);
  double yStart = ptcls[startPtcle].getCoord(1);
  double zStart = ptcls[startPtcle].getCoord(2);
  double xEnd = ptcls[endPtcle].getCoord(0);
  double yEnd = ptcls[endPtcle].getCoord(1);
  double zEnd = ptcls[endPtcle].getCoord(2);

  dbMatrix pltOrds(numOfPltPoints,dbVector(11));
  dbMatrix pltCoords(numOfPltPoints,dbVector(3));

  intMatrix suppPtcls(numOfPltPoints);

  double deltaX = xEnd-xStart;
  double deltaY = yEnd-yStart;
  double deltaZ = zEnd-zStart;

  double deltaPltPath =
    sqrt(pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ,2))/(numOfPltPoints-1);

  double plotStartValue = -deltaPltPath/2.0*(numOfPltPoints-1);

  deltaX /= (numOfPltPoints-1);
  deltaY /= (numOfPltPoints-1);
  deltaZ /= (numOfPltPoints-1);

  logFile<<"******************************************************"<<endl;
  logFile<<"determine plot points"<<endl;
  logFile<<"startPtle = "<<startPtcle<<endl;
  logFile<<"endPtcle = "<<endPtcle<<endl;
  logFile<<"xStart = "<<xStart<<" -> xEnd = "<<xEnd
	 <<"; deltaX = "<<deltaX<<endl;
  logFile<<"yStart = "<<yStart<<" -> yEnd = "<<yEnd
	 <<"; deltaY = "<<deltaY<<endl;
  logFile<<"zStart = "<<zStart<<" -> zEnd = "<<zEnd
	 <<"; deltaZ = "<<deltaZ<<endl;
  logFile<<"numOfPltPoints = "<<numOfPltPoints<<"; deltaPltPath = "
	 <<deltaPltPath<<endl;

  for(int i=0;i<numOfPltPoints;i++) {

    pltOrds[i][0] = plotStartValue + i*deltaPltPath;

    pltCoords[i][0] = xStart + i*deltaX;
    pltCoords[i][1] = yStart + i*deltaY;
    pltCoords[i][2] = zStart + i*deltaZ;


    logFile<<"Plot-POINT "<<i<<": xCoord = "<<pltCoords[i][0]<<" yCoord "
	   <<pltCoords[i][1]<<" zCoord = "<<pltCoords[i][2]<<" ord = "
	   <<pltOrds[i][0]<<endl;
  }


  /*********************************************************************/
  // set supporting particle list for all plot points

  double maxDistance = 1.0e+07;

  for(int i=0;i<numOfPltPoints;i++)

    for(int j=0;j<numOfPtcls;j++) {

      // Check if current point 'i' is supported by particle 'j'.
      if(ptcls[j].querySupported(InputData,pltCoords[i],modelData,
				 logFile))

	suppPtcls[i].push_back(j);


    }


  //#ifdef _commonDebugMode_
  logFile<<"*******************************************************"<<endl;
  logFile<<"************* plot point support list *****************"<<endl;
  for(int i=0;i<suppPtcls.size();i++) {
    logFile<<"Plot-POINT "<<i<<"("<<pltCoords[i][0]<<", "
	   <<pltCoords[i][1]<<", "<<pltCoords[i][2]<<"): ";
    for(int j=0;j<suppPtcls[i].size();j++)
      logFile<<suppPtcls[i][j]<<" ";
    logFile<<endl;
  }
  //#endif

  /*********************************************************************/
  // set custom spline and its integral

  if(radiusAlgorithm == 6 || shapefunctionType == 5) {

    for(int i=0;i<ptcls.size();i++) {

#ifdef _commonDebugMode_
      logFile<<"ptcle "<<i<<":"<<endl;
#endif

      if(symmetricRadii) {
	dbVector& radii = ptcls[i].getRadii();

	for(int k=0;k<usedDims;k++) {

	  if(radii[k] > radii[k+usedDims])
	    radii[k+usedDims] = radii[k];

	  else if(radii[k] < radii[k+usedDims])
	    radii[k] = radii[k+usedDims];

	}

      }

    }

#ifdef _commonDebugMode_
    logFile<<"**************** all spline volumes ****************"<<endl;
    for(int i=0;i<ptcls.size();i++) {
      CustomSpline& spline = *ptcls[i].getSpline();
      dbVector& radii = ptcls[i].getRadii();
      logFile<<"PTCLE "<<i<<" radii: ";
      for(int j=0;j<radii.size();j++)
	logFile<<radii[j]<<" ";
      logFile<<"volume="<<spline.getSplineIntegral()<<endl;
    }
#endif

  }

  /*********************************************************************/
  // Plot the relation of residuum to support
  int supportSize;

  if(plotSupport) {

    ofstream support,supportWeight;
    support.open("support.grf");
    supportWeight.open("supportWeight.grf");

    // Write the headers.
    support<<"# Graf: \"number of supporting particles\"\n"
	   <<"#\n"
	   <<"# X: \"plot-ordinate\" Y: \"support number\""
	   <<endl;

    supportWeight<<"# Graf: \"total support weight\"\n"
		 <<"#\n"
		 <<"# X: \"plot-ordinate\" Y: \"support-weight\""
		 <<endl;

    // Loop over all plot points and plot the number of supporting
    // particles
    double ptcleSupport,totalWeight;

    WindowFunctionSet* WFuncSet;

    for(int i=0;i<numOfPltPoints;i++) {

      ptcleSupport = suppPtcls[i].size()*scaleSupportPlot;

      // plot support
      support<<pltOrds[i][0]<<" "<<ptcleSupport<<endl;

      // support weight in total
      supportSize = suppPtcls[i].size();

      if(shapefunctionType != 5)

	WFuncSet =
	  new PrismaticWindowFunctionSet(InputData,ptcls,suppPtcls[i],
				pltCoords[i][0],pltCoords[i][1],
				pltCoords[i][2],supportSize,0,
				modelData,logFile);

      else if(shapefunctionType == 5 && windowFuncType == 4)

	WFuncSet =
	  new WindowFuncAsym(InputData,ptcls,suppPtcls[i],
			     pltCoords[i][0],pltCoords[i][1],
			     pltCoords[i][2],supportSize,0,
			     modelData,logFile,viewerSEQ);

      else {
	cerr <<"Windowfunction-shapefunction combination not supported!"
	     << endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

      dbVector& W = WFuncSet->getWindowFuncs();
      totalWeight = 0;

      for(int j=0;j<supportSize;j++)
	totalWeight += W[j];

      supportWeight<<pltOrds[i][0]<<" "<<totalWeight*scaleSupportWeightPlot
		   <<endl;

      delete WFuncSet;
    }

  }


  cout<<"finished plotting of particle support"<<endl;

  /*********************************************************************/
  // plot the support graphically

  // write the meshfile
  ofstream meshFile;
  meshFile.open("support.msh");

  // Write the header.
  meshFile<<"MESH  dimension 3  ElemType Point Nnode 1"<<endl;
  meshFile<<"Coordinates"<<endl;

  // Loop oveDOF"<<i+1<<"-displsr all particles.
  for(int i=0;i<numOfPtcls;i++) {
    meshFile<<i+1<<" ";

    // Loop over all used degrees of freedom.
    for(int j=0;j<3;j++)
      meshFile<<ptcls[i].getCoord(j)<<" ";

    meshFile<<endl;

  }

  meshFile<<"end coordinates"<<endl;
  meshFile<<endl;
  meshFile<<"Elements"<<endl;

  for(int i=0;i<numOfPtcls;i++)
    meshFile<<i+1<<" "<<i+1<<"  1"<<endl;

  meshFile<<"end elements"<<endl;

  // Write the result file marking those particles which support a
  // given point.

  // Write a header.
  ofstream resultFile;
  resultFile.open("support.res");
  resultFile<<"GiD Post Results File 1.0"<<endl;
  resultFile<<"Result \" \" \" \" 1 Scalar  OnNodes"<<endl;
  resultFile<<"ComponentNames \"ptcle-support\""<<endl;
  resultFile<<"Values"<<endl;

  // Loop over all particles to mark them if they support the given
  // point.

  for(int i=0;i<numOfPtcls;i++) {


    resultFile<<i+1<<" ";

    for(int j=0;j<suppPtcls[gidPlotPoint].size();j++)

      if(suppPtcls[gidPlotPoint][j] == i)
	resultFile<<"0 "<<endl;

      else
	resultFile<<"1 "<<endl;

  }

  resultFile<<"End Values"<<endl;

  /*********************************************************************/
  // set the particle parameter
  double xCoord,yCoord,zCoord;

  double dl;
  double norming = (double)1.0/25.000;
  int linEQSize;
  dbMatrix ptcleParams(numOfPtcls,dbVector(10));
  dbVector normFactors(3);

  double value;

  dbVector basis;
  dbVector xDerivBasis;
  dbVector yDerivBasis;
  dbVector zDerivBasis;
  dbVector xxDerivBasis;
  dbVector yyDerivBasis;
  dbVector zzDerivBasis;
  dbVector xyDerivBasis;
  dbVector yzDerivBasis;
  dbVector zxDerivBasis;

#ifdef _commonDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"************ ptcle params computation ***************"<<endl;
#endif

  // Select choosen polynom type (Pascal,Lagrangian,serendipity).
  switch(functionType) {
    
    // own creation
  case 0:

    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {
      xCoord = ptcls[i].getCoord(0)*norming;
      yCoord = ptcls[i].getCoord(1)*norming;
      zCoord = ptcls[i].getCoord(2)*norming;

      ptcleParams[i][0] = 100.0*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);

      ptcleParams[i][1] = 100.0*norming*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][2]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][3]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::dW(zCoord);

      ptcleParams[i][4]  = 100.0*pow(norming,2)*GaussWindowFunc::d2W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][5] = 100.0*pow(norming,2)*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::d2W(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][6]  = 100.0*pow(norming,2)*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::d2W(zCoord);
      ptcleParams[i][7] = 100.0*pow(norming,2)*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][8] = 100.0*pow(norming,2)*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::dW(zCoord);
      ptcleParams[i][9] = 100.0*pow(norming,2)*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::dW(zCoord);
    }

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {
      xCoord = pltCoords[i][0]*norming;
      yCoord = pltCoords[i][1]*norming;
      zCoord = pltCoords[i][2]*norming;

      pltOrds[i][1] = 100.0*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);

      pltOrds[i][2] = 100.0*norming*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][3]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][4]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::dW(zCoord);

      pltOrds[i][5] = 100.0*pow(norming,2)*GaussWindowFunc::d2W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][6] = 100.0*pow(norming,2)*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::d2W(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][7] = 100.0*pow(norming,2)*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::d2W(zCoord);
      pltOrds[i][8] = 100.0*pow(norming,2)*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][9] = 100.0*pow(norming,2)*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::dW(zCoord);
      pltOrds[i][10] = 100.0*pow(norming,2)*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::dW(zCoord);

    }

    break;

    // Pascal
  case 1:

    // Select choosen polynom order.
    switch(functionOrder) {

    case 1:
      linEQSize = 4;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);
      xxDerivBasis.resize(linEQSize);
      yyDerivBasis.resize(linEQSize);
      zzDerivBasis.resize(linEQSize);
      xyDerivBasis.resize(linEQSize);
      yzDerivBasis.resize(linEQSize);
      zxDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  PascalLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalLinear::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];

	  ptcleParams[i][4] += xxDerivBasis[j];
	  ptcleParams[i][5] += yyDerivBasis[j];
	  ptcleParams[i][6] += zzDerivBasis[j];
	  ptcleParams[i][7] += xyDerivBasis[j];
	  ptcleParams[i][8] += yzDerivBasis[j];
	  ptcleParams[i][9] += zxDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  PascalLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalLinear::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	  pltOrds[i][5] += xxDerivBasis[j];
	  pltOrds[i][6] += yyDerivBasis[j];
	  pltOrds[i][7] += zzDerivBasis[j];
	  pltOrds[i][8] += xyDerivBasis[j];
	  pltOrds[i][9] += yzDerivBasis[j];
	  pltOrds[i][10] += zxDerivBasis[j];

	}

      }

      break;

    case 2:
      linEQSize = 10;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);
      xxDerivBasis.resize(linEQSize);
      yyDerivBasis.resize(linEQSize);
      zzDerivBasis.resize(linEQSize);
      xyDerivBasis.resize(linEQSize);
      yzDerivBasis.resize(linEQSize);
      zxDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  PascalQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalQuadratic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = PascalQuadratic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = PascalQuadratic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = PascalQuadratic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = PascalQuadratic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = PascalQuadratic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = PascalQuadratic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];

	  ptcleParams[i][4] += xxDerivBasis[j];
	  ptcleParams[i][5] += yyDerivBasis[j];
	  ptcleParams[i][6] += zzDerivBasis[j];
	  ptcleParams[i][7] += xyDerivBasis[j];
	  ptcleParams[i][8] += yzDerivBasis[j];
	  ptcleParams[i][9] += zxDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  PascalQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalQuadratic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = PascalQuadratic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = PascalQuadratic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = PascalQuadratic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = PascalQuadratic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = PascalQuadratic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = PascalQuadratic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	  pltOrds[i][5] += xxDerivBasis[j];
	  pltOrds[i][6] += yyDerivBasis[j];
	  pltOrds[i][7] += zzDerivBasis[j];
	  pltOrds[i][8] += xyDerivBasis[j];
	  pltOrds[i][9] += yzDerivBasis[j];
	  pltOrds[i][10] += zxDerivBasis[j];

	}

      }

      break;

    case 3:
      linEQSize = 20;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);
      xxDerivBasis.resize(linEQSize);
      yyDerivBasis.resize(linEQSize);
      zzDerivBasis.resize(linEQSize);
      xyDerivBasis.resize(linEQSize);
      yzDerivBasis.resize(linEQSize);
      zxDerivBasis.resize(linEQSize);

      //norming = (double)1.0/1000.0; // mls-cube
      //dl = 11.0; // mls-cube
      dl = 1.0;
      norming = (double)1.0/pow(pltOrds[numOfPltPoints-1][0],3);

      for(int i=0;i<numOfPtcls;i++) {

	xCoord = (ptcls[i].getCoord(0)+dl)*norming;
	yCoord = (ptcls[i].getCoord(1)+dl)*norming;
	zCoord = (ptcls[i].getCoord(2)+dl)*norming;

	basis =  PascalCubic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalCubic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalCubic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalCubic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = PascalCubic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = PascalCubic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = PascalCubic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = PascalCubic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = PascalCubic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = PascalCubic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += norming*xDerivBasis[j];
	  ptcleParams[i][2] += norming*yDerivBasis[j];
	  ptcleParams[i][3] += norming*zDerivBasis[j];

	  ptcleParams[i][4] += pow(norming,2)*xxDerivBasis[j];
	  ptcleParams[i][5] += pow(norming,2)*yyDerivBasis[j];
	  ptcleParams[i][6] += pow(norming,2)*zzDerivBasis[j];
	  ptcleParams[i][7] += pow(norming,2)*xyDerivBasis[j];
	  ptcleParams[i][8] += pow(norming,2)*yzDerivBasis[j];
	  ptcleParams[i][9] += pow(norming,2)*zxDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = (pltCoords[i][0]+dl)*norming;
	yCoord = (pltCoords[i][1]+dl)*norming;
	zCoord = (pltCoords[i][2]+dl)*norming;

	basis =  PascalCubic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalCubic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalCubic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalCubic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = PascalCubic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = PascalCubic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = PascalCubic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = PascalCubic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = PascalCubic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = PascalCubic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += norming*xDerivBasis[j];
	  pltOrds[i][3] += norming*yDerivBasis[j];
	  pltOrds[i][4] += norming*zDerivBasis[j];

	  pltOrds[i][5] += pow(norming,2)*xxDerivBasis[j];
	  pltOrds[i][6] += pow(norming,2)*yyDerivBasis[j];
	  pltOrds[i][7] += pow(norming,2)*zzDerivBasis[j];
	  pltOrds[i][8] += pow(norming,2)*xyDerivBasis[j];
	  pltOrds[i][9] += pow(norming,2)*yzDerivBasis[j];
	  pltOrds[i][10] += pow(norming,2)*zxDerivBasis[j];

	}

      }

      break;

    default:
      logFile<<"Chosen function polynom order isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Lagrangian
  case 2:

    switch(functionOrder) {

    case 1:
      linEQSize = 8;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);
      xxDerivBasis.resize(linEQSize);
      yyDerivBasis.resize(linEQSize);
      zzDerivBasis.resize(linEQSize);
      xyDerivBasis.resize(linEQSize);
      yzDerivBasis.resize(linEQSize);
      zxDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  LagrangeLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeLinear::dPz(xCoord,yCoord,zCoord);

	xyDerivBasis = LagrangeLinear::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = LagrangeLinear::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = LagrangeLinear::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];

	  ptcleParams[i][4] += xxDerivBasis[j];
	  ptcleParams[i][5] += yyDerivBasis[j];
	  ptcleParams[i][6] += zzDerivBasis[j];
	  ptcleParams[i][7] += xyDerivBasis[j];
	  ptcleParams[i][8] += yzDerivBasis[j];
	  ptcleParams[i][9] += zxDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  LagrangeLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeLinear::dPz(xCoord,yCoord,zCoord);

	xyDerivBasis = LagrangeLinear::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = LagrangeLinear::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = LagrangeLinear::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	  pltOrds[i][5] += xxDerivBasis[j];
	  pltOrds[i][6] += yyDerivBasis[j];
	  pltOrds[i][7] += zzDerivBasis[j];
	  pltOrds[i][8] += xyDerivBasis[j];
	  pltOrds[i][9] += yzDerivBasis[j];
	  pltOrds[i][10] += zxDerivBasis[j];

	}

      }

      break;

    case 2:
      linEQSize = 27;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);
      xxDerivBasis.resize(linEQSize);
      yyDerivBasis.resize(linEQSize);
      zzDerivBasis.resize(linEQSize);
      xyDerivBasis.resize(linEQSize);
      yzDerivBasis.resize(linEQSize);
      zxDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  LagrangeQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeQuadratic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = LagrangeQuadratic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = LagrangeQuadratic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = LagrangeQuadratic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = LagrangeQuadratic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = LagrangeQuadratic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = LagrangeQuadratic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];

	  ptcleParams[i][4] += xxDerivBasis[j];
	  ptcleParams[i][5] += yyDerivBasis[j];
	  ptcleParams[i][6] += zzDerivBasis[j];
	  ptcleParams[i][7] += xyDerivBasis[j];
	  ptcleParams[i][8] += yzDerivBasis[j];
	  ptcleParams[i][9] += zxDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  LagrangeQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeQuadratic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = LagrangeQuadratic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = LagrangeQuadratic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = LagrangeQuadratic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = LagrangeQuadratic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = LagrangeQuadratic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = LagrangeQuadratic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	  pltOrds[i][5] += xxDerivBasis[j];
	  pltOrds[i][6] += yyDerivBasis[j];
	  pltOrds[i][7] += zzDerivBasis[j];
	  pltOrds[i][8] += xyDerivBasis[j];
	  pltOrds[i][9] += yzDerivBasis[j];
	  pltOrds[i][10] += zxDerivBasis[j];

	}

      }

      break;

    default:
      logFile<<"Chosen function polynom order isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Serendipity
  case 3:

    switch(functionOrder) {

    case 2:
      linEQSize = 20;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);
      xxDerivBasis.resize(linEQSize);
      yyDerivBasis.resize(linEQSize);
      zzDerivBasis.resize(linEQSize);
      xyDerivBasis.resize(linEQSize);
      yzDerivBasis.resize(linEQSize);
      zxDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  SerendipityQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = SerendipityQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = SerendipityQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = SerendipityQuadratic::dPz(xCoord,yCoord,zCoord);

	xxDerivBasis = SerendipityQuadratic::dPxx(xCoord,yCoord,zCoord);
	yyDerivBasis = SerendipityQuadratic::dPyy(xCoord,yCoord,zCoord);
	zzDerivBasis = SerendipityQuadratic::dPzz(xCoord,yCoord,zCoord);
	xyDerivBasis = SerendipityQuadratic::dPxy(xCoord,yCoord,zCoord);
	yzDerivBasis = SerendipityQuadratic::dPyz(xCoord,yCoord,zCoord);
	zxDerivBasis = SerendipityQuadratic::dPzx(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];

	  ptcleParams[i][4] += xxDerivBasis[j];
	  ptcleParams[i][5] += yyDerivBasis[j];
	  ptcleParams[i][6] += zzDerivBasis[j];
	  ptcleParams[i][7] += xyDerivBasis[j];
	  ptcleParams[i][8] += yzDerivBasis[j];
	  ptcleParams[i][9] += zxDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis.resize(linEQSize);
	xDerivBasis.resize(linEQSize);
	yDerivBasis.resize(linEQSize);
	zDerivBasis.resize(linEQSize);
	xxDerivBasis.resize(linEQSize);
	yyDerivBasis.resize(linEQSize);
	zzDerivBasis.resize(linEQSize);
	xyDerivBasis.resize(linEQSize);
	yzDerivBasis.resize(linEQSize);
	zxDerivBasis.resize(linEQSize);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	  pltOrds[i][5] += xxDerivBasis[j];
	  pltOrds[i][6] += yyDerivBasis[j];
	  pltOrds[i][7] += zzDerivBasis[j];
	  pltOrds[i][8] += xyDerivBasis[j];
	  pltOrds[i][9] += yzDerivBasis[j];
	  pltOrds[i][10] += zxDerivBasis[j];

	}

      }

      break;

    default:
      logFile<<"Choosen function polynom order isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    // tenth order spline
  case 5:

    norming = (double)1.0/35.0;

    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {
      xCoord = ptcls[i].getCoord(0)*norming;
      yCoord = ptcls[i].getCoord(1)*norming;
      zCoord = ptcls[i].getCoord(2)*norming;

      ptcleParams[i][0] = 100*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);

      ptcleParams[i][1] = 100*norming*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][2]  = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][3]  = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);

      ptcleParams[i][4] = 100*pow(norming,2)*TenthOrderSplineWinFunc::d2W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][5] = 100*pow(norming,2)*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::d2W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][6] = 100*pow(norming,2)*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::d2W(zCoord);
      ptcleParams[i][7] = 100*pow(norming,2)*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][8] = 100*pow(norming,2)*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);
      ptcleParams[i][9] = 100*pow(norming,2)*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);

      logFile<<"ptcleParams["<<i<<"]: "<<ptcleParams[i][0]<<endl;
    }

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {
      xCoord = pltCoords[i][0]*norming;
      yCoord = pltCoords[i][1]*norming;
      zCoord = pltCoords[i][2]*norming;

      pltOrds[i][1] = 100*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);

      pltOrds[i][2] = 100*norming*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][3] = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][4] = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);

      pltOrds[i][5] = 100*pow(norming,2)*TenthOrderSplineWinFunc::d2W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][6] = 100*pow(norming,2)*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::d2W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][7] = 100*pow(norming,2)*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::d2W(zCoord);
      pltOrds[i][8] = 100*pow(norming,2)*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][9] = 100*pow(norming,2)*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);
      pltOrds[i][10]= 100*pow(norming,2)*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);

    }

    break;

    // fourth order spline
  case 6:

    normFactors[0] = fabs(minCoords[0]-maxCoords[0])/2.0;
    normFactors[1] = fabs(minCoords[1]-maxCoords[1])/2.0;
    normFactors[2] = fabs(minCoords[2]-maxCoords[2])/2.0;


    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {

      xCoord = modf(fabs(ptcls[i].getCoord(0)/normFactors[0]),&value);

      if(xCoord < 0.5)

	xCoord = -1.0 + 2.0*xCoord;

      else

	xCoord = (xCoord - 0.5)*2.0;

      yCoord = modf(fabs(ptcls[i].getCoord(1)/normFactors[1]),&value);

      if(yCoord < 0.5)

	yCoord = -1.0 + 2.0*yCoord;

      else

	yCoord = (yCoord - 0.5)*2.0;

      zCoord = modf(fabs(ptcls[i].getCoord(2)/normFactors[2]),&value);

      if(zCoord < 0.5)

	zCoord = -1.0 + 2.0*zCoord;

      else

	zCoord = (zCoord - 0.5)*2.0;

      zCoord = 0;


      logFile<<"xCoord["<<i<<"]: "<<ptcls[i].getCoord(0)<<" -> "<<xCoord<<endl;
      logFile<<"yCoord["<<i<<"]: "<<ptcls[i].getCoord(1)<<" -> "<<yCoord<<endl;
      logFile<<"zCoord["<<i<<"]: "<<ptcls[i].getCoord(2)<<" -> "<<zCoord<<endl;

      ptcleParams[i][0] = QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);

      ptcleParams[i][1] = normFactors[0]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][2] = normFactors[1]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][3] = normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      ptcleParams[i][4] = pow(normFactors[0],2)*QuarticSplineWindowFunc::d2W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][5] = pow(normFactors[1],2)*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::d2W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][6] = pow(normFactors[2],2)*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::d2W(zCoord);
      ptcleParams[i][7] = normFactors[0]*normFactors[1]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][8] = normFactors[1]*normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::dW(zCoord);
      ptcleParams[i][9] = normFactors[0]*normFactors[2]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      logFile<<"ptcleParams["<<i<<"]: "<<ptcleParams[i][0]<<endl;
    }


    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      xCoord = modf(fabs(pltCoords[i][0]/normFactors[0]),&value);

      if(xCoord < 0.5)

	xCoord = -1.0 + 2.0*xCoord;

      else

	xCoord = (xCoord - 0.5)*2.0;

      yCoord = modf(fabs(pltCoords[i][1]/normFactors[1]),&value);

      if(yCoord < 0.5)

	yCoord = -1.0 + 2.0*yCoord;

      else

	yCoord = (yCoord - 0.5)*2.0;

      zCoord = modf(fabs(pltCoords[i][2]/normFactors[2]),&value);

      if(zCoord < 0.5)

	zCoord = -1.0 + 2.0*zCoord;

      else

	zCoord = (zCoord - 0.5)*2.0;


      zCoord = 0;

      logFile<<"xCoord["<<i<<"]: "<<pltCoords[i][0]<<" -> "<<xCoord<<endl;
      logFile<<"yCoord["<<i<<"]: "<<pltCoords[i][1]<<" -> "<<yCoord<<endl;
      logFile<<"zCoord["<<i<<"]: "<<pltCoords[i][2]<<" -> "<<zCoord<<endl;

      pltOrds[i][1] = QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);

      pltOrds[i][2] = normFactors[0]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][3] = normFactors[1]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][4] = normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      pltOrds[i][5] = pow(normFactors[0],2)*QuarticSplineWindowFunc::d2W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][6] = pow(normFactors[1],2)*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::d2W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][7] = pow(normFactors[2],2)*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::d2W(zCoord);
      pltOrds[i][8] = normFactors[0]*normFactors[1]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][9] = normFactors[1]*normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::dW(zCoord);
      pltOrds[i][10]= normFactors[0]*normFactors[2]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      logFile<<"func["<<i<<"]: "<<pltOrds[i][1]<<endl;
    }

    break;

  case 7:

    normFactors[0] = fabs(minCoords[0]-maxCoords[0])/2.0;
    normFactors[1] = fabs(minCoords[1]-maxCoords[1])/2.0;
    normFactors[2] = fabs(minCoords[2]-maxCoords[2])/2.0;


    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {

      xCoord = ptcls[i].getCoord(0);
      yCoord = ptcls[i].getCoord(1);
      zCoord = ptcls[i].getCoord(2);

      ptcleParams[i][0] = cos(xCoord)*cos(yCoord)*cos(zCoord);

      ptcleParams[i][1] = (-1)*sin(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][2] = (-1)*cos(xCoord)*sin(yCoord)*cos(zCoord);
      ptcleParams[i][3] = (-1)*cos(xCoord)*cos(yCoord)*sin(zCoord);

      ptcleParams[i][4] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][5] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][6] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][7] = sin(xCoord)*sin(yCoord)*cos(zCoord);
      ptcleParams[i][8] = cos(xCoord)*sin(yCoord)*sin(zCoord);
      ptcleParams[i][9] = sin(xCoord)*cos(yCoord)*sin(zCoord);

      logFile<<"ptcleParams["<<i<<"]: "<<ptcleParams[i][0]<<endl;
    }


    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      xCoord = pltCoords[i][0];
      yCoord = pltCoords[i][1];
      zCoord = pltCoords[i][2];

      pltOrds[i][1] = cos(xCoord)*cos(yCoord)*cos(zCoord);

      pltOrds[i][2] = (-1)*sin(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][3] = (-1)*cos(xCoord)*sin(yCoord)*cos(zCoord);
      pltOrds[i][4] = (-1)*cos(xCoord)*cos(yCoord)*sin(zCoord);

      pltOrds[i][5] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][6] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][7] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][8] = sin(xCoord)*sin(yCoord)*cos(zCoord);
      pltOrds[i][9] = cos(xCoord)*sin(yCoord)*sin(zCoord);
      pltOrds[i][10] = sin(xCoord)*cos(yCoord)*sin(zCoord);

    }

    break;

  default:
    logFile<<"Choosen function polynom type isn't supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  cout<<"finished calculation of the exact sample values"<<endl;

  /*********************************************************************/
  // Plot the real function.
  ofstream func,d0func,d1func,d2func,d00func,d11func,d22func,d01func,
    d12func,d20func;
  func.open("func.grf");
  d0func.open("d0func.grf");
  d1func.open("d1func.grf");
  d2func.open("d2func.grf");

  d00func.open("d00func.grf");
  d11func.open("d11func.grf");
  d22func.open("d22func.grf");

  d01func.open("d01func.grf");
  d12func.open("d12func.grf");
  d20func.open("d20func.grf");

  // Write the header of real function.
  func<<"# Graf: \"exact function\"\n"
      <<"#\n"
      <<"# X: \"plot-ordinate\" Y: \"function-value\""
      <<endl;
  d0func<<"# Graf: \"X derivation of the exact function\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: "
	<<"\"X-derivation-function-value\""
	<<endl;
  d1func<<"# Graf: \"Y derivation of the exact function\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: "
	<<"\"Y-derivation-function-value\""
	<<endl;
  d2func<<"# Graf: \"Z derivation of the exact function\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: "
	<<"\"Z-derivation-function-value\""
	<<endl;

  d00func<<"# Graf: \"XX derivation of the exact function\"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"XX-derivation-function-value\""
	 <<endl;
  d11func<<"# Graf: \"YY derivation of the exact function\"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"YY-derivation-function-value\""
	 <<endl;
  d22func<<"# Graf: \"ZZ derivation of the exact function \"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"ZZ-derivation-function-value\""
	 <<endl;

  d01func<<"# Graf: \"XY derivation of the exact function\"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"XY-derivation-function-value\""
	 <<endl;
  d12func<<"# Graf: \"YZ derivation of the exact function\"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"YZ-derivation-function-value\""
	 <<endl;
  d20func<<"# Graf: \"ZX derivation of the exact function\"\n"
	 <<"#\n"
         <<"# X: \"plot-ordinate\" Y: "
	 <<"\"ZX-derivation-function-value\""
	 <<endl;

  // Loop over all plot points and plot the true function values.

  for(int i=0;i<numOfPltPoints;i++) {

    func<<pltOrds[i][0]<<" "<<pltOrds[i][1]<<endl;

    d0func<<pltOrds[i][0]<<" "<<pltOrds[i][2]<<endl;
    d1func<<pltOrds[i][0]<<" "<<pltOrds[i][3]<<endl;
    d2func<<pltOrds[i][0]<<" "<<pltOrds[i][4]<<endl;

    if(secondOrderDerivs) {
      d00func<<pltOrds[i][0]<<" "<<pltOrds[i][5]<<endl;
      d11func<<pltOrds[i][0]<<" "<<pltOrds[i][6]<<endl;
      d22func<<pltOrds[i][0]<<" "<<pltOrds[i][7]<<endl;

      d01func<<pltOrds[i][0]<<" "<<pltOrds[i][8]<<endl;
      d12func<<pltOrds[i][0]<<" "<<pltOrds[i][9]<<endl;
      d20func<<pltOrds[i][0]<<" "<<pltOrds[i][10]<<endl;
    }
    else {
      d00func<<pltOrds[i][0]<<" 0"<<endl;
      d11func<<pltOrds[i][0]<<" 0"<<endl;
      d22func<<pltOrds[i][0]<<" 0"<<endl;

      d01func<<pltOrds[i][0]<<" 0"<<endl;
      d12func<<pltOrds[i][0]<<" 0"<<endl;
      d20func<<pltOrds[i][0]<<" 0"<<endl;
    }

  }

  cout<<"finished plotting of the real function"<<endl;


  /*********************************************************************/
  // Set shape functions and calculate the approximated function values.

  dbVector shapes;
  dbVector xDerivShapes;
  dbVector yDerivShapes;
  dbVector zDerivShapes;

  dbVector xxDerivShapes;
  dbVector yyDerivShapes;
  dbVector zzDerivShapes;

  dbVector xyDerivShapes;
  dbVector yzDerivShapes;
  dbVector zxDerivShapes;

  ShepardShapeFunc ShepardShapeSet;
  dbMatrix dShapes(3);
  dbMatrix d2Shapes(6);

  dbMatrix approxParams(numOfPltPoints,dbVector(10));

#ifdef _commonDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"********** shape function computation ***************"<<endl;
#endif

  int shapefuncType = (int)InputData->getValue("shapefunctionType");

  // calculation of the shapefunctions including their first order
  // derivation.
  for(int i=0;i<numOfPltPoints;i++) {
    supportSize = suppPtcls[i].size();

#ifdef _commonDebugMode_
    logFile<<"POINT "<<i<<": "<<endl;
#endif

    shapes.resize(supportSize);
    xDerivShapes.resize(supportSize);
    yDerivShapes.resize(supportSize);
    zDerivShapes.resize(supportSize);
    xxDerivShapes.resize(supportSize);
    yyDerivShapes.resize(supportSize);
    zzDerivShapes.resize(supportSize);
    xyDerivShapes.resize(supportSize);
    yzDerivShapes.resize(supportSize);
    zxDerivShapes.resize(supportSize);

    clearArray(shapes);
    clearArray(xDerivShapes);
    clearArray(yDerivShapes);
    clearArray(zDerivShapes);
    clearArray(xxDerivShapes);
    clearArray(yyDerivShapes);
    clearArray(zzDerivShapes);
    clearArray(xyDerivShapes);
    clearArray(yzDerivShapes);
    clearArray(zxDerivShapes);

    switch(shapefuncType) {

    case 1:

      dShapes[0] = xDerivShapes;
      dShapes[1] = yDerivShapes;
      dShapes[2] = zDerivShapes;
      d2Shapes[0] = xxDerivShapes;
      d2Shapes[1] = yyDerivShapes;
      d2Shapes[2] = zzDerivShapes;
      d2Shapes[3] = xyDerivShapes;
      d2Shapes[4] = yzDerivShapes;
      d2Shapes[5] = zxDerivShapes;

      // Shepard shapefunction set
      ShepardShapeSet.calcShapes(InputData,supportSize,suppPtcls[i],
				 ptcls,pltCoords[i][0],pltCoords[i][1],
				 pltCoords[i][2],shapes,dShapes,
				 d2Shapes,modelData,logFile,viewerSEQ);

      xDerivShapes = dShapes[0];
      yDerivShapes = dShapes[1];
      zDerivShapes = dShapes[2];
      xxDerivShapes = d2Shapes[0];
      yyDerivShapes = d2Shapes[1];
      zzDerivShapes = d2Shapes[2];
      xyDerivShapes = d2Shapes[3];
      yzDerivShapes = d2Shapes[4];
      zxDerivShapes = d2Shapes[5];


      break;

      // Calculate a EFG shapefunction set.
    case 2:

      EFGShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
			       ptcls,pltCoords[i][0],pltCoords[i][1],
			       pltCoords[i][2],shapes,
			       xDerivShapes,yDerivShapes,
			       zDerivShapes,xxDerivShapes,
			       yyDerivShapes,zzDerivShapes,
			       xyDerivShapes,yzDerivShapes,
			       zxDerivShapes,modelData,logFile,
			       viewerSEQ);

      break;

      // Calculate RKPM  shapefunction set.
    case 3:

      RKPMShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				ptcls,pltCoords[i][0],pltCoords[i][1],
				pltCoords[i][2],shapes,
				xDerivShapes,yDerivShapes,
				zDerivShapes,xxDerivShapes,
				yyDerivShapes,zzDerivShapes,
				xyDerivShapes,yzDerivShapes,
				zxDerivShapes,modelData,logFile,
				viewerSEQ);

      break;

      // Calculate orthogonalised MLS shapefunction set.
    case 4:

      OrthoShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				 ptcls,pltCoords[i][0],pltCoords[i][1],
				 pltCoords[i][2],shapes,
				 xDerivShapes,yDerivShapes,
				 zDerivShapes,xxDerivShapes,
				 yyDerivShapes,zzDerivShapes,
				 xyDerivShapes,yzDerivShapes,
				 zxDerivShapes,modelData,logFile,
				 viewerSEQ);

      break;

      // Calculate a MLS shapefunction set, where the influence zones
      // are different in negative and positive coordinate direction.
    case 5:

      AsymShapeFunc::calcShapes(InputData,supportSize,suppPtcls[i],
				ptcls,pltCoords[i][0],pltCoords[i][1],
				pltCoords[i][2],shapes,
				xDerivShapes,yDerivShapes,
				zDerivShapes,xxDerivShapes,
				yyDerivShapes,zzDerivShapes,
				xyDerivShapes,yzDerivShapes,
				zxDerivShapes,modelData,logFile,
				viewerSEQ);

      break;

    default:
      cerr<<"Choosen shape function type isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }


#ifdef _commonDebugMode_
    double PUM=0;
    logFile<<"----------------------------------------------------"<<endl;
    PUM = 0;
    logFile<<"shapes : ";
    for(int j=0;j<suppPtcls[i].size();j++) {
      PUM += shapes[j];
      logFile<<suppPtcls[i][j]<<": "<<shapes[j]<<endl;
    }
    if(fabs(1.0-PUM) > 0.000000001)
      logFile<<" PUM= "<<PUM<<endl;
#endif

    for(int j=0;j<suppPtcls[i].size();j++) {
      approxParams[i][0] += shapes[j]*ptcleParams[suppPtcls[i][j]][0];

      approxParams[i][1] += xDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][2] += yDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][3] += zDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];

      approxParams[i][4] += xxDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][5] += yyDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][6] += zzDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][7] += xyDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][8] += yzDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][9] += zxDerivShapes[j]*ptcleParams[suppPtcls[i][j]][0];
    }

  }

  cout<<"finished calculation of the approximated sample values"<<endl;

  /*********************************************************************/
  // Plot the approximation of the real function
  ofstream approx,d0approx,d1approx,d2approx,d00approx,d11approx,
    d22approx,d01approx,d12approx,d20approx;
  approx.open("approx.grf");
  d0approx.open("d0approx.grf");
  d1approx.open("d1approx.grf");
  d2approx.open("d2approx.grf");

  d00approx.open("d00approx.grf");
  d11approx.open("d11approx.grf");
  d22approx.open("d22approx.grf");

  d01approx.open("d01approx.grf");
  d12approx.open("d12approx.grf");
  d20approx.open("d20approx.grf");

  // Write the header of the approximation.
  approx<<"# Graf: \"function approximation\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: \"approximation-value\""
	<<endl;
  d0approx<<"# Graf: \"X derivation of the function approximation\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"X-derivation-approximation-value\""
	  <<endl;
  d1approx<<"# Graf: \"Y derivation of the function approximation\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"Y-derivation-approximation-value\""
	  <<endl;
  d2approx<<"# Graf: \"Z derivation of the function approximation\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"Z-derivation-approximation-value\""
	  <<endl;

  d00approx<<"# Graf: \"XX derivation of the function approximation\"\n"
	   <<"#\n"
           <<"# X: \"plot-ordinate\" Y: "
	   <<"\"XX-derivation-approximation-value\""
	   <<endl;
  d11approx<<"# Graf: \"YY derivation of the function approximation\"\n"
	   <<"#\n"
           <<"# X: \"plot-ordinate\" Y: "
	   <<"\"YY-derivation-approximation-value\""
	   <<endl;
  d22approx<<"# Graf: \"ZZ derivation of the function approximation \"\n"
	   <<"#\n"
           <<"# X: \"plot-ordinate\" Y: "
	   <<"\"ZZ-derivation-approximation-value\""
	   <<endl;

  d01approx<<"# Graf: \"XY derivation of the function approximation\"\n"
	   <<"#\n"
           <<"# X: \"plot-ordinate\" Y: "
	   <<"\"XY-derivation-approximation-value\""
	   <<endl;
  d12approx<<"# Graf: \"YZ derivation of the function approximation\"\n"
	   <<"#\n"
           <<"# X: \"plot-ordinate\" Y: "
	   <<"\"YZ-derivation-approximation-value\""
	   <<endl;
  d20approx<<"# Graf: \"ZX derivation of the function approximation\"\n"
	   <<"#\n"
           <<"# X: \"plot-ordinate\" Y: "
	   <<"\"ZX-derivation-approximation-value\""
	   <<endl;

  // Loop over all plot points and plot the approximated function
  // values.

  for(int i=0;i<numOfPltPoints;i++) {

    approx<<pltOrds[i][0]<<" "<<approxParams[i][0]<<endl;

    d0approx<<pltOrds[i][0]<<" "<<approxParams[i][1]<<endl;
    d1approx<<pltOrds[i][0]<<" "<<approxParams[i][2]<<endl;
    d2approx<<pltOrds[i][0]<<" "<<approxParams[i][3]<<endl;

    if(secondOrderDerivs) {
      d00approx<<pltOrds[i][0]<<" "<<approxParams[i][4]<<endl;
      d11approx<<pltOrds[i][0]<<" "<<approxParams[i][5]<<endl;
      d22approx<<pltOrds[i][0]<<" "<<approxParams[i][6]<<endl;

      d01approx<<pltOrds[i][0]<<" "<<approxParams[i][7]<<endl;
      d12approx<<pltOrds[i][0]<<" "<<approxParams[i][8]<<endl;
      d20approx<<pltOrds[i][0]<<" "<<approxParams[i][9]<<endl;
    }
    else {
      d00approx<<pltOrds[i][0]<<" 0"<<endl;
      d11approx<<pltOrds[i][0]<<" 0"<<endl;
      d22approx<<pltOrds[i][0]<<" 0"<<endl;

      d01approx<<pltOrds[i][0]<<" 0"<<endl;
      d12approx<<pltOrds[i][0]<<" 0"<<endl;
      d20approx<<pltOrds[i][0]<<" 0"<<endl;
    }

  }

  cout<<"finished plotting of the approximation"<<endl;

  /*********************************************************************/
  // Plot the residuum of approximation and real function
  ofstream residuum,d0residuum,d1residuum,d2residuum,d00residuum,d11residuum,
    d22residuum,d01residuum,d12residuum,d20residuum;
  residuum.open("residuum.grf");
  residuum.precision(15);
  d0residuum.open("d0residuum.grf");
  d1residuum.open("d1residuum.grf");
  d2residuum.open("d2residuum.grf");

  d00residuum.open("d00residuum.grf");
  d11residuum.open("d11residuum.grf");
  d22residuum.open("d22residuum.grf");

  d01residuum.open("d01residuum.grf");
  d12residuum.open("d12residuum.grf");
  d20residuum.open("d20residuum.grf");

  // Write the header of the residuum.
  residuum<<"# Graf: \" residuum\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: \"residuum-value[%]\""
	  <<endl;
  d0residuum<<"# Graf: \"X derivation of the residuum\"\n"
	    <<"#\n"
            <<"# X: \"plot-ordinate\" Y: "
	    <<"\"X-derivation-residuum-value[%]\""
	    <<endl;
  d1residuum<<"# Graf: \"Y derivation of the residuum\"\n"
	    <<"#\n"
            <<"# X: \"plot-ordinate\" Y: "
	    <<"\"Y-derivation-residuum-value[%]\""
	    <<endl;
  d2residuum<<"# Graf: \"Z derivation of the residuum\"\n"
	    <<"#\n"
            <<"# X: \"plot-ordinate\" Y: "
	    <<"\"Z-derivation-residuum-value[%]\""
	    <<endl;

  d00residuum<<"# Graf: \"XX derivation of the residuum\"\n"
	     <<"#\n"
             <<"# X: \"plot-ordinate\" Y: "
	     <<"\"XX-derivation-residuum-value[%]\""
	     <<endl;
  d11residuum<<"# Graf: \"YY derivation of the residuum\"\n"
	     <<"#\n"
             <<"# X: \"plot-ordinate\" Y: "
	     <<"\"YY-derivation-residuum-value[%]\""
	     <<endl;
  d22residuum<<"# Graf: \"ZZ derivation of the residuum \"\n"
	     <<"#\n"
             <<"# X: \"plot-ordinate\" Y: "
	     <<"\"ZZ-derivation-residuum-value[%]\""
	     <<endl;

  d01residuum<<"# Graf: \"XY derivation of the residuum\"\n"
	     <<"#\n"
             <<"# X: \"plot-ordinate\" Y: "
	     <<"\"XY-derivation-residuum-value[%]\""
	     <<endl;
  d12residuum<<"# Graf: \"YZ derivation of the residuum\"\n"
	     <<"#\n"
             <<"# X: \"plot-ordinate\" Y: "
	     <<"\"YZ-derivation-residuum-value[%]\""
	     <<endl;
  d20residuum<<"# Graf: \"ZX derivation of the residuum\"\n"
	     <<"#\n"
             <<"# X: \"plot-ordinate\" Y: "
	     <<"\"ZX-derivation-residuum-value[%]\""
	     <<endl;

  // Loop over all plot points and plot the residuum
  double difference,yValue;

  // plot absolute residuum
  if(absoluteResiduum) {

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      //------- f
      difference = fabs(pltOrds[i][1] - approxParams[i][0]);
      residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      //------- dfx
      difference = fabs(pltOrds[i][2] - approxParams[i][1]);
      d0residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      //------- dfy
      difference = fabs(pltOrds[i][3] - approxParams[i][2]);
      d1residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      //------- dfz
      difference = fabs(pltOrds[i][4] - approxParams[i][3]);
      d2residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      if(secondOrderDerivs) {

	//------- dfxx
	difference = fabs(pltOrds[i][5] - approxParams[i][4]);
	d00residuum<<pltOrds[i][0]<<" "<<difference<<endl;

	//------- dfyy
	difference = fabs(pltOrds[i][6] - approxParams[i][5]);
	d11residuum<<pltOrds[i][0]<<" "<<difference<<endl;

	//------- dfzz
	difference = fabs(pltOrds[i][7] - approxParams[i][6]);
	d22residuum<<pltOrds[i][0]<<" "<<difference<<endl;

	//------- dfxy
	difference = fabs(pltOrds[i][8] - approxParams[i][7]);
	d01residuum<<pltOrds[i][0]<<" "<<difference<<endl;

	//------- dfyz
	difference = fabs(pltOrds[i][9] - approxParams[i][8]);
	d12residuum<<pltOrds[i][0]<<" "<<difference<<endl;

	//------- dfzx
	difference = fabs(pltOrds[i][10] - approxParams[i][9]);
	d20residuum<<pltOrds[i][0]<<" "<<difference<<endl;
      }
      else {
	d00residuum<<pltOrds[i][0]<<" 0"<<endl;
	d11residuum<<pltOrds[i][0]<<" 0"<<endl;
	d22residuum<<pltOrds[i][0]<<" 0"<<endl;

	d01residuum<<pltOrds[i][0]<<" 0"<<endl;
	d12residuum<<pltOrds[i][0]<<" 0"<<endl;
	d20residuum<<pltOrds[i][0]<<" 0"<<endl;
      }

    }

  }
  //---------------------------------------------------------------------
  // plot relative residuum
  else {

    //dbVector storedValues(10);

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      //------- f
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][1]) - fabs(approxParams[i][0]);
	yValue = fabs(difference/pltOrds[i][1])*100.0;

	if(yValue < maxPercentage)

	  residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else
	  residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;

      }

      //------- dfx
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][2]) - fabs(approxParams[i][1]);
	yValue = fabs(difference/pltOrds[i][2])*100.0;

	if(yValue < maxPercentage)

	  d0residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else

	  d0residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
      }

      //------- dfy
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][3]) - fabs(approxParams[i][2]);
	yValue = fabs(difference/pltOrds[i][3])*100.0;

	if(yValue < maxPercentage)

	  d1residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else

	  d1residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
      }

      //------- dfz
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][4]) - fabs(approxParams[i][3]);
	yValue = fabs(difference/pltOrds[i][4])*100.0;

	if(yValue < maxPercentage)

	  d2residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else
	  d2residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;

      }

      if(secondOrderDerivs) {

	//------- dfxx
	if(fabs(pltOrds[i][1]) < minValue)
	  pltOrds[i][5] = minValue;

	difference = fabs(pltOrds[i][5]) - fabs(approxParams[i][4]);

	if(fabs(difference) < minValue && fabs(pltOrds[i][5]) < minValue)
	  d00residuum<<pltOrds[i][0]<<" 0"<<endl;
	else {
	  yValue = fabs(difference/pltOrds[i][5])*100.0;
	  if(yValue < maxPercentage)
	    d00residuum<<pltOrds[i][0]<<" "<<yValue<<endl;
	  else
	    d00residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
	}

	//------- dfyy
	if(fabs(pltOrds[i][1]) < minValue)
	  pltOrds[i][6] = minValue;

	difference = fabs(pltOrds[i][6]) - fabs(approxParams[i][5]);

	if(fabs(difference) < minValue && fabs(pltOrds[i][6]) < minValue)
	  d11residuum<<pltOrds[i][0]<<" 0"<<endl;
	else {
	  yValue = fabs(difference/pltOrds[i][6])*100.0;
	  if(yValue < maxPercentage)
	    d11residuum<<pltOrds[i][0]<<" "<<yValue<<endl;
	  else
	    d11residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
	}

	//------- dfzz
	if(fabs(pltOrds[i][1]) < minValue)
	  pltOrds[i][7] = minValue;

	difference = fabs(pltOrds[i][7]) - fabs(approxParams[i][6]);

	if(fabs(difference) < minValue && fabs(pltOrds[i][7]) < minValue)
	  d22residuum<<pltOrds[i][0]<<" 0"<<endl;
	else {
	  yValue = fabs(difference/pltOrds[i][7])*100.0;
	  if(yValue < maxPercentage)
	    d22residuum<<pltOrds[i][0]<<" "<<yValue<<endl;
	  else
	    d22residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
	}

	//------- dfxy
	if(fabs(pltOrds[i][1]) < minValue)
	  pltOrds[i][8] = minValue;

	difference = fabs(pltOrds[i][8]) - fabs(approxParams[i][7]);

	if(fabs(difference) < minValue && fabs(pltOrds[i][8]) < minValue)
	  d01residuum<<pltOrds[i][0]<<" 0"<<endl;
	else {
	  yValue = fabs(difference/pltOrds[i][8])*100.0;
	  if(yValue < maxPercentage)
	    d01residuum<<pltOrds[i][0]<<" "<<yValue<<endl;
	  else
	    d01residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
	}

	//------- dfyz
	if(fabs(pltOrds[i][1]) < minValue)
	  pltOrds[i][9] = minValue;

	difference = fabs(pltOrds[i][9]) - fabs(approxParams[i][8]);

	if(fabs(difference) < minValue && fabs(pltOrds[i][9]) < minValue)
	  d12residuum<<pltOrds[i][0]<<" 0"<<endl;
	else {
	  yValue = fabs(difference/pltOrds[i][9])*100.0;
	  if(yValue < maxPercentage)
	    d12residuum<<pltOrds[i][0]<<" "<<yValue<<endl;
	  else
	    d12residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
	}

	//------- dfzx
	if(fabs(pltOrds[i][1]) < minValue)
	  pltOrds[i][10] = minValue;

	difference = fabs(pltOrds[i][10]) - fabs(approxParams[i][9]);

	if(fabs(difference) < minValue && fabs(pltOrds[i][10]) < minValue)
	  d20residuum<<pltOrds[i][0]<<" 0"<<endl;
	else {
	  yValue = fabs(difference/pltOrds[i][10])*100.0;
	  if(yValue < maxPercentage)
	    d20residuum<<pltOrds[i][0]<<" "<<yValue<<endl;
	  else
	    d20residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
	}

      }
      else {
	d00residuum<<pltOrds[i][0]<<" 0"<<endl;
	d11residuum<<pltOrds[i][0]<<" 0"<<endl;
	d22residuum<<pltOrds[i][0]<<" 0"<<endl;

	d01residuum<<pltOrds[i][0]<<" 0"<<endl;
	d12residuum<<pltOrds[i][0]<<" 0"<<endl;
	d20residuum<<pltOrds[i][0]<<" 0"<<endl;
      }

    }

  }

  cout<<"finished plotting of the residuum"<<endl;

}

/************************************************************************/
/************************************************************************/
// Test the curve fitting property of a MaxEnt shape functions and its 
// first within the domain and on its boundary.

void testMaxEnt(InputFileData* InputData,std::vector<Particle>& ptcls,
                std::map<std::string,double>& modelData,
                std::ofstream& logFile,PetscViewer& viewerSEQ) {

  using namespace std;

  cout << " You are in testMaxEnt " << endl;

  int numOfPtcls = ptcls.size();

  int usedDims = (int)modelData["usedDimensions"];

  int windowFuncType = (int)InputData->getValue("windowFunctionType");
  int polynomialType = (int)InputData->getValue("basisPolynomType");
  int polynomialOrder = (int)InputData->getValue("basisPolynomOrder");
  int shapefunctionType = (int)InputData->getValue("shapefunctionType");

  logFile<<"******************************************************"<<endl;
  logFile<<"****************** MaxEnt test routine ******************"<<endl;
  logFile<<"DBL_EPSILON = "<<DBL_EPSILON<<endl;

  bool symmetricRadii;
  int numOfPltPoints,startPtcle,endPtcle,pltPtcle,pltPoint;
  bool absoluteResiduum;
  double minValue,maxPercentage;

  int shapeFuncPlotPosition,functionOrder,functionType;

  string dummyArray;
  ifstream inputFile("mls-test-input.dat");

  if(!inputFile) {
    logFile<<"Can't open input file mls-test-input.dat!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  // Read input file
  if(inputFile) {

    // plot path
    inputFile >> dummyArray >> numOfPltPoints;
    inputFile >> dummyArray >> startPtcle;
    inputFile >> dummyArray >> endPtcle;

    // test function
    inputFile >> dummyArray >> functionType;
    inputFile >> dummyArray >> functionOrder;

    // plot residuum

    // minimal ratio between the exact and the approximated function values
    inputFile >> dummyArray >> minValue;

    // maximal plotted percentage value for the relative residuum
    inputFile >> dummyArray >> maxPercentage;

    // plot the absolute values of the residuum
    inputFile >> dummyArray >> absoluteResiduum;

    logFile<<"numOfPltPoints = "<<numOfPltPoints<<endl;
    logFile<<"startPtcle = "<<startPtcle<<endl;
    logFile<<"endPtcle = "<<endPtcle<<endl;

    logFile<<"windowFunctionType = "<<windowFuncType<<endl;
    logFile<<"basisPolynomType = "<<polynomialType<<endl;
    logFile<<"basisPolynomOrder = "<<polynomialOrder<<endl;
    logFile<<"shapefunctionType = "<<shapefunctionType<<endl;

    logFile<<"functionType = "<<functionType<<endl;
    logFile<<"functionOrder = "<<functionOrder<<endl;

    logFile<<"minValue = "<<minValue<<endl;
    logFile<<"maxPercentage = "<<maxPercentage<<endl;
    logFile<<"absoluteResiduum = "<<absoluteResiduum<<endl;
  }

  else {
    cerr <<"Can't open input file 'mls-test-input.dat'!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }


  /*********************************************************************/
  // set a particle vector
  bool secondOrderDerivs = true;

  dbVector maxCoords(3);
  dbVector minCoords(3);

  logFile<<"******************************************************"<<endl;

  for(int i=0;i<numOfPtcls;i++) {
    double& weight = ptcls[i].getWeight();
    weight = 1.0;
    dbVector& coords = ptcls[i].getCoords();
    dbVector& radii = ptcls[i].getRadii();

    logFile<<"PARTICLE "<<i<<": "<<coords[0]<<" "<<coords[1]<<" "
	   <<coords[2]<<endl;

    for(int j=0;j<3;j++) {

      if(coords[j] > maxCoords[j])

	maxCoords[j] = coords[j];

      if(i == 0)

	minCoords[j] = coords[j];

      else if(coords[j] < minCoords[j])

	minCoords[j] = coords[j];

    }

  }

  for(int j=0;j<3;j++) {

    logFile<<"maxCoords["<<j<<"] = "<<maxCoords[j]<<endl;
    logFile<<"minCoords["<<j<<"] = "<<minCoords[j]<<endl;

  }

  /*********************************************************************/
  // Set plot points.

  if(startPtcle > ptcls.size() || endPtcle > ptcls.size()) {
    cerr <<"In testMaxEnt startPtcle or endPtcle nonexisting!"<< endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  if(startPtcle == 0 && endPtcle == 0) {
    startPtcle = 0;
    endPtcle = numOfPtcls-1;
  }
  else {
    startPtcle -= 1;
    endPtcle -= 1;
  }
  pltPtcle -= 1;

  double xStart = ptcls[startPtcle].getCoord(0);
  double yStart = ptcls[startPtcle].getCoord(1);
  double zStart = ptcls[startPtcle].getCoord(2);
  double xEnd = ptcls[endPtcle].getCoord(0);
  double yEnd = ptcls[endPtcle].getCoord(1);
  double zEnd = ptcls[endPtcle].getCoord(2);

  dbMatrix pltOrds(numOfPltPoints,dbVector(11));
  dbMatrix pltCoords(numOfPltPoints,dbVector(3));

  intMatrix suppPtcls(numOfPltPoints);

  double deltaX = xEnd-xStart;
  double deltaY = yEnd-yStart;
  double deltaZ = zEnd-zStart;

  double deltaPltPath =
    sqrt(pow(deltaX,2) + pow(deltaY,2) + pow(deltaZ,2))/(numOfPltPoints-1);

  double plotStartValue = -deltaPltPath/2.0*(numOfPltPoints-1);

  deltaX /= (numOfPltPoints-1);
  deltaY /= (numOfPltPoints-1);
  deltaZ /= (numOfPltPoints-1);

  logFile<<"******************************************************"<<endl;
  logFile<<"determine plot points"<<endl;
  logFile<<"startPtle = "<<startPtcle<<endl;
  logFile<<"endPtcle = "<<endPtcle<<endl;
  logFile<<"xStart = "<<xStart<<" -> xEnd = "<<xEnd
	 <<"; deltaX = "<<deltaX<<endl;
  logFile<<"yStart = "<<yStart<<" -> yEnd = "<<yEnd
	 <<"; deltaY = "<<deltaY<<endl;
  logFile<<"zStart = "<<zStart<<" -> zEnd = "<<zEnd
	 <<"; deltaZ = "<<deltaZ<<endl;
  logFile<<"numOfPltPoints = "<<numOfPltPoints<<"; deltaPltPath = "
	 <<deltaPltPath<<endl;

  for(int i=0;i<numOfPltPoints;i++) {

    pltOrds[i][0] = plotStartValue + i*deltaPltPath;

    pltCoords[i][0] = xStart + i*deltaX;
    pltCoords[i][1] = yStart + i*deltaY;
    pltCoords[i][2] = zStart + i*deltaZ;


    logFile<<"Plot-POINT "<<i<<": xCoord = "<<pltCoords[i][0]<<" yCoord = "
	   <<pltCoords[i][1]<<" zCoord = "<<pltCoords[i][2]<<" ord = "
	   <<pltOrds[i][0]<<endl;
  }


  /*********************************************************************/
  // set supporting particle list for all plot points

  for(int i=0;i<numOfPltPoints;i++) {
    for(int j=0;j<numOfPtcls;j++) {
      // Check if current point 'i' is supported by particle 'j'.
      if(ptcls[j].querySupported(InputData,pltCoords[i],modelData,
				 logFile)){
	suppPtcls[i].push_back(j);
      }
    }
  }

  //#ifdef _commonDebugMode_
  logFile<<"*******************************************************"<<endl;
  logFile<<"************* plot point support list *****************"<<endl;
  for(int i=0;i<suppPtcls.size();i++) {
    logFile<<"Plot-POINT "<<i<<"("<<pltCoords[i][0]<<", "
	   <<pltCoords[i][1]<<", "<<pltCoords[i][2]<<"): ";
    for(int j=0;j<suppPtcls[i].size();j++)
      logFile<<suppPtcls[i][j]<<" ";
    logFile<<endl;
  }
  //#endif


  /*********************************************************************/
  // set the particle parameter (DOF - normally unknown!)
  double xCoord,yCoord,zCoord;

  double dl;
  double norming = (double)1.0/25.000;
  int linEQSize;
  dbMatrix ptcleParams(numOfPtcls,dbVector(10));
  dbVector normFactors(3);

  double value;

  dbVector basis;
  dbVector xDerivBasis;
  dbVector yDerivBasis;
  dbVector zDerivBasis;

#ifdef _commonDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"************ ptcle params computation ***************"<<endl;
#endif

  // Select choosen polynom type (Pascal,Lagrangian,serendipity).
  switch(functionType) {
    
    // own creation
  case 0:

    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {
      xCoord = ptcls[i].getCoord(0)*norming;
      yCoord = ptcls[i].getCoord(1)*norming;
      zCoord = ptcls[i].getCoord(2)*norming;

      ptcleParams[i][0] = 100.0*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);

      ptcleParams[i][1] = 100.0*norming*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][2]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::W(zCoord);
      ptcleParams[i][3]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::dW(zCoord);
    }

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {
      xCoord = pltCoords[i][0]*norming;
      yCoord = pltCoords[i][1]*norming;
      zCoord = pltCoords[i][2]*norming;

      pltOrds[i][1] = 100.0*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);

      pltOrds[i][2] = 100.0*norming*GaussWindowFunc::dW(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][3]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::dW(yCoord)*GaussWindowFunc::W(zCoord);
      pltOrds[i][4]  = 100.0*norming*GaussWindowFunc::W(xCoord)*
	GaussWindowFunc::W(yCoord)*GaussWindowFunc::dW(zCoord);

    }

    break;

    // Pascal
  case 1:

    // Select choosen polynom order.
    switch(functionOrder) {

    case 1:
      linEQSize = 4;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  PascalLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalLinear::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  PascalLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalLinear::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	}

      }

      break;

    case 2:
      linEQSize = 10;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  PascalQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalQuadratic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  PascalQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalQuadratic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	}

      }

      break;

    case 3:
      linEQSize = 20;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);

      //norming = (double)1.0/1000.0; // mls-cube
      //dl = 11.0; // mls-cube
      dl = 1.0;
      norming = (double)1.0/pow(pltOrds[numOfPltPoints-1][0],3);

      for(int i=0;i<numOfPtcls;i++) {

	xCoord = (ptcls[i].getCoord(0)+dl)*norming;
	yCoord = (ptcls[i].getCoord(1)+dl)*norming;
	zCoord = (ptcls[i].getCoord(2)+dl)*norming;

	basis =  PascalCubic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalCubic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalCubic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalCubic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += norming*xDerivBasis[j];
	  ptcleParams[i][2] += norming*yDerivBasis[j];
	  ptcleParams[i][3] += norming*zDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = (pltCoords[i][0]+dl)*norming;
	yCoord = (pltCoords[i][1]+dl)*norming;
	zCoord = (pltCoords[i][2]+dl)*norming;

	basis =  PascalCubic::P(xCoord,yCoord,zCoord);

	xDerivBasis = PascalCubic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = PascalCubic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = PascalCubic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += norming*xDerivBasis[j];
	  pltOrds[i][3] += norming*yDerivBasis[j];
	  pltOrds[i][4] += norming*zDerivBasis[j];

	}

      }

      break;

    default:
      logFile<<"Chosen function polynom order isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Lagrangian
  case 2:

    switch(functionOrder) {

    case 1:
      linEQSize = 8;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  LagrangeLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeLinear::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  LagrangeLinear::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeLinear::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeLinear::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeLinear::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	}

      }

      break;

    case 2:
      linEQSize = 27;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  LagrangeQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeQuadratic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis =  LagrangeQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = LagrangeQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = LagrangeQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = LagrangeQuadratic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	}

      }

      break;

    default:
      logFile<<"Chosen function polynom order isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    /*******************************************************************/
    // Serendipity
  case 3:

    switch(functionOrder) {

    case 2:
      linEQSize = 20;

      basis.resize(linEQSize);
      xDerivBasis.resize(linEQSize);
      yDerivBasis.resize(linEQSize);
      zDerivBasis.resize(linEQSize);

      for(int i=0;i<numOfPtcls;i++) {
	xCoord = ptcls[i].getCoord(0);
	yCoord = ptcls[i].getCoord(1);
	zCoord = ptcls[i].getCoord(2);

	basis =  SerendipityQuadratic::P(xCoord,yCoord,zCoord);

	xDerivBasis = SerendipityQuadratic::dPx(xCoord,yCoord,zCoord);
	yDerivBasis = SerendipityQuadratic::dPy(xCoord,yCoord,zCoord);
	zDerivBasis = SerendipityQuadratic::dPz(xCoord,yCoord,zCoord);

	for(int j=0;j<linEQSize;j++) {
	  ptcleParams[i][0] += basis[j];

	  ptcleParams[i][1] += xDerivBasis[j];
	  ptcleParams[i][2] += yDerivBasis[j];
	  ptcleParams[i][3] += zDerivBasis[j];
	}

      }

      // loop over all plot points
      for(int i=0;i<numOfPltPoints;i++) {
	xCoord = pltCoords[i][0];
	yCoord = pltCoords[i][1];
	zCoord = pltCoords[i][2];

	basis.resize(linEQSize);
	xDerivBasis.resize(linEQSize);
	yDerivBasis.resize(linEQSize);
	zDerivBasis.resize(linEQSize);

	for(int j=0;j<linEQSize;j++) {
	  pltOrds[i][1] += basis[j];

	  pltOrds[i][2] += xDerivBasis[j];
	  pltOrds[i][3] += yDerivBasis[j];
	  pltOrds[i][4] += zDerivBasis[j];

	}

      }

      break;

    default:
      logFile<<"Choosen function polynom order isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
      break;
    }

    break;

    // tenth order spline
  case 5:

    norming = (double)1.0/35.0;

    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {
      xCoord = ptcls[i].getCoord(0)*norming;
      yCoord = ptcls[i].getCoord(1)*norming;
      zCoord = ptcls[i].getCoord(2)*norming;

      ptcleParams[i][0] = 100*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);

      ptcleParams[i][1] = 100*norming*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][2]  = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      ptcleParams[i][3]  = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);

      logFile<<"ptcleParams["<<i<<"]: "<<ptcleParams[i][0]<<endl;
    }

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {
      xCoord = pltCoords[i][0]*norming;
      yCoord = pltCoords[i][1]*norming;
      zCoord = pltCoords[i][2]*norming;

      pltOrds[i][1] = 100*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);

      pltOrds[i][2] = 100*norming*TenthOrderSplineWinFunc::dW(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][3] = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::dW(yCoord)*TenthOrderSplineWinFunc::W(zCoord);
      pltOrds[i][4] = 100*norming*TenthOrderSplineWinFunc::W(xCoord)*
	TenthOrderSplineWinFunc::W(yCoord)*TenthOrderSplineWinFunc::dW(zCoord);

    }

    break;

    // fourth order spline
  case 6:

    normFactors[0] = fabs(minCoords[0]-maxCoords[0])/2.0;
    normFactors[1] = fabs(minCoords[1]-maxCoords[1])/2.0;
    normFactors[2] = fabs(minCoords[2]-maxCoords[2])/2.0;


    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {

      xCoord = modf(fabs(ptcls[i].getCoord(0)/normFactors[0]),&value);

      if(xCoord < 0.5)

	xCoord = -1.0 + 2.0*xCoord;

      else

	xCoord = (xCoord - 0.5)*2.0;

      yCoord = modf(fabs(ptcls[i].getCoord(1)/normFactors[1]),&value);

      if(yCoord < 0.5)

	yCoord = -1.0 + 2.0*yCoord;

      else

	yCoord = (yCoord - 0.5)*2.0;

      zCoord = modf(fabs(ptcls[i].getCoord(2)/normFactors[2]),&value);

      if(zCoord < 0.5)

	zCoord = -1.0 + 2.0*zCoord;

      else

	zCoord = (zCoord - 0.5)*2.0;

      zCoord = 0;


      logFile<<"xCoord["<<i<<"]: "<<ptcls[i].getCoord(0)<<" -> "<<xCoord<<endl;
      logFile<<"yCoord["<<i<<"]: "<<ptcls[i].getCoord(1)<<" -> "<<yCoord<<endl;
      logFile<<"zCoord["<<i<<"]: "<<ptcls[i].getCoord(2)<<" -> "<<zCoord<<endl;

      ptcleParams[i][0] = QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);

      ptcleParams[i][1] = normFactors[0]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][2] = normFactors[1]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      ptcleParams[i][3] = normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      logFile<<"ptcleParams["<<i<<"]: "<<ptcleParams[i][0]<<endl;
    }


    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      xCoord = modf(fabs(pltCoords[i][0]/normFactors[0]),&value);

      if(xCoord < 0.5)

	xCoord = -1.0 + 2.0*xCoord;

      else

	xCoord = (xCoord - 0.5)*2.0;

      yCoord = modf(fabs(pltCoords[i][1]/normFactors[1]),&value);

      if(yCoord < 0.5)

	yCoord = -1.0 + 2.0*yCoord;

      else

	yCoord = (yCoord - 0.5)*2.0;

      zCoord = modf(fabs(pltCoords[i][2]/normFactors[2]),&value);

      if(zCoord < 0.5)

	zCoord = -1.0 + 2.0*zCoord;

      else

	zCoord = (zCoord - 0.5)*2.0;


      zCoord = 0;

      logFile<<"xCoord["<<i<<"]: "<<pltCoords[i][0]<<" -> "<<xCoord<<endl;
      logFile<<"yCoord["<<i<<"]: "<<pltCoords[i][1]<<" -> "<<yCoord<<endl;
      logFile<<"zCoord["<<i<<"]: "<<pltCoords[i][2]<<" -> "<<zCoord<<endl;

      pltOrds[i][1] = QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);

      pltOrds[i][2] = normFactors[0]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][3] = normFactors[1]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][4] = normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      pltOrds[i][5] = pow(normFactors[0],2)*QuarticSplineWindowFunc::d2W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][6] = pow(normFactors[1],2)*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::d2W(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][7] = pow(normFactors[2],2)*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::d2W(zCoord);
      pltOrds[i][8] = normFactors[0]*normFactors[1]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::W(zCoord);
      pltOrds[i][9] = normFactors[1]*normFactors[2]*QuarticSplineWindowFunc::W(xCoord)*
	QuarticSplineWindowFunc::dW(yCoord)*QuarticSplineWindowFunc::dW(zCoord);
      pltOrds[i][10]= normFactors[0]*normFactors[2]*QuarticSplineWindowFunc::dW(xCoord)*
	QuarticSplineWindowFunc::W(yCoord)*QuarticSplineWindowFunc::dW(zCoord);

      logFile<<"func["<<i<<"]: "<<pltOrds[i][1]<<endl;
    }

    break;

  case 7:

    normFactors[0] = fabs(minCoords[0]-maxCoords[0])/2.0;
    normFactors[1] = fabs(minCoords[1]-maxCoords[1])/2.0;
    normFactors[2] = fabs(minCoords[2]-maxCoords[2])/2.0;


    // loop over all particles
    for(int i=0;i<numOfPtcls;i++) {

      xCoord = ptcls[i].getCoord(0);
      yCoord = ptcls[i].getCoord(1);
      zCoord = ptcls[i].getCoord(2);

      ptcleParams[i][0] = cos(xCoord)*cos(yCoord)*cos(zCoord);

      ptcleParams[i][1] = (-1)*sin(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][2] = (-1)*cos(xCoord)*sin(yCoord)*cos(zCoord);
      ptcleParams[i][3] = (-1)*cos(xCoord)*cos(yCoord)*sin(zCoord);

      ptcleParams[i][4] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][5] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][6] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      ptcleParams[i][7] = sin(xCoord)*sin(yCoord)*cos(zCoord);
      ptcleParams[i][8] = cos(xCoord)*sin(yCoord)*sin(zCoord);
      ptcleParams[i][9] = sin(xCoord)*cos(yCoord)*sin(zCoord);

      logFile<<"ptcleParams["<<i<<"]: "<<ptcleParams[i][0]<<endl;
    }


    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      xCoord = pltCoords[i][0];
      yCoord = pltCoords[i][1];
      zCoord = pltCoords[i][2];

      pltOrds[i][1] = cos(xCoord)*cos(yCoord)*cos(zCoord);

      pltOrds[i][2] = (-1)*sin(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][3] = (-1)*cos(xCoord)*sin(yCoord)*cos(zCoord);
      pltOrds[i][4] = (-1)*cos(xCoord)*cos(yCoord)*sin(zCoord);

      pltOrds[i][5] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][6] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][7] = (-1)*cos(xCoord)*cos(yCoord)*cos(zCoord);
      pltOrds[i][8] = sin(xCoord)*sin(yCoord)*cos(zCoord);
      pltOrds[i][9] = cos(xCoord)*sin(yCoord)*sin(zCoord);
      pltOrds[i][10] = sin(xCoord)*cos(yCoord)*sin(zCoord);

    }

    break;

  default:
    logFile<<"Choosen function polynom type isn't supported!"<<endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    break;
  }

  cout<<"finished calculation of the exact sample values"<<endl;

  /*********************************************************************/
  // Plot the real function.
  ofstream func,d0func,d1func,d2func;
  func.open("func.grf");
  d0func.open("d0func.grf");
  d1func.open("d1func.grf");
  d2func.open("d2func.grf");

  // Write the header of real function.
  func<<"# Graf: \"exact function\"\n"
      <<"#\n"
      <<"# X: \"plot-ordinate\" Y: \"function-value\""
      <<endl;
  d0func<<"# Graf: \"X derivation of the exact function\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: "
	<<"\"X-derivation-function-value\""
	<<endl;
  d1func<<"# Graf: \"Y derivation of the exact function\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: "
	<<"\"Y-derivation-function-value\""
	<<endl;
  d2func<<"# Graf: \"Z derivation of the exact function\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: "
	<<"\"Z-derivation-function-value\""
	<<endl;

  // Loop over all plot points and plot the true function values.

  for(int i=0;i<numOfPltPoints;i++) {

    func<<pltOrds[i][0]<<" "<<pltOrds[i][1]<<endl;

    d0func<<pltOrds[i][0]<<" "<<pltOrds[i][2]<<endl;
    d1func<<pltOrds[i][0]<<" "<<pltOrds[i][3]<<endl;
    d2func<<pltOrds[i][0]<<" "<<pltOrds[i][4]<<endl;

  }

  cout<<"finished plotting of the real function"<<endl;


  int supportSize;

  dbVector shapes;
  dbMatrix firstDerivShapes;

  dbMatrix dShapes(3);
  dbMatrix d2Shapes(6);

  dbMatrix approxParams(numOfPltPoints,dbVector(10));

  // Calculating the Minimum distance between Particles
  double ptcleDist, minPtcleDist=9999;
  for(int i=0;i<numOfPtcls;i++) {
    for(int j=0;j<numOfPtcls;j++) {
      if (ptcls[i].querySupported(InputData,ptcls[j].getCoords(),ptcleDist,
				  modelData,logFile)){
	if (minPtcleDist>ptcleDist && i!=j)
	  minPtcleDist=ptcleDist;
      }

    }
  }

  cout << " Min PtcleDist : " << minPtcleDist << endl;


#ifdef _commonDebugMode_
  logFile<<"#####################################################"<<endl;
  logFile<<"********** shape function computation ***************"<<endl;
#endif

  int shapefuncType = (int)InputData->getValue("shapefunctionType");

  // calculation of the shapefunctions including their first order
  // derivation.
  for(int i=0;i<numOfPltPoints;i++) {
    supportSize = suppPtcls[i].size();

#ifdef _commonDebugMode_
    logFile<<"POINT "<<i<<": "<<endl;
#endif

    shapes.resize(supportSize);
    allocateArray(firstDerivShapes,usedDims,supportSize);

    clearArray(shapes);
    clearArray(firstDerivShapes);

    // Calculate a MaxEntropy shapefunction set, where the influence zones
    // are different in negative and positive coordinate direction.
    if(shapefuncType == 6) {

      MaxEntShapeFunc shapeSet(InputData,logFile);
      shapeSet.calcShapes(InputData,suppPtcls[i],
			  ptcls,pltCoords[i][0],pltCoords[i][1],
			  pltCoords[i][2],shapes,
			  firstDerivShapes,
			  modelData,logFile,viewerSEQ);

    }
    else {
      cerr<<"Choosen shape function type isn't supported!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

#ifdef _commonDebugMode_
    double PUM=0;
    logFile<<"----------------------------------------------------"<<endl;
    PUM = 0;
    logFile<<"shapes : ";
    for(int j=0;j<suppPtcls[i].size();j++) {
      PUM += shapes[j];
      logFile<<suppPtcls[i][j]<<": "<<shapes[j]<<endl;
    }
    if(fabs(1.0-PUM) > 0.000000001)
      logFile<<" PUM= "<<PUM<<endl;
#endif

    for(int j=0;j<suppPtcls[i].size();j++) {
      approxParams[i][0] += shapes[j]*ptcleParams[suppPtcls[i][j]][0];

      approxParams[i][1] += firstDerivShapes[0][j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][2] += firstDerivShapes[1][j]*ptcleParams[suppPtcls[i][j]][0];
      approxParams[i][3] += firstDerivShapes[2][j]*ptcleParams[suppPtcls[i][j]][0];
    }

  }

  cout<<"finished calculation of the approximated sample values"<<endl;

  /*********************************************************************/
  // Plot the approximation of the real function
  ofstream approx,d0approx,d1approx,d2approx,d00approx,d11approx,
    d22approx,d01approx,d12approx,d20approx;
  approx.open("approx.grf");
  d0approx.open("d0approx.grf");
  d1approx.open("d1approx.grf");
  d2approx.open("d2approx.grf");

  // Write the header of the approximation.
  approx<<"# Graf: \"function approximation\"\n"
	<<"#\n"
        <<"# X: \"plot-ordinate\" Y: \"approximation-value\""
	<<endl;
  d0approx<<"# Graf: \"X derivation of the function approximation\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"X-derivation-approximation-value\""
	  <<endl;
  d1approx<<"# Graf: \"Y derivation of the function approximation\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"Y-derivation-approximation-value\""
	  <<endl;
  d2approx<<"# Graf: \"Z derivation of the function approximation\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: "
	  <<"\"Z-derivation-approximation-value\""
	  <<endl;


  // Loop over all plot points and plot the approximated function
  // values.

  for(int i=0;i<numOfPltPoints;i++) {

    approx<<pltOrds[i][0]<<" "<<approxParams[i][0]<<endl;

    d0approx<<pltOrds[i][0]<<" "<<approxParams[i][1]<<endl;
    d1approx<<pltOrds[i][0]<<" "<<approxParams[i][2]<<endl;
    d2approx<<pltOrds[i][0]<<" "<<approxParams[i][3]<<endl;

  }

  cout<<"finished plotting of the approximation"<<endl;

  /*********************************************************************/
  // Plot the residuum of approximation and real function
  ofstream residuum,d0residuum,d1residuum,d2residuum;
  residuum.open("residuum.grf");
  residuum.precision(15);
  d0residuum.open("d0residuum.grf");
  d1residuum.open("d1residuum.grf");
  d2residuum.open("d2residuum.grf");

  // Write the header of the residuum.
  residuum<<"# Graf: \" residuum\"\n"
	  <<"#\n"
          <<"# X: \"plot-ordinate\" Y: \"residuum-value[%]\""
	  <<endl;
  d0residuum<<"# Graf: \"X derivation of the residuum\"\n"
	    <<"#\n"
            <<"# X: \"plot-ordinate\" Y: "
	    <<"\"X-derivation-residuum-value[%]\""
	    <<endl;
  d1residuum<<"# Graf: \"Y derivation of the residuum\"\n"
	    <<"#\n"
            <<"# X: \"plot-ordinate\" Y: "
	    <<"\"Y-derivation-residuum-value[%]\""
	    <<endl;
  d2residuum<<"# Graf: \"Z derivation of the residuum\"\n"
	    <<"#\n"
            <<"# X: \"plot-ordinate\" Y: "
	    <<"\"Z-derivation-residuum-value[%]\""
	    <<endl;



  // Loop over all plot points and plot the residuum
  double difference,yValue;

  // plot absolute residuum
  if(absoluteResiduum) {

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      //------- f
      difference = fabs(pltOrds[i][1] - approxParams[i][0]);
      residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      //------- dfx
      difference = fabs(pltOrds[i][2] - approxParams[i][1]);
      d0residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      //------- dfy
      difference = fabs(pltOrds[i][3] - approxParams[i][2]);
      d1residuum<<pltOrds[i][0]<<" "<<difference<<endl;

      //------- dfz
      difference = fabs(pltOrds[i][4] - approxParams[i][3]);
      d2residuum<<pltOrds[i][0]<<" "<<difference<<endl;


    }

  }
  //---------------------------------------------------------------------
  // plot relative residuum
  else {

    //dbVector storedValues(10);

    // loop over all plot points
    for(int i=0;i<numOfPltPoints;i++) {

      //------- f
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][1]) - fabs(approxParams[i][0]);
	yValue = fabs(difference/pltOrds[i][1])*100.0;

	if(yValue < maxPercentage)

	  residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else
	  residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;

      }

      //------- dfx
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][2]) - fabs(approxParams[i][1]);
	yValue = fabs(difference/pltOrds[i][2])*100.0;

	if(yValue < maxPercentage)

	  d0residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else

	  d0residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
      }

      //------- dfy
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][3]) - fabs(approxParams[i][2]);
	yValue = fabs(difference/pltOrds[i][3])*100.0;

	if(yValue < maxPercentage)

	  d1residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else

	  d1residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;
      }

      //------- dfz
      if(fabs(pltOrds[i][1]) > minValue) {

	difference = fabs(pltOrds[i][4]) - fabs(approxParams[i][3]);
	yValue = fabs(difference/pltOrds[i][4])*100.0;

	if(yValue < maxPercentage)

	  d2residuum<<pltOrds[i][0]<<" "<<yValue<<endl;

	else
	  d2residuum<<pltOrds[i][0]<<" "<<maxPercentage<<endl;

      }

    }

  }

  cout<<"finished plotting of the residuum"<<endl;

}

/************************************************************************/
/************************************************************************/
// Update all degrees of freedom.
void testUpdate(InputFileData* InputData,
                std::ofstream& logFile) {


  using namespace std;

  double pi = 3.14159265358979323846;
  double cycle = 2*pi;

  double absBeta3,absNewRotVec,absOldRotVec,absDeltaRotVec,alpha1,alpha2,
    sinusRotVec,sinusDeltaRotVec,arcTerm3;
  dbVector beta1(3),beta2(3),beta3(3);
    
  intMatrix e = getPermutations(3);

  logFile<<"#####################################################"<<endl;
  logFile<<"#################### test updating ##################"<<endl;
  logFile<<"sin(pi/2) = "<<sin(pi/2.0)<<" =? 1"<<endl;
  logFile<<"arcsin(-1) = "<<asin((double)-1.0)<<" =? "<<-pi/2.0<<endl;

  // example: product of rotations R_{i+1) = dR R_i

  // spinor R1 = alpha1 + i beta1 ( -> R_i )

  //r1 = 1.570796
  //r2 = 0
  //r3 = 0

  // spinor R2 = alpha2 + i beta2 ( -> dR )

  //w1 =  0
  //w2 = 1.570796
  //w3 = 0

  // spinor R3 = R1 R2 =
  //           = (alpha1*alpha2 - beta1*beta2) + i (alpha1 beta2 + alpha2 beta1 - beta1 x beta2) =
  //           = alpha3 + i beta3

  //r_1new = 1.209201
  //r_2new = 1.209201
  //r_3new = -1.209201

  /*********************************************************************/
  // read input
  dbVector oldDOF(3);
  oldDOF[0] = InputData->getValue("r1");
  oldDOF[1] = InputData->getValue("r2");
  oldDOF[2] = InputData->getValue("r3");

  dbVector deltaDOF(3);
  deltaDOF[0] = InputData->getValue("w1");
  deltaDOF[1] = InputData->getValue("w2");
  deltaDOF[2] = InputData->getValue("w3");

  logFile<<"oldDOF[0] = "<<oldDOF[0]<<endl;
  logFile<<"oldDOF[1] = "<<oldDOF[1]<<endl;
  logFile<<"oldDOF[2] = "<<oldDOF[2]<<endl;
  logFile<<"***"<<endl;
  logFile<<"w[0] = "<<deltaDOF[0]<<endl;
  logFile<<"w[1] = "<<deltaDOF[1]<<endl;
  logFile<<"w[2] = "<<deltaDOF[2]<<endl;
  logFile<<"***"<<endl;

  /*********************************************************************/
  // cut off full cycles of 2*pi, since the algorithm supports
  // rotation in the limit [-2*pi; 2*pi]
  int completedCycles;
  dbVector i1(3);

  absOldRotVec = sqrt(pow(oldDOF[0],2)
		      +pow(oldDOF[1],2)
		      +pow(oldDOF[2],2));

  completedCycles = (int)floor(absOldRotVec/cycle);

  logFile<<"absOldRotVec = "<<absOldRotVec<<endl;
  logFile<<"***"<<endl;
  logFile<<"completedCycles = "<<completedCycles<<endl;
  logFile<<"***"<<endl;

  // calculate normalized old rotation vector
  if(absOldRotVec > DBL_EPSILON)

    for(int j=0;j<oldDOF.size();j++)
      i1[j] = oldDOF[j]/absOldRotVec;

  else

    for(int j=0;j<oldDOF.size();j++)
      i1[j] = 0.0;

  logFile<<"i1[0] = "<<i1[0]<<endl;
  logFile<<"i1[1] = "<<i1[1]<<endl;
  logFile<<"i1[2] = "<<i1[2]<<endl;
  logFile<<"***"<<endl;

  // calculate the new absolute value of the old rotation vector
  // which in fact is the angle of the rotation about the resulting
  // axis i1
  if(completedCycles > 0) {

    absOldRotVec -= completedCycles*cycle;

    // calculate the cut-off old rotation vector
    for(int j=0;j<oldDOF.size();j++)
      oldDOF[j] = absOldRotVec*i1[j];

    logFile<<"absOldRotVec = "<<absOldRotVec<<" (cut-off)"<<endl;
    logFile<<"***"<<endl;
    logFile<<"cut-off DOF: "<<endl;
    logFile<<"oldDOF[0] = "<<oldDOF[0]<<endl;
    logFile<<"oldDOF[1] = "<<oldDOF[1]<<endl;
    logFile<<"oldDOF[2] = "<<oldDOF[2]<<endl;
    logFile<<"***"<<endl;

    absOldRotVec = sqrt(pow(oldDOF[0],2)
			+pow(oldDOF[1],2)
			+pow(oldDOF[2],2));

    logFile<<"absOldRotVec = "<<absOldRotVec<<endl;
    logFile<<"***"<<endl;
  }

  /*********************************************************************/
  // Determine the absolute value of the rotations vector and the
  // increment of the rotation vector.
  dbVector i2(3);

  // calculate the Euler parameters of the spinor which corresponds
  // to the old rotation vector

  alpha1 = cos(0.5*absOldRotVec);

  sinusRotVec = sin(0.5*absOldRotVec);

  // Loop over rotational degrees of freedom.
  for(int j=0;j<3;j++)
    beta1[j] = sinusRotVec*i1[j];

  //---------------------------------------------------------------------
  // calculate the Euler parameters of the spinor which corresponds
  // to the rotation vector increment

  absDeltaRotVec = sqrt(pow(deltaDOF[0],2)
			+pow(deltaDOF[1],2)
			+pow(deltaDOF[2],2));

  alpha2 = cos(0.5*absDeltaRotVec);


  // calculate normalized rotation vector increment
  if(absDeltaRotVec > DBL_EPSILON)

    for(int j=0;j<deltaDOF.size();j++)
      i2[j] = deltaDOF[j]/absDeltaRotVec;

  else

    for(int j=0;j<deltaDOF.size();j++)
      i2[j] = 0.0;

  logFile<<"i2[0] = "<<i2[0]<<endl;
  logFile<<"i2[1] = "<<i2[1]<<endl;
  logFile<<"i2[2] = "<<i2[2]<<endl;
  logFile<<"***"<<endl;

  sinusDeltaRotVec = sin(0.5*absDeltaRotVec);

  for(int j=0;j<3;j++)
    beta2[j] = sinusDeltaRotVec*i2[j];

  logFile<<"absOldRotVec = "<<absOldRotVec<<endl;
  logFile<<"absDeltaRotVec = "<<absDeltaRotVec<<endl;
  logFile<<"alpha1 = "<<alpha1<<endl;
  logFile<<"alpha2 = "<<alpha2<<endl;
  logFile<<"sinusRotVec = "<<sinusRotVec<<endl;
  logFile<<"sinusDeltaRotVec = "<<sinusDeltaRotVec<<endl;
  logFile<<"***"<<endl;
  logFile<<"beta1[0] = "<<beta1[0]<<endl;
  logFile<<"beta1[1] = "<<beta1[1]<<endl;
  logFile<<"beta1[2] = "<<beta1[2]<<endl;
  logFile<<"***"<<endl;
  logFile<<"beta2[0] = "<<beta2[0]<<endl;
  logFile<<"beta2[1] = "<<beta2[1]<<endl;
  logFile<<"beta2[2] = "<<beta2[2]<<endl;
  logFile<<"***"<<endl;

  /*********************************************************************/
  // calculate the new spinor
  double alpha3;

  alpha3 = alpha1*alpha2;

  for(int k=0;k<3;k++)
    alpha3 -= beta1[k]*beta2[k];

  // Loop over rotational degrees of freedom.
  for(int j=0;j<3;j++)
    beta3[j] = alpha1*beta2[j] + alpha2*beta1[j];

  // Loop over all necessary permutations.
  for(int j=0;j<6;j++)
    beta3[e[j][2]] -= e[j][3]*beta1[e[j][0]]*beta2[e[j][1]];

  absBeta3 = sqrt(pow(beta3[0],2)+pow(beta3[1],2)+
		  pow(beta3[2],2));

  logFile<<"alpha3 = "<<alpha3<<endl;
  logFile<<"***"<<endl;
  logFile<<"beta3[0] = "<<beta3[0]<<endl;
  logFile<<"beta3[1] = "<<beta3[1]<<endl;
  logFile<<"beta3[2] = "<<beta3[2]<<endl;
  logFile<<"***"<<endl;
  logFile<<"absBeta3 = "<<absBeta3<<endl;
  logFile<<"***"<<endl;

  /*********************************************************************/
  // extract the new rotation vector from its corresponding spinor
  dbVector tmpDOF(3);
  dbVector newDOF(3);

  if(fabs(alpha3 - pi/2.0) > DBL_EPSILON)
    arcTerm3 = 2*acos(alpha3)/sin(acos(alpha3));
  else
    arcTerm3 = -2.0;

  // Loop over rotational degrees of freedom.
  for(int j=0;j<3;j++)

    tmpDOF[j] = beta3[j]*arcTerm3;


  logFile<<"arcTerm3 = "<<arcTerm3<<endl;
  logFile<<"***"<<endl;
  logFile<<"tmpDOF[0] = "<<tmpDOF[0]<<endl;
  logFile<<"tmpDOF[1] = "<<tmpDOF[1]<<endl;
  logFile<<"tmpDOF[2] = "<<tmpDOF[2]<<endl;
  logFile<<"***"<<endl;

  double absRotVec = sqrt(pow(tmpDOF[0],2)
			  +pow(tmpDOF[1],2)
			  +pow(tmpDOF[2],2));

  dbVector normVec(3);

  for(int j=0;j<tmpDOF.size();j++)
    normVec[j] = tmpDOF[j]/absRotVec;

  logFile<<"i3[0] = "<<normVec[0]<<endl;
  logFile<<"i3[1] = "<<normVec[1]<<endl;
  logFile<<"i3[2] = "<<normVec[2]<<endl;
  logFile<<"***"<<endl;
  logFile<<"absTmpRotVec = "<<absRotVec<<endl;

  /*********************************************************************/
  // check if the rotation "jumps" in another 2*pi-cycle
  double absTmpRotVec;
  dbVector i3(3);

  int changeCycle = 0;
  double angle = 0;

  for(int j=0;j<tmpDOF.size();j++)
    i3[j] = beta3[j]/absBeta3;

  absTmpRotVec = sqrt(pow(tmpDOF[0],2)
		      +pow(tmpDOF[1],2)
		      +pow(tmpDOF[2],2));

  // calculate the angle between the normalized old rotation vector
  // and the normalized new rotation vector
  for(int j=0;j<tmpDOF.size();j++)
    angle += i1[j]*i3[j];

  // jump happens
  if(angle < 0) {

    // cycle increase
    if(alpha3 < 0)
      changeCycle = 1;

    // cycle decrease
    else if(alpha3 > 0)
      changeCycle = -1;

    else {
      logFile<<"In function CosseratContinuum::updateDOF jump criterion "
	     <<"can not be determined!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  }

  logFile<<"i3[0] = "<<i3[0]<<endl;
  logFile<<"i3[1] = "<<i3[1]<<endl;
  logFile<<"i3[2] = "<<i3[2]<<endl;
  logFile<<"***"<<endl;
  logFile<<"angle between i1 and i3: "<<angle<<endl;
  logFile<<"changeCycle: "<<changeCycle<<endl;
    logFile<<"***"<<endl;

    /*********************************************************************/
    // add the old cut off cycles and respect, if given, a increase or
    // decrease of the current cycle

    // still in the same cycle
    if(changeCycle == 0) {

        // calculate the actual absolute value of the new rotation vector
        // which in fact is the angle of the rotation about the resulting
        // axis i3
        absNewRotVec = completedCycles*cycle + absTmpRotVec;

        for(int j=0;j<tmpDOF.size();j++)
            newDOF[j] = absNewRotVec*i3[j];

    }

    // increase in another cycle
    else if(changeCycle == 1) {

        // reverse direction
        for(int j=0;j<i3.size();j++)
            i3[j] *= -1;


        absNewRotVec = (completedCycles+2)*cycle - absTmpRotVec;

        for(int j=0;j<tmpDOF.size();j++)
            newDOF[j] = absNewRotVec*i3[j];

    }

    // decrease in another cycle
    else {

        // reverse direction
        for(int j=0;j<i3.size();j++)
            i3[j] *= -1;

        absNewRotVec = completedCycles*cycle - absTmpRotVec;

        for(int j=0;j<tmpDOF.size();j++)
            newDOF[j] = absNewRotVec*i3[j];

    }


    logFile<<"absTmpRotVec = "<<absTmpRotVec<<endl;
    logFile<<"absNewRotVec = "<<absNewRotVec<<endl;
    logFile<<"oldCycles = "<<completedCycles<<endl;
    logFile<<"***"<<endl;

    absRotVec = sqrt(pow(newDOF[0],2)
                     +pow(newDOF[1],2)
                     +pow(newDOF[2],2));

    logFile<<"newDOF[0] = "<<newDOF[0]<<endl;
    logFile<<"newDOF[1] = "<<newDOF[1]<<endl;
    logFile<<"newDOF[2] = "<<newDOF[2]<<endl;
    logFile<<"***"<<endl;
    logFile<<"absNewRotVec = "<<absRotVec<<endl;
    logFile<<"***"<<endl;


    /*********************************************************************/
    // calculation of the new rotation tensor
    dbMatrix oldRotTens;
    dbMatrix deltaRotTens;
    dbMatrix simoRotTens;

    oldDOF[0] = InputData->getValue("r1");
    oldDOF[1] = InputData->getValue("r2");
    oldDOF[2] = InputData->getValue("r3");

    calcRotTens(oldDOF,oldRotTens,logFile);

    logFile<<"*****************************************************"<<endl;
    logFile<<"R_old"<<endl;
    for(int i=0;i<oldRotTens.size();i++)
        for(int j=0;j<oldRotTens[i].size();j++)
            logFile<<"R["<<i<<"]["<<j<<"] = "<<oldRotTens[i][j]<<endl;


    calcRotTens(deltaDOF,deltaRotTens,logFile);

    innerTensorProduct(deltaRotTens,oldRotTens,
                       simoRotTens,false,false,logFile);

    logFile<<"*****************************************************"<<endl;
    logFile<<"R_new = dR R_old (Simo)"<<endl;
    for(int i=0;i<simoRotTens.size();i++)
        for(int j=0;j<simoRotTens[i].size();j++)
            logFile<<"R["<<i<<"]["<<j<<"] = "<<simoRotTens[i][j]<<endl;

    dbMatrix spinorRotTens;

    calcRotTens(newDOF,spinorRotTens,logFile);

    logFile<<"*****************************************************"<<endl;
    logFile<<"R_new = function ( rotations_new )"<<endl;
    for(int i=0;i<spinorRotTens.size();i++)
        for(int j=0;j<spinorRotTens[i].size();j++)
            logFile<<"R["<<i<<"]["<<j<<"] = "<<spinorRotTens[i][j]<<endl;

    bool match = true;

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)

            if(fabs(simoRotTens[i][j] - spinorRotTens[i][j]) > 1.0e-06)
                match = false;

    logFile<<"*****************************************************"<<endl;
    logFile<<"R_new (Simo) = R_new (spinor): "<<match<<endl;
    cout<<"R_new (Simo) = R_new (spinor): "<<match<<endl;

}

/***********************************************************************/
/***********************************************************************/
void calcRotTens(dbVector& rotations,dbMatrix& rotationTens,
                 std::ofstream& logFile) {

  using namespace std;

  if(rotationTens.size() < 3) {
    rotationTens.resize(3);

    for(int i=0;i<3;i++)
      rotationTens[i].resize(3);

  }

  double absRotVec = sqrt(pow(rotations[0],2) + pow(rotations[1],2)
			  + pow(rotations[2],2));

  double sinusTerm,cosinusTerm;

  if (absRotVec > DBL_EPSILON) {
    sinusTerm = sin(absRotVec)/absRotVec;
    cosinusTerm = (1.0-cos(absRotVec))/pow(absRotVec,2);
  }
  else {
    sinusTerm = 1.0;
    cosinusTerm = 0.5;
  }

  logFile<<"****************************************************"<<endl;
  logFile<<"calculating rotation tensor"<<endl;
  logFile<<"absRotVec "<<absRotVec<<endl;
  logFile<<"sinusTerm "<<sinusTerm<<endl;
  logFile<<"cosinusTerm "<<cosinusTerm<<endl;

  rotationTens[0][0] = 1.0 - cosinusTerm*(pow(rotations[1],2)
					  +pow(rotations[2],2));
  rotationTens[1][0] = sinusTerm*rotations[2]
    + cosinusTerm*rotations[1]*rotations[0];

  rotationTens[2][0]= -sinusTerm*rotations[1]
    + cosinusTerm*rotations[2]*rotations[0];

  rotationTens[0][1]= -sinusTerm*rotations[2]
    + cosinusTerm*rotations[0]*rotations[1];

  rotationTens[1][1]= 1.0 - cosinusTerm*(pow(rotations[0],2)
					 +pow(rotations[2],2));
  rotationTens[2][1]= sinusTerm*rotations[0]
    + cosinusTerm*rotations[2]*rotations[1];

  rotationTens[0][2]= sinusTerm*rotations[1]
    + cosinusTerm*rotations[0]*rotations[2];

  rotationTens[1][2]= -sinusTerm*rotations[0]
    + cosinusTerm*rotations[1]*rotations[2];

  rotationTens[2][2]= 1.0 - cosinusTerm*(pow(rotations[0],2)
					 +pow(rotations[1],2));
}

/***********************************************************************/
/***********************************************************************/
// Sort the particles in 1 dimensions.
void sortParticles(std::vector<Particle>& ptcls,intVector& idx,
                   int sortingCoord,std::ofstream& logFile) {

  using namespace std;

  int particlesNum = ptcls.size();

  intVector& ptclsIdx = idx;
  dbVector ptclsVec(particlesNum);

  if(idx.size() < particlesNum) ptclsIdx.resize(particlesNum);

  // Loop over all particles and store their coordinates 1.
  for(int i=0;i<particlesNum;i++) {
    ptclsIdx[i] = i;
    ptclsVec[i] = ptcls[i].getCoord(sortingCoord);
  }

  // Sort the nodes according their coordinate 1.
  sortValuesIdx(ptclsVec,ptclsIdx,0,particlesNum-1);

  //#ifdef _commonDebugMode_
  logFile<<"************* sorted particles according "<<sortingCoord
	 <<"-direction ***********"<<endl;
  for(int i=0;i<particlesNum;i++)
    logFile<<"new: "<<i<<"; old: "<<ptclsIdx[i]
	   <<" coord = "<<ptclsVec[i]<<endl;
  //#endif

}

/***********************************************************************/
// test the modified boundary collocation method
void testBoundCollocation(InputFileData* InputData,
                          std::vector<Particle>& ptcls,
                          intMatrix& boundaryPtcls,
                          intVector& newIdx,
                          intVector& newDOFID,
                          std::map<std::string,double>& modelData,
                          std::ofstream& logFile,
                          PetscViewer& viewerMPI,
                          PetscViewer& viewerSEQ) {

  using namespace std;

  int usedDOF = (int)modelData["usedDegreesOfFreedom"];
  int particlesNum = ptcls.size();

  // choose boundary particle and a degree of freedom for testing
  int dof = 0;
  int boundPtcle = 0;

  if(boundPtcle >= particlesNum)
    boundPtcle = 0;

  boundPtcle = newIdx[boundPtcle];


  intVector matIdx(particlesNum);
  dbVector matValues(particlesNum);
  double* XCalc;

  // Set matrix containing temporary the inverted boundary particles'
  // shape functions.
  dbMatrix newShapes(particlesNum,dbVector(particlesNum));
  int m=0;

    
  logFile<<"#####################################################"<<endl;
  logFile<<"********** test boundary collocation method *********"<<endl;
  logFile<<"*****************************************************"<<endl;
  logFile<<"Initialised indices"<<endl;
  logFile<<"************ matrix A assembling ********************"<<endl;


  Mat A;
  MatCreateSeqAIJ(PETSC_COMM_SELF,particlesNum,
		  particlesNum,PETSC_NULL,PETSC_NULL,&A);
  MatSetFromOptions(A);
    
  // Loop over all particle a essential boundary conditions applied and
  // intialize Mat A row to row.
  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    intVector& sPtcls = ptcls[ptcle].getSupportPtcls();
    dbVector& shps = ptcls[ptcle].getShapeFuncs();

    // Determine the column position of matrix A entries to be stored.
    m=0;

    for(int k=0;k<sPtcls.size();k++)

      if(shps[k] != 0) {
	matIdx[m] = sPtcls[k];
	matValues[m] = shps[k];
	m++;
      }


    logFile<<j<<" Particle "<<boundaryPtcls[dof][j]<<" ";
    for(int k=0;k<sPtcls.size();k++)
      logFile<<matIdx[k]<<" ";
    logFile<<endl;


    MatSetValues(A,1,&ptcle,m,&matIdx[0],&matValues[0],
		 INSERT_VALUES);

  }

  // set the diagonal of the remaining lines
  for(int j=0;j<particlesNum;j++)

    if(newDOFID[j*usedDOF+dof] > -1)

      MatSetValue(A,j,j,1.0,INSERT_VALUES);

    
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);


  MatView(A,viewerSEQ);

  logFile<<"Matrix A"<<endl;

  for(int j=0;j<particlesNum;j++) {
    logFile<<j<<": ";

    if(newDOFID[j*usedDOF+dof] < 0) {
      dbVector& shps = ptcls[j].getShapeFuncs();

      for(int k=0;k<shps.size();k++)
	logFile<<shps[k]<<" ";
      logFile<<endl;
    }
    else
      logFile<<"non-boundary particle"<<endl;
  }
    
  // Create the solution vector and the right needed for the
  // PETSC-Solver.
  Vec X;
  Vec Y;
  VecCreateSeq(PETSC_COMM_SELF,particlesNum,&X);
  VecCreateSeq(PETSC_COMM_SELF,particlesNum,&Y);

  // Initialize Solver.
  KSP ksp;
  PC pc;
  KSPCreate(PETSC_COMM_SELF,&ksp);
  KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,"lu");
  //  PCLUSetUseInPlace(pc); // destroy the original matrix
  KSPSetType(ksp,"preonly"); // direct solving
  KSPSetFromOptions(ksp);

  // Make a loop to assemble vector Y and solve the equation system to
  // invert A.
  intVector vecIdx(2);
  dbVector oneVec(2);
  oneVec[0] = 0.0;
  oneVec[1] = 1.0;

  // Loop over all particles
  for(int j=0;j<particlesNum;j++) {

    if(j == 0) {
      VecSetValues(Y,1,&j,&oneVec[1],INSERT_VALUES);
    }
    else {
      vecIdx[0] = j-1;
      vecIdx[1] = j;
      VecSetValues(Y,2,&vecIdx[0],&oneVec[0],INSERT_VALUES);
    }

    KSPSolve(ksp,Y,X);

#ifdef _forceVecModificationDebugMode_
    //KSPView(ksp,viewerMPI);
    VecView(Y,viewerSEQ);
    VecView(X,viewerSEQ);
#endif

    VecGetArray(X,&XCalc);

    for(int k=0;k<particlesNum;k++)
      newShapes[k][j] = XCalc[k];

    VecRestoreArray(X,&XCalc);
  }
    
  // Destroy all petsc objects.
  destroyPETScSolver(ksp);
  destroyPETScMat(A);
  destroyPETScVec(X);
  destroyPETScVec(Y);

  logFile<<"******************************************************"<<endl;
  logFile<<"***** calculated modified shape functions ************"<<endl;
  for(int j=0;j<particlesNum;j++) {
    logFile<<j<<": ";
    for(int k=0;k<particlesNum;k++) {
      if(newDOFID[k*usedDOF+dof] < 0)
	logFile<<newShapes[j][k]<<" ";
    }
    logFile<<" || ";
    for(int k=0;k<particlesNum;k++) {
      if(newDOFID[k*usedDOF+dof] > -1)
	logFile<<newShapes[j][k]<<" ";
    }
    logFile<<endl;
  }

  /**********************************************************************/
  // Store for each particle a essential boundary conditions applied
  // its modified boundary neigbour shape functions.

  // Loop over all particles essential boundary particles are applied
  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    intVector& sPtcls = ptcls[ptcle].getSupportBPtcls(dof);
    dbVector& shps = ptcls[ptcle].getBShapeFuncs(dof);

    m=0;

    // Loop over all particles
    for(int k=0;k<particlesNum;k++)

      if(newShapes[ptcle][k] != 0) {

	// Store modified boundary neighbour shape functions.
	if(m < shps.size())
	  shps[m] = newShapes[ptcle][k];
	else
	  shps.push_back(newShapes[ptcle][k]);

	// Store again the supporting neighbour idenfifiers -
	// modification of shape function neighbours could changing zero
	// entries to valid values.
	if(m < sPtcls.size())
	  sPtcls[m] = k;
	else
	  sPtcls.push_back(k);

	m++;
      }

    if(m < shps.size())
      shps.resize(m);

    if(m < sPtcls.size())
      sPtcls.resize(m);

  }

  logFile<<"******************************************************"<<endl;
  logFile<<"********* stored modified shape functions ************"<<endl;

  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    intVector& sPtcls = ptcls[ptcle].getSupportBPtcls(dof);
    dbVector& shps = ptcls[ptcle].getBShapeFuncs(dof);

    logFile<<j<<".) BOUND PARTICLE "<<ptcle<<": ";

    for(int k=0;k<sPtcls.size();k++)

      if(newDOFID[sPtcls[k]*usedDOF+dof] < 0)
	logFile<<sPtcls[k]<<" ";

    logFile<<" || ";

    for(int k=0;k<sPtcls.size();k++)

      if(newDOFID[sPtcls[k]*usedDOF+dof] > -1)
	logFile<<sPtcls[k]<<" ";

    logFile<<endl;

    for(int k=0;k<sPtcls.size();k++)

      if(newDOFID[sPtcls[k]*usedDOF+dof] < 0)
	logFile<<shps[k]<<" ";


    logFile<<" || ";

    for(int k=0;k<sPtcls.size();k++)

      if(newDOFID[sPtcls[k]*usedDOF+dof] > -1)
	logFile<<shps[k]<<" ";

    logFile<<endl;
    logFile<<"*********"<<endl;

  }
    
  logFile<<"******************************************************"<<endl;
  logFile<<"********* test of boundary enforcement ***************"<<endl;

  // choose an essential boundary condition for all boundary particles
  double boundCond = 0.0;

  // set the degrees of freedom for all other particles
  double one = 1.0e+03;

  double boundDOF;

  // --------------------------------------------------------------------
  // set the modified degree of freedom vector

  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    ptcls[ptcle].setDOF(dof,boundCond);

  }

  for(int j=0;j<particlesNum;j++)

    if(newDOFID[j*usedDOF+dof] > -1)

      ptcls[j].setDOF(dof,one);

  // --------------------------------------------------------------------
  // set the real degree of freedom vector d = N-1 \tilde{d}
  // only the boundary degrees of freedom must be set - the non-boundary
  // one should be unchanged(?)

  vector<double> newDOFs(particlesNum);

  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    intVector& sPtcls = ptcls[ptcle].getSupportBPtcls(dof);
    dbVector& shps = ptcls[ptcle].getBShapeFuncs(dof);

    for(int k=0;k<sPtcls.size();k++)
      newDOFs[ptcle] += shps[k]*ptcls[sPtcls[k]].getDOF(dof);

    logFile<<"Ptcle "<<ptcle<<" newDOF="<<newDOFs[ptcle]<<endl;
  }

  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    ptcls[ptcle].setDOF(dof,newDOFs[ptcle]);
  }

  logFile<<"-----------------------------------------------------"<<endl;
  for(int j=0;j<particlesNum;j++)
    logFile<<"ptcle "<<j<<": "<<ptcls[j].getDOF(dof)<<endl;

  //---------------------------------------------------------------------
  // calculate the function value at the chosen boundary particle
  double funcValue = 0;
    
  intVector& sPtcls = ptcls[boundPtcle].getSupportPtcls();
  dbVector& shps = ptcls[boundPtcle].getShapeFuncs();
    
  for(int j=0;j<shps.size();j++) {
    funcValue += shps[j]*ptcls[sPtcls[j]].getDOF(dof);
    logFile<<j<<".) "<<sPtcls[j]<<" "<<shps[j]<<" "<<ptcls[sPtcls[j]].getDOF(dof)<<endl;
  }
    
  logFile<<"funcValue = "<<funcValue<<endl;

  for(int j=0;j<boundaryPtcls[dof].size();j++) {

    int& ptcle = boundaryPtcls[dof][j];
    intVector& sPtcls = ptcls[ptcle].getSupportBPtcls(dof);
    dbVector& shps = ptcls[ptcle].getBShapeFuncs(dof);

    for(int k=0;k<sPtcls.size();k++)

      if(newDOFID[sPtcls[k]*usedDOF+dof] > -1) {
	sPtcls.erase(sPtcls.begin()+k);
	shps.erase(shps.begin()+k);
	--k;
      }

    sPtcls.resize(0);
    shps.resize(0);
  }

}

// test Carlo's 'dexpo' fortran routine
void testDexpo(dbMatrix& IN,dbMatrix& OUT,bool seskaStorage,
               std::ofstream& logFile) {

  using namespace std;


  if(OUT.size() < 9)
    allocateArray(OUT,9,9);

  clearArray(OUT);


  dbMatrix AA;
  innerTensorProduct(IN,IN,AA,false,false,logFile);

  dbMatrix delta = getKroneckerSymbol(IN.size());
  intMatrix m2v;

  // SESKA's storage scheme
  if(seskaStorage)

    m2v = matrixToVectorFull(IN.size());

  // Carlo's storage scheme
  else {
    allocateArray(m2v,IN.size(),IN.size());

    m2v[0][0] = 0;
    m2v[1][1] = 1;
    m2v[2][2] = 2;
    m2v[0][1] = 3;
    m2v[0][2] = 5;
    m2v[1][2] = 7;
    m2v[1][0] = 4;
    m2v[2][0] = 6;
    m2v[2][1] = 8;
  }

  for(int i=0;i<IN.size();i++)

    for(int j=0;j<IN.size();j++)

      for(int r=0;r<IN.size();r++)

	for(int s=0;s<IN.size();s++)

	  OUT[m2v[i][j]][m2v[r][s]] = delta[i][r]*delta[j][s]

	    +1.0/2.0*(delta[i][r]*IN[s][j]+delta[j][s]*IN[i][r])

	    +1.0/6.0*(delta[i][r]*AA[s][j]+IN[i][r]*IN[s][j]+delta[j][s]*AA[i][r]);

}
