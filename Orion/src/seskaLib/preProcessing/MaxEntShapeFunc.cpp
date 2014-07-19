#include "MaxEntShapeFunc.h"


/************************************************************************/
/************************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions.
void MaxEntShapeFunc::calcShapes(InputFileData* InputData,
                                 intVector& sPtcls,
                                 std::vector<Particle>& Ptcls,
                                 double& x, double& y,double& z,
                                 dbVector& shapeFuncs,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile,
                                 PetscViewer& viewerSEQ){
    using namespace std;

    // Extract the number of Dimension to be used
    int usedDims = (int)modelData["usedDimensions"];

    // Define the negative Jacobian Matrix
    dbMatrix ngtv_inv_J;
    allocateArray(ngtv_inv_J,usedDims,usedDims);

    // Define the main coordinate
    dbVector mainPtcleCoord(usedDims);
    mainPtcleCoord[0] = x;
    mainPtcleCoord[1] = y;
    mainPtcleCoord[2] = z;

    // Find the number of supporting particles
    int sPtclsNum = sPtcls.size();

    // Compute the Shape Functions with Newton Iteration
    shapeFuncs = gamma_(InputData,mainPtcleCoord,Ptcls,sPtcls,
                        sPtclsNum,usedDims,ngtv_inv_J,logFile);

}

/************************************************************************/
/************************************************************************/
// Calculate at a point for all its supporting particles their 
// shape functions.
void MaxEntShapeFunc::calcShapes(InputFileData* InputData,
                                 intVector& sPtcls,
                                 std::vector<Particle>& Ptcls,
                                 double& x, double& y,double& z,
                                 dbVector& shapeFuncs,
                                 std::ofstream& logFile){
    using namespace std;

    // Extract the number of Dimension to be used
    // int usedDims = (int)modelData["usedDimensions"];
    int usedDims=3;

    // Define the negative Jacobian Matrix
    dbMatrix ngtv_inv_J;
    allocateArray(ngtv_inv_J,usedDims,usedDims);

    // Define the main coordinate
    dbVector mainPtcleCoord(usedDims);
    mainPtcleCoord[0] = x;
    mainPtcleCoord[1] = y;
    mainPtcleCoord[2] = z;

    // Find the number of supporting particles
    int sPtclsNum = sPtcls.size();

    // Compute the Shape Functions with Newton Iteration
    shapeFuncs = gamma_(InputData,mainPtcleCoord,Ptcls,sPtcls,
                        sPtclsNum,usedDims,ngtv_inv_J,logFile);

}

/************************************************************************/
/************************************************************************/
// Calculate at a point for all its supporting particles their
// shape functions.
void MaxEntShapeFunc::calcShapes(InputFileData* InputData,
                                 intVector& sPtcls,
                                 std::vector<Particle>& Ptcls,
                                 double& x, double& y,double& z,
                                 dbVector& shapeFuncs,
                                 dbMatrix& firstDerivShapes,
                                 std::map<std::string,double>& modelData,
                                 std::ofstream& logFile,
                                 PetscViewer& viewerSEQ){
    using namespace std;

    // Extract the number of Dimension to be used
    int usedDims = (int)modelData["usedDimensions"];

    // Define the negative Jacobian Matrix
    dbMatrix ngtv_inv_J;
    allocateArray(ngtv_inv_J,usedDims,usedDims);

    // Define the main coordinate
    dbVector mainPtcleCoord(usedDims);
    mainPtcleCoord[0] = x;
    mainPtcleCoord[1] = y;
    mainPtcleCoord[2] = z;

    // Find the number of supporting particles
    int sPtclsNum = sPtcls.size();

    // Compute the Shape Functions with Newton Iteration
    shapeFuncs = gamma_(InputData,mainPtcleCoord,Ptcls,sPtcls,
                        sPtclsNum,usedDims,ngtv_inv_J,logFile);

    //=======================================================================
    // Compute the First Derivatives of the Shape Functions for constant Beta
//    dbMatrix firstDerivShapes_consBeta
//            = CalFirstDerivFunc(mainPtcleCoord,Ptcls,sPtcls,
//                                sPtclsNum,shapeFuncs,ngtv_inv_J,
//                                usedDims,logFile);

     //=======================================================================
     // Calculate the First derivative of the Shape Functions for varying Beta
     dbMatrix firstDerivShapes_delBeta
             = CalFirstDerivFunc_Beta(mainPtcleCoord,Ptcls,
                                      sPtcls,sPtclsNum,shapeFuncs,
                                      ngtv_inv_J,usedDims,logFile);

    //========================================================================
    // Calculate the First derivatives using A.Arroyo Algorithm
//    dbMatrix firstDerivShapes_A
//            = CalFirstDerivFunc_A(mainPtcleCoord,Ptcls,sPtcls,
//                                  sPtclsNum,shapeFuncs,
//                                  ngtv_inv_J,usedDims,logFile);


//      firstDerivShapes = firstDerivShapes_consBeta;
      firstDerivShapes = firstDerivShapes_delBeta;
//      firstDerivShapes = firstDerivShapes_A;
}



/************************************************************************/
/************************************************************************/
// Calculate the shape function Values for a particular particles
dbVector MaxEntShapeFunc::gamma_(InputFileData* InputData,
                                 dbVector& mainPtcleCoord,
                                 std::vector<Particle>&  Ptcls,
                                 intVector& sPtcls,
                                 int& sPtclsNum, int& usedDims,
                                 dbMatrix& ngtv_inv_J,
                                 std::ofstream& logFile)
{
    using namespace std;

    // Define the vectors to be used by the Newton iteration method
    dbVector lam(usedDims), R(usedDims), dlam(usedDims);

    int niter = 0;

    double gam,Z,gamma,target_zero,TolLag,dummy ;
    dbVector shapeFuncs(sPtclsNum),temp(sPtclsNum);

    // Define the Jacobian and its inverse matrix
    dbMatrix J, inv_J;
    allocateArray (J,usedDims,usedDims);
    allocateArray (inv_J,usedDims,usedDims);

    // Extract predefined Values
    target_zero = (double)InputData->getValue("target_zero");
    TolLag = (double)InputData->getValue("TolLag");

    // ===============================================================
    string dummyArray;
    ifstream inputFile("maxent-input.dat");

    if(!inputFile) {

    }else if (inputFile){
        // Read Value from input file
        inputFile >> dummyArray >> gamma;

        // Assign value to MaxEntGamma
        InputData->setValue("MaxEntGamma",gamma);
    }
    inputFile.close();
    // ===============================================================

    // assigne predefine values in vector lam, R, dlam;
    for (int i=0; i <usedDims; i++) {
        lam[i]= 0;
        R[i]= 10;
        dlam[i]= 10;
    }

    // Newton iteration loop
    while (computeNorm(R,2,logFile)>TolLag) {

        // Loop to compute the Partition Function, Z
        Z = CalPartitionFunc(InputData,mainPtcleCoord,Ptcls,sPtcls,
                             sPtclsNum,lam,temp,usedDims,logFile);

        // Loop to compute the shape function of each supporting nodes
        for (int ii=0; ii<sPtclsNum; ii++) {
            shapeFuncs[ii]=temp[ii]/Z;
        }

        gam = log(Z);

        //Compute the 1st derivative for Netwon Method
        R = CalNewtonFirstDeriv(mainPtcleCoord,Ptcls,sPtcls,sPtclsNum,
                                shapeFuncs,usedDims,logFile);

        //Compute the 2nd derivative for Newton Method
        J = CalNewtonSecDeriv(mainPtcleCoord,Ptcls,sPtcls,sPtclsNum,
                              shapeFuncs,R,usedDims,logFile);

        // Calculate the inverse of the Jacobian(Hessian) Matrix
        calcInvDoubleDense(J,inv_J,logFile);

        // Multiply the J matrix by -1
        for (int r=0; r<usedDims;r++){
            for (int c=0; c<usedDims;c++){
                ngtv_inv_J[r][c] = inv_J[r][c] * -1;
            }
        }

        // Calculate the increment of lambda
        innerTensorProduct(ngtv_inv_J,R,dlam,false,logFile);

        // Add dlam to lam
        for (int d=0; d<usedDims; d++){
            lam[d] += dlam[d];
        }

        // Calculate the actual number of iterations that the
        // newton method has gone through
        niter++;
        if (niter==100) {
            logFile<<"Newton Failed, No convergence in 100 iterations"
                  << std::endl;
        }

    }

    return shapeFuncs;
};

/************************************************************************/
/************************************************************************/
// Calculate the Partition Function, Z
double MaxEntShapeFunc::CalPartitionFunc(InputFileData* InputData,
                                         dbVector& mainPtcleCoord,
                                         std::vector<Particle>& Ptcls,intVector& sPtcls,
                                         int& sPtclsNum, dbVector& lam,
                                         dbVector& temp,int& usedDims,
                                         std::ofstream& logFile)
{

    using namespace std;

    double beta,gamma,sum1,sum2,Z=0;

    // Extract the gamma constant
    gamma = (double)InputData->getValue("MaxEntGamma");

    for (int i=0; i<sPtclsNum; i++) {
        sum1=0;sum2=0;

        // Extract the beta value corresponding to a particular
        // supporting particle
        beta = Ptcls[sPtcls[i]].getBeta();

        for (int id=0; id<usedDims; id++) {

            // Compute the sum of the prior for all dimensions
            sum1 += pow(mainPtcleCoord[id]-Ptcls[sPtcls[i]].getCoord(id),2);

            // Calculate the sum of MaxEnt for all dimensions
            sum2 += lam[id]*(mainPtcleCoord[id]-Ptcls[sPtcls[i]].getCoord(id));

        }

        // Compute and store the Partition value of each
        // supporting particle
        temp[i] = exp((-beta*sum1) + sum2);

        // Calculate the value of the Partition Function, Z
        Z += exp((-beta*sum1) + sum2);

    }

    return Z;
};

/************************************************************************/
/************************************************************************/
// Calculate Newton First Derivatives
dbVector MaxEntShapeFunc::CalNewtonFirstDeriv(dbVector& mainPtcleCoord,
                                              std::vector<Particle>& Ptcls,
                                              intVector& sPtcls, int& sPtclsNum,
                                              dbVector& shapeFuncs,int& usedDims,
                                              std::ofstream& logFile)
{
    using namespace std;

    double sum_dgam;
    dbVector dgam(usedDims);

    for (int id=0; id<usedDims; id++) {
        sum_dgam = 0;
        for (int idg=0; idg<sPtclsNum; idg++) {
            sum_dgam += (mainPtcleCoord[id]-Ptcls[sPtcls[idg]].getCoord(id))*shapeFuncs[idg];
        }
        dgam[id]=sum_dgam;
    }

    return dgam;
};

/************************************************************************/
/************************************************************************/
// Calculate Newton Second Derivatives
dbMatrix MaxEntShapeFunc::CalNewtonSecDeriv(dbVector& mainPtcleCoord, 
                                            std::vector<Particle>& Ptcls,
                                            intVector& sPtcls,int& sPtclsNum,
                                            dbVector& shapeFuncs,dbVector dgam,
                                            int& usedDims,std::ofstream& logFile)
{

    using namespace std;

    double sum_hgam;
    dbMatrix hgam;
    allocateArray(hgam, usedDims, usedDims);

    for (int id=0; id<usedDims; id++) {
        for (int jd=0; jd<usedDims; jd++) {
            sum_hgam = 0;
            for (int j=0; j<sPtclsNum; j++) {
                sum_hgam +=
                        shapeFuncs[j]*
                        (mainPtcleCoord[id]-Ptcls[sPtcls[j]].getCoord(id))
                        *(mainPtcleCoord[jd]-Ptcls[sPtcls[j]].getCoord(jd));
            }
            hgam[id][jd] = sum_hgam-(dgam[id]*dgam[jd]);
        }
    }

    return hgam;
};


/************************************************************************/
/************************************************************************/
// Calculate the First Derivatives of the shape functions
dbMatrix MaxEntShapeFunc::CalFirstDerivFunc(dbVector& mainPtcleCoord, 
                                            std::vector<Particle>& Ptcls,
                                            intVector& sPtcls,int& sPtclsNum,
                                            dbVector& shapeFuncs,dbMatrix& ngtv_inv_J,
                                            int& usedDims,std::ofstream& logFile)
{

    using namespace std;

    dbMatrix dp;
    allocateArray(dp,usedDims,sPtclsNum);

    dbVector coordDiffVec(usedDims), firstDerivNeighbourPtclVec;

    for (int k=0; k<sPtclsNum; k++) {
        for (int kd=0; kd < usedDims; kd++) {
            coordDiffVec[kd]=(mainPtcleCoord[kd]-Ptcls[sPtcls[k]].getCoord(kd))*shapeFuncs[k];
        }

        innerTensorProduct(ngtv_inv_J, coordDiffVec, firstDerivNeighbourPtclVec, false, logFile);

        for (int kk=0; kk < usedDims; kk++) {
            dp[kk][k] =firstDerivNeighbourPtclVec[kk] ;
        }
    }

    return dp;
};

/************************************************************************/
/************************************************************************/
// Calculate the First Derivatives of the shape functions
// (including delta Beta)
dbMatrix MaxEntShapeFunc::CalFirstDerivFunc_Beta(dbVector& mainPtcleCoord,
                                                 std::vector<Particle>& Ptcls,
                                                 intVector& sPtcls,
                                                 int& sPtclsNum,
                                                 dbVector& shapeFuncs,
                                                 dbMatrix& ngtv_inv_J,
                                                 int& usedDims,
                                                 std::ofstream& logFile)
{
    using namespace std;

    // Calculate the shift coordinates
    dbMatrix shiftX;
    allocateArray(shiftX,sPtclsNum,usedDims);

    for (int sX=0; sX<sPtclsNum; sX++) {
        for(int sXd=0;sXd<usedDims;sXd++){
            shiftX[sX][sXd]=
                    mainPtcleCoord[sXd]-Ptcls[sPtcls[sX]].getCoord(sXd);
        }
    }

    // Define equation KaA : sum{P_b * |x-x_b|^2 * (x-x_b)}
    dbVector KaA(usedDims,0);
    double KaD=0,sum_sqr;

    for(int iKa=0;iKa<sPtclsNum;iKa++){
        sum_sqr=0;
        for(int iKad=0;iKad<usedDims;iKad++){
            sum_sqr += pow(shiftX[iKa][iKad],2);
        }

        for(int iKadd=0;iKadd<usedDims;iKadd++){
            KaA[iKadd] += shapeFuncs[iKa]*sum_sqr*shiftX[iKa][iKadd];
        }

        // KaD: sum{P * |x-x_b|^2}
        KaD += shapeFuncs[iKa]*sum_sqr;
    }

    // Calculate Positive Jacobian Matrix
    // Multiply the J matrix by -1
    dbMatrix inv_J;
    allocateArray(inv_J,usedDims,usedDims);
    for (int r=0; r<usedDims;r++){
        for (int c=0; c<usedDims;c++){
            inv_J[r][c] = ngtv_inv_J[r][c] * -1;
        }
    }

    // Compute the drivatives of the shape functions at each
    // supporting nodes
    dbVector shiftXvec(usedDims);
    dbVector KaVec(sPtclsNum);

    dbVector coordDiffVec(usedDims), constBetaFirstDerivVec(usedDims);

    // Define the derivative matrix
    dbMatrix dp;
    allocateArray(dp,usedDims,sPtclsNum);

    double KaC,KaE=0;

    for (int k=0; k<sPtclsNum; k++) {

        // Extract the shift coordinates of a supporting particle
        shiftXvec = shiftX[k];

        // P_a(x-x_a)
        for (int kd=0; kd < usedDims; kd++) {
            coordDiffVec[kd]=shiftXvec[kd]*shapeFuncs[k];
        }

        // -P_a(J^-1)(x-x_a)
        innerTensorProduct(ngtv_inv_J, coordDiffVec,
                           constBetaFirstDerivVec,
                           false, logFile);

        // KaE : [sum(Pa*|x-x_b|^2(x-x_b)].(J^-1)(x-x_a)
        KaE=0;
        for (int i=0;i<usedDims;i++){
            for(int j=0;j<usedDims;j++){
                KaE += KaA[i]*inv_J[i][j]*shiftXvec[j];
            }
        }

        // |x-x_a|^2
        scalarProduct(shiftXvec,shiftXvec,KaC,logFile);

        // calculating KaE - KaC + KaD
        KaVec[k] = KaE - KaC + KaD;

        // Calculate the first derivatives of the shape functions
        // for a supporting particle
        for(int bd=0;bd<usedDims;bd++){
            dp[bd][k] =  constBetaFirstDerivVec[bd] +
                shapeFuncs[k]*KaVec[k]*Ptcls[sPtcls[k]].getBetaDerivs(bd);
        }
    }

    return dp;
}


/************************************************************************/
/************************************************************************/
// Calculate the First Derivatives of the shape functions using
// A. Arroyo method
dbMatrix MaxEntShapeFunc::CalFirstDerivFunc_A(dbVector& mainPtcleCoord,
                                              std::vector<Particle>& Ptcls,
                                              intVector& sPtcls,
                                              int& sPtclsNum,
                                              dbVector& shapeFuncs,
                                              dbMatrix& ngtv_inv_J,
                                              int& usedDims,
                                              std::ofstream& logFile){

    using namespace std;

    // Calculate the shifted Coordinate: Xa-X
    dbMatrix shiftCoord;
    allocateArray(shiftCoord,sPtclsNum,usedDims);

    for (int i=0; i<sPtclsNum; i++){
        for (int id=0; id<usedDims; id++){
            shiftCoord[i][id] =
                    Ptcls[sPtcls[i]].getCoord(id) - mainPtcleCoord[id];
        }
    }

    // Compute the H (Hessian) matrix
    dbMatrix H;
    allocateArray(H,usedDims,usedDims);

    double sumH;

    for(int rH=0;rH<usedDims; rH++){
        for(int cH=0;cH<usedDims; cH++){
            sumH = 0;
            for(int sH=0;sH<sPtclsNum;sH++){
                sumH +=
                  shapeFuncs[sH]*shiftCoord[sH][rH]*shiftCoord[sH][cH];
            }
            H[rH][cH]=sumH;
        }
    }

    // Calculate Values of prior and derivative of prior
    dbVector w(sPtclsNum);
    dbMatrix wDer;
    allocateArray(wDer,sPtclsNum,usedDims);

    double sumP;

    for (int p=0; p<sPtclsNum; p++){
        sumP=0;
        for (int ip=0; ip<usedDims; ip++){
            sumP += pow(shiftCoord[p][ip],2);
        }
        w[p]=Ptcls[sPtcls[p]].getBeta()*exp(sumP);

        for (int ipd=0; ipd<usedDims; ipd++){
            wDer[p][ipd]=
                 2*w[p]*shiftCoord[p][ipd]*Ptcls[sPtcls[p]].getBeta();
        }
    }

    // Compute MA matrix
    dbMatrix MA;
    allocateArray(MA,sPtclsNum,usedDims);

    for (int m=0; m<sPtclsNum; m++){
        for (int im=0; im<usedDims; im++){
            MA[m][im]=1/w[m]*wDer[m][im];
        }
    }

    // Compute A matrix
    dbMatrix A;
    allocateArray(A,usedDims,usedDims);

    double sumA;

    for (int rA=0;rA<usedDims;rA++){
        for (int cA=0;cA<usedDims;cA++){
            sumA = 0;
            for (int sA=0;sA<sPtclsNum;sA++){
                sumA += shapeFuncs[sA]*MA[sA][cA]*shiftCoord[sA][rA];
            }
            A[rA][cA]=sumA;
        }
    }

    // Build the MC vector
    dbVector MC(usedDims);

    double sumMC;

    for(int c=0;c<usedDims;c++){
        sumMC = 0;
        for(int sc=0;sc<sPtclsNum;sc++){
            sumMC += shapeFuncs[sc]*MA[sc][c];
        }
        MC[c]=sumMC;
    }

    // --------------------------------------------------------------

    // Calculate the inverse of the Hessian Matrix, H
    dbMatrix inv_H;
    allocateArray(inv_H,usedDims,usedDims);
    calcInvDoubleDense(H, inv_H,logFile);

    // product of inv_H(Matrix X) and A(MAtrix Y)
    dbMatrix mult_invH_A;
    allocateArray(mult_invH_A,usedDims,usedDims);

    innerTensorProduct(inv_H,A,mult_invH_A,false,false,logFile);

    // invH - invH*A
    dbMatrix subMat;
    allocateArray(subMat,usedDims,usedDims);
    for (int rS=0;rS<usedDims;rS++){
        for (int cS=0;cS<usedDims;cS++){
            subMat[rS][cS] = inv_H[rS][cS]-mult_invH_A[rS][cS];
        }
    }

    // product of shiftCoord and (invH - invH*A)
    dbMatrix mult_temp;
    allocateArray(mult_temp,sPtclsNum,usedDims);

    innerTensorProduct(shiftCoord,subMat,mult_temp,false,false,logFile);

    // shapeFunc .* mult_temp
    dbMatrix mult_partA;
    allocateArray(mult_partA,sPtclsNum,usedDims);
    for (int rPA=0; rPA<sPtclsNum; rPA++){
        for (int cPA=0; cPA<usedDims; cPA++){
            mult_partA[rPA][cPA] = shapeFuncs[rPA]*mult_temp[rPA][cPA];
        }
    }

    // shapeFunc .* MA
    dbMatrix mult_partB;
    allocateArray(mult_partB,sPtclsNum,usedDims);
    for (int rPB=0; rPB<sPtclsNum; rPB++){
        for (int cPB=0; cPB<usedDims; cPB++){
            mult_partB[rPB][cPB] = shapeFuncs[rPB]*MA[rPB][cPB];
        }
    }

    // shapeFunc * MC
    dbMatrix mult_partC;
    allocateArray(mult_partC,sPtclsNum,usedDims);
    for(int rPC=0;rPC<sPtclsNum;rPC++){
        for(int cPC=0;cPC<usedDims;cPC++){
            mult_partC[rPC][cPC]=shapeFuncs[rPC]*MC[cPC];
        }
    }

    // subtraction mult_PartB and mult_PartC
    dbMatrix mult_partD;
    allocateArray(mult_partD,sPtclsNum,usedDims);
    for(int rPD=0;rPD<sPtclsNum;rPD++){
        for(int cPD=0;cPD<usedDims;cPD++){
            mult_partD[rPD][cPD]=
                    mult_partB[rPD][cPD]-mult_partC[rPD][cPD];
        }
    }

    // adding mult_PartA and mult_PartD
    dbMatrix dp_temp;
    allocateArray(dp_temp,sPtclsNum,usedDims);
    for(int rdpt=0;rdpt<sPtclsNum;rdpt++){
        for(int cdpt=0;cdpt<usedDims;cdpt++){
            dp_temp[rdpt][cdpt]=
                    mult_partA[rdpt][cdpt]+mult_partD[rdpt][cdpt];
        }
    }

    // Rearranging derivatrive matrix
    dbMatrix dp;
    allocateArray(dp,usedDims,sPtclsNum);
    for(int rdp=0;rdp<usedDims;rdp++){
        for(int cdp=0;cdp<sPtclsNum;cdp++){
            dp[rdp][cdp]=dp_temp[cdp][rdp];
        }
    }

    return dp;
}
