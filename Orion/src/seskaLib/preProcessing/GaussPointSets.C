#include "GaussPointSets.h"


/***********************************************************************/
GaussSetLine1::GaussSetLine1() {

  allocateArray(coord,1,1);
  weight = dbVector(1);

  coord[0][0] = 0.0;

  weight[0] =  2.0;

}

/***********************************************************************/
GaussSetLine2::GaussSetLine2() {

  allocateArray(coord,2,1);
  weight = dbVector(2);

  coord[0][0] = -0.577350269189626;
  coord[1][0] = 0.577350269189626;

  weight[0] =  1.0;
  weight[1] =  1.0;

}

/***********************************************************************/
GaussSetLine3::GaussSetLine3() {

  allocateArray(coord,3,1);
  weight = dbVector(3);

  coord[0][0] = -0.774596669241483;
  coord[1][0] = 0.0;
  coord[2][0] = 0.774596669241483;

  weight[0] =  0.555555555555556;
  weight[1] =  0.888888888888889;
  weight[2] =  0.555555555555556;

}

/***********************************************************************/
GaussSetLine4::GaussSetLine4() {

  allocateArray(coord,4,1);
  weight = dbVector(4);

  coord[0][0] = -0.861136311594953;
  coord[1][0] = -0.339981043584856;
  coord[2][0] = 0.339981043584856;
  coord[3][0] = 0.861136311594953;

  weight[0] =  0.347854845137454;
  weight[1] =  0.652145154862546;
  weight[2] =  0.652145154862546;
  weight[3] =  0.347854845137454;
}

/***********************************************************************/
GaussSetLine5::GaussSetLine5() {
  
  allocateArray(coord,5,1);
  weight = dbVector(5);

  coord[0][0] = 0; 
  coord[1][0] = 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)); 
  coord[2][0] = -1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)); 
  coord[3][0] = 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)); 
  coord[4][0] = -1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)); 
  
  weight[0] = 128.0/225.0;
  weight[1] = (322.0+13.0*sqrt(70.0))/900.0;
  weight[2] = (322.0+13.0*sqrt(70.0))/900.0;
  weight[3] = (322.0-13.0*sqrt(70.0))/900.0;
  weight[4] = (322.0-13.0*sqrt(70.0))/900.0;

}

/***********************************************************************/
/***********************************************************************/
GaussSetTria1::GaussSetTria1() {

  allocateArray(coord,1,3);
  weight = dbVector(1);

  coord[0][0] = 0.333333333333333;
  coord[0][1] = 0.333333333333333;
  coord[0][2] = 0.333333333333333;

  weight[0] =  1.0;

}

GaussSetTetra1::GaussSetTetra1() {

  allocateArray(coord,1,4);
  weight = dbVector(1);

  coord[0][0] = 0.25;
  coord[0][1] = 0.25;
  coord[0][2] = 0.25;
  coord[0][3] = 0.25;

  weight[0] = 1.0;

}

/***********************************************************************/
GaussSetTria3::GaussSetTria3() {

  allocateArray(coord,3,3);
  weight = dbVector(3);

  coord[0][0] = 0.5;
  coord[0][1] = 0.5;
  coord[0][2] = 0;

  coord[1][0] = 0;
  coord[1][1] = 0.5;
  coord[1][2] = 0.5;

  coord[2][0] = 0.5;
  coord[2][1] = 0;
  coord[2][2] = 0.5;

  weight[0] = 0.333333333333333;
  weight[1] = 0.333333333333333;
  weight[2] = 0.333333333333333;

}

GaussSetTetra4::GaussSetTetra4() {

  allocateArray(coord,4,4);
  weight = dbVector(4);

  // default Zienkiewcz positioning
  coord[0][0] = 0.58541020;
  coord[0][1] = 0.13819660;
  coord[0][2] = 0.13819660;
  coord[0][3] = 0.13819660;

  coord[1][0] = 0.13819660;
  coord[1][1] = 0.58541020;
  coord[1][2] = 0.13819660;
  coord[1][3] = 0.13819660;

  coord[2][0] = 0.13819660;
  coord[2][1] = 0.13819660;
  coord[2][2] = 0.58541020;
  coord[2][3] = 0.13819660;

  coord[3][0] = 0.13819660;
  coord[3][1] = 0.13819660;
  coord[3][2] = 0.13819660;
  coord[3][3] = 0.58541020;

  weight[0] = 0.25;
  weight[1] = 0.25;
  weight[2] = 0.25;
  weight[3] = 0.25;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    intVector idx(4);

    idx[0] = 1;
    idx[1] = 2;
    idx[2] = 3;
    idx[3] = 0;

    reorderVector(coord,idx); 
    reorderVector(weight,idx);
  }

}

/***********************************************************************/

GaussSetTria4::GaussSetTria4() {

  allocateArray(coord,4,3);
  weight = dbVector(4);

  coord[0][0] = 0.333333333333333;
  coord[0][1] = 0.333333333333333;
  coord[0][2] = 0.333333333333333;

  coord[1][0] = 0.6;
  coord[1][1] = 0.2;
  coord[1][2] = 0.2;

  coord[2][0] = 0.2;
  coord[2][1] = 0.6;
  coord[2][2] = 0.2;

  coord[3][0] = 0.2;
  coord[3][1] = 0.2;
  coord[3][2] = 0.6;

  weight[0] =  -0.5625;
  weight[1] =  0.520833333333333;
  weight[2] =  0.520833333333333;
  weight[3] =  0.520833333333333;

}

GaussSetTetra5::GaussSetTetra5() {


  allocateArray(coord,5,4);
  weight = dbVector(5);

  // default Zienkiewcz positioning
  coord[0][0] = 0.25;
  coord[0][1] = 0.25;
  coord[0][2] = 0.25;
  coord[0][3] = 0.25;

  coord[1][0] = 0.5;
  coord[1][1] = 0.166666666666667;
  coord[1][2] = 0.166666666666667;
  coord[1][3] = 0.166666666666667;

  coord[2][0] = 0.166666666666667;
  coord[2][1] = 0.5;
  coord[2][2] = 0.166666666666667;
  coord[2][3] = 0.166666666666667;

  coord[3][0] = 0.166666666666667;
  coord[3][1] = 0.166666666666667;
  coord[3][2] = 0.5;
  coord[3][3] = 0.166666666666667;

  coord[4][0] = 0.166666666666667;
  coord[4][1] = 0.166666666666667;
  coord[4][2] = 0.166666666666667;
  coord[4][3] = 0.5;

  weight[0] = -0.8;
  weight[1] = 0.45;
  weight[2] = 0.45;
  weight[3] = 0.45;
  weight[4] = 0.45;

  // GiD conform nodal numbering
  if(Gauss_Point_Positioning == 1) {

    std::cerr<<"In GaussSetTetra5::GaussSetTetra5 GiD Gauss point plotting not\n" 
	     <<"supported!"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/***********************************************************************/
GaussSetRect1::GaussSetRect1() {

  allocateArray(coord,1,2);
  weight = dbVector(1);

  coord[0][0] = 0.0;
  coord[0][1] = 0.0;
  
  weight[0] =  4;
  
}

GaussSetCube1::GaussSetCube1() {

  allocateArray(coord,1,3);
  weight = dbVector(1);
  
  coord[0][0] = 0.0;
  coord[0][1] = 0.0;
  coord[0][2] = 0.0;

  weight[0] = 8;
  
}

/***********************************************************************/
GaussSetRect4::GaussSetRect4() {

  allocateArray(coord,4,2);
  weight = dbVector(4);

  coord[0][0] = 0.577350269189626;
  coord[0][1] = 0.577350269189626;

  coord[1][0] = -0.577350269189626;
  coord[1][1] = 0.577350269189626;

  coord[2][0] = -0.577350269189626;
  coord[2][1] = -0.577350269189626;

  coord[3][0] = 0.577350269189626;
  coord[3][1] = -0.577350269189626;
  
  weight[0] =  1;
  weight[1] =  1;
  weight[2] =  1;
  weight[3] =  1;
  
}

GaussSetCube8::GaussSetCube8() {

  allocateArray(coord,8,3);
  allocateArray(weight,8);

  // default Zienkiewcz positioning
  coord[0][0] = 0.577350269189626; // 6
  coord[0][1] = 0.577350269189626;
  coord[0][2] = 0.577350269189626;

  coord[1][0] = -0.577350269189626;  // 7
  coord[1][1] = 0.577350269189626;
  coord[1][2] = 0.577350269189626;

  coord[2][0] = -0.577350269189626; // 4
  coord[2][1] = -0.577350269189626;
  coord[2][2] = 0.577350269189626;

  coord[3][0] = 0.577350269189626; // 5
  coord[3][1] = -0.577350269189626;
  coord[3][2] = 0.577350269189626;

  coord[4][0] = 0.577350269189626; //2
  coord[4][1] = 0.577350269189626;
  coord[4][2] = -0.577350269189626;

  coord[5][0] = -0.577350269189626; // 3
  coord[5][1] = 0.577350269189626;
  coord[5][2] = -0.577350269189626;

  coord[6][0] = -0.577350269189626; // 0
  coord[6][1] = -0.577350269189626;
  coord[6][2] = -0.577350269189626;

  coord[7][0] = 0.577350269189626; // 1
  coord[7][1] = -0.577350269189626;
  coord[7][2] = -0.577350269189626;

  weight[0] = 1;
  weight[1] = 1;
  weight[2] = 1;
  weight[3] = 1;
  weight[4] = 1;
  weight[5] = 1;
  weight[6] = 1;
  weight[7] = 1;

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

    reorderVector(coord,idx); 
    reorderVector(weight,idx);
  }

  
}

/***********************************************************************/
GaussSetRect9::GaussSetRect9() {
  
  allocateArray(coord,9,2);
  weight = dbVector(9);
  
  coord[0][0] = 0.774596669241483;
  coord[0][1] = 0.000000000000000;

  coord[1][0] = 0.774596669241483;
  coord[1][1] = 0.774596669241483;

  coord[2][0] = 0.000000000000000;
  coord[2][1] = 0.774596669241483;

  coord[3][0] = -0.774596669241483;
  coord[3][1] = 0.774596669241483;

  coord[4][0] = -0.774596669241483;
  coord[4][1] = 0.000000000000000;

  coord[5][0] = -0.774596669241483;
  coord[5][1] = -0.774596669241483;

  coord[6][0] = 0.0000000000000000;
  coord[6][1] = -0.774596669241483;

  coord[7][0] = 0.774596669241483;
  coord[7][1] = -0.774596669241483;

  coord[8][0] = 0.000000000000000;
  coord[8][1] = 0.000000000000000;
  
  weight[0] = 25.0/81.0;
  weight[1] = 40.0/81.0;
  weight[2] = 25.0/81.0;
  weight[3] = 40.0/81.0;
  weight[4] = 25.0/81.0;
  weight[5] = 40.0/81.0;
  weight[6] = 25.0/81.0;
  weight[7] = 40.0/81.0;	
  weight[8] = 64.0/81.0;

}

GaussSetCube27::GaussSetCube27() {

  allocateArray(coord,27,3);
  allocateArray(weight,27);

  
  coord[0][0] = 0.774596669241483; // 17
  coord[0][1] = 0.000000000000000;
  coord[0][2] = 0.774596669241483;

  coord[1][0] = 0.774596669241483; // 6
  coord[1][1] = 0.774596669241483;
  coord[1][2] = 0.774596669241483;

  coord[2][0] = 0.000000000000000;
  coord[2][1] = 0.774596669241483; // 18
  coord[2][2] = 0.774596669241483;

  coord[3][0] = -0.774596669241483; // 7
  coord[3][1] = 0.774596669241483;
  coord[3][2] = 0.774596669241483;

  coord[4][0] = -0.774596669241483;
  coord[4][1] = 0.000000000000000; // 19
  coord[4][2] = 0.774596669241483;

  coord[5][0] = -0.774596669241483; // 4
  coord[5][1] = -0.774596669241483;
  coord[5][2] = 0.774596669241483;

  coord[6][0] = 0.0000000000000000;
  coord[6][1] = -0.774596669241483; // 16
  coord[6][2] = 0.774596669241483;

  coord[7][0] = 0.774596669241483; // 5
  coord[7][1] = -0.774596669241483;
  coord[7][2] = 0.774596669241483;

  coord[8][0] = 0.000000000000000;
  coord[8][1] = 0.000000000000000; // 25
  coord[8][2] = 0.774596669241483;


  coord[9][0] = 0.774596669241483;
  coord[9][1] = 0.000000000000000; // 22
  coord[9][2] = 0.000000000000000;

  coord[10][0] = 0.774596669241483;
  coord[10][1] = 0.774596669241483; // 14
  coord[10][2] = 0.000000000000000;

  coord[11][0] = 0.000000000000000;
  coord[11][1] = 0.774596669241483; // 23
  coord[11][2] = 0.000000000000000;

  coord[12][0] = -0.774596669241483;
  coord[12][1] = 0.774596669241483; // 15
  coord[12][2] = 0.000000000000000;

  coord[13][0] = -0.774596669241483;
  coord[13][1] = 0.000000000000000; // 24
  coord[13][2] = 0.000000000000000;

  coord[14][0] = -0.774596669241483; // 12
  coord[14][1] = -0.774596669241483;
  coord[14][2] = 0.000000000000000;

  coord[15][0] = 0.0000000000000000;
  coord[15][1] = -0.774596669241483; // 21
  coord[15][2] = 0.000000000000000;

  coord[16][0] = 0.774596669241483;
  coord[16][1] = -0.774596669241483; // 13
  coord[16][2] = 0.000000000000000;

  coord[17][0] = 0.000000000000000;
  coord[17][1] = 0.000000000000000; // 26
  coord[17][2] = 0.000000000000000;


  coord[18][0] = 0.774596669241483;
  coord[18][1] = 0.000000000000000; // 9
  coord[18][2] = -0.774596669241483;

  coord[19][0] = 0.774596669241483; // 2
  coord[19][1] = 0.774596669241483;
  coord[19][2] = -0.774596669241483;

  coord[20][0] = 0.000000000000000;
  coord[20][1] = 0.774596669241483; // 10
  coord[20][2] = -0.774596669241483;

  coord[21][0] = -0.774596669241483;
  coord[21][1] = 0.774596669241483; // 3
  coord[21][2] = -0.774596669241483;

  coord[22][0] = -0.774596669241483;
  coord[22][1] = 0.000000000000000; // 11
  coord[22][2] = -0.774596669241483;

  coord[23][0] = -0.774596669241483; // 0
  coord[23][1] = -0.774596669241483;
  coord[23][2] = -0.774596669241483;

  coord[24][0] = 0.0000000000000000; // 8
  coord[24][1] = -0.774596669241483;
  coord[24][2] = -0.774596669241483;

  coord[25][0] = 0.774596669241483; // 1
  coord[25][1] = -0.774596669241483;
  coord[25][2] = -0.774596669241483;

  coord[26][0] = 0.000000000000000; 
  coord[26][1] = 0.000000000000000; // 20
  coord[26][2] = -0.774596669241483;

  
  weight[0] = 200.0/729.0;
  weight[1] = 125.0/729.0;
  weight[2] = 200.0/729.0;
  weight[3] = 125.0/729.0;
  weight[4] = 200.0/729.0;
  weight[5] = 125.0/729.0;
  weight[6] = 200.0/729.0;
  weight[7] = 125.0/729.0;
  weight[8] = 320.0/729.0;
  weight[9] = 320.0/729.0;
  weight[10] = 200.0/729.0;
  weight[11] = 320.0/729.0;
  weight[12] = 200.0/729.0;
  weight[13] = 320.0/729.0;
  weight[14] = 200.0/729.0;
  weight[15] = 320.0/729.0;
  weight[16] = 200.0/729.0;
  weight[17] = 512.0/729.0;
  weight[18] = 200.0/729.0;
  weight[19] = 125.0/729.0;
  weight[20] = 200.0/729.0;
  weight[21] = 125.0/729.0;
  weight[22] = 200.0/729.0;
  weight[23] = 125.0/729.0;
  weight[24] = 200.0/729.0;
  weight[25] = 125.0/729.0;
  weight[26] = 320.0/729.0;


  if(Gauss_Point_Positioning == 1) {

    intVector idx(27);

    // idx[new] = old(Zienkiewicz)
    idx[0] = 17;
    idx[1] = 6;
    idx[2] = 18;
    idx[3] = 7;
    idx[4] = 19;
    idx[5] = 4;
    idx[6] = 16;
    idx[7] = 5;  
    idx[8] = 25;
    idx[9] = 22;
    idx[10] = 14;
    idx[11] = 23;
    idx[12] = 15;
    idx[13] = 24;
    idx[14] = 12;
    idx[15] = 21;
    idx[16] = 13;
    idx[17] = 26;
    idx[18] = 9;
    idx[19] = 2;
    idx[20] = 10;
    idx[21] = 3;
    idx[22] = 11;
    idx[23] = 0;
    idx[24] = 8;
    idx[25] = 1;
    idx[26] = 20;

    reorderVector(coord,idx); 
    reorderVector(weight,idx);
  }

}

/***********************************************************************/
GaussSetRect16::GaussSetRect16() {
  
  allocateArray(coord,16,2);
  weight = dbVector(16);

  dbVector iCoords(4);
  dbVector iWeights(4);

  iCoords[0] = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0); 	
  iCoords[1] = -sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0); 	
  iCoords[2] = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0); 
  iCoords[3] = -sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0); 	

  iWeights[0] = (18.0+sqrt(30.0))/36.0;
  iWeights[1] = (18.0+sqrt(30.0))/36.0;
  iWeights[2] = (18.0-sqrt(30.0))/36.0;
  iWeights[3] = (18.0-sqrt(30.0))/36.0;

  for(int i=0;i<4;i++) {

    for(int j=0;j<4;j++) {

      weight[i*4+j] = iWeights[i]*iWeights[j];

      coord[i*4+j][0] = iCoords[i];
      coord[i*4+j][1] = iCoords[j];
    }

  }

}

GaussSetCube64::GaussSetCube64() {

  allocateArray(coord,64,3);
  weight = dbVector(64);

  dbVector iCoords(4);
  dbVector iWeights(4);

  iCoords[0] = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0); 	
  iCoords[1] = -sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0); 	
  iCoords[2] = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0); 
  iCoords[3] = -sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0); 	

  iWeights[0] = (18.0+sqrt(30.0))/36.0;
  iWeights[1] = (18.0+sqrt(30.0))/36.0;
  iWeights[2] = (18.0-sqrt(30.0))/36.0;
  iWeights[3] = (18.0-sqrt(30.0))/36.0;


  for(int i=0;i<4;i++) {

    for(int j=0;j<4;j++) {

      for(int k=0;k<4;k++) {

	weight[((i*4)+j)*4+k] = iWeights[i]*iWeights[j]*iWeights[k];

	coord[((i*4)+j)*4+k][0] = iCoords[i];
	coord[((i*4)+j)*4+k][1] = iCoords[j];
	coord[((i*4)+j)*4+k][2] = iCoords[k];
      }

    }

  }

}

/***********************************************************************/
GaussSetRect25::GaussSetRect25() {
  
  allocateArray(coord,25,2);
  weight = dbVector(25);

  dbVector iCoords(5);
  dbVector iWeights(5);

  iCoords[0] = 0; 
  iCoords[1] = 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)); 
  iCoords[2] = -1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)); 
  iCoords[3] = 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)); 
  iCoords[4] = -1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)); 
  
  iWeights[0] = 128.0/225.0;
  iWeights[1] = (322.0+13.0*sqrt(70.0))/900.0;
  iWeights[2] = (322.0+13.0*sqrt(70.0))/900.0;
  iWeights[3] = (322.0-13.0*sqrt(70.0))/900.0;
  iWeights[4] = (322.0-13.0*sqrt(70.0))/900.0;

  for(int i=0;i<5;i++) {

    for(int j=0;j<5;j++) {

      weight[i*5+j] = iWeights[i]*iWeights[j];

      coord[i*5+j][0] = iCoords[i];
      coord[i*5+j][1] = iCoords[j];
    }

  }

}

GaussSetCube125::GaussSetCube125() {
  
  allocateArray(coord,125,2);
  weight = dbVector(125);

  dbVector iCoords(5);
  dbVector iWeights(5);

  iCoords[0] = 0; 
  iCoords[1] = 1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)); 
  iCoords[2] = -1.0/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0)); 
  iCoords[3] = 1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)); 
  iCoords[4] = -1.0/3.0*sqrt(5.0+2.0*sqrt(10.0/7.0)); 
  
  iWeights[0] = 128.0/225.0;
  iWeights[1] = (322.0+13.0*sqrt(70.0))/900.0;
  iWeights[2] = (322.0+13.0*sqrt(70.0))/900.0;
  iWeights[3] = (322.0-13.0*sqrt(70.0))/900.0;
  iWeights[4] = (322.0-13.0*sqrt(70.0))/900.0;

  for(int i=0;i<5;i++) {

    for(int j=0;j<5;j++) {

      for(int k=0;k<5;k++) {

	weight[((i*5)+j)*5+k] = iWeights[i]*iWeights[j]*iWeights[k];

	coord[((i*5)+j)*5+k][0] = iCoords[i];
	coord[((i*5)+j)*5+k][1] = iCoords[j];
	coord[((i*5)+j)*5+k][2] = iCoords[k];
      }

    }

  }

}
