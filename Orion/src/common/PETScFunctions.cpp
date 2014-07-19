#include "PETScFunctions.h"



/************************************************************************/
// Create a dense sequential PETSc coefficient matrix.
void createDenseSequentialPETScMat(int rowNum,int colNum,
				   std::map<std::string,bool>& mOptions,
				   Mat& mat,std::ofstream& logFile) {

  using namespace std;

  MatCreate(PETSC_COMM_SELF,&mat);
  MatSetSizes(mat,rowNum,colNum,rowNum,colNum);
  MatSetType(mat,MATSEQDENSE);

  //   for(map<string,bool>::iterator p = mOptions.begin();
  //       p!=mOptions.end();++p) {

  //     if(p->first == "MAT_SYMMETRIC")
  //       MatSetOption(mat,MAT_SYMMETRIC,PETSC_TRUE);

  // #ifdef  _calculationDebugMode_
  //     logFile<<"createSparseParallelPETScMat"<<endl;
  //     logFile<<"Kmat option: "<<p->first<<" = "<<p->second<<endl;
  // #endif

  //   }


  //MatSeqDenseSetPreallocation(mat,PETSC_NULL);
  MatSetUp(mat);

#ifdef  _KmatAssemblingDebugMode_
  logFile<<"****************** matrix allocation *****************"<<endl;
  logFile<<"rowNum="<<rowNum<<endl;
  logFile<<"colNum="<<colNum<<endl;
#endif


}

/************************************************************************/
// Create a dense parallel PETSc coefficient matrix.
void createDenseParallelPETScMat(int localRowNum,int localColNum,
				 int globalRowNum,int globalColNum,
				 std::map<std::string,bool>& mOptions,
				 Mat& mat,std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  MatCreate(PETSC_COMM_WORLD,&mat);
  MatSetSizes(mat,localRowNum,localColNum,globalRowNum,
	      globalColNum);
  MatSetType(mat,MATMPIDENSE);

  //   for(map<string,bool>::iterator p = mOptions.begin();
  //       p!=mOptions.end();++p) {

  //     if(p->first == "MAT_SYMMETRIC")
  //       MatSetOption(mat,MAT_SYMMETRIC,PETSC_TRUE);

  // #ifdef  _calculationDebugMode_
  //     logFile<<"createSparseParallelPETScMat"<<endl;
  //     logFile<<"Kmat option: "<<p->first<<" = "<<p->second<<endl;
  // #endif

  //   }

  //MatMPIDenseSetPreallocation(mat,PETSC_NULL);
  MatSetUp(mat);

#ifdef  _KmatAssemblingDebugMode_
  logFile<<"****************** matrix allocation *****************"<<endl;
  logFile<<"localRowNum="<<localRowNum<<endl;
  logFile<<"localColNum="<<localColNum<<endl;
  logFile<<"globalRowNum="<<globalRowNum<<endl;
  logFile<<"globalColNum="<<globalColNum<<endl;
#endif


}

/************************************************************************/
// Create a sparse parallel PETSc coefficient matrix.
void createSparseParallelPETScMat(int localRowNum,int localColNum,
				  int globalRowNum,int globalColNum,
				  int maxDiagRowNonzeros,
				  intVector& diagEntries,
				  int maxOffDiagRowNonzeros,
				  intVector& offDiagEntries,
				  std::map<std::string,bool>& mOptions,
				  Mat& mat,std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);


  MatCreate(PETSC_COMM_WORLD,&mat);
  MatSetSizes(mat,localRowNum,localColNum,globalRowNum,
	      globalColNum);
  MatSetType(mat,MATAIJ);

  // nonzeros for each local row are specified
  if(diagEntries.size() == localRowNum &&
     offDiagEntries.size() == localRowNum) {

    MatSeqAIJSetPreallocation(mat,0,&diagEntries[0]);
    MatMPIAIJSetPreallocation(mat,0,&diagEntries[0],0,&offDiagEntries[0]);
  }

  // default number of nonzeros applicable for all rows specified
  else if(maxDiagRowNonzeros > 0 && maxOffDiagRowNonzeros >= 0) {

    MatSeqAIJSetPreallocation(mat,maxDiagRowNonzeros,PETSC_NULL);
    MatMPIAIJSetPreallocation(mat,maxDiagRowNonzeros,PETSC_NULL,
 			      maxOffDiagRowNonzeros,PETSC_NULL);
  }

  // PETSc decides of a default allocation
  else
    MatSetUp(mat);


  for(map<string,bool>::iterator p = mOptions.begin();
      p!=mOptions.end();++p) {

    if(p->first == "MAT_SYMMETRIC" && (bool)p->second)
      MatSetOption(mat,MAT_SYMMETRIC,PETSC_TRUE);

    else if(p->first == "MAT_SYMMETRIC" && !(bool)p->second)
      MatSetOption(mat,MAT_SYMMETRIC,PETSC_FALSE);

    logFile<<"Kmat option: "<<p->first<<" = "<<p->second<<endl;

  }


}

/************************************************************************/
// Create a parallel PETSc vector.
void createParallelPETScVec(int localRowNum,int globalRowNum,Vec& vec) {


  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  VecCreate(PETSC_COMM_WORLD,&vec);

  if(size == 1)
    VecSetType(vec,VECSEQ);
  else
    VecSetType(vec,VECMPI);

  VecSetSizes(vec,localRowNum,globalRowNum);

}



/************************************************************************/
// Create a parallel PETSc solver object.
void createParallelPETScSolver(Mat& mat,
			       std::map<std::string,double>& problemData,
			       std::vector<std::string> solvingMethod,
			       KSP& ksp,PC& pc,
			       std::map<std::string,double>& calcData,
			       std::map<std::string,double>& modelData,
			       std::ofstream& logFile) {


  using namespace std;

  int precond = (int)problemData["preconditioner"];
  int solver = (int)problemData["solver"];
  double absTol = problemData["absoluteResidualTolerance"];
  double relDecTol = problemData["relativeResidualDecreaseTolerance"];
  double relIncTol = problemData["relativeResidualIncreaseTolerance"];
  int numOfMaxIts = (int)problemData["numOfMaxIterations"];
  int directCalc = (int)problemData["directSolving"];

  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);


  if((int)problemData["initializeKSP"] == 1) {

//     // already done in zeroiseSolverEQSet
//     MatSetOption(mat,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
//     MatSetOption(mat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);

    MatInfo info;
    MatGetInfo(mat,MAT_LOCAL,&info);

    if(rank == 0) {
      logFile<<"------------------------"<<endl;
      logFile<<"PETSc matrix allocation:"<<endl;
      logFile<<"allocated nonzeros="<<info.nz_allocated<<endl;
      logFile<<"number of mallocs="<<info.mallocs<<endl;
      logFile<<"used nonzeros="<<info.nz_used<<endl;
      logFile<<"not needed nonzeros="<<info.nz_unneeded<<endl;
    }

    KSPCreate(PETSC_COMM_WORLD,&ksp);
  }


  // needed to compute matrix conditioning
  if((bool)problemData["deallocateParallelKmat"] &&
     calcData["calculationStep"] == 1 && directCalc != 1)

    KSPSetComputeSingularValues(ksp,PETSC_TRUE);


  /*********************************************************************/

  if(solvingMethod.size() < 2)
    solvingMethod.resize(2);


  int m,n;
  int localBlocks,blocks;
  int mglevels;


  PetscInt nlocal,first;
  PC subpc;
  KSP* subksp;
  Mat F;



  // Select preconditioning method for parallel solving using Krylov
  // subspace methods etc.

  if(size == 1)
    directCalc = 1;

  if(directCalc != 1) {

    KSPGetPC(ksp,&pc);

    //     if(size == 1) {
    //       PCSetType(pc,"lu");
    //       solvingMethod[0] = "lu-factorizing";
    //     }
    //     else {
    //       logFile<<"LU-factorization is only for sequential calculation "
    // 	     <<"available!"<<endl;
    //       MPI_Abort(PETSC_COMM_WORLD,1);
    //     }

    switch(precond) {

    case 0:
      PCSetType(pc,"none");
      solvingMethod[0] = "none";
      break;

    case 1:

      if(size == 1) {
	PCSetType(pc,"lu");
	solvingMethod[0] = "lu-factorizing";
      }
      else {
	logFile<<"LU-factorization is only for sequential calculation "
	       <<"available!"<<endl;
	MPI_Abort(PETSC_COMM_WORLD,1);
      }

      break;

    case 2:

      if(size == 1) {
	PCSetType(pc,"ilu");
	solvingMethod[0] = "Incomplete LU";
      }
      else {
	PCSetType(pc,"asm");

	//       MatGetLocalSize(mat,&m,&n);
	//       localBlocks = (int)(m/problemData["PCSubMatSizeLimit"]);
	//       MPI_Allreduce(&localBlocks,&blocks,1,MPI_INT,MPI_MAX,
	// 		    PETSC_COMM_WORLD);

	//       if(blocks > 0)
	// 	PCASMSetOverlap(pc,blocks);
	//       else
	// 	blocks = 1;

	solvingMethod[0] = "Additive Schwarz with with local Incomplete LU";
      }

      break;

      // Cholesky factorization double as efficient than LU but matrix
      // needs to be symmetric definite
    case 3:

      if((bool)problemData["tangentSymmetric"]) {

	PCSetType(pc,"cholesky");
	solvingMethod[0] = "Cholesky";
	break;
      }

      else {
	logFile<<"Cholesky-factorization supports symmetric definite "
	       <<"matrices."<<endl;
	MPI_Abort(PETSC_COMM_WORLD,1);
      }

      break;

    case 4:
      PCSetType(pc,"icc");
      solvingMethod[0] = "Incomplete Cholesky";
      break;
    case 5:
      PCSetType(pc,"jacobi");
      solvingMethod[0] = "Jacobi";
      break;
    case 6:
      PCSetType(pc,"bjacobi");
      solvingMethod[0] = "Block Jacobi";
      break;
    case 7:
      PCSetType(pc,"pbjacobi");
      solvingMethod[0] = "pbjacobi";
      break;
    case 8:
      PCSetType(pc,"asm");

      //     PCASMSetType(pc,PC_ASM_BASIC);

      //     MatGetLocalSize(mat,&m,&n);
      //     localBlocks = (int)(m/problemData["PCSubMatSizeLimit"]);
      //     MPI_Allreduce(&localBlocks,&blocks,1,MPI_INT,MPI_MAX,
      // 		  PETSC_COMM_WORLD);

      //     if(blocks > 0)
      //       PCASMSetOverlap(pc,blocks);
      //     else
      //       blocks = 1;

      //     solvingMethod[0] = "Additive Schwarz  (" +
      //       convertIntToString((int)blocks) + " sub-blocks)";

      solvingMethod[0] = "Additive Schwarz";
      break;

    case 9:
      PCSetType(pc,"gasm");

      //     MatGetLocalSize(mat,&m,&n);
      //     localBlocks = (int)(m/problemData["PCSubMatSizeLimit"]);
      //     MPI_Allreduce(&localBlocks,&blocks,1,MPI_INT,MPI_MAX,
      // 		  PETSC_COMM_WORLD);

      //     if(blocks > 0)
      //       PCGASMSetOverlap(pc,blocks);
      //     else
      //       blocks = 1;

      //     solvingMethod[0] = "Gram-Schmidt orthogonalization  (" +
      //       convertIntToString((int)blocks) + " sub-blocks)";

      solvingMethod[0] = "Gram-Schmidt orthogonalization";
      break;

    case 10:
      PCSetType(pc,"sor");
      solvingMethod[0] = "SOR (and SSOR)";
      break;
    case 11:
      PCSetType(pc,"eisenstat");
      solvingMethod[0] = "SOR with Eisenstat trick";
      break;
    case 12:
      PCSetType(pc,"galerkin");
      solvingMethod[0] = "galerkin";
    case 13:
      PCSetType(pc,"redundant");
      solvingMethod[0] = "redundant";
      break;
    case 14:
      PCSetType(pc,"composite");
      solvingMethod[0] = "combination of preconditioners";
      break;
    case 15:
      PCSetType(pc,"mg");
      solvingMethod[0] = "MG";
      break;
    case 16:
      PCSetType(pc,"ksp");
      solvingMethod[0] = "linear solver";
      break;
    case 17:
      PCSetType(pc,"spai");
      solvingMethod[0] = "spai";
      break;
    case 18:
      PCSetType(pc,"nn");
      solvingMethod[0] = "nn";
      break;
    case 19:
      PCSetType(pc,"samg");
      solvingMethod[0] = "samg";
      break;
    case 20:
      PCSetType(pc,"mat");
      solvingMethod[0] = "mat";
      break;
    case 21:
      PCSetType(pc,"hypre");
      solvingMethod[0] = "hypre";
      break;
    case 22:
      PCSetType(pc,"fieldsplit");
      solvingMethod[0] = "fieldsplit";
      break;
    case 23:
      PCSetType(pc,"tfs");
      solvingMethod[0] = "tfs";
      break;
    case 24:
      PCSetType(pc,"ml");
      solvingMethod[0] = "ml";
      break;
    case 25:
      PCSetType(pc,"prometheus");
      solvingMethod[0] = "prometheus";
      break;
    case 26:
      PCSetType(pc,"shell");
      solvingMethod[0] = "shell";
      break;
    default:
      logFile<<"In function PETScFunctions::createParallelPETScSolver\n"
	     <<"preconditioner '"<<precond<<"' is not available!"<<endl;
      MPI_Abort(PETSC_COMM_WORLD,1);
    }

    //PCLUSetUseInPlace(pc); // destroy the original matrix

    // avoid zero pivots in tangent matrix
    if(problemData["tangentZeroPivotShift"] > 0) {
      //PCFactorSetShiftType(pc,PETSC_DECIDE);
      //PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE);
      PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);
    }

  }

  // --------------------------------------------------------------------
  // select the solver method

  // preconditioning for direct solving using parallelized PETSc routines
  if(directCalc == 1 && problemData["directSolver"] == 1) {

    if(size > 1) {
      cout<<"Standard PETSc direct solver does not support parallel runs.\n"
	  <<"Change input data 'directSolver' to <2,3>"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    solvingMethod[0] = "lu-factorizing";
    solvingMethod[1] = "direct solver (PETSc)";

    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);

    //     PCFactorSetMatSolverPackage(pc,MATSOLVERPETSC);
    //     PCFactorSetUpMatSolverPackage(pc); /* call MatGetFactor() to create F */
    //     PCFactorGetMatrix(pc,&F);

    KSPSetUp(ksp);
  }

  else if(directCalc == 1 && problemData["directSolver"] == 2) {

    solvingMethod[0] = "lu-factorizing";
    solvingMethod[1] = "direct solver (MUMPS)";

    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);

    //     if((bool)problemData["deallocateParallelKmat"]) {
    //       logFile<<"MUMPS does not support yet intermediate equation system\n"
    // 	     <<"deallocation."<<endl;
    //       MPI_Abort(PETSC_COMM_WORLD,1);
    //     }

    //    if(!(bool)problemData["tangentSymmetric"]) {

    //    PCSetType(pc,PCLU);
    //    }
    //    else {

    //      PCSetType(pc,PCCHOLESKY);
    //      solvingMethod[0] = "Cholesky";
    //    }

    PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverPackage(pc);
    PCFactorGetMatrix(pc,&F);

    int icntl=7; int ival = 2;
    MatMumpsSetIcntl(F,icntl,ival);

    KSPSetUp(ksp);
  }
  // preconditioning for direct solving using parallelized SuperLU
  else if(directCalc == 1 && problemData["directSolver"] == 3) {

    if(size == 1) {
      cout<<"SuperLU_Dist does not support sequential runs.\n"
	  <<"Change input data 'directSolver' to <1,2>"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

    solvingMethod[0] = "lu-factorizing";
    solvingMethod[1] = "direct solver (SuperLU)";

    KSPSetType(ksp,KSPPREONLY);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCLU);

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);

    PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
    PCFactorSetUpMatSolverPackage(pc);
    PCFactorGetMatrix(pc,&F);
    //MatSuperluSetILUDropTol(F,1.e-8);

    KSPSetUp(ksp);
  }

  // --------------------------------------------------------------------
  // parallel solving using ILU pre-conditioning and GMRES on sub-blocks

  else if(size > 1 && precond == 2 &&
	  (bool)problemData["subBlockKSPDefinition"]) {
    int localBlockNum;
    KSP *subksp;
    PC  subpc;

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);
    KSPSetUp(ksp);

    PCASMGetSubKSP(pc,&localBlockNum,PETSC_NULL,&subksp);

    // Loop over the local blocks, setting various KSP options
    // for each block.
    for (int i=0;i<localBlockNum;i++) {
      KSPGetPC(subksp[i],&subpc);
      PCSetType(subpc,"lu");
      KSPSetType(subksp[i],"preonly");
      KSPSetTolerances(subksp[i],relDecTol,absTol,relIncTol,numOfMaxIts);

    }

    KSPSetFromOptions(ksp);

    solvingMethod[1] =
      "Generalized Minimal Residual with ILU and direct LU solves on "
      + convertIntToString((int)localBlockNum) + " local sub-blocks";
  }

  // --------------------------------------------------------------------
  // parallel solving using Additive Schwarz combined with
  // Jacobi pre-conditioning and GMRES or direct solves on sub-blocks
  else if(size > 1 && precond == 8 &&
	  (bool)problemData["subBlockKSPDefinition"]) {

    int localBlockNum,first;
    KSP *subksp;
    PC  subpc;

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);
    KSPSetUp(ksp);

    PCASMGetSubKSP(pc,&localBlockNum,&first,&subksp);

    solvingMethod[1] =
      "Generalized Minimal Residual with Additive Schwarz and ";


    for (int i=0;i<localBlockNum;i++) {
      KSPGetPC(subksp[i],&subpc);

      if((bool)problemData["subBlockSolverPreconditioner"] > 0) {
	PCSetType(subpc,"jacobi");
	solvingMethod[1] = solvingMethod[1] + "Jacobi";
	KSPSetType(subksp[i],"gmres");
      }
      else {
	PCSetType(subpc,"lu");
	KSPSetType(subksp[i],"preonly");
	solvingMethod[1] = solvingMethod[1] + "direct LU solves";
      }

      KSPSetTolerances(subksp[i],relDecTol,absTol,relIncTol,numOfMaxIts);
    }

    KSPSetFromOptions(ksp);

    solvingMethod[1] = solvingMethod[1] + " on "
      + convertIntToString((int)localBlockNum) + " local sub-blocks";
  }

  // --------------------------------------------------------------------
  // parallel solving using Gram-Schmidt orthogonalization combined with
  // Jacobi pre-conditioning and GMRES direct solves on sub-blocks
  else if(size > 1 && precond == 9 &&
	  (bool)problemData["subBlockKSPDefinition"]) {

    int localBlockNum,first;
    KSP *subksp;
    PC  subpc;

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);
    KSPSetUp(ksp);

    PCGASMGetSubKSP(pc,&localBlockNum,&first,&subksp);

    solvingMethod[1] =
      "Generalized Minimal Residual with Gram-Schmid and ";


    for (int i=0;i<localBlockNum;i++) {
      KSPGetPC(subksp[i],&subpc);

      if((bool)problemData["subBlockSolverPreconditioner"] > 0) {
	PCSetType(subpc,"jacobi");
	solvingMethod[1] = solvingMethod[1] + "Jacobi";
	KSPSetType(subksp[i],"gmres");
      }
      else {
	PCSetType(subpc,"lu");
	KSPSetType(subksp[i],"preonly");
	solvingMethod[1] = solvingMethod[1] + "direct LU solves";
      }

      KSPSetTolerances(subksp[i],relDecTol,absTol,relIncTol,numOfMaxIts);
    }

    KSPSetFromOptions(ksp);

    solvingMethod[1] = solvingMethod[1] + " on "
      + convertIntToString((int)localBlockNum) + " local sub-blocks";
  }

  // --------------------------------------------------------------------
  // other iterative solving
  else {

    // Select a iterative solver method.
    switch(solver) {
    case 0:
      KSPSetType(ksp,"gmres");
      //KSPGMRESSetRestart(ksp,(int)problemData["gmresRestartIterationNum"]);
      //KSPGMRESSetCGSRefinementType(ksp,(KSPGMRESCGSRefinementType)
      //			     KSP_GMRES_CGS_REFINEMENT_IFNEEDED);

      solvingMethod[1] = "Generalized Minimal Residual";
      break;
    case 1:
      KSPSetType(ksp,"fgmres");
      solvingMethod[1] = "fgmres";
      break;
    case 2:
      KSPSetType(ksp,"lgmres");
      solvingMethod[1] = "lgmres";
      break;
    case 3:
      if((bool)problemData["tangentSymmetric"]) {
	KSPSetType(ksp,"cg");
	solvingMethod[1] = "Conjugate Gradient";
      }
      else {
	KSPSetType(ksp,"bcgs");
	solvingMethod[1] = "Biconjugate Gradient Stabilized (BiCGSTAB)";
      }
      break;
    case 4:
      KSPSetType(ksp,"bicg");
      solvingMethod[1] = "Bijugate Gradient";
      break;
    case 5:
      KSPSetType(ksp,"bcgs");
      solvingMethod[1] = "Biconjugate Gradient Stabilized (BiCGSTAB)";
      break;
    case 6:
      KSPSetType(ksp,"cgs");
      solvingMethod[1] = "Conjugate Gradient Squared";
      break;
    case 7:
      KSPSetType(ksp,"tfqmr");
      solvingMethod[1] = "Transpose-Free Quasi-Minimal Residual (1)";
      break;
    case 8:
      KSPSetType(ksp,"tcqmr");
      solvingMethod[1] = "Transpose-Free Quasi-Minimal Residual (2)";
      break;
    case 9:
      if((bool)problemData["tangentSymmetric"]) {
	KSPSetType(ksp,"cr");
	solvingMethod[1] = "Conjugate Residual";
      }
      else {
	KSPSetType(ksp,"bcgs");
	solvingMethod[1] = "Biconjugate Gradient Stabilized (BiCGSTAB)";
      }
      break;
    case 10:
      KSPSetType(ksp,"lsqr");
      solvingMethod[1] = "Least Square Method";
      break;
    case 11:
      KSPSetType(ksp,"richardson");
      solvingMethod[1] = "Richardson";
      break;
    case 12:
      KSPSetType(ksp,"chebychev");
      solvingMethod[1] = "Chebychev";
      break;
    case 13:
      KSPSetType(ksp,"qcg");
      solvingMethod[1] = "";
      break;
    case 14:
      KSPSetType(ksp,"cgne");
      solvingMethod[1] = "cgne";
      break;
    case 15:
      KSPSetType(ksp,"stcg");
      solvingMethod[1] = "stcg";
      break;
    case 16:
      KSPSetType(ksp,"bcgs");
      solvingMethod[1] = "bcgs";
      break;
    case 17:
      KSPSetType(ksp,"bcgsl");
      solvingMethod[1] = "bcgsl";
      break;
    case 18:
      KSPSetType(ksp,"minres");
      solvingMethod[1] = "minres";
      break;
    case 19:
      KSPSetType(ksp,"symmlq");
      solvingMethod[1] = "symmlq";
      break;
    case 20:
      KSPSetType(ksp,"lcd");
      solvingMethod[1] = "lcd";
      break;
    default:
      logFile<<"In function PETScFunctions::createParallelPETScSolver\n"
	     <<"solver '"<<solver<<"' is not available!"<<endl;
      MPI_Abort(PETSC_COMM_WORLD,1);
    }


    KSPSetTolerances(ksp,relDecTol,absTol,relIncTol,numOfMaxIts);
    KSPSetFromOptions(ksp);

    KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);
    KSPSetUp(ksp);
  }


  //   if (!rank) {
  //     if (i%2) {
  //       PCSetType(subpc,PCILU);
  //     } else {
  //       PCSetType(subpc,PCNONE);
  //       KSPSetType(subksp[i],KSPBCGS);
  //       KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  //     }
  //   } else {
  //     PCSetType(subpc,PCJACOBI);
  //     KSPSetType(subksp[i],KSPGMRES);
  //     KSPSetTolerances(subksp[i],1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  //   }

  logFile<<"Preconditioner: "<<solvingMethod[0]<<endl;
  logFile<<"Solver: "<<solvingMethod[1]<<endl;
}

/************************************************************************/
// specify the reason for the iterative solver breakdown
void writeKSPConvergenceReason(KSPConvergedReason reason,
			       int solverIterations,
			       std::ofstream& logFile) {

  using namespace std;

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  if(reason < 1 && reason != -3) {

    logFile<<"solving equation system has diverged (reason="
	   <<reason<<") after "
	   <<solverIterations<<" solver iterations."<<endl;

    if(rank == 0)
      cout<<"solving equation system has diverged after "
	  <<solverIterations<<" solver iterations."<<endl;

  }
  else if(reason == -3) {

    logFile<<"solving equation system hasn't converged (reason="
	   <<reason<<") after "
	   <<solverIterations<<" solver iterations."<<endl;

    if(rank == 0)
      cout<<"solving equation system hasn't converged after "
	  <<solverIterations<<" solver iterations."<<endl;
  }

  else {

    logFile<<"solving equation system has converged after "
	   <<solverIterations<<" solver iterations."<<endl;

    if(rank == 0)
      cout<<"solving equation system has converged after "
	  <<solverIterations<<" solver iterations."<<endl;

  }

}

/************************************************************************/
// check whether PETSc equation solve was successful
bool checkPETScConvergence(KSP& ksp,int& solverIterations,
			   std::ofstream& logFile) {

  using namespace std;

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  KSPGetIterationNumber(ksp,&solverIterations);

  KSPConvergedReason globalReason;
  MPI_Allreduce(&reason,&globalReason,1,MPI_INT,MPI_MIN,
		MPI_COMM_WORLD);

  int globalIterations;
  MPI_Allreduce(&solverIterations,&globalIterations,1,MPI_INT,MPI_MAX,
		MPI_COMM_WORLD);
  solverIterations = globalIterations;

  //
  writeKSPConvergenceReason(globalReason,solverIterations,logFile);

  bool flag;

  // definitely divergenced
  if(globalReason < 1 && globalReason != -3)

    flag = false;

  // iteration limit reached - not converged yet
  else if(globalReason == -3)

    flag = true;

  else

    flag = true;


  return flag;
}

/************************************************************************/
// Assemble a vector containing all global entries of a parallel PETSc
// vector.
void getGlobalParallelPETScVec(int localSize,int globalSize,Vec& vec,
			       dbVector& globalVec) {

  using namespace std;

  if(globalVec.size() < globalSize)
    globalVec.resize(globalSize);

  clearArray(globalVec);

  double* localVec;
  VecGetArray(vec,&localVec);

  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  intVector recvCounts(size);
  intVector displs(size);

  MPI_Allgather(&localSize,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		PETSC_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];

  MPI_Allgatherv(localVec,recvCounts[rank],MPI_DOUBLE,&globalVec[0],
		 &recvCounts[0],&displs[0],MPI_DOUBLE,PETSC_COMM_WORLD);

  VecRestoreArray(vec,&localVec);
  //delete[] localVec;
}

/************************************************************************/
void destroyPETScSolver(KSP& ksp) {


  KSPDestroy(&ksp);
}

/************************************************************************/
void destroyPETScMat(Mat& mat) {

  //PetscBool flag;
  //MatValid(mat,&flag);
  //if(flag == PETSC_TRUE)
  MatDestroy(&mat);

}

/************************************************************************/
void destroyPETScVec(Vec& vec) {

  //PetscBool flag;
  //VecValid(vec,&flag);
  //if(flag == PETSC_TRUE)
  VecDestroy(&vec);

}

/************************************************************************/
// compute the conditioning number of a sparse parallel matrix
void checkMatrixConditioning(Mat& mat,
			     std::map<std::string,double>& problemData,
			     std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  double absTol = problemData["absoluteResidualTolerance"];
  double relDecTol = problemData["relativeResidualDecreaseTolerance"];
  double relIncTol = problemData["relativeResidualIncreaseTolerance"];
  int numOfMaxIts = (int)problemData["numOfMaxIterations"];

  int m,n,M,N;
  double emax,emin;

  MatGetLocalSize(mat,&m,&n);
  MatGetSize(mat,&M,&N);

  KSP ksp;
  PC pc;

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPGetPC(ksp,&pc);

  PCSetType(pc,"none");
  //KSPSetType(ksp,"cg");
  KSPSetType(ksp,"gmres");
  //KSPGMRESSetRestart(ksp,200);
  KSPSetTolerances(ksp,relDecTol,absTol,relIncTol,numOfMaxIts);

  KSPSetComputeSingularValues(ksp,PETSC_TRUE);
  KSPSetFromOptions(ksp);

  KSPSetOperators(ksp,mat,mat,DIFFERENT_NONZERO_PATTERN);
  KSPSetUp(ksp);

  // ---------------------------------------------------------------------

  intVector dummyVecIdx(n);
  dbVector dummyVec(n);
  intVector recvCounts(size);
  intVector displs(size);

  Vec X,Y;
  createParallelPETScVec(n,N,X);
  createParallelPETScVec(n,N,Y);

  MPI_Allgather(&n,1,MPI_INT,&recvCounts[0],1,MPI_INT,
		PETSC_COMM_WORLD);

  displs[0] = 0;

  for(int i=1;i<size;i++)
    displs[i] = displs[i-1]+recvCounts[i-1];


  for(int i=displs[rank],j=0;j<n;i++,j++) {
    dummyVec[j] = 1.0;
    dummyVecIdx[j] = i;
  }

  VecSetValues(Y,n,&dummyVecIdx[0],&dummyVec[0],INSERT_VALUES);
  VecAssemblyBegin(Y);
  VecAssemblyEnd(Y);

  // ---------------------------------------------------------------------

  int solverIterations;
  KSPConvergedReason reason;

  KSPSolve(ksp,Y,X);
  KSPGetIterationNumber(ksp,&solverIterations);
  KSPGetConvergedReason(ksp,&reason);

#ifdef _commonDebugMode_
  writeKSPConvergenceReason(reason,solverIterations,logFile);
#endif

  KSPComputeExtremeSingularValues(ksp,&emax,&emin);

  long double condNum = emax/emin;

  int digits = (int)log10(condNum);

  int accuracy = 16 - digits;

  logFile<<"******************************************************"<<endl;
  logFile<<"********* NOT preconditioned matrix ******************"<<endl;
  logFile<<"max eigenvalue = "<<emax<<endl;
  logFile<<"min eigenvalue = "<<emin<<endl;
  logFile<<"matrix conditioning number: "<<condNum<<endl;
  logFile<<"Solving of linear equation system accurate to "<<accuracy
 	 <<" decimal digits."<<endl;
  if(rank == 0) {
    cout<<"un-conditioned matrix conditioning number: "<<condNum<<endl;
    cout<<"Solving of linear equation system accurate to "<<accuracy
	<<" decimal digits."<<endl;
  }

  destroyPETScSolver(ksp);
  destroyPETScVec(X);
  destroyPETScVec(Y);

}

/************************************************************************/
// compute the conditioning number of a sparse parallel matrix
void checkMatrixConditioning(KSP& ksp,
			     std::map<std::string,double>& problemData,
			     std::ofstream& logFile) {

  using namespace std;

  int rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  double emax,emin;
  KSPComputeExtremeSingularValues(ksp,&emax,&emin);

  long double condNum = emax/emin;

  int digits = (int)log10(condNum);

  int accuracy = 16 - digits;

  logFile<<"******************************************************"<<endl;
  logFile<<"************* preconditioned matrix ******************"<<endl;
  logFile<<"max eigenvalue = "<<emax<<endl;
  logFile<<"min eigenvalue = "<<emin<<endl;
  logFile<<"matrix conditioning number: "<<condNum<<endl;
  logFile<<"Solving of linear equation system accurate to "<<accuracy
 	 <<" decimal digits."<<endl;
  if(rank == 0) {
    cout<<"pre-conditioned matrix conditioning number: "<<condNum<<endl;
    cout<<"Solving of linear equation system accurate to "<<accuracy
	<<" decimal digits."<<endl;
  }

}


