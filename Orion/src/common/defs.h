/***********************************************************************/
/***********************************************************************/
// Definitions used to provide debugging informations
// add comment
// Common define statements
// Add comment symbol "//" in-front of #undef line to activate debugging algorithm

#ifndef _commonDebugMode_
#define _commonDebugMode_
#endif 
#undef _commonDebugMode_

#ifndef _commonDebugMode_SVD_
#define _commonDebugMode_SVD_
#endif
#undef _commonDebugMode_SVD_

// ----------------------------------------------------------------------------
// ROM Calculation
#ifndef _ROMCalcDebugMode_
#define _ROMCalcDebugMode_
#endif
#undef _ROMCalcDebugMode_

// ----------------------------------------------------------------------------
// Database debug statements
#ifndef _DatabaseDebugMode_
#define _DatabaseDebugMode_
#endif
#undef _DatabaseDebugMode_

#ifndef _DataDebugMode_
#define _DataDebugMode_
#endif
#undef _DataDebugMode_

// ----------------------------------------------------------------------------
// GridNodes
#ifndef _GridNodesDebugMode_
#define _GridNodesDebugMode_
#endif
#undef _GridNodesDebugMode_

#ifndef _NodesDebugMode_
#define _NodesDebugMode_
#endif
#undef _NodesDebugMode_

// ----------------------------------------------------------------------------
// POD calculation
#ifndef _PODCalcDebugMode_
#define _PODCalcDebugMode_
#endif
#undef _PODCalcDebugMode_

// POD calculation checks
#ifndef _PODCalcCheckMode_
#define _PODCalcCheckMode_
#endif
#undef _PODCalcCheckMode_

// ----------------------------------------------------------------------------
// PODI calculation
#ifndef _PODICalcDebugMode_
#define _PODICalcDebugMode_
#endif
#undef _PODICalcDebugMode_

// PODI calculation checks
#ifndef _PODICalcCheckMode_
#define _PODICalcCheckMode_
#endif
#undef _PODICalcCheckMode_

// ----------------------------------------------------------------------------
// MLS debug mode
#ifndef _MLSDebugMode_
#define _MLSDebugMode_
#endif
#undef _MLSDebugMode_

// ----------------------------------------------------------------------------
// Activate checking algorithms
#ifndef _checkMode_
#define _checkMode_
#endif
#undef _checkMode_

// ----------------------------------------------------------------------------
// Geometry define statements
#ifndef _FEdebugMode_
#define _FEdebugMode_
#endif
#undef _FEdebugMode_

#ifndef _FEdebugMode_PointInPoly_
#define _FEdebugMode_PointInPoly_
#endif
#undef _FEdebugMode_PointInPoly_

// ****************************************************************************
// ****************************************************************************
// Seska stuff

#ifndef _geometryDebugMode_
#define _geometryDebugMode_
#endif
#undef _geometryDebugMode_

#ifndef _anisotropyInterpolationDebugMode_
#define _anisotropyInterpolationDebugMode_
#endif
#undef _anisotropyInterpolationDebugMode_

#ifndef _influenceRadiusDebugMode_
#define _influenceRadiusDebugMode_
#endif
#undef _influenceRadiusDebugMode_

#ifndef _supportDebugMode_
#define _supportDebugMode_
#endif
#undef _supportDebugMode_

//#ifndef _influenceSphereDebugMode_
//#define _influenceSphereDebugMode_
//#endif
//#undef _influenceSphereDebugMode_

//#ifndef _partitioningDebugMode_
//#define _partitioningDebugMode_
//#endif
//#undef _partitioningDebugMode_

// kinematics define statements
//#ifndef _kinematicsDebugMode_
//#define _kinematicsDebugMode_
//#endif
//#undef _kinematicsDebugMode_

// multiplicative rotation update
//#ifndef _rotationUpdateDebugMode_
//#define _rotationUpdateDebugMode_
//#endif
//#undef _rotationUpdateDebugMode_

// constitutive laws define statements
//#ifndef _constitutiveDebugMode_
//#define _constitutiveDebugMode_
//#endif
//#undef _constitutiveDebugMode_

// Assembling discrete equation system define statements
//#ifndef _modelDebugMode_
//#define _modelDebugMode_
//#endif
//#undef _modelDebugMode_

//#ifndef _KmatModificationDebugMode_
//#define _KmatModificationDebugMode_
//#endif
//#undef _KmatModificationDebugMode_

//#ifndef _forceVecModificationDebugMode_
//#define _forceVecModificationDebugMode_
//#endif
//#undef _forceVecModificationDebugMode_

//#ifndef _KmatAssemblingDebugMode_
//#define _KmatAssemblingDebugMode_
//#endif
//#undef _KmatAssemblingDebugMode_

//#ifndef _forceVecAssemblingDebugMode_
//#define _forceVecAssemblingDebugMode_
//#endif
//#undef _forceVecAssemblingDebugMode_

//#ifdef _forceVecModificationDebugMode_
//#define _KmatAssemblingDebugMode_
//#define _forceVecAssemblingDebugMode_
//#elif defined _KmatModificationDebugMode_
//#define _KmatAssemblingDebugMode_
//#define _forceVecAssemblingDebugMode_
//#endif

//#ifndef _checkParallelKmatDebugMode_
//#define _checkParallelKmatDebugMode_
//#endif
//#undef _checkParallelKmatDebugMode_

//#ifndef _checkEQSystemDebugMode_
//#define _checkEQSystemDebugMode_
//#endif
//#undef _checkEQSystemDebugMode_

// Calculation define statements
//#ifndef _calculationDebugMode_
//#define _calculationDebugMode_
//#endif
//#undef _calculationDebugMode_

// Postprocessing define statements
//#ifndef _postProcDebugMode_
//#define _postProcDebugMode_
//#endif
//#undef _postProcDebugMode_

// Shape functions debug mode
//#ifndef _shapeFunctionsDebugging_
//#define _shapeFunctionsDebugging_
//#endif
//#undef _shapeFunctionsDebugging_

//#ifndef _MLSDebugMode_
//#define _MLSDebugMode_
//#endif
//#undef _MLSDebugMode_

//#ifndef _testBoundCollocation_
//#define _testBoundCollocation_
//#endif
//#undef _testBoundCollocation_

//#ifndef _debugMode_
//#define _debugMode_
//#endif
//#undef _debugMode_

//#ifndef _calculateSeries_
//#define _calculateSeries_
//#endif
//#undef _calculateSeries_

//#ifndef _checkMemory_
//#define _checkMemory_
//#endif
//#undef _checkMemory_

//#ifndef _checkMatrixConditioning_
//#define _checkMatrixConditioning_
//#endif
//#undef _checkMatrixConditioning_

//#ifndef _cuda_
//#define _cuda_
//#endif
//#undef _cuda_


//#ifndef _cudaDebug_
//#define _cudaDebug_
//#endif
//#undef _cudaDebug_

//inverse calcs
//#ifndef _inverseDebugMode_
//#define _inverseDebugMode_
//#endif
//#undef _inverseDebugMode_

//#ifndef _inverseDisplayDebugMode_
//#define _inverseDisplayDebugMode_
//#endif
//#undef _inverseDisplayDebugMode_

// determining the unloaded configuration
//#ifndef _loadedConfigDebugMode_
//#define _loadedConfigDebugMode_
//#endif
//#undef _loadedConfigDebugMode_
 
