#ifndef commonTypedefs_h_
#define commonTypedefs_h_

#include <valarray>
#include <vector>

typedef std::vector<bool> blVector;
typedef std::vector<int> intVector;
typedef std::vector<long int> longIntVector;
typedef std::vector<long unsigned> longUnsignVector;
typedef std::vector<double> dbVector;
typedef std::vector<std::vector<bool> > blMatrix;
typedef std::vector<std::vector<int> > intMatrix;
typedef std::vector<std::vector<long unsigned> > longUnsignMatrix;
typedef std::vector<std::vector<double> > dbMatrix;
typedef std::vector<std::vector<int> > intMatrix;
typedef std::vector<std::vector<std::vector<int> > > intMatrix3;
typedef std::vector<std::vector<std::vector<double> > > dbMatrix3;
typedef std::vector<std::vector<std::vector<std::vector<int> > > > intMatrix4;
typedef std::vector<std::vector<std::vector<std::vector<double> > > > dbMatrix4;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > dbMatrix5;


//typedef std::valarray<double> dbVector;
//typedef std::valarray<std::valarray<double> > dbMatrix;
typedef std::valarray<long double> ldbVector;
typedef std::valarray<std::valarray<long double> > ldbMatrix;


// =============================================================================
// Orion Precision Type definitions
// (Allows to switch between single and double precision)

typedef double oPType;
//typedef float oPType;

typedef std::vector<oPType> oPTVector;
typedef std::vector<std::vector<oPType> > oPTMatrix;
typedef std::vector<std::vector<std::vector<oPType> > > oPTMatrix3;
typedef std::vector<std::vector<std::vector<std::vector<oPType> > > > oPTMatrix4;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<oPType> > > > > oPTMatrix5;

#endif
