// Calculate the prismatic window function and its derivations at a given
// point.

#ifndef PrismaticWindowFunction_h_
#define PrismaticWindowFunction_h_

#include <fstream>
#include <iostream>
#include <vector>

#include "commonTypedefs.h"
#include "mpi.h"
#include "Particle.h"
#include "petscsys.h"

namespace CubicSplineWindowFunc {
      
  // The window function.
  double W1(double& coord);
  double W2(double& coord);

  // First order derivation of the cubic spline.
  double dW1(double& coord);
  double dW2(double& coord);

  // Second order derivation of the cubic spline.
  double d2W1(double& rad);
  double d2W2(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// Quartic spline window function.
namespace QuarticSplineWindowFunc {
      
  // The window function.
  double W(double& coord);

  // First order derivation of the quartic spline.
  double dW(double& coord);

  // Second order derivation of the quartic spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// Tenth order spline window function.
namespace TenthOrderSplineWinFunc {
      
  // The window function.
  double W(double& coord);

  // First order derivation of the quartic spline.
  double dW(double& coord);

  // Second order derivation of the quartic spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// Exponential window function.
namespace ExponentialWindowFunc {
      
  // The window function.
  double W(double& coord);

  // First order derivation of the exponential spline.
  double dW(double& coord);

  // Second order derivation of the exponential spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// Gauss spline window function.
namespace GaussWindowFunc {
      
  // The window function.
  double W(double& coord);

  // First order derivation of the exponential spline.
  double dW(double& coord);

  // Second order derivation of the exponential spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// test Gauss spline window function (C^\infty).
namespace GaussWindowFunc2 {
      
  // The window function.
  double W(double& r,double& radius,double& c,double& k);

  // First order derivation of the exponential spline.
  double dW(double& r,double& radius,double& c,double& k);

  // Second order derivation of the exponential spline.
  double d2W(double& r,double& radius,double& c,double& k);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& c,
		   Particle& ptcle,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,
			    double& c,Particle& ptcle,
			    double& xDWFunc,double& yDWFunc,
			    double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& c,Particle& ptcle,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// C^3 spline consisting of 6 quartic polynomials
namespace QuarticSplineWindowFunc2 {

  // The window function.
  double W(double& coord);

  // First order derivation of the cubic spline.
  double dW(double& coord);

  // Second order derivation of the cubic spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// C^4 spline consisting of 5 quintic polynomials
namespace QuinticSplineWindowFunc {

  // The window function.
  double W(double& coord);

  // First order derivation of the cubic spline.
  double dW(double& coord);

  // Second order derivation of the cubic spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// Fifth order spline window function.
namespace FifthOrderSplineWinFunc {
      
  // The window function.
  double W(double& coord);

  // First order derivation of the quartic spline.
  double dW(double& coord);

  // Second order derivation of the quartic spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

/***********************************************************************/
// Eighth order spline window function.
namespace EighthOrderSplineWinFunc {
      
  // The window function.
  double W(double& coord);

  // First order derivation of the quartic spline.
  double dW(double& coord);

  // Second order derivation of the quartic spline.
  double d2W(double& rad);

  // Calculate the window function in 3-D.
  void calcWinFunc(double& x,double& y,double& z,double& wFunc);

  // Calculate the first order derivations the window function in 3-D.
  void calcWinFunc1stDerivs(double& x,double& y,double& z,double& xDWFunc,
			    double& yDWFunc,double& zDWFunc);

  // Calculate the second order derivations of the window function 
  // in 3-D.
  void calcWinFunc2ndDerivs(double& x,double& y,double& z,
			    double& xxDWFunc,double& yyDWFunc,
			    double& zzDWFunc,double& xyDWFunc,
			    double& yzDWFunc,double& zxDWFunc);

};

#endif
