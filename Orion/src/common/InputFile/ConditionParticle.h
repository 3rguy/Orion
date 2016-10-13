/// Stores all properties of a condition particle.

#ifndef ConditionParticle_h_
#define ConditionParticle_h_

#include <iostream>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "mpi.h"


class ConditionParticle {

  public:

    ConditionParticle() : ID(0), materialID(0), materialDirectionAngle(0),
    materialRotationAngle(0), surfaceID(0), surfaceTangentDirection(1),
    restSarcomereLength(0) {};
    ~ConditionParticle() {};

    /// Particle identity number.
    int& getID() { return ID; };
    int& getMaterialID() { return materialID; };

    /// Coordinates manipulating
    dbVector& getCoords() { return coords; };

    /// surfaces and elements
    dbVector& getSurfaceNormal(int ID);
    dbMatrix& getAllSurfaceNormals() { return surfaceNormals; };
    intVector& getElems() { return elems; };

    /*******************************************************************/
    /// anisotropy
    int& getSurfaceID() { return surfaceID; };
    int& getSurfaceTangentDirection() { return surfaceTangentDirection; };
    double& getMaterialDirectionAngle() { return materialDirectionAngle; };
    double& getMaterialRotationAngle() {return materialRotationAngle;}
    dbMatrix3& getMaterialDirections() { return materialDirections; };



    /*******************************************************************/
    /// cardiac mechanics
    ///
    /// relaxed sarcomere length
    double& getRestSarcomereLength() { return restSarcomereLength; };


    /********************************************************************/
    /// Clear all class arrays which are not locally needed.
    void clearArrays();

  private:

    int ID,materialID;
    int surfaceID;
    int surfaceTangentDirection;
    double materialDirectionAngle;
    double materialRotationAngle;
    dbVector coords;
    dbMatrix surfaceNormals;
    intVector elems;

    /********************************************************************/
    /// anisotropy
    dbMatrix3 materialDirections;

    /********************************************************************/
    /// cardiac mechanics
    double restSarcomereLength;


};

#endif
