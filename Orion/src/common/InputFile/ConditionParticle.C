// Stores all properties of a condition particle.

#include "ConditionParticle.h"


/**********************************************************************/
/**********************************************************************/
// Clear all class arrays which are not locally needed.
void ConditionParticle::clearArrays() {

  using namespace std;

  resizeArray(coords,0);
  resizeArray(surfaceNormals,0);
  resizeArray(elems,0);
  resizeArray(materialDirections,0);
}

