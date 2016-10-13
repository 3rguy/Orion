/// Stores all properties of an anisotropy condition

#include "AnisotropyCondition.h"

AnisotropyCondition::AnisotropyCondition() :
    rotationAxisType(1), rotationAngle(0) {

  using namespace std;

  // Anisotropy conditions
  pushBackVector(admissibleAnisotropyTypes,
                 (string) "Surface-Sarcomere-Length-Condition");
  pushBackVector(admissibleAnisotropyTypes,
                 (string) "Surface-Material-Direction-Angle-Condition");

}

/***********************************************************************/
/***********************************************************************/
// check whether the condition is an anisotropy condition.
bool AnisotropyCondition::isAnisotropyType(std::string& t) {

  using namespace std;

  bool flag = false;

  if(find(admissibleAnisotropyTypes.begin(),
          admissibleAnisotropyTypes.end(),t)
    != admissibleAnisotropyTypes.end()) flag = true;

  return flag;
}
