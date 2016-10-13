// Stores all properties of a resultant reaction on a
// set of elements or nodes.

#include "ResultantReactionCondition.h"

ResultantReactionCondition::ResultantReactionCondition() {

  using namespace std;

  // resultant reactions
  pushBackVector(admissibleResultantReactionsTypes,
                 (string) "Resultant-Internal-Traction");
  pushBackVector(admissibleSurfaceResultantReactionTypes,
                 (string) "Resultant-Internal-Traction");


}

/***********************************************************************/
/***********************************************************************/
// check whether specified resultant reaction type is admissible
bool ResultantReactionCondition::isResultantReactionType(std::string& t) {

  using namespace std;

  bool flag = false;

  if(find(admissibleResultantReactionsTypes.begin(),
          admissibleResultantReactionsTypes.end(),t)
    != admissibleResultantReactionsTypes.end()) flag = true;

  return flag;
}

// check whether condition type is a surface reaction
bool ResultantReactionCondition::isSurfaceResultantReaction() {

  using namespace std;

  return (contains(admissibleSurfaceResultantReactionTypes,type));
}
