/// Stores all properties of a resultant reaction on a
/// set of elements or nodes.

#ifndef ResultantReactionCondition_h_
#define ResultantReactionCondition_h_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "Condition.h"
#include "ConditionElement.h"
#include "ConditionParticle.h"

class ResultantReactionCondition : public virtual Condition {

  public:

    ResultantReactionCondition();
    ~ResultantReactionCondition() {}

    /// check whether specified resultant reaction type is admissible
    bool isResultantReactionType(std::string& t);
    bool isSurfaceResultantReaction();

  protected:

    std::vector<std::string> admissibleResultantReactionsTypes;
    std::vector<std::string> admissibleSurfaceResultantReactionTypes;

};

#endif
