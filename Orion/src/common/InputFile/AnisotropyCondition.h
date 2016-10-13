/// Stores all properties of an anisotropy condition

#ifndef AnisotropyCondition_h_
#define AnisotropyCondition_h_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "commonFunctions.h"
#include "commonTypedefs.h"
#include "Condition.h"

class AnisotropyCondition : public virtual Condition {

  public:

    AnisotropyCondition();
    ~AnisotropyCondition() {}

    bool setParam(std::string name,double value) {
      using namespace std;
      bool info = false;

//      if(name == "rotationAngle") {
//        rotationAngle = value;
//        info = true;
//      }
//      else if(name == "rotationAxisType") {
      if(name == "rotationAxisType") {
        rotationAxisType = value;
        info = true;
      }
      else if(Condition::setParam(name,value))
        info = true;
      else {
        std::cerr << "In AnisotropyCondition::setParam parameter "<< name << " does not exist." << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

      return info;

    };

    /// check whether the condition is an anisotropy condition.
    bool isAnisotropyType(std::string& t);
    int& getRotationAxisType() { return rotationAxisType; };
    dbVector& getRotationAngle() { return rotationAngle; };

  private:

    std::vector<std::string> admissibleAnisotropyTypes;
    int rotationAxisType; // 0: user-specified, 1,2,3: local axes
    dbVector rotationAngle;

};

#endif
