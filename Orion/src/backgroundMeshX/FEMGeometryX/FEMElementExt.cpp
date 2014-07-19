#include "FEMElementExt.h"

FEMElementExt::FEMElementExt(int usedDOF) : FEMElement(usedDOF) {}

FEMElementExt::FEMElementExt(FEMElement& element) : FEMElement(element) {}


