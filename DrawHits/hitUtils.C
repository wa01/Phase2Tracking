#include "Math/DisplacementVector3D.h"
#include "Math/Vector3Dfwd.h"
#include "TObject.h"
#include <iostream>

float dot2D(float v1x, float v1z, float v2x, float v2z) {
  float inProd = v1x*v2x + v1z*v2z;
  float rho1 = sqrt(v1x*v1x + v1z*v1z);
  float rho2 = sqrt(v2x*v2x + v2z*v2z);
  return inProd/rho1/rho2;
};

