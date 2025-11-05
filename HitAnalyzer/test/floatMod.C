#include <iostream>

float floatMod(float a, float b) {
  //
  // change fmod behaviour for negative numbers
  //   (fmod will return -fmod(|a|,b) for a<0)
  //   
  float result = fmod(a,b);
  if ( result < 0. )  result += b;
  return result;
}
