#include "vec.h"
#include "mat.h"
#include "quaternion.h"

#include <iostream>
using namespace std;

int main()
{
  Qtrni q1(1,2,3,4), q2(3,4,2,1);
  cout<< q1 <<endl;
  cout<< conjugate(q1) <<endl;
  cout<< prod(q1, q2) <<endl;
  cout<< prod(q1, conjugate(q1)) <<endl;
  return 0;
}
