//////////////////////////////////////////////////////////////////////
/// 3D Curve Library 
///
/// File: curve3d.h
/// Author: Philip Yang
//////////////////////////////////////////////////////////////////////

#include <vector>
#include <list>
#include "../math/vec.h"
#include "../math/mat.h"
#include "../math/quaternion.h"
#include "database.h"


class Curve
{
 public:
  /// initialization
  Curve(){};
  Curve(Gesture* gesture); /// from the curve data type used in wiimote  

  /// perform polyline simplification 
  void simplify();

  /// sample n points from the curve data which makes it smaller
  /// should choose from uniform sampling and non-uniform sampling
  /// the latter case is based on curvature
  void sample(int n);

  /// compute the minkovski distance
  /// including Manhattan, Euclidean and L_{\infty} norm
  float minkovski_dist(int l);

  /// compute the frechet distance
  float frechet_dist();

  /// compute the hausdorff distance
  float hausdorff_dist();
  
  /// transform the data and give it to the svm learner
  void transform(char* file_name);

  /// print the data to a outside 
  friend inline ostream& operator << ( ostream& out, Curve& curve );

 private:
  /// a list of 3d points as data
  list<Vec3f> _data;
};
