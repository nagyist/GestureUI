#include <cstdio>
#include <cstdlib>
#include "curve3d.h"

/// obtain the curve from a raw gesture data
Curve::Curve(Gesture* gesture)
{
  /// get the data from an integrated gesture type
  //_data = gesture->integrate();  
}


/// perform polyline simplification 
void Curve::simplify()
{

}

/// sample n points from the curve data which makes it smaller
/// should choose from uniform sampling and non-uniform sampling
/// the latter case is based on curvature
void Curve::sample(int n)
{
  int num_pts = _data.size();
  srand(time(NULL));
  
  //list<Vec3f>::iterator iter = _data.begin();
  int idx = 0;
  float thresh = 0.5;
  // while ( num_pts > n )
  //   {
  //     Vec3f curr;// = _data[idx];
  //     bool flag = false;
  //     if ( rand()/RAND_MAX > thresh )
  // 	{
  // 	  flag = true;
  // 	  _data.remove(idx);
  // 	  if (--num_pts == n)
  // 	    break;
  // 	}
  //     else
  // 	++idx;     
  //   }
}

/// compute the minkovski distance
/// including Manhattan, Euclidean and L_{\infty} norm
float minkovski_dist(int l);

/// compute the frechet distance
float frechet_dist();

/// compute the hausdorff distance
float hausdorff_dist();




/// print the data to a outside 
inline ostream& operator << ( ostream& out, Curve& curve );

