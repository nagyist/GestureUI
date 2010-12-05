///////////////////////////////////////////////////////////
/// 
/// Quaternion Math Library
///
///////////////////////////////////////////////////////////

#ifndef __QUATERNION_HEADER__
#define __QUATERNION_HEADER__

#include <iostream>
#include <cstring>
#include "mat.h"
#include "vec.h"


#include <algorithm>
#ifndef _WIN32
// This is actually here to correct for
// a bug in the MSVC STL.
using std::min;
using std::max;
#endif

#include <iostream>
#include <cstring>
#include <cmath>
#include <assert.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include <GL/gl.h>

//==========[ Forward References ]=========================

template <class T> class Vec;
template <class T> class Vec3;
template <class T> class Vec4;
template <class T> class Mat3;
template <class T> class Mat4;

// quaternion and double quaternion
template <class T> class Qtrn;
//template <class T> class DQtrn;

template <class T> T operator*( const Vec<T>& a, const Vec<T>& b ); 
template <class T> Vec<T> operator-( const Vec<T>& v );
template <class T> Vec<T> operator*( const Vec<T>& a, const double d );
template <class T> Vec<T> operator*( const double d, const Vec<T>& a );
template <class T> Vec<T> operator/( const Vec<T>& a, const double d );
template <class T> Vec<T> operator^( const Vec<T>& a, const Vec<T>& b );
template <class T> bool operator==( const Vec<T>& a, const Vec<T>& xb );
template <class T> bool operator!=( const Vec<T>& a, const Vec<T>& b );
template <class T> std::ostream& operator<<( std::ostream& os, const Vec<T>& a );
template <class T> std::istream& operator>>( std::istream& is, const Vec<T>& a );
template <class T> Vec<T> minimum( const Vec<T>& a, const Vec<T>& b );
template <class T> Vec<T> maximum( const Vec<T>& a, const Vec<T>& b );
template <class T> Vec<T> prod( const Vec<T>& a, const Vec<T>& b );

template <class T> inline T operator *( const Vec3<T>& a, const Qtrn<T>& b );
template <class T> inline T operator *( const Qtrn<T>& b, const Vec3<T>& a );
template <class T> inline Vec3<T> operator -(const Vec3<T>& v);
template <class T> inline Vec3<T> operator *(const Vec3<T>& a, const double d );
template <class T> inline Vec3<T> operator *(const double d, const Vec3<T>& a);
template <class T> inline Vec3<T> operator *(const Mat4<T>& a, const Vec3<T>& v);
template <class T> inline Vec3<T> operator *(const Vec3<T>& v, Mat4<T>& a);
template <class T> inline T operator *(const Vec3<T>& a, const Vec3<T>& b);
template <class T> inline Vec3<T> operator *( const Mat3<T>& a, const Vec3<T>& v );
template <class T> inline Vec3<T> operator *( const Vec3<T>& v, const Mat3<T>& a );
template <class T> inline Vec3<T> operator /(const Vec3<T>& a, const double d);
template <class T> inline Vec3<T> operator ^(const Vec3<T>& a, const Vec3<T>& b);
template <class T> inline bool operator ==(const Vec3<T>& a, const Vec3<T>& b);
template <class T> inline bool operator !=(const Vec3<T>& a, const Vec3<T>& b);
template <class T> inline std::ostream& operator <<( std::ostream& os, const Vec3<T>& v );
template <class T> inline std::istream& operator >>( std::istream& is, Vec3<T>& v );
template <class T> inline Vec3<T> minimum( const Vec3<T>& a, const Vec3<T>& b );
template <class T> inline Vec3<T> maximum(const Vec3<T>& a, const Vec3<T>& b);
template <class T> inline Vec3<T> prod(const Vec3<T>& a, const Vec3<T>& b );

template <class T> inline Qtrn<T> operator -( const Qtrn<T>& v );
template <class T> inline Qtrn<T> operator *(const Qtrn<T>& a, const double d );
template <class T> inline Qtrn<T> operator *(const double d, const Qtrn<T>& a);
template <class T> inline T operator *(const Qtrn<T>& a, const Qtrn<T>& b);
template <class T> inline Qtrn<T> operator *(const Mat4<T>& a, const Qtrn<T>& v);
template <class T> inline Qtrn<T> operator *( const Qtrn<T>& v, Mat4<T>& a );
template <class T> inline Qtrn<T> operator /(const Qtrn<T>& a, const double d);
template <class T> inline bool operator ==(const Qtrn<T>& a, const Qtrn<T>& b);
template <class T> inline bool operator !=(const Qtrn<T>& a, const Qtrn<T>& b);
template <class T> inline std::ostream& operator <<( std::ostream& os, const Qtrn<T>& v );
template <class T> inline std::istream& operator >>( std::istream& is, Qtrn<T>& v );
template <class T> inline Qtrn<T> minimum( const Qtrn<T>& a, const Qtrn<T>& b );
template <class T> inline Qtrn<T> maximum( const Qtrn<T>& a, const Qtrn<T>& b);
template <class T> inline Qtrn<T> prod(const Qtrn<T>& a, const Qtrn<T>& b );
template <class T> inline Qtrn<T> prod(const Vec3<T>& a, const Qtrn<T>& b );
template <class T> inline Qtrn<T> prod(const Qtrn<T>& a, const Vec3<T>& b );
template <class T> inline Qtrn<T> conjugate( const Qtrn<T>& q );
template <class T> Vec3<T> Qtrn2Vec3(Qtrn<T> &v);

//==========[ Exception Classes ]==========================

//class VectorSizeMismatch {};

//==========[ class Qtrn (Quaternion) ]=================================

template <class T>
class Qtrn {

  //---[ Private Variable Declarations ]-------

  /// q = w + xi + yj + zk
  T n[4];

 public:
	
  //---[ Constructors ]------------------------

  Qtrn() { n[0] = 0.0; n[1] = 0.0; n[2] = 0.0; n[3] = 0.0; }
  Qtrn( const T x, const T y, const T z, const T w )
    { n[0] = x; n[1] = y; n[2] = z; n[3] = w; }
  Qtrn( const Qtrn& v )
    { n[0] = v.n[0]; n[1] = v.n[1]; n[2] = v.n[2]; n[3] = v.n[3]; }

  //---[ Equal Operators ]---------------------

  Qtrn<T>& operator =( const Qtrn<T>& v )
    { n[0] = v.n[0]; n[1] = v.n[1]; n[2] = v.n[2]; n[3] = v.n[3];
      return *this; }
  Qtrn<T>& operator +=( const Qtrn<T>& v )
    { n[0] += v.n[0]; n[1] += v.n[1]; n[2] += v.n[2]; n[3] += v.n[3];
      return *this; }
  Qtrn<T>& operator -= ( const Qtrn<T>& v )
    { n[0] -= v.n[0]; n[1] -= v.n[1]; n[2] -= v.n[2]; n[3] -= v.n[3];
      return *this; }

  // scalar operation
  Qtrn<T>& operator *= ( const T d )
    { n[0] *= d; n[1] *= d; n[2] *= d; n[3] *= d; return *this; }
  Qtrn<T>& operator /= ( const T d )
    { n[0] /= d; n[1] /= d; n[2] /= d; n[3] /= d; return *this; }

  //---[ Access Operators ]--------------------

  T& operator []( int i )
    { return n[i]; }
  T operator []( int i ) const 
  { return n[i]; }

  //---[ Arithmetic Operators ]----------------

  Qtrn<T> operator-( const Qtrn<T>& a ) { return Qtrn<T>(n[0]-a.n[0],n[1]-a.n[1],n[2]-a.n[2],n[3]-a.n[3]); }
  Qtrn<T> operator+( const Qtrn<T>& a ) { return Qtrn<T>(a.n[0]+n[0],a.n[1]+n[1],a.n[2]+n[2],a.n[3]+n[3]); }

  //---[ Length Methods ]----------------------

  double length2() const
  { return n[0]*n[0] + n[1]*n[1] + n[2]*n[2] + n[3]*n[3]; }
  double length() const
  { return sqrt( length2() ); }

  //---[ Zero Test ]---------------------------

  bool isZero() const { return n[0]==0&&n[1]==0&&n[2]==0&&n[3]==0; }
  void zeroElements() { memset(n,0,4*sizeof(T)); }

  //---[ Normalization ]-----------------------

  void normalize() {
    double len = length();
    assert(len != 0);
    n[0] /= len; n[1] /= len; n[2] /= len; n[3] /= len;
  }
       

  //---[ Friend Methods ]----------------------

#ifdef WIN32
  // VCC is non-standard
  template <class U> friend T operator *( const Vec3<T>& a, const Qtrn<T>& b );
  template <class U> friend T operator *( const Qtrn<T>& b, const Vec3<T>& a );
  //	template <class U> friend Qtrn<T> operator -( const Qtrn<T>& v );
  template <class U> friend Qtrn<T> operator *( const Qtrn<T>& a, const double d );
  template <class U> friend Qtrn<T> operator *( const double d, const Qtrn<T>& a );
  template <class U> friend T operator *( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend Qtrn<T> operator *( const Mat4<T>& a, const Qtrn<T>& v );
  //	template <class U> friend Qtrn<T> operator *( const Qtrn<T>& v, const Mat4<T>& a );
  template <class U> friend Qtrn<T> operator /( const Qtrn<T>& a, const double d );
  //	template <class U> friend Qtrn<T> operator ^( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend bool operator ==( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend bool operator !=( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend std::ostream& operator <<( std::ostream& os, const Qtrn<T>& v );
  template <class U> friend std::istream& operator >>( std::istream& is, Qtrn<T>& v );
  template <class U> friend Qtrn<T> minimum( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend Qtrn<T> maximum( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend Qtrn<T> prod( const Qtrn<T>& a, const Qtrn<T>& b );
  template <class U> friend Qtrn<T> conjugate( const Qtrn<T>& q );

#else

  // 3d vector and quaternion product
  friend T operator * <>( const Vec3<T>& a, const Qtrn<T>& b );
  friend T operator * <>( const Qtrn<T>& b, const Vec3<T>& a );

  // scalar quaternion product
  //	friend Qtrn<T> operator - <>( const Qtrn<T>& v );
  friend Qtrn<T> operator * <>( const Qtrn<T>& a, const double d );
  friend Qtrn<T> operator * <>( const double d, const Qtrn<T>& a );

  friend T operator * <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend Qtrn<T> operator * <>( const Mat4<T>& a, const Qtrn<T>& v );
  //	friend Qtrn<T> operator * <>( const Qtrn<T>& v, const Mat4<T>& a );
  friend Qtrn<T> operator / <>( const Qtrn<T>& a, const double d );
  //	friend Qtrn<T> operator ^ <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend bool operator == <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend bool operator != <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend std::ostream& operator << <>( std::ostream& os, const Qtrn<T>& v );
  friend std::istream& operator >> <>( std::istream& is, Qtrn<T>& v );
  friend Qtrn<T> minimum <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend Qtrn<T> maximum <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend Qtrn<T> prod <>( const Qtrn<T>& a, const Qtrn<T>& b );
  friend Qtrn<T> prod <>( const Vec3<T>& a, const Qtrn<T>& b );
  friend Qtrn<T> prod <>( const Qtrn<T>& a, const Vec3<T>& b );
  friend Qtrn<T> conjugate <>( const Qtrn<T>& q );

#endif
};

typedef Qtrn<int> Qtrni;
typedef Qtrn<float> Qtrnf;
typedef Qtrn<double> Qtrnd;

//==========[ Vec Methods ]================================


template <class T>
inline Qtrn<T> operator -( const Qtrn<T>& v ) 
{
  return Qtrn<T>( -v.n[0], -v.n[1], -v.n[2], -v.n[3] );
}

template <class T>
inline Qtrn<T> operator *(const Qtrn<T>& a, const double d ) 
{
  return Qtrn<T>( a.n[0] * d, a.n[1] * d, a.n[2] * d, a.n[3] * d );
}

template <class T>
inline Qtrn<T> operator *(const double d, const Qtrn<T>& a) 
{
  return a * d;
}

template <class T>
inline T operator *(const Qtrn<T>& a, const Qtrn<T>& b) {
  return a.n[0]*b.n[0] + a.n[1]*b.n[1] + a.n[2]*b.n[2] + a.n[3]*b.n[3];
}

template <class T>
inline Qtrn<T> operator *(const Mat4<T>& a, const Qtrn<T>& v) 
{
  return Qtrn<T>( a.n[0]*v.n[0]+a.n[1]*v.n[1]+a.n[2]*v.n[2]+a.n[3]*v.n[3],
		  a.n[4]*v.n[0]+a.n[5]*v.n[1]+a.n[6]*v.n[2]+a.n[7]*v.n[3],
		  a.n[8]*v.n[0]+a.n[9]*v.n[1]+a.n[10]*v.n[2]+a.n[11]*v.n[3],
		  a.n[12]*v.n[0]+a.n[13]*v.n[1]+a.n[14]*v.n[2]+a.n[15]*v.n[3]);
}

template <class T>
inline Qtrn<T> operator *( const Qtrn<T>& v, Mat4<T>& a )
{
  return a.transpose() * v;
}

template <class T>
inline Qtrn<T> operator /(const Qtrn<T>& a, const double d) 
{
  return Qtrn<T>( a.n[0] / d, a.n[1] / d, a.n[2] / d, a.n[3] / d );
}

template <class T>
inline bool operator ==(const Qtrn<T>& a, const Qtrn<T>& b) 
{
  return a.n[0] == b.n[0] && a.n[1] == b.n[1] && a.n[2] == b.n[2] 
    && a.n[3] == b.n[3];
}

template <class T>
inline bool operator !=(const Qtrn<T>& a, const Qtrn<T>& b) 
{
  return !( a == b );
}

template <class T>
inline std::ostream& operator <<( std::ostream& os, const Qtrn<T>& v ) 
{
  return os << v.n[0] << " " << v.n[1] << " " << v.n[2] << " " << v.n[3];
}

template <class T>
inline std::istream& operator >>( std::istream& is, Qtrn<T>& v ) 
{
  return is >> v.n[0] >> v.n[1] >> v.n[2] >> v.n[3];
}

template <class T>
inline Qtrn<T> minimum( const Qtrn<T>& a, const Qtrn<T>& b ) 
{
  return Qtrn<T>( min(a.n[0],b.n[0]), min(a.n[1],b.n[1]), min(a.n[2],b.n[2]),
		  min(a.n[3],b.n[3]) );
}

template <class T>
inline Qtrn<T> maximum( const Qtrn<T>& a, const Qtrn<T>& b) 
{
  return Qtrn<T>( max(a.n[0],b.n[0]), max(a.n[1],b.n[1]), max(a.n[2],b.n[2]),
		  max(a.n[3],b.n[3]) );
}


//---[ Conjugate ]---------------------------
// q* = w - xi - yj - zk
template <class T>
inline Qtrn<T> conjugate( const Qtrn<T>& q )
{
  return Qtrn<T>(q.n[0], -q.n[1], -q.n[2], -q.n[3]);
}


/// Quaternion product
template <class T>
inline Qtrn<T> prod(const Qtrn<T>& a, const Qtrn<T>& b ) 
{
  //return Qtrn<T>( a.n[0]*b.n[0], a.n[1]*b.n[1], a.n[2]*b.n[2], a.n[3]*b.n[3] );
  return Qtrn<T>( 
		 -a.n[0]*b.n[0] - a.n[1]*b.n[1] - a.n[2]*b.n[2] + a.n[3]*b.n[3],
		 a.n[0]*b.n[3] + a.n[1]*b.n[2] - a.n[2]*b.n[1] + a.n[3]*b.n[0],
		 -a.n[0]*b.n[2] + a.n[1]*b.n[3] + a.n[2]*b.n[0] + a.n[3]*b.n[1],
		 a.n[0]*b.n[1] - a.n[1]*b.n[0] + a.n[2]*b.n[3] + a.n[3]*b.n[2]  );
}

/// Quaternion product
template <class T>
inline Qtrn<T> prod(const Qtrn<T>& a, const Vec3<T>& b ) 
{
  //return Qtrn<T>( a.n[0]*b.n[0], a.n[1]*b.n[1], a.n[2]*b.n[2], a.n[3]*b.n[3] );
  return Qtrn<T>( 
		 -a.n[1]*b.n[0] - a.n[2]*b.n[1] + a.n[3]*b.n[2],
		 a.n[0]*b.n[2] + a.n[1]*b.n[1] - a.n[2]*b.n[0],
		 -a.n[0]*b.n[1] + a.n[1]*b.n[2] + a.n[3]*b.n[0],
		 a.n[0]*b.n[0] + a.n[2]*b.n[2] + a.n[3]*b.n[1]  );
}

/// Quaternion product
template <class T>
inline Qtrn<T> prod(const Vec3<T>& a, const Vec3<T>& b ) 
{
  //return Qtrn<T>( a.n[0]*b.n[0], a.n[1]*b.n[1], a.n[2]*b.n[2], a.n[3]*b.n[3] );
  return Qtrn<T>( 
		 - a.n[0]*b.n[1] - a.n[1]*b.n[2] + a.n[2]*b.n[3],
		 a.n[0]*b.n[2] - a.n[1]*b.n[1] + a.n[2]*b.n[0],
		 a.n[0]*b.n[3] + a.n[1]*b.n[0] + a.n[2]*b.n[1],
		 a.n[0]*b.n[0] + a.n[1]*b.n[3] + a.n[2]*b.n[2]  );
}


/// Convert the quaternion to 3d vector
/// return the vector portion
template <class T>
Vec3<T> Qtrn2Vec3(Qtrn<T> &v) 
{
  return Vec3<T>(v.n[1], v.n[2], v.n[3]);
}


#endif //__QUATERNION_HEADER__
