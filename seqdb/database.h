#ifndef __GESTURE_DATABASE_H__
#define __GESTURE_DATABASE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include <cstring>


using namespace std;

class Gesture;
class GestureDB;

// define the datatype that represents a gesture
struct Point
{
  Point(){}
  Point(float x, float y, float z, float acc_x, float acc_y, float acc_z, float pitch, float roll, float yaw)
    :x(x), y(y), z(z), acc_x(acc_x), acc_y(acc_y), acc_z(acc_z), pitch(pitch), roll(roll), yaw(yaw){}
  Point(Point* p)
  {
    x = p->x;
    y = p->y;
    z = p->z;
    acc_x = p->acc_x;
    acc_y = p->acc_y;
    acc_z = p->acc_z;
    pitch = p->pitch;
    roll = p->roll;
    yaw = p->yaw;
  }
  

  // define io related functions for debugging
  friend istream& operator>> (istream& in, Point& p);
  friend ostream& operator<< (ostream& out, Point& p);

  float x;
  float y;
  float z;
  float acc_x;
  float acc_y;
  float acc_z;
  float pitch;
  float roll;
  float yaw;
};


/////////////////////////////////////////////////
/// compute the distance of two gestures
/////////////////////////////////////////////////

inline float minkovski_l1(Point* p1, Point* p2);
inline float minkovski_l2(Point* p1, Point* p2);
inline float minkovski_inf(Point* p1, Point* p2);
inline float sigmoid(Point* p1, Point* p2);



typedef enum
  {
    PUSH,
    WAVE,
    SQUARE,
    CIRCLE,
    FOUR,
    NINE
  } GestureType;
#define NUM_GESTURETYPE 6

GestureType idx2ges(int idx);
int ges2idx(GestureType ges_t);
const char* ges2str(GestureType ges_t);


/**
 * \class Gesture
 *
 * Each gesture is composed of a sequence of motion points
 * and a label depicting its semantic
 */
class Gesture 
{
public:

  /// visualize the gesture
  void draw();

  /// append a new point at the end of the gesture sequence
  void append(Point* point);
  
  void setType(GestureType type){type_gesture = type;}
  GestureType getType(){return type_gesture;}

  int numPoint(){return data.size();}

  /// sample n points from the curve data which makes it smaller
  /// should choose from uniform sampling and non-uniform sampling
  /// the latter case is based on curvature
  void sample(int n);
  
  /// vectorize the points using VQ
  void vectorize();
  
  /// apply the filters described in the German paper
  void filter();
  
  /// transform the data and give it to the svm learner
  void format_svm(char* file_name);

  Point* operator[](unsigned int idx){return idx<data.size() ? data[idx]: NULL;} 

  friend istream& operator>> (istream& in, Gesture& g);
  friend ostream& operator<< (ostream& out, Gesture& g); 

private:
  vector<Point*> data;
  GestureType type_gesture;
};



/////////////////////////////////////////////////
/// compute the distance of two gestures
/////////////////////////////////////////////////
inline float score(Point* p1, Point* p2);
inline float dist_l1(Gesture* g1, Gesture* g2);
inline float dist_l2(Gesture* g1, Gesture* g2);
inline float frechet(Gesture* g1, Gesture* g2);
//float manhattan(Gesture* g1, Gesture* g2);
float smith_waterman(Gesture* g1, Gesture* g2);



// define a database holding gestures
typedef vector<Gesture*> GestureSeq;
class GestureDB 
{
public:
  GestureDB(){}
  GestureDB(const char* file_name):file_name(file_name){}
  void load();
  void save();
  void append(Gesture* g);
  void remove();  
  size_t size(){ return data.size(); };
  Gesture* findNext();  

  Gesture* operator[](size_t idx){ return data[idx]; }

  friend ostream& operator<< (ostream& out, GestureDB& db);

  /// Construct a HMM model for different type of gestures
  void GenHMM();

  // find the k-nearest neighbors
  //void findNearest(int k);  
  GestureType kNN(Gesture* query, int k);
  GestureType predict(Gesture *g);

private:
  const char* file_name;
  GestureSeq data;
};


#endif
