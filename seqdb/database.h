#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>

using namespace std;

// define the datatype that represents a gesture
struct Point
{
  Point(){}
  Point(float x, float y, float z, float pitch, float roll, float yaw)
    :x(x), y(y), z(z), pitch(pitch), roll(roll), yaw(yaw){}

  // define io related functions for debugging
  friend istream& operator>> (istream& in, Point& p);
  friend ostream& operator<< (ostream& out, Point& p);

  float x;
  float y;
  float z;
  float pitch;
  float roll;
  float yaw;
};

inline float score_l1(Point* p1, Point* p2);
inline float score_l2(Point* p1, Point* p2);


typedef enum
  {
    PUSH,
    WAVE
  } GestureType;
#define NUM_GESTURETYPE 2



/**
 * \class Gesture
 *
 * Each gesture is composed of a sequence of motion points
 * and a label depicting its semantic
 */
class Gesture 
{
public:

  // visualize the gesture
  void draw();

  // append a new point at the end of the gesture sequence
  void append(Point* point);
  
  void setType(GestureType type){type_gesture = type;}
  GestureType getType(){return type_gesture;}

  int numPoint(){return data.size();}

  /// sample n points from the curve data which makes it smaller
  /// should choose from uniform sampling and non-uniform sampling
  /// the latter case is based on curvature
  void sample(int n);
  
  /// apply the filters described in the German paper
  void filter();
  
  /// transform the data and give it to the svm learner
  void transform(char* file_name);

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
  Gesture* findNext();  

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


