#include <cmath>
#include <cstring>
#include <iostream>
#include <list>

#include "database.h"

using namespace std;


GestureType idx2ges(int idx)
{
  return (GestureType)idx;
}


int ges2idx(GestureType ges_t)
{
  return (int)ges_t;
}

const char* ges2str(GestureType ges_t)
{
  switch(ges_t)
    {
    case PUSH:
      return "PUSH";
      break;
    case WAVE:
      return "WAVE";
      break;
    case CIRCLE:
      return "CIRCLE";
      break;
    case SQUARE:
      return "SQUARE";
      break;
    case FOUR:
      return "FOUR";
      break;
    case NINE:
      return "NINE";
      break;      
    default:
      return "UNKNOWN";
      break;
    }
}


inline ostream& operator<< (ostream& out, GestureDB& db)
{
  Gesture *ptr_g;
  int count = 0;
  
  // we also want to add meta data to the record
  while ( (ptr_g = db.findNext()) && (++count) )
    out<<"# gesture #"<<count<<endl<<*ptr_g<<endl;
  return out;
}


inline istream& operator>> (istream& in, Point& p)
{
  return in>>p.x>>p.y>>p.z
	   >>p.acc_x>>p.acc_y>>p.acc_z
	   >>p.pitch>>p.roll>>p.yaw;
}


inline ostream& operator<< (ostream& out, Point& p)
{
  return out<<p.x<<" "
	    <<p.y<<" "
	    <<p.z<<" "
	    <<p.acc_x<<" "
	    <<p.acc_y<<" "
	    <<p.acc_z<<" "
	    <<p.pitch<<" "
	    <<p.roll<<" "
	    <<p.yaw;  
}


inline float minkovski_l1(Point* p1, Point* p2)
{
  return ( abs(p1->x - p2->x) + 
	   abs(p1->y - p2->y) +
	   abs(p1->z - p2->z) +
	   abs(p1->acc_x - p2->acc_x) + 
	   abs(p1->acc_y - p2->acc_y) +
	   abs(p1->acc_z - p2->acc_z) +
	   abs(p1->pitch - p2->pitch) + 
	   abs(p1->roll - p2->roll) + 
	   abs(p1->yaw - p2->yaw) );
}


inline float minkovski_l2(Point* p1, Point* p2)
{
  return sqrt( pow(p1->x - p2->x, 2) + 
	       pow(p1->y - p2->y, 2) +
	       pow(p1->z - p2->z, 2) +
	       pow(p1->acc_x - p2->acc_x, 2) + 
	       pow(p1->acc_y - p2->acc_y, 2) +
	       pow(p2->acc_z - p2->acc_z, 2) +
	       pow(p1->pitch - p2->pitch, 2) + 
	       pow(p1->roll - p2->roll, 2) + 
	       pow(p1->yaw - p2->yaw, 2) );
}


/// we can't find a good way to determine weights
inline float sigmoid(Point* p1, Point* p2)
{
  return 1/(1 + exp(sqrt( pow(p1->x - p2->x, 2) + 
			  pow(p1->y - p2->y, 2) +
			  pow(p1->z - p2->z, 2) +
			  pow(p1->acc_x - p2->acc_x, 2) + 
			  pow(p1->acc_y - p2->acc_y, 2) +
			  pow(p2->acc_z - p2->acc_z, 2) +
			  pow(p1->pitch - p2->pitch, 2) + 
			  pow(p1->roll - p2->roll, 2) + 
			  pow(p1->yaw - p2->yaw, 2) )));
}


inline istream& operator>> (istream& in, Gesture& g)
{
  return in;
}


inline ostream& operator<< (ostream& out, Gesture& g)
{
  switch (g.type_gesture)
    {
    case PUSH:
      out<<"@ PUSH"<<endl;
      break;
    case WAVE:
      out<<"@ WAVE"<<endl;
      break;

    case CIRCLE:
      out<<"@ CIRCLE"<<endl;
      break;
   
    case SQUARE:
      out<<"@ SQUARE"<<endl;
      break;
      
    case FOUR:
      out<<"@ FOUR"<<endl;
      break;

    case NINE:
      out<<"@ NINE"<<endl;
      break;
      
    default:
      out<<"@ UNKNOWN"<<endl;
      break;
    }
  
  for (vector<Point*>::iterator iter = g.data.begin(); iter != g.data.end(); ++iter)
    {
      out<<**iter<<endl;
    }
  return out<<endl;
}


// visualize the gesture
void Gesture::draw(){}


// append a new point at the end of the gesture sequence
void Gesture::append(Point* point)
{
  data.push_back(point);
}

/////////////////////////////////////////////////
/// Distance
/////////////////////////////////////////////////

inline float score(Point* p1, Point* p2)
{
  return (p1 == NULL || p2 == NULL)? -2 : sigmoid(p1, p2);
  //return (p1 == NULL || p2 == NULL)? -1 : 1/(1+minkovski_l2(p1, p2));
}


inline float dist_l1(Gesture* g1, Gesture* g2)
{
  float res = 0.0f;
  for (int i=0; i<min(g1->numPoint(), g2->numPoint()); ++i)
    {
      //res += abs((*g1)[i] - (*g2)[i]);
      res += score((*g1)[i], (*g2)[i]);
    }

  return res;  
}


inline float dist_l2(Gesture* g1, Gesture* g2)
{
  float res = 0.0f;
  for (int i=0; i<min(g1->numPoint(), g2->numPoint()); ++i)
    {
      //res += pow((*g1)[i] - (*g2)[i], 2);
      res += pow(score((*g1)[i], (*g2)[i]), 2);
    }
  
  //return res;
  return sqrt(res);  
}
  


float smith_waterman(Gesture* g1, Gesture* g2)
{
  int len1 = g1->numPoint();
  int len2 = g2->numPoint();
  float *score_buf = new float[len1];
  memset((void*)score_buf, 0, len1*sizeof(float)); 
  //for (int i=0; i<len1; score_buf[i]=4e7, ++i);

  float prev_l = 0;
  float best = 0;
  for ( int j=1; j<len2; ++j )
    for ( int i=1; i<len1; ++i )
      {	
  	float tmp = max( max( score( (*g1)[i], (*g2)[j-1] ) + score_buf[i],
  			      score( (*g1)[i-1], (*g2)[j] ) + prev_l ),
  			 score( (*g1)[i-1], (*g1)[j-1] ) + score_buf[i-1] );
	
  	score_buf[i-1] = prev_l;
  	prev_l = tmp;
  	best = max(best, tmp);
      }

  return best;
}



inline string chomp(string& str)
{                                                                                
  int begin_idx = str.find_first_not_of(" \t");                                 
  int length = str.find_last_not_of(" \t") - begin_idx + 1;                     
  str = str.substr(begin_idx, length);                                          
  return str;
}


void GestureDB::load()
{
  ifstream fin;
  fin.open(file_name);

  if ( !fin )
    {
      cerr<<"Error: cannot open the database file "<<file_name<<endl;
      return;
    }
  
  /// check if the file is corrupted
  /// check if the file is empty 
  
  string line;
  bool first_time = 1;
  bool is_reading = 0;
  stringstream ss ( stringstream::in | stringstream::out );
  Gesture* ptr_g = new Gesture();
  Point* ptr_p = new Point();	 

  while (getline(fin, line))
    {
      // an empty line marks the end of a gesture record
      if ( line.size() == 0 || chomp(line).size() == 0 )
	{
	  if ( first_time )
	    first_time = 0;      	  
	  else if ( is_reading )
	    {
	      data.push_back(ptr_g);
	      ptr_g = new Gesture();	  
	      is_reading = false;
	    }
	 	  
	  continue;
	}

      // # for comment and meta data for human reader
      if ( line[0] == '#' )
	continue; 

      // @ for class label
      if ( line[0] == '@' )
	{
	  string label = line.substr(1, line.size()-1);
	  chomp(label);
	  if ( label == "PUSH" )
	    ptr_g->setType(PUSH);
	  else if ( label == "WAVE" )
	    ptr_g->setType(WAVE);
	  else if ( label == "CIRCLE" )
	    ptr_g->setType(CIRCLE);
	  else if ( label == "SQUARE" )
	    ptr_g->setType(SQUARE);
 	  else if ( label == "FOUR" )
	    ptr_g->setType(FOUR);
	  else if ( label == "NINE" )
	    ptr_g->setType(NINE);

	  continue;
	}

      is_reading = true;     

      // if we are in the middle of a record
      ss << line;
      ss >> *ptr_p;
      ss.clear(); // this is vital...!!!
      //cout<<*ptr_p<<endl;
      ptr_g->append(ptr_p);	  
      ptr_p = new Point();	  
    }  


  //cout<<*this<<endl;

  fin.close();
}
  

void GestureDB::save()
{
  ofstream fout;
  // all the previously loaded examples will also be written back
  fout.open(file_name);
  fout<<*this<<endl;
  fout.close();
}


void GestureDB::append(Gesture* g)
{
  data.push_back(g);
}
 

void GestureDB::remove()
{
  /// do nothing
}


/// find the next element in the sequence
Gesture* GestureDB::findNext()
{
  static size_t curr_idx = 0;
  if ( curr_idx >= data.size() )
    {
      curr_idx = 0;
      return NULL;
    }
  return data[curr_idx++];
}


/// compute the distance of two gestures 
/// calls one of the distance metrices 
/// in the gesture library
inline float distance(Gesture* g1, Gesture* g2)
{
  //return dist_l2(g1, g2);
  return smith_waterman(g1,g2);
}


/// find the k-nearest neighbors
GestureType GestureDB::kNN(Gesture* query, int k)
{
  //map<float, Gesture*> knn_map;
  list<Gesture*> knn;
  list<float> weight;
  
  printf("\t\t database size: %d\n", data.size());
  Gesture* curr_ges = data[0]; //findNext();
  knn.push_front(curr_ges);
  weight.push_front( distance(query, curr_ges) );     

  //while ( curr_ges = findNext() )
  for ( int idx=1; idx<data.size(); ++idx )
    {
      curr_ges = data[idx];

      // insert the gesture to the appropriate position
      float curr_dist = distance(query, curr_ges);
      list<float>::iterator weight_iter = weight.begin();
      list<Gesture*>::iterator knn_iter = knn.begin();	  
      
      bool is_added = false;
      for (; knn_iter != knn.end(); ++knn_iter, ++weight_iter )	  
	if ( *weight_iter < curr_dist )
	  {
	    weight.insert(weight_iter, curr_dist);
	    knn.insert(knn_iter, curr_ges);
	    is_added = true;
	    break;
	  }

      if ( (knn.size() < k) && !is_added )
	{
	  weight.push_back(curr_dist);
	  knn.push_back(curr_ges);	  
	}

      if ( knn.size() > k )
	{
	  weight.pop_back();
	  knn.pop_back();
	}
    }
  
  int count[NUM_GESTURETYPE] = {0};
  list<float>::iterator weight_iter = weight.begin(); 
  list<Gesture*>::iterator knn_iter = knn.begin();
  for (; knn_iter != knn.end(); ++knn_iter, ++weight_iter )	  
    {
      printf("\t\t kNN: gesture type %s, smith_waterman %f\n", 
	     ges2str((*knn_iter)->getType()), *weight_iter);
      
      //++count[(int)(*knn_iter)->getType()];
      ++count[ges2idx((*knn_iter)->getType())];
    }

  GestureType majority = (GestureType) 0;
  int max_count = 0;
  for ( int i=0; i<NUM_GESTURETYPE; ++i )
    {
      if ( max_count < count[i] )
	{
	  max_count = count[i];
	  majority = (GestureType) i;
	}
    }
  return majority;
}


// // // test the result
// int main()
// {
//   Point* ptr_p = new Point();//(1,2,3,4,54,6);
//   Gesture* ptr_g = new Gesture();
//   GestureDB* ptr_db = new GestureDB("tmp_gesture.db");
//   ptr_db->load();  

//   // cin>>*ptr_p;
//   // cout<<*ptr_p<<endl;
//   // delete ptr_p;
  
//   for ( int k=0; k<2; ++k )
//     {
//       for ( int i=0; i<3; ++i)
// 	{
// 	  ptr_p = new Point(1,2,3,4,54,6);
// 	  ptr_g->append(ptr_p);
// 	}
//       ptr_db->append(ptr_g);
//       ptr_g = new Gesture();
//     }  

//   //cout<<*ptr_db;
//   ptr_db->save();
//   //cout<<*ptr_g;
// }


