#include <iostream>
#include "database.h"

int main()
{
  GestureDB db("gesture_data.db");
  db.load();

  //cout<<*db<<endl;
  /// fetch the 7th gesture
  Gesture* gesture = db[7];
  cout<<*gesture<<endl;

  gesture->sample(13);
  cout<<*gesture<<endl;

  //db.save();

  return 0;
}