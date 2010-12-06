/*
 *	wiiuse
 *
 *	Written By:
 *		Michael Laforest	< para >
 *		Email: < thepara (--AT--) g m a i l [--DOT--] com >
 *
 *	Copyright 2006-2007
 *
 *	This file is part of wiiuse.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	$Header$
 *
 */

/**
 *	@file
 *
 *	@brief Example using the wiiuse API.
 *
 *	This file is an example of how to use the wiiuse library.
 */

#include <stdio.h>
#include <stdlib.h>

// we plan to support all the platforms
// including windows and MacOS
#ifndef WIN32
#include <unistd.h>
#endif

#include "wiiuse.h"
#include "database.h"

// change the number of wiimotes here
#define MAX_WIIMOTES  2

typedef enum 
  {
    GESTURE_RECORDING,
    GESTURE_RECOGNIZING,
    IDLE
  } State;

GestureDB* database;

/**
 *	@brief Callback that handles an event.
 *
 *	@param wm		Pointer to a wiimote_t structure.
 *
 *	This function is called automatically by the wiiuse library when an
 *	event occurs on the specified wiimote.
 */
void handle_event(struct wiimote_t* wm) 
{
  // the event handler should be implemented as a state machine
  //static bool is_gesture_sensing = false;
  static State state = IDLE;
  static bool is_training = true; /// default being training
  static Gesture* gesture = new Gesture();

  wiiuse_motion_sensing(wm, 1);
  wiiuse_set_ir(wm, 1);
      
  // switch (state)
  //   {
  //   case IDLE:
  //     break;
  //   case GESTURE_RECORDING:
  //     printf("STATE: recording\n");
  //     break;
  //   case GESTURE_RECOGNIZING:
  //     printf("STATE: detecting\n");
  //     break;
  //   }

  // decide if we are training or not
  if ( IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_A) )
    {
      printf("\t\t Learning Stage Changed\t");
      if ( is_training )
	  printf("from Trainig --> Testing\n");
      else
	printf("from Testing --> Training\n");

      is_training = !is_training;
    }

  // save the current database to the system
  else if ( IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_HOME) )
    {
      printf("\t\t saving the database to disk");
      database->save();
    }
  
  // perform state transition
  else if ( IS_JUST_PRESSED(wm, WIIMOTE_BUTTON_B) )
    {
      if ( state == IDLE )
	{
	  // toggle changes the state of wumbling
	  wiiuse_toggle_rumble(wm);
	  printf("gesture recording started:\n");

	  state = GESTURE_RECORDING;	  
	}

      else if ( state == GESTURE_RECORDING )  
	{
	  wiiuse_toggle_rumble(wm);

	  printf("gesture recording ended\n");
	  printf("Which gesture is that? (1: Push/Pull 2:Waving \n");
	  state = GESTURE_RECOGNIZING;	  
	}
    }
  
  if ( state == GESTURE_RECORDING )  
    {
      printf("Wiimote data recorded");
    
      // rotation
      printf("Roll  = %f [%f]\n", wm->orient.roll, wm->orient.a_roll);
      printf("Pitch = %f [%f]\n", wm->orient.pitch, wm->orient.a_pitch);
      printf("Yaw   = %f\n", wm->orient.yaw);
      
      // acceleration
      printf("Graviational Force x = %f\n", wm->gforce.x);
      printf("Graviational Force y = %f\n", wm->gforce.y);
      printf("Graviational Force z = %f\n", wm->gforce.z);

      // ir tracking
      // at most four could be tracked by a single wiimote
      for (int i=0; i < 4; ++i) 
      	{
      	  if (wm->ir.dot[i].visible)
      	    printf("IR source %i: (%u, %u)\n", i, wm->ir.dot[i].x, wm->ir.dot[i].y);
      	}
      
      printf("IR cursor: (%u, %u)\n", wm->ir.x, wm->ir.y);
      printf("IR z distance: %f\n", wm->ir.z);
      printf("\n\n");
      
      /// add a new point to the gesture sequence
      gesture->append(new Point(wm->ir.x,
				wm->ir.y,
				wm->ir.z,
				wm->gforce.x,
				wm->gforce.y,
				wm->gforce.z,
				wm->orient.pitch,
				wm->orient.roll,
				wm->orient.yaw));     
    }

  /// recognize the gesture just recorded
  else if ( state == GESTURE_RECOGNIZING )
    {                  
      // in traing mode, user will be prompted to recognize the gesture
      if ( is_training )
	{	  
	  if ( IS_PRESSED(wm, WIIMOTE_BUTTON_ONE) )
	    {
	      printf("//////////////  Training ////////////\n");
	      //printf("\t\t Gesture Push/Pull Recorded\n");
	      printf("\t\t Gesture FOUR Recorded\n");
	      // get the gesture back to the gesture database
	      //gesture->setType(PUSH);
	      gesture->setType(FOUR);
	      database->append(gesture);	      
	      gesture = new Gesture();

	      state = IDLE;
	    }
	  else if ( IS_PRESSED(wm, WIIMOTE_BUTTON_TWO) )
	    {
	      printf("//////////////  Training ////////////\n");
	      //printf("\t\t Gesture Wave Recorded\n");	      
	      printf("\t\t Gesture NINE Recorded\n");	      
	      // get the gesture back to the gesture database
	      gesture->setType(NINE);
	      database->append(gesture);
	      gesture = new Gesture();

	      state = IDLE;
	    }	  
	}

      
      // in testing mode, the program will give the prediction
      else 
	{
	  printf("//////////////  Testing ////////////\n");
	  //TODO: perform simple testing
	  //printf("\t\t I can't decide which gesture that is\n");
	  printf("\t\t Prediction: \n");
	  GestureType pred_t = database->kNN(gesture, 7);
	  printf("\t\t\t %s\n", ges2str(pred_t));
	  
	  // switch ( database->kNN(gesture, 7) )
	  //   {
	  //   case PUSH:
	  //     printf("\t\t\t PUSH\n");
	  //     break;
	  //   case WAVE:
	  //     printf("\t\t\t WAVE\n");
	  //     break;
	  //   case CIRCLE:
	  //     printf("\t\t\t CIRCLE\n");
	  //     break;
	  //   case SQUARE:
	  //     printf("\t\t\t SQUARE\n");
	  //     break;
	  //   case FOUR:
	  //     printf("\t\t\t FOUR\n");
	  //     break;
	  //   case NINE:
	  //     printf("\t\t\t NINE\n");
	  //     break;
	  //   default:
	  //     printf("\t\t\t UNDECIDED\n");
	  //     break;
	  //   }
	  
	  delete gesture;
	  gesture = new Gesture();
	  state = IDLE;
	}
    }      
}


/**
 *	@brief Callback that handles a read event.
 *
 *	@param wm		Pointer to a wiimote_t structure.
 *	@param data		Pointer to the filled data block.
 *	@param len		Length in bytes of the data block.
 *
 *	This function is called automatically by the wiiuse library when
 *	the wiimote has returned the full data requested by a previous
 *	call to wiiuse_read_data().
 *
 *	You can read data on the wiimote, such as Mii data, if
 *	you know the offset address and the length.
 *
 *	The \a data pointer was specified on the call to wiiuse_read_data().
 *	At the time of this function being called, it is not safe to deallocate
 *	this buffer.
 */
void handle_read(struct wiimote_t* wm, byte* data, unsigned short len) 
{
  int i = 0;
  
  printf("\n\n--- DATA READ [wiimote id %i] ---\n", wm->unid);
  printf("finished read of size %i\n", len);
  for (; i < len; ++i) 
    {
      if (!(i%16))
	printf("\n");
      printf("%x ", data[i]);
    }
  printf("\n\n");
}


/**
 *	@brief Callback that handles a controller status event.
 *
 *	@param wm				Pointer to a wiimote_t structure.
 *	@param attachment		Is there an attachment? (1 for yes, 0 for no)
 *	@param speaker			Is the speaker enabled? (1 for yes, 0 for no)
 *	@param ir				Is the IR support enabled? (1 for yes, 0 for no)
 *	@param led				What LEDs are lit.
 *	@param battery_level	Battery level, between 0.0 (0%) and 1.0 (100%).
 *
 *	This occurs when either the controller status changed
 *	or the controller status was requested explicitly by
 *	wiiuse_status().
 *
 *	One reason the status can change is if the nunchuk was
 *	inserted or removed from the expansion port.
 */
void handle_ctrl_status(struct wiimote_t* wm) 
{
  printf("\n\n--- CONTROLLER STATUS [wiimote id %i] ---\n", wm->unid);
  
  printf("attachment:      %i\n", wm->exp.type);
  printf("speaker:         %i\n", WIIUSE_USING_SPEAKER(wm));
  printf("ir:              %i\n", WIIUSE_USING_IR(wm));
  printf("leds:            %i %i %i %i\n", WIIUSE_IS_LED_SET(wm, 1), WIIUSE_IS_LED_SET(wm, 2), WIIUSE_IS_LED_SET(wm, 3), WIIUSE_IS_LED_SET(wm, 4));
  printf("battery:         %f %%\n", wm->battery_level);
}


/**
 *	@brief Callback that handles a disconnection event.
 *
 *	@param wm				Pointer to a wiimote_t structure.
 *
 *	This can happen if the POWER button is pressed, or
 *	if the connection is interrupted.
 */
void handle_disconnect(wiimote* wm) 
{
  printf("\n\n--- DISCONNECTED [wiimote id %i] ---\n", wm->unid);
}


/** 
 * 
 */
void test(struct wiimote_t* wm, byte* data, unsigned short len) 
{
  printf("test: %i [%x %x %x %x]\n", len, data[0], data[1], data[2], data[3]);
}


/**
 *	@brief main()
 *
 *	Connect to up to two wiimotes and print any events
 *	that occur on either device.
 */
int main(int argc, char** argv) 
{
  wiimote** wiimotes;
  int num_found, num_connected;

  wiimotes =  wiiuse_init(MAX_WIIMOTES);

  num_found = wiiuse_find(wiimotes, MAX_WIIMOTES, 5);
  if (!num_found) 
    {
      printf ("No wiimotes found.");
      return 0;
    }

  num_connected = wiiuse_connect(wiimotes, MAX_WIIMOTES);
  if (num_connected)
    printf("Num_Connected to %i wiimotes (of %i found).\n", num_connected, num_found);
  else 
    {
      printf("Failed to connect to any wiimote.\n");
      return 0;
    }  
 

  for (int i=0; i<num_connected; ++i)
    {
      //Philip: TODO: add support for muliple wiimote synergy
      wiiuse_set_leds(wiimotes[i], WIIMOTE_LED_1);
      wiiuse_rumble(wiimotes[i], 1);
    }
  // wiiuse_set_leds(wiimotes[1], WIIMOTE_LED_2);
  // wiiuse_set_leds(wiimotes[2], WIIMOTE_LED_3);
  // wiiuse_set_leds(wiimotes[3], WIIMOTE_LED_4);

  //wiiuse_rumble(wiimotes[1], 1);

#ifndef WIN32
  usleep(200000);
#else
  Sleep(200);
#endif

  wiiuse_rumble(wiimotes[0], 0);
  wiiuse_status(wiimotes[0]);

  ///////////////////////////////////////////////////////
  /// User Input
  ///////////////////////////////////////////////////////
 
  // load the gesture database
  database = new GestureDB("data/gesture_data.db");
  database->load();  
  cout<<"\tGesture Data Library loaded: "<<database->size()<<" transactions in total"<<endl;

  if ( database->size() == 0 )
    {
      cerr<<"Error: database contains no data, please collect some training data first"<<endl;
    }
  //cout<<*database;
  
  printf("Press B to record a gesture\n");
  while (1) 
    {
      if (wiiuse_poll(wiimotes, num_connected)) 
	{
	  /*
	   *	This happens if something happened on any wiimote.
	   *	So go through each one and check if anything happened.
	   */
	  for (int i=0; i < num_connected; ++i) 
	    {
	      switch (wiimotes[i]->event) 
		{
		case WIIUSE_EVENT:
		  /* a generic event occured */
		  handle_event(wiimotes[i]);
		  break;
		  
		case WIIUSE_STATUS:
		  /* a status event occured */
		  handle_ctrl_status(wiimotes[i]);
		  break;

		case WIIUSE_DISCONNECT:
		case WIIUSE_UNEXPECTED_DISCONNECT:
		  /* the wiimote disconnected */
		  handle_disconnect(wiimotes[i]);
		  database->save();
		  break;

		case WIIUSE_READ_DATA:
		  break;

		default:
		  break;
		}
	    }
	}
    }

  /*
   *	Disconnect the wiimotes
   */
  database->save();
  //cout<<*database;
  wiiuse_cleanup(wiimotes, num_connected);

  return 0;
}
