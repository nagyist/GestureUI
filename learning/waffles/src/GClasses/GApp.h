/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GAPP_H__
#define __GAPP_H__

#include <stddef.h>
#ifdef WIN32
#	include <windows.h>
#else
#	include <termios.h>
#endif

namespace GClasses {

typedef void (*DaemonMainFunc)(void* pArg);

/// Contains some generally useful functions for launching applications
class GApp
{
public:
	/// returns the process ID of the daemon to the parent process.
	/// The daemon process will exit (rather than return) when it finishes.
	/// Note: this doesn't freopen stdout or stderr. You should do that.
	static int launchDaemon(DaemonMainFunc pDaemonMain, void* pArg);

	/// Returns the full name (including path) of the executing application.
	/// Returns the length of the returned string, or -1 if it failed.
	/// If clipFilename is true, it will omit the filename and just return
	/// the folder in which the application resides.
	static int appPath(char* pBuf, size_t len, bool clipFilename);

	/// Executes the specified system command. Output is directed to
	/// the console. (If you want to hide the output or direct it to a
	/// file, use SystemExecute.) The child app runs inside a security
	/// sandbox (at least on Windows, I'm not totally sure about Linux),
	/// so you can't necessarily access everything the parent process
	/// could access. Throws if there is an
	/// error executing the command. Otherwise, returns the exit code
	/// of the child process. "show" does nothing on Linux.
	static int systemCall(const char* szCommand, bool wait, bool show);

	/// Executes the specified system command. (The application inherits
	/// the full environment and privileges from this app. Uses the 
	/// szCommand specifies both the application and the parameters.
	/// blocks until the process exits if (and only if) wait is true.
	/// output of the application is written to the files specified by
	/// szStdOutFilename and szStdErrFilename. These filenames can be
	/// NULL, in which case the output is absorbed (/dev/null).
	/// (If you want the output to go to the console, use SystemCall.)
	/// Throws if there is an error executing the command. Otherwise,
	/// returns the exit code of the application.
	static int systemExecute(const char* szCommand, bool wait, const char* szStdOutFilename, const char* szStdErrFilename);

	/// If you're having trouble with bogus floating point values
	/// (like NAN or INF), this method will cause an exception to
	/// be thrown when overflow, divide-by-zero, or a NAN condition
	/// occurs.
	static void enableFloatingPointExceptions();
};

typedef void (*sighandler_t)(int);

/// Temporarily handles certain signals. (When this object is destroyed, it puts all the signal
/// handlers back the way they were.) Periodically call "check" to see if a signal has occurred.
class GSignalHandler
{
public:
#ifndef WIN32
	sighandler_t m_prevSigInt;
	sighandler_t m_prevSigTerm;
	sighandler_t m_prevSigPipe;
	sighandler_t m_prevSigSegV;
#endif
	int m_gotSignal;

	GSignalHandler();
	~GSignalHandler();

	/// Call this periodically. Returns 0 if no signal has occurred. Otherwise, returns the number of the signal.
	int check();

	/// You can call this to simulate a signal.
	void onSignal(int sig);
};


/// This class provides a non-blocking method for reading characters from
/// stdin. (If there are no characters ready in stdin, it immediately
/// returns '\0'.) The constructor sets flags on the console so that it
/// passes characters to the stream immediately (instead of when Enter is
/// pressed), and so that it doesn't echo the keys, and it makes stdin
/// non-blocking. The destructor puts all those things back the way they were.
class GPassiveConsole
{
#	ifdef WIN32
protected:
	DWORD m_oldMode;
	HANDLE m_hStdin;
	INPUT_RECORD m_inputRecord;
public:
	GPassiveConsole(bool echo);
	~GPassiveConsole();
	char getChar();
#	else
protected:
	struct termios m_old;
	int m_stdin;
	int m_oldStreamFlags;
public:
	GPassiveConsole(bool echo);
	~GPassiveConsole();
	char getChar();
#	endif
};




/// Parses command-line args and provides methods
/// to conveniently process them.
class GArgReader
{
	int m_argc;
	char** m_argv;
	int m_argPos;

public:
	/// Pass the args that are passed in to main
	GArgReader(int argc, char* argv[]);

	/// Returns the current position--that is, the argument number.
	int get_pos();

	/// Sets the current position
	void set_pos(int n);

	/// Returns the current arg (without advancing)
	const char* peek();

	/// Returns the current arg as a string, and advances.
	/// Throws an exception if the end of the args was already
	/// reached before this call.
	const char* pop_string();

	/// Returns the current arg as a uint, and advances.
	/// Throws an exception if the end of the args was already
	/// reached before this call.
	int pop_uint();

	/// Returns the current arg as a double, and advances.
	/// Throws an exception if the end of the args was already
	/// reached before this call.
	double pop_double();

	/// If the current arg matches flagName, advances and returns true.
	/// Otherwise, returns false
	bool if_pop(const char* flagName);

	/// Returns the number of remaining args
	int size();

	/// Returns true of the next arg begins with '-'
	bool next_is_flag();
};

} // namespace GClasses

#endif // __GAPP_H__
