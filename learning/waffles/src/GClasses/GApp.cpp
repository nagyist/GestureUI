/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GApp.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef WIN32
#	include <windows.h>
#	include <direct.h>
#	include <float.h>
#else // WIN32
#	define _GNU_SOURCE 1
#	include <fenv.h>
#	include <unistd.h>
#	include <signal.h>
#	include <sys/wait.h>
//#	include <termios.h>
//#	include <fcntl.h>
#	ifndef __linux__
#		include <mach-o/dyld.h>
#	endif
#endif // !WIN32
#include "GMacros.h"
#include "GHolders.h"
#include "GFile.h"
#include <errno.h>
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>

using namespace GClasses;
using std::cout;
using std::string;


/*static*/ int GApp::launchDaemon(DaemonMainFunc pDaemonMain, void* pArg)
{
#ifdef WIN32
	// Windows isn't POSIX compliant and it has its own process system that
	// isn't really friendly to launching daemons.  You're supposed to create
	// a "service", but I don't know how to do that (and I'm too lazy to learn
	// something that can't generalize off a proprietary platform) so let's
	// just launch it like a normal app and be happy with that.
	pDaemonMain(pArg);
	return 0;
#else // WIN32

	// Fork the first time
	int firstPid = fork();
	if(firstPid < 0)
		throw "Error forking (the first time) in GApp::LaunchDaemon";
	if(firstPid)
		return firstPid;

	// Fork the second time
	int secondPid = fork();
	if(secondPid < 0)
		throw "Error forking (the second time) in GApp::LaunchDaemon";
	if(secondPid)
		return secondPid;

	// Drop my process group leader and become my own process group leader
	// (so the process isn't terminated when the group leader is killed)
	setsid();

	// Set the file creation mask. (I don't know why we do this.)
	umask(0);

	// Get off any mounted drives so that they can be unmounted without
	// killing the daemon
	if(chdir("/") != 0)
	{
	}

	// Launch the daemon
	pDaemonMain(pArg);

	exit(0);
#endif // !WIN32
}

/*static*/ int GApp::appPath(char* pBuf, size_t len, bool clipFilename)
{
#ifdef WIN32
	int bytes = GetModuleFileName(NULL, pBuf, (unsigned int)len);
	if(bytes == 0)
		return -1;
	else
	{
		if(clipFilename)
		{
			while(bytes > 0 && pBuf[bytes - 1] != '\\' && pBuf[bytes - 1] != '/')
				pBuf[--bytes] = '\0';
		}
		return bytes;
	}
#else
#	ifdef __linux__
	std::ostringstream os;
	os << "/proc/" << getpid() << "/exe";
	string tmp = os.str();
	int bytes = MIN((int)readlink(tmp.c_str(), pBuf, len), (int)len - 1);
#	else
	// Darwin
	uint32_t l = len;
	_NSGetExecutablePath(pBuf, &l); // todo: this returns the path to the symlink, not the actual app.
	int bytes = strlen(pBuf);
#	endif
	if(bytes >= 0)
	{
		pBuf[bytes] = '\0';
		if(clipFilename)
		{
			while(bytes > 0 && pBuf[bytes - 1] != '/')
				pBuf[--bytes] = '\0';
		}
	}
	return bytes;
#endif
}

int GApp_measureParamLen(const char* sz)
{
	int len = 0;
	while(true)
	{
		if(*sz == '"')
		{
			len++;
			sz++;
			while(*sz != '"' && *sz != '\0')
			{
				len++;
				sz++;
			}
			if(*sz == '"')
			{
				len++;
				sz++;
			}
			continue;
		}
		else if(*sz == '\'')
		{
			len++;
			sz++;
			while(*sz != '\'' && *sz != '\0')
			{
				len++;
				sz++;
			}
			if(*sz == '\'')
			{
				len++;
				sz++;
			}
			continue;
		}
		else if(*sz <= ' ')
			break;
		len++;
		sz++;
	}
	return len;
}

int GApp_measureWhitespaceLen(const char* sz)
{
	int len = 0;
	while(*sz <= ' ' && *sz != '\0')
	{
		len++;
		sz++;
	}
	return len;
}

int GApp_CountArgs(const char* sz)
{
	int count = 0;
	while(true)
	{
		sz += GApp_measureWhitespaceLen(sz);
		if(*sz == '\0')
			break;
		count++;
		sz += GApp_measureParamLen(sz);
	}
	return count;
}

void GApp_ParseArgs(char* sz, char* argv[], int cap)
{
	int count = 0;
	while(true)
	{
		if(*sz != '\0' && *sz <= ' ')
		{
			*sz = '\0';
			sz++;
		}
		sz += GApp_measureWhitespaceLen(sz);
		if(*sz == '\0')
			break;
		argv[count++] = sz;
		if(count >= cap)
			break;
		sz += GApp_measureParamLen(sz);
	}
	argv[count] = NULL;
}

int GApp::systemCall(const char* szCommand, bool wait, bool show)
{
#ifdef WIN32
	// Parse the args
	GTEMPBUF(char, szCopy, (int)strlen(szCommand) + 1);
	strcpy(szCopy, szCommand);
	int argc = GApp_CountArgs(szCopy);
	if(argc == 0)
		return 0;
	char* argv[3];
	GApp_ParseArgs(szCopy, argv, 2);

	// Call it
	SHELLEXECUTEINFO sei;
	memset(&sei, '\0', sizeof(SHELLEXECUTEINFO));
	sei.cbSize = sizeof(SHELLEXECUTEINFO);
	sei.fMask = SEE_MASK_NOCLOSEPROCESS | SEE_MASK_FLAG_DDEWAIT;
	sei.hwnd = NULL;
	sei.lpVerb = NULL;
	sei.lpFile = argv[0];
	sei.lpParameters = argv[1];
	sei.lpDirectory = NULL;
	sei.nShow = show ? SW_SHOW : SW_HIDE;
	if(!ShellExecuteEx(&sei))
		ThrowError("An error occurred while executing the command \"", argv[0], " ", argv[1], "\"");
	DWORD ret = 0;
	if(wait)
	{
		WaitForSingleObject(sei.hProcess, INFINITE);
		if(!GetExitCodeProcess(sei.hProcess, &ret))
		{
			CloseHandle(sei.hProcess);
			ThrowError("Failed to obtain exit code");
		}
		CloseHandle(sei.hProcess);
	}
	return ret;
#else
	string s = szCommand;
	if(!wait)
		s += " &";
	int status = system(s.c_str());
	if(status == -1)
		ThrowError("Failed to execute command");
	return WEXITSTATUS(status);
#endif
}

int GApp::systemExecute(const char* szCommand, bool wait, const char* szStdOutFilename, const char* szStdErrFilename)
{
#ifdef WIN32
	// Set the bInheritHandle flag so pipe handles are inherited.
	SECURITY_ATTRIBUTES saAttr;
	saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	saAttr.bInheritHandle = TRUE; // pipe handles are inherited
	saAttr.lpSecurityDescriptor = NULL;

	// Create pipes for the child process's stdout and stdin
	HANDLE hChildStdinRd, hChildStdinWr, hChildStdoutRd, hChildStdoutWr, hChildStderrRd, hChildStderrWr;
	if(!CreatePipe(&hChildStderrRd, &hChildStderrWr, &saAttr, 0))
		ThrowError("Failed to create a pipe for the child processes stderr");
	SetHandleInformation(hChildStderrRd, HANDLE_FLAG_INHERIT, 0); // Ensure that the read handle to the child process's pipe for stdout is not inherited
	if(!CreatePipe(&hChildStdoutRd, &hChildStdoutWr, &saAttr, 0))
		ThrowError("Failed to create a pipe for the child processes stdout");
	SetHandleInformation(hChildStdoutRd, HANDLE_FLAG_INHERIT, 0); // Ensure that the read handle to the child process's pipe for stdout is not inherited
	if(!CreatePipe(&hChildStdinRd, &hChildStdinWr, &saAttr, 0))
		ThrowError("Failed to create a pipe for child processes stdin");
	SetHandleInformation( hChildStdinWr, HANDLE_FLAG_INHERIT, 0); // Ensure that the write handle to the child process's pipe for stdin is not inherited.

	// Create the child process.
	STARTUPINFO siStartInfo;
	ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
	siStartInfo.cb = sizeof(STARTUPINFO); 
	siStartInfo.hStdError = hChildStderrWr;
	siStartInfo.hStdOutput = hChildStdoutWr;
	siStartInfo.hStdInput = hChildStdinRd;
	siStartInfo.dwFlags |= STARTF_USESTDHANDLES;
	PROCESS_INFORMATION piProcInfo;
	ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));
	if(!CreateProcess(NULL,
		(LPSTR)szCommand,     // command line
		NULL,          // process security attributes
		NULL,          // primary thread security attributes
		TRUE,          // handles are inherited
		0, /*CREATE_NO_WINDOW*/             // creation flags
		NULL,          // NULL means use parent's environment
		NULL,          // NULL means use parent's current directory
		&siStartInfo,  // STARTUPINFO pointer
		&piProcInfo)) // receives PROCESS_INFORMATION
		ThrowError("Failed to create process");

	DWORD ret = 0;
	if(wait)
	{
		WaitForSingleObject(piProcInfo.hProcess, INFINITE);
		if(!GetExitCodeProcess(piProcInfo.hProcess, &ret))
		{
			CloseHandle(piProcInfo.hProcess);
			CloseHandle(piProcInfo.hThread);
			ThrowError("Failed to obtain exit code");
		}
	}
	CloseHandle(piProcInfo.hProcess);
	CloseHandle(piProcInfo.hThread);

	// Close the child processes' stdin pipe, since we don't really need it
	CloseHandle(hChildStdoutWr);
	CloseHandle(hChildStderrWr);

	// Read the output from the child process and save to file
	if(szStdErrFilename)
	{
		DWORD dwRead;
		CHAR chBuf[4096];
		FILE* pFileOut = fopen(szStdErrFilename, "wb");
		FileHolder hFileOut(pFileOut);
		while(true)
		{
			if(!ReadFile(hChildStderrRd, chBuf, 4096, &dwRead, NULL) || dwRead == 0)
				break;
			if(fwrite(chBuf, dwRead, 1, pFileOut) != 1)
				ThrowError("Failed to write to the stderr file");
		}
	}
	if(szStdOutFilename)
	{
		DWORD dwRead;
		CHAR chBuf[4096];
		FILE* pFileOut = fopen(szStdOutFilename, "wb");
		FileHolder hFileOut(pFileOut);
		while(true)
		{
			if(!ReadFile(hChildStdoutRd, chBuf, 4096, &dwRead, NULL) || dwRead == 0)
				break;
			if(fwrite(chBuf, dwRead, 1, pFileOut) != 1)
				ThrowError("Failed to write to the stdout file");
		}
	}

	// Clean up
	CloseHandle(hChildStdinWr);
	CloseHandle(hChildStdinRd);
	CloseHandle(hChildStdoutRd);
	CloseHandle(hChildStderrRd);

	return ret;
#else
	// Parse the args
	GTEMPBUF(char, szCopy, strlen(szCommand) + 1);
	strcpy(szCopy, szCommand);
	int argc = GApp_CountArgs(szCopy);
	if(argc == 0)
		return 0;
	char* argv[argc + 1];
	GApp_ParseArgs(szCopy, argv, 2147483647);

	// Call it
	pid_t pid = fork();
	if(pid == 0)
	{
		// Redirect the output streams
		FILE* pStdOut = NULL;
		FILE* pStdErr = NULL;
		if(szStdOutFilename)
			pStdOut = freopen(szStdOutFilename, "wb", stdout);
		else
			pStdOut = freopen("/dev/null", "wb", stdout);
		if(szStdErrFilename)
			pStdErr = freopen(szStdErrFilename, "wb", stderr);
		else
			pStdErr = freopen("/dev/null", "wb", stderr);

		// Call the child process
		execvp(argv[0], argv);

		// execvp only returns if there is an error. Otherwise, it replaces the current process.
		cout << "Error calling execvp. errno=" << errno << "\n";
		ThrowError("Error calling execvp. errno=", gformat(errno));
	}
	else if(pid > 0)
	{
		// The calling process
		if(wait)
		{
			int status;
			waitpid(pid, &status, 0);
			if(WIFEXITED(status))
			{
				int ret = WEXITSTATUS(status);
				return ret;
			}
			else if(WIFSIGNALED(status))
				ThrowError("The process was interruped with signal ", gformat(WSTOPSIG(status)));
			else
				ThrowError("The process stopped without exiting, and it cannot be restarted.");
		}
	}
	else
		ThrowError("There was an error forking the process");
	return 0;
#endif
}

// static
void GApp::enableFloatingPointExceptions()
{
#ifdef WIN32
	unsigned int cw = _control87(0, 0) & MCW_EM; // should we use _controlfp instead?
	cw &= ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	_control87(cw,MCW_EM);
#else
#	ifdef DARWIN
	// todo: Anyone know how to do this on Darwin?
#	else
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#	endif
#endif
}


GSignalHandler* g_pSignalHandler = NULL;
#ifndef WIN32
void GApp_onSigInt(int n)
{
	g_pSignalHandler->onSignal(SIGINT);
}

void GApp_onSigTerm(int n)
{
	g_pSignalHandler->onSignal(SIGTERM);
}

void GApp_onSigPipe(int n)
{
	g_pSignalHandler->onSignal(SIGPIPE);
}

void GApp_onSigSegV(int n)
{
	g_pSignalHandler->onSignal(SIGSEGV);
}
#endif // !WIN32

GSignalHandler::GSignalHandler()
{
	if(g_pSignalHandler)
		ThrowError("GSignalHandler is not reentrant, so it cannot be nested");
	g_pSignalHandler = this;
	m_gotSignal = 0;
#ifndef WIN32
	m_prevSigInt = signal(SIGINT, GApp_onSigInt); if(m_prevSigInt == SIG_ERR) m_prevSigInt = SIG_DFL;
	m_prevSigTerm = signal(SIGTERM, GApp_onSigTerm); if(m_prevSigTerm == SIG_ERR) m_prevSigInt = SIG_DFL;
	m_prevSigPipe = signal(SIGPIPE, GApp_onSigPipe); if(m_prevSigPipe == SIG_ERR) m_prevSigInt = SIG_DFL;
	m_prevSigSegV = signal(SIGSEGV, GApp_onSigSegV); if(m_prevSigSegV == SIG_ERR) m_prevSigInt = SIG_DFL;
#endif // !WIN32
}

GSignalHandler::~GSignalHandler()
{
#ifndef WIN32
	signal(SIGINT, m_prevSigInt);
	signal(SIGTERM, m_prevSigTerm);
	signal(SIGPIPE, m_prevSigPipe);
	signal(SIGSEGV, m_prevSigSegV);
#endif // !WIN32
	g_pSignalHandler = NULL;
}

void GSignalHandler::onSignal(int sig)
{
	g_pSignalHandler->m_gotSignal = sig;
}

int GSignalHandler::check()
{
	return m_gotSignal;
}





#ifdef WIN32
GPassiveConsole::GPassiveConsole(bool echo)
{
	m_hStdin = GetStdHandle(STD_INPUT_HANDLE);
	GetConsoleMode(m_hStdin, &m_oldMode);
	DWORD newMode = m_oldMode & (~ENABLE_LINE_INPUT);
	if(!echo)
		newMode &= (~ENABLE_ECHO_INPUT);
	SetConsoleMode(m_hStdin, newMode);
}

GPassiveConsole::~GPassiveConsole()
{
	SetConsoleMode(m_hStdin, m_oldMode);
}

char GPassiveConsole::getChar()
{
	DWORD n;
	while(true)
	{
		if(!PeekConsoleInput(m_hStdin, &m_inputRecord, 1, &n))
			ThrowError("PeekConsoleInput failed");
		if(n == 0)
			return '\0';
		if(m_inputRecord.EventType == KEY_EVENT && m_inputRecord.Event.KeyEvent.bKeyDown && m_inputRecord.Event.KeyEvent.uChar.AsciiChar != 0)
			break;
		if(!ReadConsoleInput(m_hStdin, &m_inputRecord, 1, &n))
			ThrowError("ReadConsoleInput failed");
	}
	char c;
	if(ReadConsole(m_hStdin/*hConsoleInput*/, &c, 1, &n, NULL))
	{
		GAssert(n > 0);
		return c;
	}
	else
	{
		ThrowError("ReadConsole failed");
		return '\0';
	}
}
#else
GPassiveConsole::GPassiveConsole(bool echo)
{
	if(tcgetattr(0, &m_old) < 0)
		ThrowError("Error getting terminal settings");
	struct termios tmp;
	memcpy(&tmp, &m_old, sizeof(struct termios));
	tmp.c_lflag &= ~ICANON;
	if(!echo)
		tmp.c_lflag &= ~ECHO;
	tmp.c_cc[VMIN] = 1;
	tmp.c_cc[VTIME] = 0;
	if(tcsetattr(0, TCSANOW, &tmp) < 0)
		ThrowError("Error setting terminal settings");
	m_stdin = fileno(stdin);
	m_oldStreamFlags = fcntl(m_stdin, F_GETFL, 0);
	if(fcntl(m_stdin, F_SETFL, m_oldStreamFlags | O_NONBLOCK) == -1)
		ThrowError("Error setting stdin to non-blocking");
}

GPassiveConsole::~GPassiveConsole()
{
	if(tcsetattr(0, TCSANOW, &m_old) < 0)
		ThrowError("Error restoring terminal settings");
	if(fcntl(m_stdin, F_SETFL, m_oldStreamFlags) == -1)
		ThrowError("Error restoring stdin flags");
}

char GPassiveConsole::getChar()
{
	char c;
	if(read(m_stdin, &c, 1) > 0)
		return c;
	else
		return '\0';
}
#endif







GArgReader::GArgReader(int argc, char* argv[])
: m_argc(argc), m_argv(argv), m_argPos(0)
{
}

int GArgReader::get_pos()
{
	return m_argPos;
}

void GArgReader::set_pos(int n)
{
	m_argPos = n;
}

const char* GArgReader::peek()
{
	return m_argv[m_argPos];
}

const char* GArgReader::pop_string()
{
	if(m_argPos >= m_argc)
		ThrowError("Unexpected end of arguments");
	return m_argv[m_argPos++];
}

int GArgReader::pop_uint()
{
	const char* str = pop_string();
	for(int i = 0; str[i] != '\0'; i++)
	{
		if(str[i] < '0' || str[i] > '9')
			ThrowError("Expected an unsigned integer value for parameter ", gformat(m_argPos), ". (Got \"", str, "\".)");
	}
	return atoi(str);
}

double GArgReader::pop_double()
{
	const char* str = pop_string();
	return atof(str);
}

bool GArgReader::if_pop(const char* flagName)
{
	if(m_argPos >= m_argc)
		return false;
	if(_stricmp(peek(), flagName) == 0)
	{
		m_argPos++;
		return true;
	}
	else
		return false;
}

int GArgReader::size()
{
	return m_argc - m_argPos;
}

bool GArgReader::next_is_flag()
{
	if(size() == 0)
		return false;	
	return (peek()[0] == '-');
}
