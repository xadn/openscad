/*
 Enable easy piping under Windows(TM) command line.

 We use the 'devenv'(TM) method, which means we have two binary files:

  openscad.com, with IMAGE_SUBSYSTEM_WINDOWS_CUI flag set
  openscad.exe, with IMAGE_SUBSYSTEM_WINDOWS_GUI flag set

 The .com version is a 'wrapper' for the .exe version. If you call
 'openscad' with no extension from a script or shell, the .com version
 is prioritized by the OS and feeds the GUI stdout to the console. We use
 pure C to minimize binary size when cross-compiling (~10kbytes). See Also:

 http://stackoverflow.com/questions/493536/can-one-executable-be-both-a-console-and-gui-app
 http://blogs.msdn.com/b/oldnewthing/archive/2009/01/01/9259142.aspx
 http://blogs.msdn.com/b/junfeng/archive/2004/02/06/68531.aspx
 http://msdn.microsoft.com/en-us/library/aa298534%28v=vs.60%29.aspx
 http://cournape.wordpress.com/2008/07/29/redirecting-stderrstdout-in-cmdexe/
 Open Group popen() documentation
 inkscapec by Jos Hirth work at http://kaioa.com
 Nop Head's OpenSCAD_cl at github.com

 TODO:
 Work with unicode: http://www.i18nguy.com/unicode/c-unicode.html

 TODO: rewrite to use strsafe.h, replace mingw32 with
  mingw-w64's 32-bit system in the build
*/

#define _WIN32_WINNT 0x0500
#include <windows.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXCMDLEN 64000
#define BUFFSIZE 42

void guard_wcscat( wchar_t * s, wchar_t * s2 )
{
	if (lstrlenW(s)+lstrlenW(s2)+2 >= MAXCMDLEN) return;
	wcsncat( s, s2, lstrlenW(s2) );
}

int wmain( int argc, wchar_t * wargv[] )
{
	FILE *cmd_stdout;
	wchar_t cmd[MAXCMDLEN];
	wchar_t buffer[BUFFSIZE];
	wchar_t *fgets_result;
	int eof = 0;
	int pclose_result;
	int i;
	int result = 0;

	guard_wcscat( cmd, L"openscad.exe" );
	for ( i = 1 ; i < argc ; ++i ) {
		guard_wcscat( cmd, L" " );
		guard_wcscat( cmd, wargv[i] );
	}
	guard_wcscat( cmd, L" " );
	guard_wcscat( cmd, L" 2>&1" ); // capture stderr and stdout

	cmd_stdout = _wpopen( cmd, L"rt" );
	if ( cmd_stdout == NULL ) {
		wprintf( L"Error opening _popen for command: %s", cmd );
		_wperror( L"Error message:" );
		return 1;
	}

	while ( !eof )
	{
		fgets_result = fgetws( buffer, BUFFSIZE, cmd_stdout );
		if ( fgets_result == NULL ) {
			if ( ferror( cmd_stdout ) ) {
				wprintf(L"Error reading from stdout of %s\n", cmd);
				result = 1;
			}
			if ( feof( cmd_stdout ) ) {
				eof = 1;
			}
		} else {
			fputws( buffer, stdout );
		}
	}

	pclose_result = _pclose( cmd_stdout );
	if ( pclose_result < 0 ) {
		_wperror(L"Error while closing stdout for command:");
		result = 1;
	}

	return result;
}

#if defined( __MINGW32__ ) || defined( __MINGW64__ )
int main( int argc, char * argv[] )
{
	(void)argc;
	(void)argv;
	int wargc;
	wchar_t * wcmdline = GetCommandLineW();
	wchar_t ** wargv;
	wargv = CommandLineToArgvW( wcmdline, &wargc );
	return wmain(wargc, wargv);
}
#endif

