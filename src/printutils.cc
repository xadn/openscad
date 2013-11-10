#include "printutils.h"
#include <sstream>
#include <stdio.h>

std::list<std::string> print_messages_stack;
OutputHandlerFunc *outputhandler = NULL;
void *outputhandler_data = NULL;
bool OpenSCAD::debug = false;

void set_output_handler(OutputHandlerFunc *newhandler, void *userdata)
{
	outputhandler = newhandler;
	outputhandler_data = userdata;
}

void print_messages_push()
{
	print_messages_stack.push_back(std::string());
}

void print_messages_pop()
{
	std::string msg = print_messages_stack.back();
	print_messages_stack.pop_back();
	if (print_messages_stack.size() > 0 && !msg.empty()) {
		if (!print_messages_stack.back().empty()) {
			print_messages_stack.back() += "\n";
		}
		print_messages_stack.back() += msg;
	}
}

void PRINT(const std::string &msg)
{
	if (msg.empty()) return;
	if (print_messages_stack.size() > 0) {
		if (!print_messages_stack.back().empty()) {
			print_messages_stack.back() += "\n";
		}
		print_messages_stack.back() += msg;
	}
	PRINT_NOCACHE(msg);
}

void PRINT_NOCACHE(const std::string &msg)
{
	if (msg.empty()) return;
	if (!outputhandler) {
		fprintf(stderr, "%s\n", msg.c_str());
		fflush(stderr);
	} else {
		outputhandler(msg, outputhandler_data);
	}
}

void PRINTDEBUG(const std::string &msg)
{
	// what if.. we could at cmdline say --debug=CGAL_Nef_polyhedron.cc
	// and it would only dump DEBUG from that file?
	// is that possible? even more, what about regex case-insense a str
	if (OpenSCAD::debug) {
		fprintf(stderr,("DEBUG: " + msg + "\n").c_str() );
		//PRINT( "DEBUG: " + msg );
	}
}

std::string two_digit_exp_format( std::string doublestr )
{
#ifdef _WIN32
	size_t exppos = doublestr.find('e');
	if ( exppos != std::string::npos) {
		exppos += 2;
		if ( doublestr[exppos] == '0' ) doublestr.erase(exppos,1);
	}
#endif
	return doublestr;
}

std::string two_digit_exp_format( double x )
{
	std::stringstream s;
	s << x;
	return two_digit_exp_format( s.str() );
}
