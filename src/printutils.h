#ifndef PRINTUTILS_H_
#define PRINTUTILS_H_

#include <string>
#include <list>
#include <iostream>
#include <boost/format.hpp>

typedef void (OutputHandlerFunc)(const std::string &msg, void *userdata);
extern OutputHandlerFunc *outputhandler;
extern void *outputhandler_data;
namespace OpenSCAD { extern std::string debug; }

void set_output_handler(OutputHandlerFunc *newhandler, void *userdata);

extern std::list<std::string> print_messages_stack;
void print_messages_push();
void print_messages_pop();

void PRINT(const std::string &msg);
#define PRINTB(_fmt, _arg) do { PRINT(str(boost::format(_fmt) % _arg)); } while (0)

void PRINT_NOCACHE(const std::string &msg);
#define PRINTB_NOCACHE(_fmt, _arg) do { PRINT_NOCACHE(str(boost::format(_fmt) % _arg)); } while (0)

// usage:
// PRINTD(" Outputting 3 points: ");
// PRINTDB("point0, point1, point2: %s %s %s", p0 % p1 % p2 );

void PRINTDEBUG(const std::string &filename,const std::string &msg);
#define PRINTD(_arg) do { PRINTDEBUG(std::string(__FILE__),_arg); } while (0)
#define PRINTDB(_fmt, _arg) do { PRINTDEBUG(std::string(__FILE__),str(boost::format(_fmt) % _arg)); } while (0)

void PRINT_CONTEXT(const class Context *ctx, const class Module *mod, const class ModuleInstantiation *inst);

std::string two_digit_exp_format( std::string doublestr );
std::string two_digit_exp_format( double x );


#endif
