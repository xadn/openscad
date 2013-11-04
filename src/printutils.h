#ifndef PRINTUTILS_H_
#define PRINTUTILS_H_

#include <string>
#include <list>
#include <iostream>
#include <boost/format.hpp>

typedef void (OutputHandlerFunc)(const std::string &msg, void *userdata);
extern OutputHandlerFunc *outputhandler;
extern void *outputhandler_data;
extern bool debug;

void set_output_handler(OutputHandlerFunc *newhandler, void *userdata);

extern std::list<std::string> print_messages_stack;
void print_messages_push();
void print_messages_pop();

void PRINT(const std::string &msg);
#define PRINTB(_fmt, _arg) do { PRINT(str(boost::format(_fmt) % _arg)); } while (0)

void PRINT_NOCACHE(const std::string &msg);
#define PRINTB_NOCACHE(_fmt, _arg) do { PRINT_NOCACHE(str(boost::format(_fmt) % _arg)); } while (0)

void PRINTDEBUG(const std::string &msg);
#define PRINTD(_arg) do { PRINTDEBUG(_arg); } while (0)
#define PRINTDB(_fmt, _arg) do { PRINTDEBUG(str(boost::format(_fmt) % _arg)); } while (0)

void PRINT_CONTEXT(const class Context *ctx, const class Module *mod, const class ModuleInstantiation *inst);

std::string two_digit_exp_format( std::string doublestr );
std::string two_digit_exp_format( double x );


#endif
