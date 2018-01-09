#ifndef __COMMAND__
#define __COMMAND__
#include "global.h"
void Read_INI(char* Config_File,unsigned & MAXCOUNT,FMFILES & F,BATPARAMETERS & BP);
void Parse_Command_line(int argc, char* argv[],char* & GENOME,unsigned & MAXCOUNT,FMFILES & F,BATPARAMETERS & BP);
#endif
