#ifndef __FASTSW__
#define __FASTSW__
#include "limits.h"
#include "assert.h"
#include "string.h"
#include "common-l.h"
#include "global-l.h"
#include "math.h"
#include "filters-l.h"

Alignment Fast_SW(char *Ref,char *Read,char *Quality,int OFF,char Sign,char* rawRef,READ & RawR);
Alignment Fast_SWX(char Sign,READ & RawR,char* rawRef,char *Ref,char *Read,int OFF);
#endif
