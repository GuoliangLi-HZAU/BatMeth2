/*

   Timing.c		Measuring Program running time

   This module contains functions for measuring program running time.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef RUSAGE
#include <sys/resource.h>
#else
#ifdef TIME_BY_CLOCK
#include <time.h>
#else
#include <sys/time.h>
#endif
#endif

#include "Timing.h"


double setStartTime() {

#ifdef RUSAGE

	double usertime, systime;
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);

	usertime = (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1000000.0;
	systime = (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1000000.0;

	return(usertime + systime);

#else
#ifdef TIME_BY_CLOCK

	return (double)clock() / (double)CLOCKS_PER_SEC;

#else 

	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double)tp.tv_sec + (double)tp.tv_usec / (double)1000000;

#endif
#endif

}

double getElapsedTime(double startTime) {

#ifdef RUSAGE

	double usertime, systime;
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);

	usertime = (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1000000.0;
	systime = (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1000000.0;

	return (usertime + systime) - startTime;

#else
#ifdef TIME_BY_CLOCK

	return (double)clock() / (double)CLOCKS_PER_SEC - startTime;

#else 

	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double)tp.tv_sec + (double)tp.tv_usec / (double)1000000 - startTime;

#endif
#endif

}

void printElapsedTime(FILE *file, const int printHour, const int printMin, const int printSec, 
					  const int secNumberOfDecimal, const double seconds) {

	printElapsedTimeNoNewLine(file, printHour, printMin, printSec, 0, secNumberOfDecimal, seconds);
	fprintf(file, "\n");

}

void printElapsedTimeNoNewLine(FILE *file, const int printHour, const int printMin, const int printSec, 
					  const int secMinPrintLength, const int secNumberOfDecimal, const double seconds) {

	int hour, min;
	double sec;
	char secondDisplay[8] = "%0.0f s";
	
	#ifdef DEBUG
	if (printHour && !printMin && printSec) {
		fprintf(stderr, "printElapsedTime(): Cannot skip minute only!\n");
		exit(1);
	}
	if (secNumberOfDecimal > 9) {
		fprintf(stderr, "printElapsedTime(): secNumberOfDecimal > 9!\n");
		exit(1);
	}
	#endif

	secondDisplay[1] = secondDisplay[1] + (char)secMinPrintLength;
	secondDisplay[3] = secondDisplay[3] + (char)secNumberOfDecimal;

	sec = seconds;
	min = (int)(seconds / 60);
	if (!printSec && printMin) {
		if (seconds - min * 60 >= 30) {
			min++;
		}
	}
	if (printMin) {
		sec -= min * 60;
	}

	hour = min / 60;
	if (!printMin) {
		min = hour * 60;
		if (min >= 30) {
			hour++;
		}
	}
	if (printHour) {
		min -= hour * 60;
	}

	if (printHour) {
        fprintf(file, "%d h ", hour);
	}
	if (printMin) {
		fprintf(file, "%d m ", min);
	}
	if (printSec) {
		fprintf(file, secondDisplay, sec);
	}
	
}

