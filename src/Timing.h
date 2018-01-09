/*

   Timing.h		Measuring Program running time

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

#ifndef __TIMING_H__
#define __TIMING_H__


double setStartTime();
double getElapsedTime(double startTime);
void printElapsedTime(FILE *file, const int printHour, const int printMin, const int printSec, 
					  const int secNumberOfDecimal, const double seconds);
void printElapsedTimeNoNewLine(FILE *file, const int printHour, const int printMin, const int printSec, 
					  const int secMinPrintLength, const int secNumberOfDecimal, const double seconds);

#endif

