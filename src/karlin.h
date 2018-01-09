/*
  see karlin.c for details.
*/

#ifndef _karlin_
#define _karlin_

/*
wss:
The following definitions are added to remove dependency on blast.h
*/
#include <stdint.h>
#define int4 int32_t
#define uint4 uint32_t
#define int8 int64_t
#define uint8 uint64_t
#define int2 int16_t
#define uint2 uint16_t

double BlastKarlin_lambda;
double BlastKarlin_K;
double BlastKarlin_H;

void BlastKarlinBlkCalc(double* scoreProbabilities, int4 min, int4 max);

int4 BlastComputeLengthAdjustment(double K,
                             double logK,
                             double alpha_d_lambda,
                             double beta,
                             int4 query_length,
                             uint4 db_length,
                             int4 db_num_seqs,
                             int4 *length_adjustment);
#endif

