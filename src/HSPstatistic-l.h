/* ==========================================================================
 *
 * HSPstatistic.h
 *
 * author: Wong Swee Seong
 * date: 30-03-2006
 *
 * gapped and ungapped extension for string matching
 * 
 * ==========================================================================*/


#ifndef _HSP_statistic_
#define _HSP_statistic_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define int4 int32_t
#define uint4 uint32_t
#define int8 int64_t
#define uint8 uint64_t
#define int2 int16_t
#define uint2 uint16_t


typedef struct stat__score_block
{
  int match, mismatch, gapOpen, gapExt;
} STAT_scoreBlock;

typedef struct stat__Xdropoffs
{
  double ungapXdropoff, gapXdropoff, gapXdropoffFinal;
} STAT_Xdropoffs;

void initializeHSPstatistic(uint8 databaseSize, uint8 numOfSequences, uint8 minSeqLength, double *dbCharProb,
							uint8 querySize, uint8 numOfContext, double *queryCharProb, int scoringMatrix[16][16]);
void printHSPstatistic(FILE * outFile);
int4 calcUngapCutoffScore();
int4 calcGapCutoffScore();
int4 getUngapXdropoff();
int4 getGapXdropoff();
int4 getGapXdropoffFinal();

int4 stat_ungapNormalized2nominal(double normalizedScore);
double stat_ungapNominal2normalized(int4 nominalScore);
double stat_ungapCalcEvalue(double normalizedScore);
int4 stat_ungapEvalue2nominal(double evalue);
double stat_gapNominal2normalized(int4 nominalScore);
int4 stat_gapNormalized2nominal(double normalizedScore);
double stat_gapCalcEvalue(double normalizedScore);
int4 stat_gapEvalue2nominal(double evalue);
 
#ifdef __cplusplus
}
#endif

#endif
