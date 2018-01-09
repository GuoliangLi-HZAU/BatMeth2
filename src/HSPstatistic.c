/* ==========================================================================
 *
 * HSPstatistic.c
 *
 * author: Wong Swee Seong
 * date: 30-03-2006
 *
 * gapped and ungapped extensions for string matching
 *
 * There are many references to BLAST statistic, a good development reading
 * would from the NCBI BLAST pc version, file doc/scoring.pdf  
 *
 * This files contains fragmented codes extracted and modified from FSA-BLAST 
 * by MIcheal Cameron, and as such the following statement is provided.

   Copyright (c) 2005, Michael Cameron
   Permission to use this code is freely granted under the BSD license agreement,
   provided that this statement is retained.
   
   Statistical functions related to BLAST (ie. converting nominal scores
   to normalized scores and to e-values)
 * ==========================================================================*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "karlin.h"
#include "HSPstatistic.h"

// default cutoff evalue for blastn, used in ungapped extension
// as define in NCBI BLAST, file: include/algo/blast/core/blast_parameters.h
#define CUTOFF_E_BLASTN 0.05 

double stat_log2;
int4 stat_ungapXdropoff;
int4 stat_gapXdropoff;
int4 stat_gapXdropoffFinal;
double stat_expectationValue = 10;
STAT_scoreBlock stat_scoreBlock = {1,-3,5,2};
STAT_Xdropoffs stat_XdropoffsInBits = {20,30,50};

/* ungap parameters */
double stat_ungapK;
double stat_ungapLogK;
double stat_ungapH;
double stat_ungapLambda;

/* gap parameters */
double stat_gapK;
double stat_gapLogK;
double stat_gapH;
double stat_gapLambda;

/* search space parameters */
uint8 stat_databaseSize;
uint8 stat_numOfSequences;
uint8 stat_averageSubjectSize; // average sequence size in database
uint8 stat_minimumSubjectSize; // minimum sequence size in database
int8 stat_adjustmentLength;
uint8 stat_querySize;
uint8 stat_numOfContext;
uint8 stat_effectiveQuerySize;
uint8 stat_effectiveDatabaseSize;
uint8 stat_effectiveSearchSpaceSize;


//enum karlinParam_ENUM {LAMBDA,parmK,parmH};
//int HSPscoreBlock[4] = {1, -3, 5, 2};
/*  HSPstatBlk -
 *  1. ungapXdropOff = 10, default 20 in bits
 *  2. gapXdropOff = 15, default 30 in bits
 *  3. gapXdropOffFinal = 25, default 50 in bits
 *  4. ungapCutoffScore = 15, to be calculated
 *  5. expectationValue = 10
 * */
//float HSPXdropoffInBits[3] = {10, 15, 25};

/* Karlin Parameters extracted from NCBI tool-kit, 
 * file: blast_stat.c
 * Will add more in future
 * Only storing (lambda,k) pair for each pre-defined score value set 
 * (match,miss), we will the need other parameter if
 * we want to calculate score for gapped alignment
 */
/*
int KarlinParamSelect = 0;
const int KarlinParamSet = 3;
const int scoreList[][2] = {{1,-3}, 
                            {1,-4},  
                            {1,-2}};  
const double KarlinParam[][2] = {{1.374, 0.711},
                                 {1.383, 0.738},
                                 {1.28, 0.46}};
*/

/*  forward declaration */
void stat_calcUngapKarlinParameters(double *dbCharProb, double *queryCharProb, int scoringMatrix[16][16]);
double stat_calcLengthAdjust();

/*
  Initialize the statistical paramters for extension 
  Need only to call once, with the database size and number of sequences in 
  the database.
*/
void initializeHSPstatistic(uint8 databaseSize, uint8 numOfSequences, uint8 minSeqLength, double *dbCharProb,
							uint8 querySize, uint8 numOfContext, double *queryCharProb, int scoringMatrix[16][16])
{
 
  stat_databaseSize = databaseSize;
  stat_numOfSequences = numOfSequences;
  stat_minimumSubjectSize = minSeqLength;
  stat_averageSubjectSize = (uint4)(databaseSize /numOfSequences);  

  stat_log2 = log(2.0);

  // Compute for the Karlin-Atschul parameters: ungap K, H and Lambda
  // For DNA, assuming that each character A,T,G,C are equally probable at 0.25
  stat_calcUngapKarlinParameters(dbCharProb, queryCharProb, scoringMatrix); 

  // for blastn, they are the same
  stat_gapK = stat_ungapK;
  stat_gapH = stat_ungapH;
  stat_gapLambda = stat_ungapLambda;
  stat_ungapLogK = log(stat_ungapK);
  stat_gapLogK = stat_ungapLogK;
 
  stat_ungapXdropoff = (int4)ceil(stat_XdropoffsInBits.ungapXdropoff * stat_log2 /
                                  stat_ungapLambda); 
  stat_gapXdropoff = (int4)floor(stat_XdropoffsInBits.gapXdropoff * stat_log2 /
                                 stat_gapLambda); 
  stat_gapXdropoffFinal = (int4)floor(stat_XdropoffsInBits.gapXdropoffFinal * 
                                      stat_log2 / stat_gapLambda); 

  stat_querySize = querySize;
  stat_numOfContext = numOfContext;
  // compute the adjustment length
  // stat_calcLengthAdjust is a down-sized calc for nucleotide only
  stat_adjustmentLength = (int4)floor(stat_calcLengthAdjust());
  
  stat_effectiveQuerySize = stat_querySize - stat_adjustmentLength;
  stat_effectiveDatabaseSize = stat_databaseSize - 
                               stat_numOfSequences * stat_adjustmentLength;
  stat_effectiveSearchSpaceSize = stat_effectiveDatabaseSize * stat_effectiveQuerySize;
}

void stat_calcUngapKarlinParameters(double *dbCharProb, double *queryCharProb, int scoringMatrix[16][16])
{
  int4 scoreProbSize, highest, lowest;
  double * scoreProb;
  int i, j;

  lowest = stat_scoreBlock.mismatch;
  highest = stat_scoreBlock.match;
  scoreProbSize = highest-lowest+1;
  scoreProb = (double *) calloc(scoreProbSize, sizeof(double)); 

  for (i=0; i<16; i++) {
	  for (j=0; j<16; j++) {
		  scoreProb[scoringMatrix[i][j] - lowest] += dbCharProb[i] * queryCharProb[j];
	  }
  }

  // calculate the Karlin parameters, refer to karlin.c
  BlastKarlinBlkCalc(scoreProb-lowest, lowest, highest);
  stat_ungapH = BlastKarlin_H;
  stat_ungapK = BlastKarlin_K;
  stat_ungapLambda = BlastKarlin_lambda;
  
  if (stat_ungapK  <= 0.001)
  {
    fprintf(stderr,"ERROR: unable to calculate the KARLIN statistic for the given scoring scheme.\n");
    exit(-1);
  }  

  free(scoreProb);
}


/*
   Codes extracted and modified from statistic.c in FSA BLAST by Michael Cameron 
*/
// Calculate length adjustment (aka. "effective HSP length") by iteratively
// applying equation from "Local Alignment Statistics" Altschul,Gish
// using ungapped K & H values
double stat_calcLengthAdjust()
{
  int4 count;
  double lengthAdjust;
  double minimumQueryLength;

  // Query can not be shorter than 1/K
  minimumQueryLength = (double)1 / stat_ungapK;

  count = 0;
  lengthAdjust = 0;
  // Calculate length adjustment value
  while (count < 5)
  {
    lengthAdjust = (stat_ungapLogK +
      log(((double)stat_querySize - lengthAdjust) *
                 ((double)stat_databaseSize - (double)stat_numOfSequences *
                  lengthAdjust))) / stat_ungapH;

    // Length adjust can not shorten query to less than minimum query length
    if (lengthAdjust > stat_querySize - minimumQueryLength)
    {
      lengthAdjust = stat_querySize - minimumQueryLength;
    }

    count++;
  }

  return lengthAdjust;
}

int4 calcUngapCutoffScore()
{
  uint8 minSize = stat_minimumSubjectSize;
  int4 gappedScore, ungappedScore;

  if (minSize > stat_querySize)
     minSize = stat_querySize;

  ungappedScore = (int4)ceil(
           (log((double)minSize * (double)stat_numOfContext * ((double)stat_averageSubjectSize / (double)CUTOFF_E_BLASTN)) 
           + stat_ungapLogK) / stat_ungapLambda );
  gappedScore = stat_gapEvalue2nominal(stat_expectationValue);
  
  if (gappedScore > ungappedScore)
     return ungappedScore;
  else
     return gappedScore;
}

int4 calcGapCutoffScore()
{
    return stat_gapEvalue2nominal(stat_expectationValue);
}

int4 getUngapXdropoff()
{
  return stat_ungapXdropoff;
}

int4 getGapXdropoff()
{
  return stat_gapXdropoff;
}

int4 getGapXdropoffFinal()
{
  return stat_gapXdropoffFinal;
}

void printHSPstatistic(FILE * outFile)
{
  int vtmp;
  fprintf(outFile,"Lambda        K        H\n");
  fprintf(outFile,"%10.2f  %10.2f    %10.2f\n", stat_ungapLambda, stat_ungapK, stat_ungapH);
  fprintf(outFile,"Gapped\n");
  fprintf(outFile,"Lambda        K        H\n");
  fprintf(outFile,"%10.2f  %10.2f    %10.2f\n\n", stat_gapLambda, stat_gapK, stat_gapH);
  fprintf(outFile,"Matrix: blastn matrix: %d %d\n", stat_scoreBlock.match, stat_scoreBlock.mismatch );
  fprintf(outFile,"Gap Penalties: Existence: %d, Extension: %d\n", stat_scoreBlock.gapOpen, stat_scoreBlock.gapExt);
  fprintf(outFile,"Length of query: %d\n", stat_querySize);
  fprintf(outFile,"Length of database: %llu\n", stat_databaseSize);
  fprintf(outFile,"Effective HSP length: %d\n", stat_adjustmentLength);
  fprintf(outFile,"Effective length of query: %d\n", stat_effectiveQuerySize);
  fprintf(outFile,"Effective length of database: %llu\n", stat_effectiveDatabaseSize);
  fprintf(outFile,"Effective search space: %llu\n",stat_effectiveSearchSpaceSize);
  fprintf(outFile,"X1: %d (%2.1lf bits)\n",stat_ungapXdropoff, 
       stat_ungapXdropoff * stat_ungapLambda / stat_log2);
  fprintf(outFile,"X2: %d (%2.1lf bits)\n",stat_gapXdropoff,
       stat_gapXdropoff * stat_gapLambda / stat_log2);
  vtmp = calcUngapCutoffScore();
  fprintf(outFile,"S1: %d (%2.1f bits)\n",vtmp,
               stat_ungapNominal2normalized(vtmp));
  vtmp = calcGapCutoffScore();
  fprintf(outFile,"S2: %d (%2.1f bits)\n",vtmp,
               stat_gapNominal2normalized(vtmp));
}

/*--------------------------------------------------------------------------
 * Routines for conversion between evalue and score, as well as from 
   normalise (in bits) to nominal.
 */

// Convert a normalized score to a nominal score
int4 stat_ungapNormalized2nominal(double normalizedScore)
{
  double nominalScore;

  nominalScore = (normalizedScore * stat_log2 + stat_ungapLogK) / stat_ungapLambda;

  return (int4)floor(nominalScore);
}

// Convert an ungapped nominal score to a normalized score
double stat_ungapNominal2normalized(int4 nominalScore)
{
  double normalizedScore;

  normalizedScore = (((double)nominalScore * stat_ungapLambda) - stat_ungapLogK)
                    / stat_log2;

  return normalizedScore;
}

// Calculate the evalue for a given ungapped normalizedScore
double stat_ungapCalcEvalue(double normalizedScore)
{
  return stat_effectiveSearchSpaceSize / pow(2, normalizedScore);
}

// Given an evalue (such as a cutoff) calculate the minimum ungapped nominal score needed to attain it
int4 stat_ungapEvalue2nominal(double evalue)
{
  double normalizedScore;

  normalizedScore = log(stat_effectiveSearchSpaceSize / evalue) / (double)stat_log2;
  return (int4)ceil((stat_log2 * normalizedScore + stat_ungapLogK)
         / stat_ungapLambda);
}

// Convert a gapped nominal score to a normalized score
double stat_gapNominal2normalized(int4 nominalScore)
{
  double normalizedScore;

  normalizedScore = (((double)nominalScore * stat_gapLambda)
                    - stat_gapLogK) / stat_log2;

  return normalizedScore;
}

// Given a normalized score calculate the lowest gapped nominal score needed to achieve it
int4 stat_gapNormalized2nominal(double normalizedScore)
{
  return (int4)ceil(((stat_log2 * normalizedScore + stat_gapLogK)
         / stat_gapLambda) - 0.5);
}
// Calculate the evalue for a given gapped normalizedScore
double stat_gapCalcEvalue(double normalizedScore)
{
  return stat_effectiveSearchSpaceSize / pow(2, normalizedScore);
}

// Given an evalue (such as a cutoff) calculate the minimum gapped nominal score needed to attain it
int4 stat_gapEvalue2nominal(double evalue)
{
  return (int4)ceil(log(stat_gapK * stat_effectiveSearchSpaceSize / evalue) / stat_gapLambda);
}


//---------------------------------------------------------------------------

