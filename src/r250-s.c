/* r250.c	the r250 uniform random number algorithm

		Kirkpatrick, S., and E. Stoll, 1981; "A Very Fast
		Shift-Register Sequence Random Number Generator",
		Journal of Computational Physics, V.40

		also:

		see W.L. Maier, DDJ May 1991



*/

#include <limits.h>
#include "r250-s.h"

// static functions
static unsigned int randlcg();


#define MSB          0x80000000L
#define ALL_BITS     0xffffffffL
#define HALF_RANGE   0x40000000L
#define STEP         7
#define BITS         32

static unsigned int r250_buffer[ 250 ];
static int r250_index;

static unsigned int randlcg(int sd)       /* returns a random unsigned integer */
{
		static int quotient1  = LONG_MAX / 16807L;
		static int remainder1 = LONG_MAX % 16807L;

        if ( sd <= quotient1 )
                sd = (sd * 16807L) % LONG_MAX;
        else
        {
                int high_part = sd / quotient1;
                int low_part  = sd % quotient1;

                int test = 16807L * low_part - remainder1 * high_part;

                if ( test > 0 )
                        sd = test;
                else
                        sd = test + LONG_MAX;

        }

        return sd;
}

void r250_init(int sd)
{
	int j, k;
	unsigned int mask, msb;

	r250_index = 0;
	for (j = 0; j < 250; j++)      /* fill r250 buffer with BITS-1 bit values */
		sd = r250_buffer[j] = randlcg(sd);


	for (j = 0; j < 250; j++)	/* set some MSBs to 1 */
		if ( (sd=randlcg(sd)) > HALF_RANGE )
			r250_buffer[j] |= MSB;


	msb = MSB;	        /* turn on diagonal bit */
	mask = ALL_BITS;	/* turn off the leftmost bits */

	for (j = 0; j < BITS; j++)
	{
		k = STEP * j + 3;	/* select a word to operate on */
		r250_buffer[k] &= mask; /* turn off bits left of the diagonal */
		r250_buffer[k] |= msb;	/* turn on the diagonal bit */
		mask >>= 1;
		msb  >>= 1;
	}

}

unsigned int r250()		/* returns a random unsigned integer */
{
	register int	j;
	register unsigned int new_rand;

	if ( r250_index >= 147 )
		j = r250_index - 147;	/* wrap pointer around */
	else
		j = r250_index + 103;

	new_rand = r250_buffer[ r250_index ] ^ r250_buffer[ j ];
	r250_buffer[ r250_index ] = new_rand;

	if ( r250_index >= 249 )	/* increment pointer for next time */
		r250_index = 0;
	else
		r250_index++;

	return new_rand;

}


double dr250()		/* returns a random double in range 0..1 */
{
	register int	j;
	register unsigned int new_rand;

	if ( r250_index >= 147 )
		j = r250_index - 147;	/* wrap pointer around */
	else
		j = r250_index + 103;

	new_rand = r250_buffer[ r250_index ] ^ r250_buffer[ j ];
	r250_buffer[ r250_index ] = new_rand;

	if ( r250_index >= 249 )	/* increment pointer for next time */
		r250_index = 0;
	else
		r250_index++;

	return (double)new_rand / ALL_BITS;

}

