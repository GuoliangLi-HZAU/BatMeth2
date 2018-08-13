/*
 *======================================================================================================================
 * File: bfix.c
 *
 * Description:
 *     These routines are used to insert and extract bit fields from an array of characters pointed to by an unsigned
 *     char* pointer.
 *
 *======================================================================================================================
 */

/*
 * Notes:
 *     a. in the following notes any annotation of the form n/m means n for 32 bit systems and m for 64 bit systems.
 *        operation on 32 or 64 bit systems should be transparent.
 *
 *     b. normally, no checking for reasonable argument values is done.  if BFIX_DEBUG is defined then reasonable
 *        values for bit_offset and bit_len are checked.  a reasonable value for bit_len is dependent on the starting
 *        bit_offset value as explained in note b.  a check if both BFIX_BIG_ENDIAN and BFIX_LITTLE_ENDIAN are
 *        defined at the same time is made.  see the code for BFIX_DEBUG inclusions.
 *
 *     c. bit_len should be <=32/64 to 25/57. it depends on the value of bit_offset.  the method always uses a memmove()
 *        of 4/8 bytes to a long temporary storage in which logical operations can be performed.  this means that the
 *        bit_len can be at most 4/8 bytes, but in a case in which the start bit is not at the beginning of a byte,
 *        then the bit_len can only extend until the end of the 4/8'th byte.  if the start bit is the last bit of a
 *        byte this will limit bit_len to 25/57 bits - the last bit of the first byte plus the next 3/7 bytes.
 *        a bit_len of zero is ok - with bfx() returning 0 and bfi() inserting nothing.
 *
 *     d. and bit_offset+bit_len should not overrun the array.
 *
 *     e. value should not be too long to fit into the bit field.  if it is, the high order bits in front of the low
 *        order bit_len bits will be truncated.
 *
 *     f. all bit_len bits will be set and no bit outside the bit field will be changed.
 *
 *     g. value may be negative and prefix 2,s complement sign bits are truncated to fit into bit_len bits.
 *
 *     h. 4/8 bytes are always read from the unsigned char array, modified and then written back.  this means that if
 *        you set the last bit of the array, then the next 3/7 bytes will be read and written back, thus seemingly
 *        overrunning the array.  this may or may not be a problem.  if the 4/8 bytes does not overrun the array
 *        then no bits beyond the end of the array will be changed.  if the 4/8 bytes does overrun the array
 *        some provision must be made to deal with this possibility.  the array could be padded by 3/7 extra bytes.
 *
 *     i. there are three ways to handle endianness.  by default endianness is determined at run time.  this, however,
 *        introduces some run time penalty.  if BFIX_BIG_ENDIAN or BFIX_LITTLE_ENDIAN(but not both) is defined in
 *        bfix.h then that particular behavior will be compiled in.
 */

/*
 *----------------------------------------------------------------------------------------------------------------------
 * include files
 *----------------------------------------------------------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>

#include "bfix.h"

/*
 *----------------------------------------------------------------------------------------------------------------------
 * constants
 *----------------------------------------------------------------------------------------------------------------------
 */

/* comment in one of these for compiled in endianness - default to run time endian check */
/* #define BFIX_BIG_ENDIAN */
/* #define BFIX_LITTLE_ENDIAN */

/*
 *----------------------------------------------------------------------------------------------------------------------
 * data types
 *----------------------------------------------------------------------------------------------------------------------
 */

/*
 *----------------------------------------------------------------------------------------------------------------------
 * global variables
 *----------------------------------------------------------------------------------------------------------------------
 */

/*
 *----------------------------------------------------------------------------------------------------------------------
 * local function prototypes
 *----------------------------------------------------------------------------------------------------------------------
 */

/*
 *----------------------------------------------------------------------------------------------------------------------
 * functions
 *----------------------------------------------------------------------------------------------------------------------
 */




/*
 *----------------------------------------------------------------------------------------------------------------------
 * bfi:
 *     bit field insertion
 *
 * Arguments:
 *     Bit field insert:
 *         bfi(cptr, bit_offset, bit_len, value);
 *
 *             cptr       - pointer to unsigned char array
 *             bit_offset - bit offset(starting from 1) from the start of
 *                          the char array
 *             bit_len    - bit length of field to insert
 *             value      - value to be inserted
 *
 * Returns:
 *     None if !BFIX_DEBUG
 *     int 1 if fail and 0 on success
 *----------------------------------------------------------------------------------------------------------------------
 */

#ifdef BFIX_DEBUG
    int
#else
    void
#endif
bfi(
    unsigned char *cptr,
    unsigned long bit_offset,
    unsigned long bit_len,
    long value)
{
    /* machine dependencies */
    const unsigned int BITS_PER_BYTE = 8;
    unsigned int BYTES_PER_LONG;
    unsigned int BITS_PER_LONG;

    unsigned long l;

    #if !defined(BFIX_BIG_ENDIAN) && !defined(BFIX_LITTLE_ENDIAN)
    unsigned long m;
    #endif

    #if !defined(BFIX_BIG_ENDIAN)
    unsigned int i, j, size;
    unsigned char* c = (unsigned char*)&l;
    unsigned char tmp;
    #endif

    unsigned long mask;
    unsigned long byte_offset;
    unsigned long pre_shift;
    unsigned long post_shift;


    BYTES_PER_LONG = sizeof(unsigned long);
    BITS_PER_LONG = BYTES_PER_LONG * BITS_PER_BYTE;
    for (i=0 ; i<BYTES_PER_LONG ; i++)
    {
       ((unsigned char *)&mask)[i] = 0xff;
    }

    #ifdef BFIX_DEBUG
        #if defined(BFIX_BIG_ENDIAN) && defined(BFIX_LITTLE_ENDIAN)
            #error bfi(): Both BFIX_BIG_ENDIAN and BFIX_LITTLE_ENDIAN defined at the same time.
        #endif

        if (bit_offset < 1)
        {
            fprintf(stderr, "bfi: arg #2, bit_offset = %ld is < 1.\n", bit_offset);
            return 1;
        }

        if (bit_len < 0)
        {
            fprintf(stderr, "bfi: arg #3, bit_len = %ld is < 0.\n", bit_len);
            return 1;
        }
    #endif

    /*
     * calculate byte offset(first byte=0) of start of
     * BYTES_PER_LONG bytes containing bit field
     */
    byte_offset = (bit_offset-1)/BITS_PER_BYTE;

    /*
     * calculate how many bits to shift bit field left
     * to clear bits above bit field
     */
    pre_shift = bit_offset - byte_offset*BITS_PER_BYTE - 1;

    #ifdef BFIX_DEBUG
        if (bit_len > (BITS_PER_LONG-pre_shift))
        {
            fprintf(stderr, "bfi: arg #3, bit_len = %ld to long.\n", bit_len);
            return 1;
        }
    #endif

    /*
     * calculate how many bits to shift bit field left
     * to clear bits below bit field
     */
    post_shift = BITS_PER_LONG - bit_len - pre_shift;

    /*
     * move bit field into position over bits to set in l
     * corresponding to correct bits in cptr
     */
    value <<= post_shift;

    /* zero out mask bits after bit field */
    mask <<= post_shift;

    /* zero out mask bits before bit field */
    mask <<= pre_shift;
    mask >>= pre_shift;

    /* zero out value bits before and after bit field */
    value &= mask;

    /* move BYTES_PER_LONG bytes to tmp storage */
    memmove((unsigned char *)&l, &cptr[byte_offset], BYTES_PER_LONG);

    #if defined(BFIX_BIG_ENDIAN)
    #elif defined(BFIX_LITTLE_ENDIAN)
        size = sizeof(l);
        for (i=0;i < size/2; i++)
        {
            j = size - i - 1;
            tmp = c[i];
            c[i] = c[j];
            c[j] = tmp;
        }
    #else
        m = l;
        l = 1;
        if (c[0])
        {
            /* little endian */
            l = m;
            size = sizeof(l);
            for (i=0;i < size/2; i++)
            {
                j = size - i - 1;
                tmp = c[i];
                c[i] = c[j];
                c[j] = tmp;
            }
            for (i=0;i < size; i++)
            {
            }
        }
        else
        {
            /* big endian */
            l = m;
        }
    #endif

    /* zero out bit field bits and then or value bits into them */
    l = (l & (~mask)) | value;

    #if defined(BFIX_BIG_ENDIAN)
    #elif defined(BFIX_LITTLE_ENDIAN)
        size = sizeof(l);
        for (i=0;i < size/2; i++)
        {
            j = size - i - 1;
            tmp = c[i];
            c[i] = c[j];
            c[j] = tmp;
        }
    #else
        m = l;
        l = 1;
        if (c[0])
        {
            /* little endian */
            l = m;
            size = sizeof(l);
            for (i=0;i < size/2; i++)
            {
                j = size - i - 1;
                tmp = c[i];
                c[i] = c[j];
                c[j] = tmp;
            }
            for (i=0;i < size; i++)
            {
            }
        }
        else
        {
            /* big endian */
            l = m;
        }
    #endif

    /* move tmp storage back to cptr array */
    memmove(&cptr[byte_offset], (unsigned char *)&l, BYTES_PER_LONG);

    #ifdef BFIX_DEBUG
        return 0;
    #endif
}



/*
 *------------------------------------------------------------------------------
 * bfx:
 *     bit field extraction
 *
 * Arguments:
 *     Bit field extract:
 *         l = bfx(cptr, bit_offset, bit_len);
 *
 *             cptr       - pointer to unsigned char array
 *             bit_offset - bit offset(starting from 1) from the start of
 *                          the char array
 *             bit_len    - bit length of field to extract
 *
 * Returns:
 *     on success unsigned long extracted bit field, 0 on failure
 *------------------------------------------------------------------------------
 */

    unsigned long
bfx(
    const unsigned char *cptr,
    unsigned long bit_offset,
    unsigned long bit_len)
{
    /* machine dependencies */
    const unsigned int BITS_PER_BYTE = 8;
    unsigned int BYTES_PER_LONG;
    unsigned int BITS_PER_LONG;

    unsigned long l;

    #if !defined(BFIX_BIG_ENDIAN) && !defined(BFIX_LITTLE_ENDIAN)
    unsigned long m;
    #endif

    #if !defined(BFIX_BIG_ENDIAN)
    unsigned int i, j, size;
    unsigned char* c = (unsigned char*)&l;
    unsigned char tmp;
    #endif

    unsigned long byte_offset;
    unsigned long left_shift;
    unsigned long right_shift;


    BYTES_PER_LONG = sizeof(unsigned long);
    BITS_PER_LONG = BYTES_PER_LONG * BITS_PER_BYTE;

    #ifdef BFIX_DEBUG
        #if defined(BFIX_BIG_ENDIAN) && defined(BFIX_LITTLE_ENDIAN)
            #error bfx(): Both BFIX_BIG_ENDIAN and BFIX_LITTLE_ENDIAN defined at the same time.
        #endif

        if (bit_offset < 1)
        {
            fprintf(stderr, "bfx: arg #2, bit_offset < 1.\n");
            return 0;
        }

        if (bit_len < 0)
        {
            fprintf(stderr, "bfx: arg #3, bit_len < 0.\n");
            return 0;
        }
    #endif

    if (bit_len == 0)
    {
       return 0;
    }

    /*
     * calculate byte offset(first byte=0) of start of
     * BYTES_PER_LONG bytes containing bit field
     */
    byte_offset = (bit_offset-1)/BITS_PER_BYTE;

    /*
     * calculate how many bits to shift bit field left
     * to clear bits above bit field
     */
    left_shift = bit_offset - byte_offset*BITS_PER_BYTE - 1;

    #ifdef BFIX_DEBUG
        if (bit_len > (BITS_PER_LONG-left_shift))
        {
            fprintf(stderr, "bfx: arg #3, bit_len to long.\n");
            return 0;
        }
    #endif

    /*
     * calculate how many bits to shift bit field right
     * to right justify bit field
     */
    right_shift = BITS_PER_LONG - bit_len;

    /* move BYTES_PER_LONG bytes to tmp storage */
    #if defined(BFIX_BIG_ENDIAN)
        memmove((unsigned char *)&l, &cptr[byte_offset], BYTES_PER_LONG);
    #elif defined(BFIX_LITTLE_ENDIAN)
        memmove((unsigned char *)&l, &cptr[byte_offset], BYTES_PER_LONG);
        size = sizeof(l);
        for (i=0;i < size/2; i++)
        {
            j = size - i - 1;
            tmp = c[i];
            c[i] = c[j];
            c[j] = tmp;
        }
    #else
        m = l;
        l = 1;
        if (c[0])
        {
            /* little endian */
            l = m;
            memmove((unsigned char *)&l, &cptr[byte_offset], BYTES_PER_LONG);
            size = sizeof(l);
            for (i=0;i < size/2; i++)
            {
                j = size - i - 1;
                tmp = c[i];
                c[i] = c[j];
                c[j] = tmp;
            }
        }
        else
        {
            /* big endian */
            l = m;
            memmove((unsigned char *)&l, &cptr[byte_offset], BYTES_PER_LONG);
        }
    #endif

    /*
     * clear bits above bit field, right justify bit
     * field, and return this value
     */
    return (l << left_shift) >> right_shift;
}

