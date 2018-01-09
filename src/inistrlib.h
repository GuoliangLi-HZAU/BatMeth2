
/*-------------------------------------------------------------------------*/
/**
  @file     inistrlib.h
  @author   N. Devillard
  @date     Jan 2001
  @version  $Revision: 1.3 $
  @brief    Various string handling routines to complement the C lib.

  This modules adds a few complementary string routines usually missing
  in the standard C library.
*/
/*--------------------------------------------------------------------------*/

/*
	$Id: inistrlib.h,v 1.3 2001/10/19 08:31:41 ndevilla Exp $
	$Author: ndevilla $
	$Date: 2001/10/19 08:31:41 $
	$Revision: 1.3 $
*/

#ifndef _INISTRLIB_H_
#define _INISTRLIB_H_

/*---------------------------------------------------------------------------
   								Includes
 ---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------
  							Function codes
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
  @brief    Convert a string to lowercase.
  @param    s   String to convert.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing a lowercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */
/*--------------------------------------------------------------------------*/
char * inistrlwc(char * s);

/*-------------------------------------------------------------------------*/
/**
  @brief    Convert a string to uppercase.
  @param    s   String to convert.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string
  containing an uppercased version of the input string. Do not free
  or modify the returned string! Since the returned string is statically
  allocated, it will be modified at each function call (not re-entrant).
 */
/*--------------------------------------------------------------------------*/
char * inistrupc(char * s);

/*-------------------------------------------------------------------------*/
/**
  @brief    Skip blanks until the first non-blank character.
  @param    s   String to parse.
  @return   Pointer to char inside given string.

  This function returns a pointer to the first non-blank character in the
  given string.
 */
/*--------------------------------------------------------------------------*/
char * inistrskp(char * s);

/*-------------------------------------------------------------------------*/
/**
  @brief    Remove blanks at the end of a string.
  @param    s   String to parse.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end of the string have been removed.
  Do not free or modify the returned string! Since the returned string
  is statically allocated, it will be modified at each function call
  (not re-entrant).
 */
/*--------------------------------------------------------------------------*/
char * inistrcrop(char * s);

/*-------------------------------------------------------------------------*/
/**
  @brief    Remove blanks at the beginning and the end of a string.
  @param    s   String to parse.
  @return   ptr to statically allocated string.

  This function returns a pointer to a statically allocated string,
  which is identical to the input string, except that all blank
  characters at the end and the beg. of the string have been removed.
  Do not free or modify the returned string! Since the returned string
  is statically allocated, it will be modified at each function call
  (not re-entrant).
 */
/*--------------------------------------------------------------------------*/
char * inistrstrip(char * s) ;

/*-------------------------------------------------------------------------*/
/**
  @brief    Duplicate a string by allocating memory and copy the contents.
  @param    s   String to duplicate.
  @return   ptr to dynamically allocated string.

  This function returns a pointer to a dynamically allocated string,
  which is identical to the input string.
  Free the returned string when it is no longer used.
 */
/*--------------------------------------------------------------------------*/
char *inistrdup(const char *s);

#endif
