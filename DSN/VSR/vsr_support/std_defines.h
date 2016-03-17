/********************************************************************************
* FILE NAME: std_defines.h
* CHANGE HISTORY: $Log: std_defines.h,v $
* CHANGE HISTORY: Revision 1.2  2009-06-26 16:11:13  ceg
* CHANGE HISTORY: more 64bit compiler upgrade stuff
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2006-04-20 00:20:43  jerry
* CHANGE HISTORY: Initial version.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2005/09/09 17:51:53  ceg
* CHANGE HISTORY: restore removed typedefs, wrapped with "ifndef linux"
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2005/08/31 19:41:21  soriano
* CHANGE HISTORY: remove typedefs which conflict with Linux typedefs
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2004/09/08 23:14:01  ceg
* CHANGE HISTORY: initial add from vsr 2.5
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1.1.1  2001/10/24 00:29:06  andre
* CHANGE HISTORY: create VSR repository
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.9  2000/06/12 23:44:57  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.8  2000/03/08 23:00:05  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.7  2000/02/16 22:19:14  andre
* CHANGE HISTORY: clean up rs/doc, adjust makefiles, /rs/rsr now compiles
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.6  2000/02/15 19:09:57  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.5  2000/01/28 01:35:55  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.4  2000/01/21 23:16:14  andre
* CHANGE HISTORY: remove dependence on fsp.h fsp_com.h fsx_monitor.h
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2000/01/20 19:27:40  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
********************************************************************************/

#ifndef STD_DEFINES_H
#define STD_DEFINES_H

/* This file defines some constants that are used throughout the 
 * code.  If the constant has not been previously defined it is
 * defined.  If it has been previously defined for the same 
 * value as we want nothing is done.  If the constant has been
 * previously defined to value different than we want then the
 * constant is redefined without using an #undef so that a compiler
 * warning is generated to notify us of a potenial problem.  */

#if !defined(TRUE)
#define TRUE 1
#endif
#if TRUE != 1
#define TRUE 1
#endif

#if !defined(FALSE)
#define FALSE 0
#endif 
#if FALSE != 0
#define FALSE 0
#endif 

#if !defined(YES)
#define YES 1
#endif 
#if YES != 1
#define YES 1
#endif 

#if !defined(NO)
#define NO 0
#endif 
#if NO != 0
#define NO 0
#endif 

#if !defined(ON)
#define ON 1
#endif 
#if ON != 1
#define ON 1
#endif 

#if !defined(OFF)
#define OFF 0
#endif 
#if OFF != 0
#define OFF 0
#endif 

#if !defined(OK)
#define OK 0
#endif 
#if OK != 0
#define OK 0
#endif 

#if !defined(ERROR)
#define ERROR -1
#endif 
#if ERROR != -1
#define ERROR -1
#endif 


/* gcc complains if we test NULL so it better be right */
#if !defined(NULL)
#define NULL ((void *) 0)
#endif 

#if !defined(STR_SZ)
#define STR_SZ 256
#endif 
#if STR_SZ != 256
#define STR_SZ 256
#endif 


/* some handy typedefs */

typedef int BOOLEAN; 
#ifdef linux
typedef __int8_t i8, I8, sint8;
typedef __uint8_t u8, U8, uint8;
typedef __int16_t i16, I16, sint16;
typedef __uint16_t u16, U16, uint16;
typedef __int32_t i32, I32, sint32;
typedef __uint32_t u32, U32, uint32;
#endif


/* some handy macros */

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) (((a) >= 0) ? (a) : -(a))
#endif

#ifndef SQR
#define SQR(a) (a * a)
#endif

#ifndef DIM
#define DIM(a) (sizeof(a)/sizeof(a[0]))
#endif

#ifndef IN_RANGE
#define IN_RANGE(n,lo,hi) ((lo) <= (n) && (n) <= (hi))
#endif

#ifndef SWAP
#define SWAP(a,b) ({ long tmp; tmp = (a), (a) = (b), (b) = tmp})
#endif



#endif /* STD_DEFINES_H */  

/*       1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

