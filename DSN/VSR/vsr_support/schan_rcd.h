/********************************************************************************
* FILE NAME: schan_rcd.h
* CHANGE HISTORY: $Log: schan_rcd.h,v $
* CHANGE HISTORY: Revision 1.1  2006/04/26 22:23:06  spr
* CHANGE HISTORY: support for fft.c
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2006/04/20 18:22:35  jerry
* CHANGE HISTORY: Initial version.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2005/09/12 21:29:40  ceg
* CHANGE HISTORY: undo last change
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2005/09/08 20:15:37  soriano
* CHANGE HISTORY: add #include <std_types.h>
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2004/09/08 23:14:39  ceg
* CHANGE HISTORY: initial add from vsr 2.5
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2002/05/29 22:10:46  andre
* CHANGE HISTORY: use SCHAN_RCD_HDR_SIZE macro
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2002/05/24 23:02:22  andre
* CHANGE HISTORY: remove if_source from rcd - its in the vsr log file
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2002/05/23 21:34:38  andre
* CHANGE HISTORY: switch to new vsr data record format
* CHANGE HISTORY:
********************************************************************************/
#if !defined(SCHAN_RCD_H)
#define SCHAN_RCD_H


/***  SYSTEM INCLUDE FILES  ***/

#include <stddef.h>


/***  PROJECT INCLUDE FILES  ***/

#include "std_defines.h"
#include "fsr_types.h"
/* #include <date_time.h> */


/***  PROJECT DEFINES  ***/

/*******************************************************************************
* SET_SCHAN_RCD_VERSION_NUMBER sets the version number of the schan_rcd data
* structure.  The RCS Revision string ($Revision: 1.1 $) is parsed for the
* {V}.{R} numbers at the end.  The version number is set to 1000*(V-1) + R,
* which simplifies to R for V=1.
*******************************************************************************/
#define SET_SCHAN_RCD_VERSION_NUMBER(version)	do		\
{							\
   char *s = "$Revision: 1.1 $";			\
							\
   /* skip leading non-digits */			\
   while ( *s != '\0'  &&  !isdigit((int)*s) ) s++;	\
							\
   version = 1000*(strtol(s, &s, 10) - 1);		\
   version += strtol(s+1, &s, 10);			\
} while (0);


/*******************************************************************************
* The following macros, in particular VME_SAMPLE_I and VME_SAMPLE_Q,
* are used to extract the I and Q components of the samples contained
* in the "vme_data" array.
*******************************************************************************/

/* number of samples in a single vme_data */
#define	SAMPLES_PER_DATUM(bits_per_sample)		(16/(bits_per_sample))

/* array index of vme_data holding sample "n" */
#define	SAMPLE_TO_DATA_INDEX(n, bits_per_sample)	\
   ((n)/SAMPLES_PER_DATUM(bits_per_sample))

/* bit offset of sample "n" within a single vme_data */
#define	SAMPLE_TO_BIT_OFFSET(n, bits_per_sample)	\
   (((n)%SAMPLES_PER_DATUM(bits_per_sample))*(bits_per_sample))

#define	_BIT_FIELD_(data, start, width)			\
   ((((unsigned long)(data))>>(start))&((1<<(width))-1))

#define	_VME_BITS_(offset, vme_data, n, bits_per_sample)	\
   _BIT_FIELD_(							\
      (vme_data)[SAMPLE_TO_DATA_INDEX((n), (bits_per_sample))],	\
      SAMPLE_TO_BIT_OFFSET((n), (bits_per_sample)) + (offset),	\
      (bits_per_sample))

#define	VME_BITS_I(vme_data, n, bits_per_sample)	\
   _VME_BITS_(0, (vme_data), (n), (bits_per_sample))

#define	VME_BITS_Q(vme_data, n, bits_per_sample)	\
   _VME_BITS_(16, (vme_data), (n), (bits_per_sample))

#define	BITS_TO_SAMPLE(bits, bits_per_sample)		\
   ((bits) -						\
      (((bits) < (1<<((bits_per_sample)-1))) ?		\
      0 : (1<<(bits_per_sample))))

#define	VME_SAMPLE_I(vme_data, n, bits_per_sample)			\
   BITS_TO_SAMPLE(VME_BITS_I((vme_data), (n), (bits_per_sample)),	\
      (bits_per_sample))

#define	VME_SAMPLE_Q(vme_data, n, bits_per_sample)			\
   BITS_TO_SAMPLE(VME_BITS_Q((vme_data), (n), (bits_per_sample)),	\
      (bits_per_sample))


#define SCHAN_RCD_HDR_SIZE offsetof(schan_rcd_t, data)
#define SCHAN_RCD_DATA_SIZE(sample_rate, bits_per_sample) \
   ((sample_rate)*1000*(bits_per_sample)/4)
#define SCHAN_RCD_NUM_DATA(sample_rate, bits_per_sample) \
   ((sample_rate)*1000*(bits_per_sample)/16)


#define ULONG_ASCII(A,B,C,D)	(((A)<<24)|((B)<<16)|((C)<<8)|(D))
#define SCHAN_RCD_LABEL ULONG_ASCII('V','S','R','D')


/***  DEFINE PROJECT TYPES  ***/

/*******************************************************************************
* Important Notice:  This interface was developed using SUN Ultra 	
* Sparc II CPUs.  Multibyte data objects, including floating point	
* values, are defined to have the big endian byte order which is	
* native to that processor.			
*******************************************************************************/

/*******************************************************************************
* To use schan_rcd_t allocate a block of memory large enough
* to hold an anc_data_t and the vme_data[] array.  Then use 
* a schan_rcd_t pointer to read/write into the allocated memory.
* The vme_data[] array contains one second of data, its size
* depends on the schan bandwidth and numbits per sample.
*******************************************************************************/

typedef struct
{
   U32 label;                       /* four byte label always set to 'VSRD' */
   U32 record_length;               /* the length of the record in bytes */
   U16 record_version;              /* version number (increment with changes to this file */
   U16 software_version;            /* version number (increment with changes to this file */
   U16 spc_id;                      /* station id - 10, 40, 60, 21 */
   U16 vsr_id;                      /* vsr1a, vsr1b ... from enum */                          
   U16 schan_id;                    /* subchannel id 0,1,2,3 */
   U16 bits_per_sample;             /* number of bits per sample - 1, 2, 4, 8, or 16 */         
   U16 sample_rate;                 /* number of samples per second in kilo-samples per second */           
   U16 sdplr;                       /* algoithm used to apply doppler correction to sfro */
   U16 prdx_dss_id;                 /* DSS number */                                          
   U16 prdx_sc_id;                  /* spacecraft number listed in predicts file */                                   
   U16 prdx_pass_num;               /* pass number listed in predicts file */                                     
   char prdx_uplink_band[2];        /* uplink RF band */                                      
   char prdx_downlink_band[2];      /* downlink RF band */                                    
   U16 trk_mode;                    /* tracking mode */                                       
   U16 uplink_dss_id;              /* uplink dss id when in 3 way */                         
   U16 ddc_lo;                      /* DDC LO frequency in MHz */                             
   U16 rf_to_if_lo;                 /* RF to IF LO frequency in MHz */                        
   U16 data_error;                  /* hw error flag, dma error or num_samples error, 0 ==> no errors */
   U16 year;                        /* time tag - year */
   U16 doy;                         /* time tag - day of year */
   I32 sec;                         /* time tag - second of day */
   I32 data_time_offset;            /* in nano-seconds */
   double frov;                     /* in Hz */
   double fro;                      /* in Hz */
   double frr;                      /* in Hz/sec */
   double sfro;                     /* in Hz */
   double rf_freq;                  /* in Hz */
   double schan_accum_phase;        /* number of accumulated full turns */
   double schan_phase_poly[4];      /* coefficients for baseband phase polynomial */
   char schan_label[16];            /* ascii label - indicates source of data - carrier, subcarrier, etc */
}
schan_rcd_hdr_t;


typedef struct
{
   schan_rcd_hdr_t hdr;             /* record header */
   U32 data[1];                 /* first entry in variably sized array */
}
schan_rcd_t;

#endif /* SCHAN_RCD_H */


/*       1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

/********************************************************************************
********************************************************************************/

