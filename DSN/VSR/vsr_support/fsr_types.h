/********************************************************************************
* FILE NAME: fsr_types.h
* CHANGE HISTORY: $Log: fsr_types.h,v $
* CHANGE HISTORY: Revision 1.1  2006/04/20 00:20:43  jerry
* CHANGE HISTORY: Initial version.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2004/09/08 23:14:00  ceg
* CHANGE HISTORY: initial add from vsr 2.5
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.4  2002/07/02 17:44:59  andre
* CHANGE HISTORY: cleanup
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2002/02/15 22:32:49  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2002/01/14 16:27:14  andre
* CHANGE HISTORY: mon data is now a union of fsr and schan segments
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1.1.1  2001/10/24 00:29:06  andre
* CHANGE HISTORY: create VSR repository
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.4  2001/08/17 00:31:36  andre
* CHANGE HISTORY: reorder schan_mode_t enum so that enum position == num_schans
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2001/08/08 20:40:53  andre
* CHANGE HISTORY: -adjust Makefile to create mdspecs.h from mdspecs.dat
* CHANGE HISTORY: -modified fsr_monitor to work better with MDS stuff
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2000/09/20 16:22:07  andre
* CHANGE HISTORY: mv fsr_id_t to fsr_types.h
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2000/01/28 22:46:49  miyatake
* CHANGE HISTORY: initial add: created from remnants of cmd_prdx.h
* CHANGE HISTORY:
********************************************************************************/

#ifndef FSR_TYPES_H      
#define FSR_TYPES_H      


/***  DEFINE CONSTANTS  ***/
   
#define NUM_FSR_CHAN 4
#define SCHAN_PER_CHAN 4
#define NUM_FSR_SCHAN (NUM_FSR_CHAN * SCHAN_PER_CHAN)

    
/* enum for rsr id */
typedef enum
{
  FSR_ID_1A = 0,
  FSR_ID_1B,
  FSR_ID_2A,
  FSR_ID_2B,
  FSR_ID_3A,
  FSR_ID_3B,
  FSR_ID_4A,
  FSR_ID_4B,
  FSR_ID_NA
}
fsr_id_t;


/* enum for band in rmt and tlm predicts files */
typedef enum
{
   RF_BAND_NA,
   RF_BAND_L,
   RF_BAND_S,
   RF_BAND_X,
   RF_BAND_KA
} rf_band_t;


/* enum for tracking mode in rmt and tlm predicts files */
typedef enum
{
   TRK_MODE_NA,
   TRK_MODE_1_WAY,
   TRK_MODE_2_WAY,
   TRK_MODE_3_WAY  
} trk_mode_t;


/* enum for rsr state */
typedef enum
{
   IDLE = 0,
   CNF,
   RUN
}
state_t;


/* enum for rsr schan dplr */
typedef enum
{
   SCHAN_DPLR_NA = 0,
   SCHAN_DPLR_NONE = 0,  /* no doppler correction on sfro */
   SCHAN_DPLR_1WAY,      /* doppler from 1way prdx applied to sfro */
   SCHAN_DPLR_CARR,      /* doppler from current track mode applied to sfro */
}
schan_dplr_t;                              

                              
#endif /* FSR_TYPES_H */  
