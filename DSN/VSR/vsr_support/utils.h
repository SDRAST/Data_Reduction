/********************************************************************************
* FILE NAME: utils.h
* CHANGE HISTORY: $Log: utils.h,v $
* CHANGE HISTORY: Revision 1.1  2006/04/20 15:55:45  jerry
* CHANGE HISTORY: Initial version.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.5  2006/01/27 16:42:02  ceg
* CHANGE HISTORY: add "daemon_socket" and "daemon_accept", which allow an INETD spawned
* CHANGE HISTORY: process to merge with a previously spawned instance, producing a
* CHANGE HISTORY: single process with multiple connections.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.4  2005/12/09 18:34:48  ceg
* CHANGE HISTORY: add logging to command handler.
* CHANGE HISTORY: accomodate rsp_read -> dpc_read move.
* CHANGE HISTORY: add terminate_on_close.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2005/08/08 17:30:24  andre
* CHANGE HISTORY: add get_xil_id function
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2005/03/15 15:46:39  ceg
* CHANGE HISTORY: strip out VSRisms
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1  2005/03/15 15:38:57  ceg
* CHANGE HISTORY: from VSR 2.5
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.4  2002/06/05 23:16:57  andre
* CHANGE HISTORY: add tpgm command
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2002/02/15 22:35:39  andre
* CHANGE HISTORY: RSR to VSR changes
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2001/12/04 19:21:27  andre
* CHANGE HISTORY: add set_schan get_schan functions
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.1.1.1  2001/10/24 00:29:14  andre
* CHANGE HISTORY: create VSR repository
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.8  2001/09/13 23:19:10  andre
* CHANGE HISTORY: add parse_fa func for use by clients
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.7  2001/08/08 20:42:55  andre
* CHANGE HISTORY: removed some apparently useless functions
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.6  2000/11/29 21:47:56  andre
* CHANGE HISTORY: add get_short_rf_band
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.5  2000/11/20 19:05:55  miyatake
* CHANGE HISTORY: add monitor data publishing conversions.
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.4  2000/07/12 21:40:36  andre
* CHANGE HISTORY: *** empty log message ***
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.3  2000/03/27 21:20:34  andre
* CHANGE HISTORY: added schan config stuff
* CHANGE HISTORY:
* CHANGE HISTORY: Revision 1.2  2000/02/23 20:59:15  andre
* CHANGE HISTORY: fsr_cmnd and rt_cmnd can now process cnf, vfy, run, halt, uvfy, hi,
* CHANGE HISTORY: set, att and way commands
* CHANGE HISTORY:
********************************************************************************/
#if !defined(UTILS_H)
#define UTILS_H

/***  SYSTEM INCLUDE FILES  ***/

/***  PROJECT INCLUDE FILES  ***/

#include "std_defines.h"


/***  DEFINE CONSTANTS  ***/

/***  FUNCTION PROTOTYPES  ***/

int readn(int fd, char * ptr, int nbytes);
int writen(int fd, char * ptr, int nbytes);

int set_bool(int *bool, char *string);
char *get_bool(int bool, char *string);

/* from "terminate_on_close.c" */
void terminate_on_close( int fd );

/* from "log.c" */
void init_log( int argc, char **argv );

/* from "daemon.c" */
int daemon_socket( char *name );
int daemon_accept( int sock );

#endif /* UTILS_H */

/*       1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

/********************************************************************************
********************************************************************************/




