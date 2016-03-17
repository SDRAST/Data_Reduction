/* $Id: bls.h,v 1.2 2009-06-26 14:17:11 ceg Exp $ */
#if !defined(BLS_H)
#define BLS_H

#include <sys/types.h>

#define	BLS_PORT	9999

/* I/O to the BLS file system is in multiples of BLS_BLOCK_SIZE bytes */
/* "bls_lseek" offsets will be rounded away from zero.  "bls_write" will */
/* be padded.  "bls_read" will skip forward before next read. */
#define	BLS_BLOCK_SIZE		512

/* round length N up to the next BLS_BLOCK_SIZE boundary */
#define	BLS_LENGTH(N)	(((N) + (BLS_BLOCK_SIZE-1)) & ~(BLS_BLOCK_SIZE-1))

/* how many blocks are required to hold N bytes */
#define	BLS_BLOCKS(N)	(((N) + (BLS_BLOCK_SIZE-1)) / BLS_BLOCK_SIZE)

/* Data blocks in the BLS file system are allocated in fixed size clusters */
#define	BLS_BLOCKS_PER_CLUSTER	128
#define	BLS_CLUSTER_SIZE	(BLS_BLOCKS_PER_CLUSTER*BLS_BLOCK_SIZE)

/* A BLS cluster is identified by its disk number (0,1,2,...) and */
/* per disk cluster number.  These two values are merged into a single */
/* __uint32_t value. */
typedef __uint32_t bls_cid_t;
#define	BLS_CID_CLUSTER_BITS	27
#define	BLS_CID_MAX_NCLUSTER	(1<<BLS_CID_CLUSTER_BITS)
#define	BLS_CID_DISK_BITS	5
#define	BLS_CID_MAX_NDISK	(1<<BLS_CID_DISK_BITS)
#define	BLS_CID_DISK_MASK	(BLS_CID_MAX_NDISK-1)
#define	BLS_DC2CID(DISK, CLUSTER)	\
	(((CLUSTER)<<BLS_CID_DISK_BITS) | (DISK) )
#define	BLS_CID2DISK(CID)	((CID) & BLS_CID_DISK_MASK)
#define	BLS_CID2CLUSTER(CID)	((CID) >> BLS_CID_DISK_BITS)

/* macros for long long <-> hi / lo __uint32_t conversion */
#define	BLS_LL2HI(LL)	((__uint32_t)((LL)>>32))
#define	BLS_LL2LO(LL)	((__uint32_t)((LL)&((1LL<<32)-1)))
#define	BLS_HI2LL(HI)	(((long long)(HI))<<32)
#define	BLS_LO2LL(LO)	((long long)(LO))

/* A BLS metafile consists of an eight byte file length, followed by */
/* one or more bls_cid_t cluster identifiers. */

typedef struct
{
   __uint32_t token;
   __uint32_t arg[3];
} bls_packet_t;

#define BLS_MAGIC       	(('B'<<24)|('L'<<16)|('S'<<8))
#define BLS_RESPONSE(X) 	((X)^0xFF)
#define	BLS_RESPONSE_ERRNO		0
#define	BLS_RESULT_HI			1
#define	BLS_RESULT_LO			2
#define BLS_OPEN        	(BLS_MAGIC + 'o')
#define	BLS_OPEN_SIZE			0
#define	BLS_OPEN_OFLAG			1
#define	BLS_QUERY		(BLS_MAGIC + 'q')
#define	BLS_QUERY_MINIMUM		0
#define	BLS_QUERY_EOFP			1
#define	BLS_SEEK		(BLS_MAGIC + 's')
#define	BLS_SEEK_OFFSET_HI		0
#define	BLS_SEEK_OFFSET_LO		1
#define	BLS_SEEK_WHENCE			2
#define	BLS_READ		(BLS_MAGIC + 'r')
#define	BLS_READ_COUNT			0
#define	BLS_WRITE		(BLS_MAGIC + 'w')
#define	BLS_WRITE_COUNT			0
#define	BLS_LS			(BLS_MAGIC + 'l')
#define	BLS_DISK_TOTAL		(BLS_MAGIC + 't')
#define	BLS_DISK_FREE		(BLS_MAGIC + 'f')
#define	BLS_MTIME		(BLS_MAGIC + 'm')
#define	BLS_INFO		(BLS_MAGIC + 'i')

/* In response to a BLS_INFO request, BLSD sends a bls_packet_t with */
/* more than three args.  The total number of args is sent as arg[0] */
#define	BLSI_ARGS		0
#define	BLSI_TIMESTAMP		1
#define	BLSI_BITS		2
#define	BLSI_BITS_CFG_OK	(1<<0)
#define	BLSI_BITS_INIT_OK	(1<<1)
#define	BLSI_DISK_TOTAL_HI	3
#define	BLSI_DISK_TOTAL_LO	4
#define	BLSI_DISK_FREE_HI	5
#define	BLSI_DISK_FREE_LO	6
#define	BLSI_LOCAL_IN_HI	7
#define	BLSI_LOCAL_IN_LO	8
#define	BLSI_LOCAL_OUT_HI	9
#define	BLSI_LOCAL_OUT_LO	10
#define	BLSI_LOCAL_N_RO_WO	11	/* two unsigned shorts */
#define	BLSI_LOCAL_N_RW_NOIO	12	/* two unsigned shorts */
#define	BLSI_REMOTE_IN_HI	13
#define	BLSI_REMOTE_IN_LO	14
#define	BLSI_REMOTE_OUT_HI	15
#define	BLSI_REMOTE_OUT_LO	16
#define	BLSI_REMOTE_N_RO_WO	17	/* two unsigned shorts */
#define	BLSI_REMOTE_N_RW_NOIO	18	/* two unsigned shorts */
#define	BLSI_MAX_ARG		18
#define	SIZEOF_BLS_INFO(N)	\
	(offsetof(bls_packet_t, arg) + (N)*sizeof(__uint32_t))

/* The results of a BLS_INFO request are presented to the user */
/* via the following structure */
typedef struct
{
   __uint32_t timestamp;	/* milliseconds since epoch (wraps) */
   unsigned cfg_ok : 1;
   unsigned init_ok : 1;
   struct
   {
      long long total;		/* total disk space, in bytes */
      long long free;		/* free disk space, in bytes */
   } disk;
   struct
   {
      long long in;		/* input byte count (since daemon startup) */
      long long out;		/* output byte count (since daemon startup) */
      unsigned short n_ro;	/* current number of read-only connections */
      unsigned short n_wo;	/* current number of write-only connections */
      unsigned short n_rw;	/* current number of read-write connections */
      unsigned short n_noio;	/* current number of no I/O connections */
   } local, remote;
} bls_info_t;

/* for clients with different "errno.h" than the BLSD (running on Solaris) */
/* errno from the BLSD is shifted out of the normal errno range */
#ifdef VXWORKS
#define	BLSD_ERRNO(n)		((n)<<16)
#else
#define	BLSD_ERRNO(n)		(n)
#endif

/* values for "bls_open()" "oflag" argument */
#define	BLSO_RDOK			(1<<16)
#define	BLSO_WROK			(1<<17)
#define	BLSO_CREAT			(1<<18)
#define	BLSO_TRUNC			(1<<19)
#define	BLSO_EXCL			(1<<20)
#define	BLSO_MWROK			(1<<21)
#define	BLSO_RDONLY			(BLSO_RDOK)
#define	BLSO_WRONLY			(BLSO_WROK)
#define	BLSO_RDWR			(BLSO_RDOK|BLSO_WROK)
#define	BLSO_ALLOFLAGS		\
	(BLSO_RDOK|BLSO_WROK|BLSO_CREAT|BLSO_TRUNC|BLSO_EXCL|BLSO_MWROK)

/* values for "bls_lseek()" "whence" argument */
#define	BLSS_SET			0
#define	BLSS_CUR			1
#define	BLSS_END			2
#define	BLSS_ALLOC			3
#define	BLSS_EXTEND			4

#define	BLS_MAX_CID	256
typedef struct
{
   int sock;
   void *fp;	/* opaque FILE pointer for buffered output to "sock" */
   __uint32_t oflag;
   __uint32_t in_progress;
   int nowait;
   int n_cid;
   bls_cid_t cid[BLS_MAX_CID];
#ifdef VXWORKS
   long long filepos;
   struct { int disk, sector, count; } scsiWrtSecs_error;
#endif
} BLS;

/* macro to set the "nowait" flag and return the BLS pointer */
/* this allows the usage e.g. "bls_query(BLS_NOWAIT(bp), 0)" */
#define	BLS_NOWAIT(bp)	( (bp)->nowait = 1, (bp) )

/* prototypes */
BLS *bls_open( char *name, unsigned oflag );
int bls_delete( char *name );
int bls_rename( char *old, char *new );
unsigned bls_query( BLS *bp, unsigned minimum, long long *eofp );
long long bls_seek( BLS *bp, long long offset, int whence );
void bls_close( BLS *bp );
int bls_read( BLS *bp, void *buffer, int count );
int bls_write( BLS *bp, void *buffer, int count );
long long bls_ls( BLS *bp, char *buffer, int bufsize, long long *timep );
long long bls_disk_total( BLS *bp );
long long bls_disk_free( BLS *bp );
time_t bls_mtime( BLS *bp );
int bls_info( BLS *bp, bls_info_t *infop );
int init_bls( int n );
int bls_label( char *block );

#endif /* BLS_H */
