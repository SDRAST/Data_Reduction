/******************************************************************************
This takes 'num_32bit_word' words from one subchannel record (schan_rcd_t *rcd)
and parses each for 'samp_per_32bit_word' samples/word and turns them into
complex data (FFTW_COMPLEX *data).

Consider the case of 8 bit samples.
   The function takes the 32 bit word and shifts it left 24 bits, leaving the
rightmost eight.  It then shifts that back 24 bits leaving the real part of
the first sample.
   It then shifts the word left 8 bits and right 24 bits. This gives the
imaginary part of the first sample.
   Then it takes the word and shifts it left 16 bits and right 24 bits for
the real part of the second sample.
   Finally, it shifts the word left 0 bits and right 24 bits for the imaginary
part of the second sample.
   So, the packing is [I2][I1][R2][R1] for 8-bit samples.
   For 4-bit words it is:
R1: left 28, right 28
I1: left 12, right 28
R2: left 24, right 28
I2: left  8, right 28, and so on
So the packing is [I4][I3][I2][I1][R4][R3][R2][R1].

******************************************************************************/
/* system includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fftw.h>

/* project includes */
#include "utils.h"
#include "bls.h"
#include "schan_rcd.h"

int get_vdr_data_samples(int          num_32bit_word,
			 int          samp_per_32bit_word, 
                         schan_rcd_t  *rcd,
			 FFTW_COMPLEX *data)
{
   int i, j;
   U16 bits_per_sample;
   I32 isamp, qsamp;
   int time_point_count = 0;

   /* endian conversion */
   bits_per_sample = ntohs(rcd->hdr.bits_per_sample);

   for(i = 0; i < num_32bit_word; i++)
   {
      for(j = 0; j < samp_per_32bit_word; j++)
      {
         /* left shift & right shift data to obtain isamp and qsamp, this
          * isolates and sign extends the current sample in the vme_data word */
         isamp = ntohl(rcd->data[i]) << (32 - (bits_per_sample * (j + 1)));
         isamp >>= (32 - bits_per_sample);
         
         qsamp = ntohl(rcd->data[i]) << (16 - (bits_per_sample * (j + 1)));
         qsamp >>= (32 - bits_per_sample);
      
         /* we need to add 0.5 to remove dc offset caused by the
          * truncation of samples in last stage of ASD filters */
         c_re(data[time_point_count]) = (float)isamp + 0.5;
         c_im(data[time_point_count]) = (float)qsamp + 0.5;
         
         time_point_count++;
      }
   }

   return OK;
}
