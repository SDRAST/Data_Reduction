/******************************************************************************
This takes 'num_words' words from one subchannel record (schan_rcd_t *rcd)
and parses each for 'cmplx_samp_per_32bit_word' samples/word and turns them into
complex data (FFTW_COMPLEX *data).

Consider the case of 8 bit samples on a 32-bit machine. The packing on disk for
8-bit samples after Endian byte-swapping with ntohl() is [I2][I1][R2][R1].
   The function takes the 32 bit word and shifts it left 24 bits, leaving the
rightmost eight.  It then shifts that back 24 bits leaving the real part of
the first sample.
   It then shifts the word left 8 bits and right 24 bits. This gives the
imaginary part of the first sample.
   Then it takes the word and shifts it left 16 bits and right 24 bits for
the real part of the second sample.
   Finally, it shifts the word left 0 bits and right 24 bits for the imaginary
part of the second sample.

On a 64-bit machine, a 64-bit word is packed, after Endian byte-swapping with
nthol() is [I5][I3][R4][R3][I2][I1][R2][R1].

So the algorithm for extracting samples from one unsigned integer is:
for (j = 0; j < samp_per_word; j++)
{
     // left shift & right shift data to obtain isamp and qsamp, this
     // isolates and sign extends the current sample in the vme_data word
    isamp = ntohl(rcd->data[i]) << (bits_per_word - (bits_per_sample * (j + 1)));
    isamp >>= (bits_per_word - bits_per_sample);
    qsamp = ntohl(rcd->data[i]) << (bits_per_word - (bits_per_sample * (j + 3)));
    qsamp >>= (bits_per_word – bits_per_sample);
}
******************************************************************************/
/* system includes */
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

/* project includes */
#include "utils.h"
#include "bls.h"
#include "schan_rcd.h"

int diag = 0;

static PyObject * set_diag_flag(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "i", &diag))
    return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject * unpack_vdr_data(PyObject *self, PyObject *args)
{
  /* Note that samp_per_32bit_word is computed from
   * (BITS_PER_WORD/2)/ntohs(rcd->hdr.bits_per_sample)
   * where BITS_PER_WORD depends on the machine architecture */

  int           num_32bit_words;           /* number of 32-bit words in the buffer */
  int           integer_samp_per_32bit_word;
  int           cmplx_samp_per_32bit_word; /* complex samples */
  const U32 *   buffer;
  int           buflen;
  
  int           i, j;                            /* loop counters */
  Py_complex    cmplx;
  PyObject * Cmplx;
  int           L2;
  PyObject * Clist;
  PyObject * result;
  U16        bits_per_sample;
  I32        word_32bit, isamp, qsamp;
  Py_ssize_t time_point_count = 0;

  if (!PyArg_ParseTuple(args, "iis#", &num_32bit_words,
                                      &cmplx_samp_per_32bit_word,
                                      &buffer,
                                      &buflen))
    return NULL;
  if (diag) {
    printf("unpack_vdr_data: Buffer length is %d 32-bit words\n",num_32bit_words);
    fflush(stdout);
  }
  if (diag) {
    printf("unpack_vdr_data: %d complex samples per 32-bit word\n",
           cmplx_samp_per_32bit_word);
  }
  if (diag) {
    printf("unpack_vdr_data: Buffer size %d bytes\n", buflen);
    fflush(stdout);
  }
  integer_samp_per_32bit_word = cmplx_samp_per_32bit_word*2;
  bits_per_sample = 32/integer_samp_per_32bit_word;
  if (diag) {
    printf("unpack_vdr_data: Bits per sample is %d\n", bits_per_sample);
    fflush(stdout);
  }
  
  L2 = num_32bit_words*cmplx_samp_per_32bit_word;
  if (diag) {
    printf("unpack_vdr_data: List will have %d complex numbers\n",L2);
  }
  Clist = PyList_New((Py_ssize_t) L2);
  
  for(i = 0; i < num_32bit_words; i++)
  {
    word_32bit = ntohl(buffer[i]);
    // word_32bit = buffer[i];
    if (diag) {
      printf("raw 32-bit word %d = %08x\n",i,buffer[i]);
      printf("word %d = %08x\n",i,word_32bit);
      fflush(stdout);
    }
    for(j = 0; j < cmplx_samp_per_32bit_word; j++)
    {
      /* left shift & right shift data to obtain isamp and qsamp, this
       * isolates and sign extends the current sample in the vme_data
       * word */
      isamp = word_32bit << (32 - (bits_per_sample * (j + 1)));
      if (diag) {
        printf("sample %d, real - step 1: %08x\n",j,isamp);
        fflush(stdout);
      }
      isamp >>= (32 - bits_per_sample);
      if (diag) {
        printf("sample %d, real - step 2: %08x\n",j,isamp);
        fflush(stdout);
      }
         
      qsamp = word_32bit << (16 - (bits_per_sample * (j + 1)));
      if (diag) {
        printf("sample %d, imag - step 1: %08x\n",j,qsamp);
        fflush(stdout);
      }
      qsamp >>= (32 - bits_per_sample);
      if (diag) {
       printf("sample %d, imag - step 2: %08x\n",j,qsamp);
       fflush(stdout);
      }
      /* we need to add 0.5 to remove dc offset caused by the
      * truncation of samples in last stage of ASD filters */
      cmplx.real = (float)isamp + 0.5;
      cmplx.imag = (float)qsamp + 0.5;
      if (diag) {
        printf("complex pair = %f,%f\n",cmplx.real,cmplx.imag);
        fflush(stdout);
      }
      Cmplx = PyComplex_FromCComplex(cmplx);
      if (diag) {
        printf("Created Python complex object %d\n", time_point_count);
      }
      if (PyList_SetItem(Clist, time_point_count, Cmplx)) {
        fprintf(stderr,"Could not set list item %d\n",time_point_count);
      }
      time_point_count++;
    }
  }
  result = Py_BuildValue("O", Clist);
  Py_CLEAR(Cmplx);
  Py_CLEAR(Clist);
  return result;
}

/* This is the Method table.  The name MUST be the same as the Python
 module name.  The last entry in the table is NULLs*/

static PyMethodDef
vdr_support[] = {
  {"unpack_vdr_data",  unpack_vdr_data, METH_VARARGS,
  "Unpack a buffer with VDR data"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* This initializes the module
 The name MUST be 'init' followed directly by the Python module = method
definition name, no underscore! */

PyMODINIT_FUNC
initvdr_support(void)
{
  (void) Py_InitModule("vdr_support", vdr_support);
}
