"""
from https://documen.tician.de/pycuda/tutorial.html
"""
import pycuda
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy

# create 32-bit float random data
a = numpy.random.randn(4,4)
a = a.astype(numpy.float32)
print a
a_gpu = cuda.mem_alloc(a.nbytes)
cuda.memcpy_htod(a_gpu, a)

# define the function to be executed in the GPU
mod = SourceModule("""
   __global__ void doublify(float *a)
   {
     int idx = threadIdx.x + threadIdx.y*4;
     a[idx] *= 2;
   }""")

# invoke the function
func = mod.get_function("doublify")
func(a_gpu, block=(4,4,1))

#retrieve and display the result
a_doubled = numpy.empty_like(a)
cuda.memcpy_dtoh(a_doubled, a_gpu)
print a_doubled
