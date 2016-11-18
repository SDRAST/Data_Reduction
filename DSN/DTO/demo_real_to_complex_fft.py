"""
This example illustrates how to attach a transformation to an FFT computation 
object that will make it operate on real-valued inputs.

The test is to compare the GPU result with a numpy.fft result.

From::
  https://github.com/fjarri/reikna/blob/develop/examples/demo_real_to_complex_fft.py
"""
import numpy
from reikna.cluda import dtypes, any_api
from reikna.fft import FFT
from reikna.core import Annotation, Type, Transformation, Parameter

def get_complex_trf(arr):
    """
    A transformation that transforms a real array to a complex one by adding a
    zero imaginary part
    """
    complex_dtype = dtypes.complex_for(arr.dtype)
    return Transformation(
        [Parameter('output', Annotation(Type(complex_dtype, arr.shape), 'o')),
        Parameter('input', Annotation(arr, 'i'))],
        """
        ${output.store_same}(
            COMPLEX_CTR(${output.ctype})(
                ${input.load_same},
                0));
        """)


# make the test data
arr = numpy.random.normal(size=3000).astype(numpy.float32)
trf = get_complex_trf(arr)

# Pick the first available GPGPU API and make a Thread on it.
api = any_api()
thr = api.Thread.create()

# Create the FFT computation and attach the transformation above to its input.
fft = FFT(trf.output) # (A shortcut: using the array type saved in the transformation)
fft.parameter.input.connect(trf, trf.output, new_input=trf.input)
cfft = fft.compile(thr)

# Run the computation
arr_dev = thr.to_device(arr)
res_dev = thr.array(arr.shape, numpy.complex64)
cfft(res_dev, arr_dev)
result = res_dev.get()

# Compare to numpy
reference = numpy.fft.fft(arr)
assert numpy.linalg.norm(result - reference) / numpy.linalg.norm(reference) < 1e-6
