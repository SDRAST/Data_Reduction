import logging
import os
import datetime
try:
    import pickle as pickle
except ImportError as err:
    import pickle

import numpy as np
import scipy.optimize as op
import Pyro5.serializers

#from ..configuration import config

module_logger = logging.getLogger(__name__)

__all__ = [
    "gauss_function",
    "line",
    "gauss_fit",
    "register_cls_Pyro5_serializer",
    "ModelProxy"
]


def gauss_function(x, a, x0, sigma, slope=0, intercept=0):
    return a*np.exp(-(x - x0)**2/(2*sigma**2)) + slope * x + intercept


def line(x, a, b, c, slope, intercept):
    val = slope * x + intercept
    return val


def gauss_fit(offsets, tsys, mean=0.0, sigma=2.1):
    slope = (tsys[-1] - tsys[0])/(offsets[-1] - offsets[0])
    intercept = (tsys[-1] + tsys[0]) / 2.
    try:
        popt, pcov = op.curve_fit(
            gauss_function, offsets, tsys, p0=(
                1, mean, sigma, slope, intercept))
    except Exception as err:
        # module_logger.warning(
        #   "gauss_fit: Optimal results could not be obtained for gaussian fit.")
        # module_logger.warning(
        #   "gauss_fit: Error: {}".format(err))
        popt = np.zeros(5)
        pcov = np.zeros((5, 5))
    return {
        "popt": popt,
        "pcov": pcov
    }


def gaussian(x, amp, mu, sigma):
    return amp*np.exp(-(x - mu)**2/(2*sigma**2))


def n_gaussian(x, *params):
    n_params = 3
    res = np.zeros(x.shape[0])
    for i in range(0, len(params), n_params):
        param = params[i:i+n_params]
        res += gaussian(x, *param)
    return res


def find_most_likely_value(samples, **histogram_kwargs):
    hist, bin_edges = np.histogram(samples, **histogram_kwargs)
    highest_prob_idx = np.argmax(hist)
    return bin_edges[highest_prob_idx]


def quadratic(x, a, x0, y0):
    return (a*(x-x0)**2)+y0


def register_cls_Pyro5_serializer(cls):
    if hasattr(cls, "to_dict"):
        module_logger.debug("registering {} to SerializerBase".format(cls))
        Pyro5.serializers.SerializerBase.register_class_to_dict(
            cls, getattr(cls, "to_dict")
        )
        Pyro5.serializers.SerializerBase.register_dict_to_class(
            "{}.{}".format(cls.__module__, cls.__name__),
            getattr(cls, "from_dict")
        )


class ModelProxy(object):

    def __init__(self, model_file_path=None):
        module_logger.debug("ModelProxy.__init__: path is %s", model_file_path)
        self._model_file_path = model_file_path
        self.model = None
        self.params = None
        try:
            self.load()
            module_logger.debug(
                "ModelProxy.__init__: loaded model file")
        except Exception as err:
            module_logger.error(
                "ModelProxy.__init__: couldn't load model: {}".format(err)
            )

    @property
    def model_file_path(self):
        return self._model_file_path

    @model_file_path.setter
    def model_file_path(self, val):
        self._model_file_path = val
        try:
            self.load()
            module_logger.debug(
                "ModelProxy.model_file_path: loaded model file")
        except Exception as err:
            module_logger.error(
                "ModelProxy.model_file_path: couldn't load model: {}".format(err)
            )

    def load(self):
        """
        """
        with open(self.model_file_path, "rb") as f:
            module_logger.debug("ModelProxy.load: file opened as %s", f)
            model_dict = pickle.load(f, encoding="latin1")
            module_logger.debug("ModelProxy.load: loaded: %s", model_dict)
            self.model = model_dict["clf"]
            self.params = model_dict["params"]

    def predict(self, fits):
        if self.model is None:
            return
        if hasattr(fits, "keys"):
            fits = [fits[k] for k in self.params]
        fits = np.array(fits)
        if fits.ndim == 1:
            fits = fits.reshape(1, -1)
        fits[np.isposinf(fits)] = 0.0
        fits[np.isneginf(fits)] = 0.0
        fits[np.isnan(fits)] = 0.0
        module_logger.debug("ModelProxy.predict: fits: {}".format(fits))
        return self.model.predict(fits)
