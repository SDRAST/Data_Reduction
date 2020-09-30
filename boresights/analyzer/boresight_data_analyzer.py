"""
boresight_data_analyzer.py
"""
import operator
import logging

import numpy as np
import scipy.stats as stats

import Data_Reduction.boresights.util as DRbu
import MonitorControl.DSS_server_cfg as config

module_logger = logging.getLogger(__name__)

__all__ = [
    "SingleAxisBoresightDataAnalyzer",
    "ChannelDataAnalyzer"
]

_model = DRbu.ModelProxy(
                   model_file_path=config.tams_config.boresight_model_file_path)


def _set_model_file_path(*args):
    _model.model_file_path = config.tams_config.boresight_model_file_path


#config.tams_config.on("boresight_model_file_path", _set_model_file_path)
#config.tams_config.on("project_dir", _set_model_file_path)
#config.tams_config.on("model_dir", _set_model_file_path)


class ChannelDataAnalyzer(object):
    operator_method_name_mapping = {
        "+": "add",
        "-": "sub",
        "*": "mul",
        "/": "div"
    }

    def __init__(self, offset_data, tsys_data, channel="0"):
        """
        """
        self.offset_data = np.array(offset_data)
        self.tsys_data = np.array(tsys_data)
        self.channel = channel
        self.fit = {
            "chi_sqr": None,
            "chi_sqr_dof": None,
            "popt": None,
            "pcov": None,
            "amplitude": None,
            "offset": None,
            "sigma": None,
            "amplitude_err": None,
            "offset_err": None,
            "sigma_err": None,
            "score": None
        }

    def calculate_fit(self, **fit_kwargs):
        """
        Fit tsys_data to a gaussian plus a line.
        Args:
            fit_kwargs (dict): passed to util.gauss_fit
        """
        if "mean" not in fit_kwargs:
            fit_kwargs["mean"] = self.offset_data[np.argmax(self.tsys_data)]

        def chi_sqr(x, y, *model_param):
            model = DRbu.gauss_function(x, *model_param)
            return np.sum(np.sqrt((model-y)**2))

        fit = DRbu.gauss_fit(self.offset_data, self.tsys_data, **fit_kwargs)
        fit_chi_sqr = chi_sqr(self.offset_data, self.tsys_data, *fit["popt"])

        self.fit["chi_sqr"] = fit_chi_sqr
        self.fit["chi_sqr_dof"] = fit_chi_sqr / float(len(self.offset_data) - len(fit["popt"]))
        self.fit["popt"] = fit["popt"]
        self.fit["pcov"] = fit["pcov"]

        self.fit["amplitude"] = fit["popt"][0]
        self.fit["offset"] = fit["popt"][1]
        self.fit["sigma"] = fit["popt"][2]

        self.fit["amplitude_err"] = np.sqrt(fit["pcov"][0, 0])
        self.fit["offset_err"] = np.sqrt(fit["pcov"][1, 1])
        self.fit["sigma_err"] = np.sqrt(fit["pcov"][2, 2])
        return self

    def calculate_score(self):
        score = _model.predict(self.fit)
        self.fit["score"] = score
        return self

    def __sub__(self, other):
        new_data = (self.tsys_data - other.tsys_data)
        new_channel = "{} - {}".format(self.channel, other.channel)
        new_offset_data = np.copy(self.offset_data)
        return ChannelDataAnalyzer(
            new_offset_data, new_data, channel=new_channel)

    def __add__(self, other):
        new_data = self.tsys_data + other.tsys_data
        new_channel = "{} + {}".format(self.channel, other.channel)
        new_offset_data = np.copy(self.offset_data)
        return ChannelDataAnalyzer(
            new_offset_data, new_data, channel=new_channel)

    def __getitem__(self, item):
        if item in self.fit:
            return self.fit[item]

    def to_dict(self):
        obj_dict = {}
        obj_dict["__class__"] = "{}.{}".format(
            self.__module__, self.__class__.__name__)
        obj_dict["offset_data"] = self.offset_data.tolist()
        obj_dict["tsys_data"] = self.tsys_data.tolist()
        obj_dict["channel"] = self.channel
        obj_dict["fit"] = self.fit
        for param in obj_dict["fit"]:
            if obj_dict["fit"][param] is not None:
                vals = np.array(obj_dict["fit"][param])
                vals[np.logical_not(np.isfinite(vals))] = 1e99
                obj_dict["fit"][param] = vals.tolist()
            # if isinstance(obj_dict["fit"][param], np.ndarray):
                # obj_dict["fit"][param] = obj_dict["fit"][param].tolist()

        return obj_dict

    @classmethod
    def from_dict(cls, cls_name, obj_dict):
        obj = cls(obj_dict["offset_data"], obj_dict["tsys_data"],
                  channel=obj_dict["channel"])
        obj.fit = obj_dict["fit"].copy()
        return obj


class SingleAxisBoresightDataAnalyzer(object):

    def __init__(self,
                 offset_data=None,
                 tsys_data=None,
                 axis="el",
                 direction="right"):

        offset_data = np.array(offset_data)
        tsys_data = np.array(tsys_data)
        self.channels = {
            str(i): ChannelDataAnalyzer(
                offset_data, tsys_data[:, i], channel=str(i))
            for i in range(tsys_data.shape[1])
        }
        self.axis = axis
        self.direction = direction

        module_logger.debug("{}.__init__: n channels: {}".format(
            self.__class__.__name__, len(self.channels)))
        module_logger.debug("{}.__init__: axis: {}".format(
            self.__class__.__name__, self.axis))
        module_logger.debug("{}.__init__: direction: {}".format(
            self.__class__.__name__, self.direction))
        # module_logger.debug("{}.__init__: offset data: {}...{}".format(
        #     self.__class__.__name__,
        #     offset_data[0],
        #     offset_data[-1])
        # )
        # module_logger.debug("{}.__init__: tsys data: {}...{}".format(
        #     self.__class__.__name__, tsys_data[0], tsys_data[-1]))
        self.offset_data = offset_data
        self.tsys_data = tsys_data

    def get_fits(self):
        """
        Turn numpy arrays in fits attribute into normal Python lists and return
        """
        fits = {key: {} for key in self.channels}
        for channel in self.channels:
            for param in channels[channel].fit:
                if isinstance(self.channels[channel].fit[param], np.ndarray):
                    fits[channel][param] = \
                        channels[channel].fit[param].tolist()
                else:
                    fits[channel][param] = channels[channel].fit[param]
        return fits

    def calculate_offsets(self, channels=None, **fit_kwargs):
        """
        Calculate fit parameters.
        """
        def _calculate_offsets(channel):

            if hasattr(channel, "format"):  # is it a string?
                if channel in self.channels:
                    self.channels[channel].calculate_fit().calculate_score()

            elif hasattr(channel, "append"):  # is this a list?
                chan1, operator_method_name, chan2 = channel
                if operator_method_name in ChannelDataAnalyzer.operator_method_name_mapping:
                    operator_method_name = ChannelDataAnalyzer.operator_method_name_mapping[operator_method_name]
                operator_method = getattr(operator, operator_method_name)
                if chan1 in self.channels and chan2 in self.channels:
                    new_channel = operator_method(
                        self.channels[chan1],
                        self.channels[chan2]
                    ).calculate_fit(**fit_kwargs).calculate_score()
                    self.channels[new_channel.channel] = new_channel

        if channels is None:
            channels = list(self.channels.keys())
        for c in channels:
            _calculate_offsets(c)

    def calculate_offset_distributions(self,
                                       channel=None,
                                       callback=DRbu.gauss_function,
                                       n_walkers=100,
                                       iterations=500):
        """
        Instead of simplly doing a fit, do MCMC in order to determine best
        fit paramter distributions. Use as starting points the results of
        scipy.optimize.curve_fit.
        """
        import emcee

        def lnprior_gaussian(init_params):
            beam_fwhm = 13.5
            amp, off, sig, slope, intercept = init_params
            prob_dist = [
                stats.norm(loc=amp, scale=1.0),
                stats.norm(loc=off, scale=beam_fwhm),
                stats.norm(loc=sig, scale=beam_fwhm),
                stats.norm(loc=slope, scale=0.01),
                stats.norm(loc=intercept, scale=0.01),
            ]
            def _lnprior(fit_params):
                res = 1.0
                for i in range(len(fit_params)):
                    dist = prob_dist[i]
                    p = fit_params[i]
                    res *= dist.pdf(p)
                return np.log(res)
            return _lnprior

        def lnprior_uniform(init_params):
            amp, off, sig, slope, intercept = init_params
            beam_fwhm = 13.5
            quarter_beam = beam_fwhm/4.0
            quarter_beam = 1.0
            # quarter_beam =50.0

            param_ranges = [
                [amp - (0.5*amp), amp + (0.5*amp)],
                [off-beam_fwhm, off+beam_fwhm],
                [sig-quarter_beam, sig+quarter_beam],
                [slope-(abs(slope)), slope+(abs(slope))],
                [intercept-0.01, intercept+0.01]
            ]
            def _lnprior(fit_params):
                """
                Assume a uniform distribution over model parameters
                """
                for i in range(len(fit_params)):
                    p = fit_params[i]
                    p_range = param_ranges[i]
                    if p > p_range[1] or p < p_range[0]:
                        return -np.inf
                return 0.0
            return _lnprior

        def lnlikelihood(fit_params, x, y):
            """
            Log likelihood function. We don't have error on our y values
            (allthought we definitely should), so we're not including that
            in the calculation
            Args:
                fit_params (tuple/list):
                x (np.ndarray):
                y (np.ndarray):
            returns
                np.ndarray
            """
            model = callback(x, *fit_params)
            return -0.5 * np.sum((y-model)**2)

        def lnprob(prior_fn, likelihood_fn):
            def _lnprob(fit_params, x, y):
                prior = prior_fn(fit_params)
                if not np.isfinite(prior):
                    return -np.inf
                return prior + likelihood_fn(fit_params, x, y)
            return _lnprob

        def create_sampler(n_walkers, init_params, x, y):
            likelihood_fn = lnlikelihood
            prior_fn = lnprior_gaussian(init_params)
            prob_fn = lnprob(prior_fn, likelihood_fn)
            module_logger.debug("calculate_offset_distributions.create_sampler: n_walkers: {}, init_params: {}".format(n_walkers, init_params))
            n_dim = len(init_params)
            init_pos = [init_params + 1e-4*np.random.rand(n_dim) for i in range(n_walkers)]
            module_logger.debug("calculate_offset_distributions.create_sampler: init_pos[0] {}".format(init_pos[0]))
            sampler = emcee.EnsembleSampler(n_walkers, n_dim, prob_fn, args=(x, y))
            return sampler, init_pos

        def get_samples(channel):
            if self.channels[channel]["popt"] is None:
                self.calculate_offsets(channel=channel)
            best_fit = self.channels[channel]["popt"]
            sampler, init_pos = create_sampler(n_walkers,best_fit,self.offset_data, self.tsys_data[:,channel])
            pos,prob,state = sampler.run_mcmc(init_pos,iterations)
            samples = sampler.chain[:, int(iterations/5):, :].reshape((-1, len(best_fit)))
            return samples

        samples = []
        if channel is None:
            samples = [get_samples(i) for i in range(self.n_channels)]
        else:
            samples = get_samples(channel)
        return samples

    def generate_report(self, delimiter=" "):
        fit_params = ["amplitude", "offset", "sigma"]
        report = []
        for channel in self.channels:
            for fit_param in fit_params:
                report.append(
                    "{:.3f}+/-{:.3f}".format(
                        self.channels[channel][fit_param],
                        self.channels[channel][fit_param+"_err"]
                    )
                )
        return delimiter.join(report)

    def _plot_matplotlib(self,
                         channel="0",
                         ax=None,
                         x_label=None,
                         y_label=None,
                         title=None):
        module_logger.debug("_plot_matplotlib: called")
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        if x_label is None:
            x_label = "{} offset, {} direction".format(
                self.axis, self.direction)
        if y_label is None:
            y_label = "Power meter reading"
        if title is None:
            title = "Channel {}\n".format(channel)
            title += "Calculated offset: {:.3f}+/-{:.3f}, ".format(
                self.channels[channel]["offset"],
                self.channels[channel]["offset_err"]
            )
            title += "Calculated Amplitude {:.3f}+/-{:.3f}".format(
                self.channels[channel]["amplitude"],
                self.channels[channel]["amplitude_err"]
            )
        offset_data = self.channels[channel].offset_data
        fit_x = np.linspace(offset_data[0], offset_data[-1], 200)
        fit_y = DRbu.gauss_function(fit_x, *self.channels[channel]["popt"])

        res_scatter = ax.scatter(offset_data,
                                 self.channels[channel].tsys_data,
                                 color="orange",
                                 alpha=0.8)
        res_plot = ax.plot(fit_x, fit_y)
        ax.set_xlabel(x_label)
        ax.set_title(title)
        ax.set_ylabel(y_label)
        module_logger.debug("_plot_matplotlib: ax.lines: {}".format(ax.lines))
        return ax

    def _plot_gnuplotlib(self, channel="0", title=None, **kwargs):
        import gnuplotlib as gp
        if title is None:
            title = "Channel {} Tsys (K) vs Offset in {}".format(
                channel, self.axis)
        offset_data = self.channels[channel].offset_data
        fit_x = np.linspace(offset_data[0], offset_data[-1],200)
        gp.plot(
            (offset_data, self.channels[channel].tsys_data),
            (fit_x, DRbu.gauss_function(fit_x, *self.channels[channel]["popt"])),
            _with='lines', terminal='dumb 80,40', unset='grid'
        )

    def plot(self, **kwargs):
        if config.backend == "matplotlib":
            return self._plot_matplotlib(**kwargs)
        elif config.backend == "gnuplotlib":
            return self._plot_gnuplotlib(**kwargs)

    def __getitem__(self, item):
        if item in self.channels:
            return self.channels[item].fit
        else:
            raise IndexError(
                "Can't find item {} in self.channels.keys: {}".format(
                    item, list(self.channel.keys())))

    def __iter__(self):
        for channel in self.channels:
            yield self.channels[channel].fit

    def __contains__(self, item):
        if item in self.channels:
            return True
        else:
            return False

    def to_dict(self):
        obj_dict = {}
        obj_dict["__class__"] = "{}.{}".format(
            self.__module__, self.__class__.__name__)
        obj_dict["offset_data"] = self.offset_data.tolist()
        obj_dict["tsys_data"] = self.tsys_data.tolist()
        obj_dict["axis"] = self.axis
        obj_dict["direction"] = self.direction
        obj_dict["channels"] = {}
        for chan in self.channels:
            obj_dict["channels"][chan] = self.channels[chan].to_dict()
        return obj_dict

    @classmethod
    def from_dict(cls, cls_name, obj_dict):
        # channels = {chan_name: ChannelDataAnalyzer.from_dict(
        #     "", obj_dict["channels"][chan_name]
        # ) for chan_name in obj_dict["channels"]}
        # obj = cls(channels=channels,
        #           axis=obj_dict["axis"],
        #           direction=obj_dict["direction"])
        obj = cls(offset_data=obj_dict["offset_data"],
                  tsys_data=obj_dict["tsys_data"],
                  axis=obj_dict["axis"],
                  direction=obj_dict["direction"])
        # obj.channels = {chan_name: ChannelDataAnalyzer.from_dict(
        #     "",
        #     obj_dict["channels"][chan_name]
        # ) for chan_name in obj_dict["channels"]}
        return obj

