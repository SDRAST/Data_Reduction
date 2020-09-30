
import logging

# import matplotlib.pyplot as plt # only required by
# FluxCalibrationDataAnalyzer method plot()


module_logger = logging.getLogger(__name__)


class FluxCalibrationDataAnalyzer(object):

    def __init__(self, calibration_data):
        """
        """
        self.data = calibration_data
        self.fits = {
            "gains": None,
            "linear": None,
            "quadratic": None,
            "nd_temp": None,
            "non_linearity": None,
            "tsys_factors": None
        }

    def calculate_tsys_factors(self):
        """
        """
        from ...server.dss_server import DSSServer
        res = DSSServer.process_minical_calib(self.data)
        for key in res:
            self.fits[key] = res[key]

    def generate_report(self):
        pass

    # I don't think this is used anywhere
    def plot(self,
             channel="all",
             ax=None,
             x_label=None,
             y_label=None,
             title=None):
        if channel is None or channel == "all":
            channel = list(range(4))
        else:
            channel = [channel]
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        if x_label is None:
            x_label = "Power meter readings (W)"
        if y_label is None:
            y_label = "Power meter readings (K)"
        if title is None:
            title = "Flux calibration fit results"
            if len(channel) == 1:
                title += " for channel {}".format(channel[0])
            else:
                title += " for all channels"
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        for c in channel:
            # need to plot linear and quadratic fits
            x = [self.data[key][c] for key
                 in ["sky", "sky+ND", "load", "load+ND"]]
            # linear = [self.fits["linear"][i][c] for i in xrange(len(self.fits["linear"]))]
            quadratic = [self.fits["quadratic"][i][c]
                         for i in range(len(self.fits["quadratic"]))]
            # ax.plot(x, linear, label="Linear, channel {}".format(c), marker="o")
            ax.plot(
                x, quadratic,
                label="Quadratic, channel {}\ntsys factor: {:.3E}".format(
                    c, self.fits["tsys_factors"][c]),
                marker="o"
            )

        ax.legend()
        return ax
