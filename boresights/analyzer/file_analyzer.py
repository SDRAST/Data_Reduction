"""
file_analyzer.py
"""
import logging
import datetime
import os
import time


__all__ = [
    "FileAnalyzer",
]

module_logger = logging.getLogger(__name__)


class FileAnalyzer(object):
    """
    FileAnalyzer superclass
    """
    def __init__(self, file_path=None):
        file_name = None
        if hasattr(file_path, "format"):
            file_name = os.path.basename(file_path)
        self._file_path = file_path
        self._file_name = file_name
        self.data = {}
        self.meta_data = {
            "file_path": self._file_path,
            "file_name": self._file_name
        }

    @property
    def file_path(self):
        return self._file_path

    @file_path.setter
    def file_path(self, new_file_path):
        new_file_name = os.path.basename(new_file_path)
        self._file_path = new_file_path
        self._file_name = new_file_name

    @property
    def file_name(self):
        return self._file_name

    def _check_existing(self, file_path, overwrite):
        if os.path.exists(file_path) and not overwrite:
            return False
        else:
            return True

    def _plot_save_path(self, save_dir=None):
        if save_dir is None:
            save_dir = os.path.dirname(self.file_path)
        save_file_name = self.file_name.replace(
            ".hdf5", "") + ".png"
        save_file_path = os.path.join(save_dir, save_file_name)
        return save_file_path, save_file_name

    def load(self, **kwargs):
        if self.file_path.endswith(".hdf5"):
            self._load_hdf5(**kwargs)
        elif self.file_path.endswith(".json"):
            self._load_json(**kwargs)
        elif self.file_path.endswith(".dat"):
            self._load_pickle(**kwargs)

    def _load_json(self):
        raise NotImplementedError()

    def _load_hdf5(self):
        raise NotImplementedError()

    def _load_pickle(self):
        raise NotImplementedError()

    def dump(self):
        if self.file_path.endswith(".hdf5"):
            self._dump_hdf5()
        elif self.file_path.endswith(".json"):
            self._dump_json()
        elif self.file_path.endswith(".dat"):
            self._dump_pickle()

    def _dump_json(self):
        raise NotImplementedError()

    def _dump_hdf5(self):
        raise NotImplementedError()

    def _dump_pickle(self):
        raise NotImplementedError()

    def plot(self):
        raise NotImplementedError()

    def to_dict(self):
        raise NotImplementedError()

    @classmethod
    def from_dict(cls):
        raise NotImplementedError()





# class TwoDirectionScanningBoresightFileAnalyzer(BoresightFileAnalyzer):
#
#     def __init__(self, *args):
#         super(TwoDirectionScanningBoresightFileAnalyzer, self).__init__(
#                 *args, boresight_type="scanning")
#
#         self.data = {
#             "el": {
#                 "left": None,
#                 "right": None,
#             },
#             "xel": {
#                 "left": None,
#                 "right": None
#             }
#         }
#
#     def load_data(self, **kwargs):
#         file_obj = h5py.File(self.file_path, "r")
#         self.data["el"]["right"] = self._load_data_hdf5(
#             file_obj, "EL", "right", **kwargs)
#         self.data["el"]["left"] = self._load_data_hdf5(
#             file_obj, "EL", "left", **kwargs)
#         self.data["xel"]["right"] = self._load_data_hdf5(
#             file_obj, "XEL", "right", **kwargs)
#         self.data["xel"]["left"] = self._load_data_hdf5(
#             file_obj, "XEL", "left", **kwargs)
#         self.load_meta_data()
#         file_obj.close()
#
#     def report_meta_data(self,
#                          delimiter="\n", line_delimiter=" | ", header=True):
#         report = []
#         report.append(
#             super(
#                 TwoDirectionScanningBoresightFileAnalyzer,
#                 self
#             ).report_meta_data(
#                     delimiter=delimiter,
#                     line_delimiter=line_delimiter,
#                     header=header
#             )
#         )
#
#         if "rate" in self.meta_data:
#             line_header = ["rate", "sample_rate", "limit"]
#             line = ["{}={}".format(i, str(self.meta_data[i])) for
#                     i in ["rate", "sample_rate", "limit"]]
#             if header:
#                 report.append(line_delimiter.join(line_header))
#             report.append(line_delimiter.join(line))
#         return delimiter.join(report)
#
#     def plot(self, channel="0", save_dir=None, overwrite=True, axs=None):
#         save_file_path, save_file_name = self._plot_save_path(save_dir)
#         if not self._check_existing(save_file_path, overwrite):
#             return
#         fig = None
#         if axs is None:
#             fig, axs = plt.subplots(2, 2, figsize=(10, 10/1.3))
#         self.data["el"]["right"].plot(ax=axs[0][0], channel=channel)
#         self.data["el"]["left"].plot(ax=axs[0][1], channel=channel)
#         self.data["xel"]["right"].plot(ax=axs[1][0], channel=channel)
#         self.data["xel"]["left"].plot(ax=axs[1][1], channel=channel)
#         suptitle = save_file_name + "\n" +\
#             self.report_meta_data(
#                 header=False,
#                 line_delimiter=", ",
#                 delimiter="\n")
#
#         if fig is not None:
#             fig.suptitle(suptitle)
#             top = 0.98-(0.03*float(len(suptitle.split('\n'))))
#             fig.tight_layout(rect=[0, 0.03, 1, top])
#             fig.savefig(save_file_path)
#             plt.close(fig)
#
#
# class SingleDirectionScanningBoresightFileAnalyzer(BoresightFileAnalyzer):
#
#     def __init__(self, *args):
#         super(
#             SingleDirectionScanningBoresightFileAnalyzer, self).__init__(
#                 *args, boresight_type="scanning")
#         self.data = {
#             "el": {
#                 "right": None,
#             },
#             "xel": {
#                 "right": None
#             }
#         }
#
#     def load_data(self, **kwargs):
#
#         file_obj = h5py.File(self.file_path, "r")
#         self.data["el"]["right"] = self._load_data_hdf5(
#             file_obj, "EL", "right", **kwargs)
#         self.data["xel"]["right"] = self._load_data_hdf5(
#             file_obj, "XEL", "right", **kwargs)
#         self.load_meta_data()
#         file_obj.close()
#
#     def report_meta_data(self,
#                          delimiter="\n",
#                          line_delimiter=" | ",
#                          header=True):
#         report = []
#         report.append(super(SingleDirectionScanningBoresightFileAnalyzer, self).report_meta_data(
#             delimiter=delimiter, line_delimiter=line_delimiter, header=header
#         ))
#         module_logger.debug(report)
#         if "rate" in self.meta_data:
#             line_header = ["rate", "sample_rate", "limit"]
#             line = ["{}={}".format(i, str(self.meta_data[i]))
#                     for i in ["rate", "sample_rate", "limit"]]
#             if header:
#                 report.append(line_delimiter.join(line_header))
#             report.append(line_delimiter.join(line))
#         return delimiter.join(report)
#
#     def plot(self, channel="0", save_dir=None, overwrite=True, axs=None):
#         # module_logger.debug(
#         #   "SingleDirectionScanningBoresightFileAnalyzer.plot: called")
#         save_file_path, save_file_name = self._plot_save_path(save_dir)
#         if not self._check_existing(save_file_path, overwrite):
#             return
#         fig = None
#         if axs is None:
#             fig, axs = plt.subplots(2, 1, figsize=(10, 10/1.3))
#         self.data["el"]["right"].plot(ax=axs[0], channel=channel)
#         self.data["xel"]["right"].plot(ax=axs[1], channel=channel)
#         suptitle = save_file_name + "\n" +\
#             self.report_meta_data(header=False,
#                                   line_delimiter=", ",
#                                   delimiter="\n")
#         if fig is not None:
#             fig.suptitle(suptitle)
#             top = 0.98-(0.03*float(len(suptitle.split('\n'))))
#             fig.tight_layout(rect=[0, 0.03, 1, top])
#             fig.savefig(save_file_path)
#             plt.close(fig)
#         return axs
#
#
# class SingleDirectionSteppingBoresightFileAnalyzer(BoresightFileAnalyzer):
#
#     def __init__(self, *args):
#         super(
#             SingleDirectionSteppingBoresightFileAnalyzer, self).__init__(
#                 *args, boresight_type="stepping")
#         self.data = {
#             "el": {
#                 "right": None,
#             },
#             "xel": {
#                 "right": None
#             }
#         }
#
#     def load_data(self, **kwargs):
#         file_obj = h5py.File(self.file_path, "r")
#         self.data["el"]["right"] = self._load_data_hdf5(
#             file_obj, "EL", "right", **kwargs)
#         self.data["xel"]["right"] = self._load_data_hdf5(
#             file_obj, "XEL", "right", **kwargs)
#         self.load_meta_data()
#         file_obj.close()
#
#     def report_meta_data(self, delimiter="\n", line_delimiter=" | ", header=True):
#         report = []
#         report.append(super(SingleDirectionSteppingBoresightFileAnalyzer, self).report_meta_data(
#             delimiter=delimiter, line_delimiter=line_delimiter, header=header
#         ))
#         if "settle_time" in self.meta_data:
#             line_header = ["settle_time", "integration_time"]
#             line = ["{}={}".format(i, str(self.meta_data[i]))
#                     for i in ["settle_time", "integration_time"]]
#             if header:
#                 report.append(line_delimiter.join(line_header))
#             report.append(line_delimiter.join(line))
#
#         return delimiter.join(report)
#
#     def plot(self, channel="0", save_dir=None, overwrite=False, axs=None):
#         save_file_path, save_file_name = self._plot_save_path(save_dir)
#         if not self._check_existing(save_file_path, overwrite):
#             return
#         fig = None
#         if axs is None:
#             fig, axs = plt.subplots(2, 1, figsize=(10, 10/1.3))
#         self.data["el"]["right"].plot(ax=axs[0], channel=channel)
#         self.data["xel"]["right"].plot(ax=axs[1], channel=channel)
#         if fig is not None:
#             suptitle = save_file_name + "\n" +\
#                 self.report_meta_data(header=False,
#                                       line_delimiter=", ",
#                                       delimiter="\n")
#             fig.suptitle(suptitle)
#             top = 0.98-(0.03*float(len(suptitle.split('\n'))))
#             fig.tight_layout(rect=[0, 0.03, 1, top])
#             fig.savefig(save_file_path)
#             plt.close(fig)
#
#     def to_dict(self):
#         module_logger.debug("SingleDirectionSteppingBoresightFileAnalyzer.to_dict: called")
#         return super(SingleDirectionSteppingBoresightFileAnalyzer, self).to_dict()
#
#
# class TwoDirectionSteppingBoresightFileAnalyzer(BoresightFileAnalyzer):
#
#     def __init__(self, *args):
#         super(
#             TwoDirectionSteppingBoresightFileAnalyzer, self).__init__(
#                 *args, boresight_type="stepping")
#         self.data = {
#             "el": {
#                 "right": None,
#                 "left": None
#             },
#             "xel": {
#                 "right": None,
#                 "left": None
#             }
#         }
#
#     def load_data(self, **kwargs):
#         file_obj = h5py.File(self.file_path, "r")
#         self.data["el"]["right"] = self._load_data_hdf5(
#             file_obj, "EL", direction="right", **kwargs)
#         self.data["xel"]["right"] = self._load_data_hdf5(
#             file_obj, "XEL", direction="right", **kwargs)
#         self.data["el"]["left"] = self._load_data_hdf5(
#             file_obj, "EL", direction="left", **kwargs)
#         self.data["xel"]["left"] = self._load_data_hdf5(
#             file_obj, "XEL", direction="left", **kwargs)
#         self.load_meta_data()
#         file_obj.close()
#
#     def report_meta_data(self, delimiter="\n", line_delimiter=" | ", header=True):
#         report = []
#         report.append(super(TwoDirectionSteppingBoresightFileAnalyzer, self).report_meta_data(
#             delimiter=delimiter, line_delimiter=line_delimiter, header=header
#         ))
#         if "settle_time" in self.meta_data:
#             line_header = ["settle_time", "integration_time"]
#             line = ["{}={}".format(i, str(self.meta_data[i]))
#                     for i in ["settle_time", "integration_time"]]
#             if header:
#                 report.append(line_delimiter.join(line_header))
#             report.append(line_delimiter.join(line))
#
#         return delimiter.join(report)
#
#     def plot(self, channel="0", save_dir=None, overwrite=False, axs=None):
#         save_file_path, save_file_name = self._plot_save_path(save_dir)
#         if not self._check_existing(save_file_path, overwrite):
#             return
#         fig = None
#         if axs is None:
#             fig, axs = plt.subplots(2, 2, figsize=(10, 10/1.3))
#         self.data["el"]["right"].plot(ax=axs[0][0], channel=channel)
#         self.data["el"]["left"].plot(ax=axs[0][1], channel=channel)
#         self.data["xel"]["right"].plot(ax=axs[1][0], channel=channel)
#         self.data["xel"]["left"].plot(ax=axs[1][1], channel=channel)
#         if fig is not None:
#             suptitle = save_file_name + "\n" +\
#                 self.report_meta_data(header=False,
#                                       line_delimiter=", ",
#                                       delimiter="\n")
#             fig.suptitle(suptitle)
#             top = 0.98-(0.03*float(len(suptitle.split('\n'))))
#             fig.tight_layout(rect=[0, 0.03, 1, top])
#             fig.savefig(save_file_path)
#             plt.close(fig)
#

# for cls in [FluxCalibrationFileAnalyzer,
#             BoresightFileAnalyzer]:
#             # SingleDirectionScanningBoresightFileAnalyzer,
#             # SingleDirectionSteppingBoresightFileAnalyzer,
#             # TwoDirectionSteppingBoresightFileAnalyzer,
#             # TwoDirectionScanningBoresightFileAnalyzer]:
#     if hasattr(cls, "to_dict"):
#         module_logger.debug("registering {} to SerializerBase".format(cls))
#         Pyro4.util.SerializerBase.register_class_to_dict(
#             cls, getattr(cls, "to_dict")
#         )
#         Pyro4.util.SerializerBase.register_dict_to_class(
#             "{}.{}".format(cls.__module__, cls.__name__),
#             getattr(cls, "from_dict")
#         )
