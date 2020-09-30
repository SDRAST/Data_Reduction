
import datetime
import logging

import h5py

from .file_analyzer import FileAnalyzer
from .flux_calibration_data_analyzer import FluxCalibrationDataAnalyzer

module_logger = logging.getLogger(__name__)


class FluxCalibrationFileAnalyzer(FileAnalyzer):

    def __init__(self, file_path):
        super(FluxCalibrationAnalyzer, self).__init__()
        self.file_path = file_path
        self.file_name = os.path.basename(self.file_path)
        self.timestamp = "2018-" + \
            self.file_name.replace(".hdf5", "").split("_")[-1]
        self.timestamp_obj = datetime.datetime.strptime(
            self.timestamp, "%Y-%j-%Hh%Mm%Ss")
        self.meta_data["timestamp"] = self.timestamp
        self.meta_data["file_path"] = self.file_path
        self.meta_data["file_name"] = self.file_name

    def load_data(self):
        calib_data = {}
        with h5py.File(self.file_path, "r") as f:
            for key in list(f.keys()):
                calib_data[key] = f[key][...]
        self.data = FluxCalibrationDataAnalyzer(calib_data)

    def load_meta_data(self):
        pass

    def report_meta_data(self, header=False,
                         line_delimiter=", ", delimiter="\n"):
        return ""

    # this is not used here
    def plot(self, save_dir=None, overwrite=True, ax=None):
        save_file_path, save_file_name = self._plot_save_path(save_dir)
        if not self._check_existing(save_file_path, overwrite):
            return
        fig = None
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10/1.3))
        self.data.plot(ax=ax)
        suptitle = save_file_name + "\n" + self.report_meta_data(
            header=False, line_delimiter=", ", delimiter="\n")
        if fig is not None:
            fig.suptitle(suptitle)
            top = 0.98-(0.03*float(len(suptitle.split('\n'))))
            fig.tight_layout(rect=[0, 0.03, 1, top])
            fig.savefig(save_file_path)
            plt.close(fig)
