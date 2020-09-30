
import logging
import datetime
import json
try:
    import pickle as pickle
except ImportError as err:
    import pickle

import h5py
import numpy as np
#import matplotlib.pyplot as plt # only used in BoresightFileAnalyzer.plot()

from support import hdf5_util
from .file_analyzer import FileAnalyzer
from .boresight_data_analyzer import (
    ChannelDataAnalyzer, SingleAxisBoresightDataAnalyzer
)

module_logger = logging.getLogger(__name__)

__all__ = ["BoresightFileAnalyzer"]


class BoresightFileAnalyzer(FileAnalyzer):
    """
    """
    date_formatter = "%Y-%j-%Hh%Mm%Ss"
    algorithm_params = {
        "scanning": ["rate", "sample_rate", "limit"],
        "stepping": ["integration_time", "settle_time"]
    }

    def __init__(self, file_path=None, boresight_type=None):
        super(BoresightFileAnalyzer, self).__init__(file_path=file_path)
        self.meta_data.update({
            "boresight_type": boresight_type,
        })
        self.data.update({
            "el": {
                "right": None,
                "left": None,
            },
            "xel": {
                "right": None,
                "left": None,
            }
        })

    @property
    def directions(self):
        """
        Determine whether this boresight run was single direction or
        two direction. If this returns 0, then it means that there is no
        associated with this object.

        Returns:
            int: number of directions (0, 1 or 2)
        """
        el_data = self.data["el"]
        if el_data["right"] is None and el_data["left"] is None:
            return 0
        if el_data["right"] is None or el_data["left"] is None:
            return 1
        else:
            return 2

#   This is not used here or in dss_server.py
    def plot(self, channel="0", save_dir=None, overwrite=True, axs=None):
        """

        Args:
            channel (str): default "0"
            save_dir (str): default None
            overwrite (bool): default True
            axs (list/np.ndarray): default None
        Returns:
            axs
        """
        if self.directions == 0:
            module_logger.error(
                ("BoresightFileAnalyzer.plot: "
                 "no data associated with this object")
            )
        save_file_path, save_file_name = self._plot_save_path(save_dir)
        if not self._check_existing(save_file_path, overwrite):
            return
        fig = None
        if axs is None:
            if self.directions == 1:
                fig, axs = plt.subplots(2, 1, figsize=(10, 10/1.3))
            else:
                fig, axs = plt.subplots(2, 2, figsize=(10, 10/1.3))
        if self.directions == 1:
            self.data["el"]["right"].plot(ax=axs[0], channel=channel)
            self.data["xel"]["right"].plot(ax=axs[1], channel=channel)
        elif self.directions == 2:
            self.data["el"]["right"].plot(ax=axs[0][0], channel=channel)
            self.data["el"]["left"].plot(ax=axs[0][1], channel=channel)
            self.data["xel"]["right"].plot(ax=axs[1][0], channel=channel)
            self.data["xel"]["left"].plot(ax=axs[1][1], channel=channel)

        suptitle = save_file_name + "\n" +\
            self.report_meta_data(header=False,
                                  line_delimiter=", ",
                                  delimiter="\n")
        if fig is not None:
            fig.suptitle(suptitle)
            top = 0.98-(0.03*float(len(suptitle.split('\n'))))
            fig.tight_layout(rect=[0, 0.03, 1, top])
            fig.savefig(save_file_path)
            plt.close(fig)
        return axs

    def report_meta_data(self,
                         delimiter="\n",
                         line_delimiter=" | ",
                         header=True):
        """
        Create a formatted string with information from ``meta_data`` attribute.

        Examples:

        .. code-block::python

            >>> print(a.meta_data)
            {u'el': 0.994504988193512, u'DECJ2000': -0.6363220818328413,
             u'name': '0521-365', 'boresight_type': None,
             'timestamp': '2018-125-07h19m31s', u'RAJ2000': 1.4092068106703093,
             u'initial_xel': 22.9, u'flux': u'3.50', u'rate': 2.0,
             ...}
            >>> print(a.report_meta_data())
            name | az | el
            0521-365 | AZ=255.881 | EL=56.981
            Initial El offset: 4.700 | Initial xEL offset: 22.900
            xt: 0.0 | yt: 0.0 ...

        Args:
            delimiter (str, optional): Delimiter between lines of report.
                Defaults to "\n".
            line_delimiter (str, optional): Delimite inside individual lines
                of report. Defaults to "|".
            header (boolean, optional): Defaults to True.
        Returns:
            str: Nicely formatted string containing useful bits of meta
                data information.
        """

        report = []

        if "name" in self.meta_data:  # means we have source meta data
            convert = 180./np.pi
            line_header = ["name", "az", "el"]
            line = [self.meta_data["name"]]
            line.append(
                "AZ={:.3f}".format(convert * float(self.meta_data["az"])))
            line.append(
                "EL={:.3f}".format(convert * float(self.meta_data["el"])))
            if header:
                report.append(line_delimiter.join(line_header))
            report.append(line_delimiter.join(line))

        if "initial_el" in self.meta_data:
            line = [
                "Initial El offset: {:.3f}".format(
                    self.meta_data["initial_el"]),
                "Initial xEL offset: {:.3f}".format(
                    self.meta_data["initial_xel"])
            ]
            report.append(line_delimiter.join(line))

        if self.meta_data["boresight_type"] in self.algorithm_params:
            alg_line = []
            for param in self.algorithm_params[self.meta_data["boresight_type"]]:
                alg_line.append(
                    "{}: {}".format(
                        param, self.meta_data[param]
                    )
                )
            report.append(line_delimiter.join(alg_line))

        sr_line = []
        for tilt_param in ["xt", "yt", "z", "x", "y"]:
            if tilt_param in self.meta_data:
                sr_line.append(
                    "{}: {}".format(
                        tilt_param, self.meta_data[tilt_param]))
            if tilt_param.upper() in self.meta_data:
                sr_line.append(
                    "{}: {}".format(
                        tilt_param, self.meta_data[tilt_param.upper()]))
        report.append(line_delimiter.join(sr_line))
        return delimiter.join(report)

    def _traverse_data(self, data, callback):
        """
        Apply some callback to every data analyzer object present in ``data``

        Examples:

        .. code-block::python

            def callback(obj):
                print(obj.axis, obj.direction)

            file_analyzer_obj._traverse_data(a.data, callback)

        Might output something like the following:

        .. code-block::none

            EL right
            XEL right
            EL left
            XEL left

        Args:
            data (dict): data through which to move
            callback (callable): callback to apply to each data analyzer object.
                Must take as argument an object of SingleAxisBoresightDataAnalyzer.
        """
        for key in data:
            sub_data = data[key]
            if isinstance(sub_data, dict):
                self._traverse_data(sub_data, callback)
            elif isinstance(sub_data, SingleAxisBoresightDataAnalyzer):
                callback(sub_data)

    def generate_report(self,
                        delimiter="\n",
                        line_delimiter=" | ",
                        header=True):
        """
        """
        report = []

        def generate_report_callback(single_axis_obj):
            sub_report = []
            sub_report.append(
                single_axis_obj.axis)
            sub_report.append(
                single_axis_obj.direction)
            sub_report.append(
                single_axis_obj.generate_report(delimiter=line_delimiter))
            sub_report_str = line_delimiter.join(sub_report)
            report.append(sub_report_str)
        report.append(
            self.report_meta_data(
                delimiter=delimiter,
                line_delimiter=line_delimiter,
                header=header
            )
        )
        self._traverse_data(self.data, generate_report_callback)
        return delimiter.join(report)

    def calculate_offsets(self, **kwargs):
        """
        """
        def calculate_offsets_callback(single_axis_obj):
            single_axis_obj.calculate_offsets(**kwargs)
        self._traverse_data(self.data, calculate_offsets_callback)

    # def __getattr__(self, attr):
    #     if attr in self.meta_data:
    #         return self.meta_data[attr]
    #     else:
    #         raise AttributeError(
    #             "{} doesn't have attr {}".format(
    #                 self.__class__.__name__, attr))

    def __iter__(self):
        for axis in self.data:
            for direction in self.data[axis]:
                yield axis, direction, self.data[axis][direction]

    def __getitem__(self, item):

        attrs = ["data", "meta_data"]
        if item in attrs:
            return getattr(self, item)
        raise KeyError("{} not in {}".format(item, attrs))

    def __contains__(self, item):
        if item in self.meta_data:
            return True
        else:
            return False

    # file loading/dumping
    def _dump_json(self):
        with open(self.file_path, "w") as f:
            json.dump(self.to_dict(), f)

    def _dump_hdf5(self):
        data_dict = {}
        # data_dict = self.to_dict()["data"]
        # for axis in data_dict.keys():
        #     data_dict_axis = data_dict[axis]
        #     for direction in data_dict_axis.keys():
        #         if not data_dict[axis][direction]:
        #             del data_dict[axis][direction]

                # for channel in data_dict[axis][direction]["channels"]:
                #     del data_dict[axis][direction]["channels"][channel]["fit"]
        for axis in self.data:
            data_dict[axis] = {}
            for direction in self.data[axis]:
                if self.data[axis][direction] is not None:
                    single_axis_data = {
                        "offset": self.data[axis][direction].offset_data,
                        "tsys": self.data[axis][direction].tsys_data
                    }
                    data_dict[axis][direction] = single_axis_data

        with h5py.File(self.file_path, "w") as f:
            for attr_name in self.meta_data:
                attr = self.meta_data[attr_name]
                module_logger.debug("_dump_hdf5: attr_name: {}, attr: {}, type(attr): {}".format(
                    attr_name, attr, type(attr)
                ))
                if hasattr(attr, "strftime"):
                    attr = attr.strftime(self.date_formatter)
                f.attrs[attr_name] = attr
            hdf5_util.dump_dict(f, data_dict)

    def _dump_pickle(self):
        with open(self.file_path, "w") as f:
            pickle.dump(self.to_dict(), f)

    def _get_timestamp_hdf5(self):

        def create_timestamp_obj(timestamp_str):
            try:
                return datetime.datetime.strptime(
                    timestamp_str, self.date_formatter
                )
            except ValueError as err:
                return None

        timestamp_str = self.file_name.strip(".hdf5").split("_")[-1]
        timestamp_obj = create_timestamp_obj(timestamp_str)
        if timestamp_obj is None:
            timestamp_str = "2018-" + timestamp_str
            timestamp_obj = create_timestamp_obj(timestamp_str)

        return {
            "timestamp": timestamp_str,
            "timestamp_obj": timestamp_obj
        }

    def _load_data_hdf5(self, file_obj, axis,
                        direction="right", offset_key=None, tsys_key=None):
        """
        Get system temperature and offset data for a given offset axis and
        direction.

        Args:
            file_obj (h5py.File): The file object containing data.
            axis (str): position offset axis, either "EL" or "XEL".
            direction (str, optional): "left" or "right" boresight direction.
                Defaults to "right".
            offset_key (list, optional): Key in ``file_obj`` for offset data.
                Defaults to ["points", "offset"]
            tsys_key (list, optional): Key in ``file_obj`` for
                system temperature data. Defaults to ["tsys"]
        Returns:
            SingleAxisBoresightDataAnalyzer or None
        """
        if offset_key is None:
            offset_key = ["points", "offset"]
        if hasattr(offset_key, "format"):
            offset_key = [offset_key]

        if tsys_key is None:
            tsys_key = ["tsys"]
        if hasattr(tsys_key, "format"):
            tsys_key = [tsys_key]

        if axis not in list(file_obj.keys()):
            axis = axis.lower()

        tsys = []
        points = []
        single_axis_obj = None

        try:
            for key in offset_key:
                if key in file_obj[axis][direction]:
                    points.extend(file_obj[axis][direction][key][...])
                    break

            for key in tsys_key:
                if key in file_obj[axis][direction]:
                    tsys.extend(file_obj[axis][direction][key][...])
                    break
        except KeyError as err:
            return None

        if len(tsys) != 0 and len(points) != 0:
            single_axis_obj = SingleAxisBoresightDataAnalyzer(
                points, tsys,
                axis=axis.lower(),
                direction=direction
            )

        return single_axis_obj

    def _load_meta_data_hdf5(self, file_obj):
        meta_data = {}
        for attr in file_obj.attrs:
            meta_data[attr] = file_obj.attrs[attr]
        return meta_data

    def _load_hdf5(self, **kwargs):

        with h5py.File(self.file_path, "r") as file_obj:

            self.meta_data.update(self._load_meta_data_hdf5(file_obj))

            if "timestamp" not in self.meta_data:
                self.meta_data.update(self._get_timestamp_hdf5())

            self.data["el"]["right"] = self._load_data_hdf5(
                file_obj, "EL", "right", **kwargs)
            self.data["el"]["left"] = self._load_data_hdf5(
                file_obj, "EL", "left", **kwargs)
            self.data["xel"]["right"] = self._load_data_hdf5(
                file_obj, "XEL", "right", **kwargs)
            self.data["xel"]["left"] = self._load_data_hdf5(
                file_obj, "XEL", "left", **kwargs)

        self.calculate_offsets()

    def _load_json(self, **kwargs):
        with open(self.file_path, "r") as f:
            obj_dict = json.load(f)
        self = self.from_dict("", obj_dict)

    def _load_pickle(self, **kwargs):
        with open(self.file_path, "r") as f:
            obj_dict = pickle.load(f)
        self = self.from_dict("", obj_dict)

    # Serialization
    def to_dict(self):
        obj_dict = {}
        obj_dict["__module__"] = self.__module__
        obj_dict["__class__"] = "{}.{}".format(
            self.__module__, self.__class__.__name__)
        obj_dict["file_path"] = self.file_path
        obj_dict["file_name"] = self.file_name
        obj_dict["meta_data"] = self.meta_data.copy()
        if "timestamp_obj" in obj_dict["meta_data"]:
            timestamp_obj = obj_dict["meta_data"]["timestamp_obj"]
            timestamp_str = None
            if hasattr(timestamp_obj, "strftime"):
                timestamp_str = timestamp_obj.strftime(self.date_formatter)
            obj_dict["meta_data"]["timestamp_obj"] = timestamp_str

        obj_dict["data"] = {}
        for axis in self.data:
            obj_dict["data"][axis] = {}
            for direction in self.data[axis]:
                if hasattr(self.data[axis][direction], "to_dict"):
                    obj_dict["data"][axis][direction] = \
                        self.data[axis][direction].to_dict()
                else:
                    obj_dict["data"][axis][direction] = None
        return obj_dict

    @classmethod
    def from_dict(cls, cls_name, obj_dict):
        obj = cls(obj_dict["file_path"])
        obj.meta_data = obj_dict["meta_data"].copy()
        for axis in obj_dict["data"]:
            if axis not in obj.data:
                obj.data[axis] = {}
            for direction in obj_dict["data"][axis]:
                if obj_dict["data"][axis][direction] is not None:
                    obj.data[axis][direction] = \
                        SingleAxisBoresightDataAnalyzer.from_dict(
                            "",
                            obj_dict["data"][axis][direction])
                else:
                    obj.data[axis][direction] = None
        return obj
