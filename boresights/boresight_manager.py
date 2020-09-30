import os
import json
import logging
import multiprocessing
import time
import functools
import datetime
import Pyro5

from local_dirs import *
from support.pyro.pyro5_server import Pyro5Server
from MonitorControl.DSS_server_cfg import tams_config
from .analyzer import (BoresightFileAnalyzer,
                       FluxCalibrationFileAnalyzer)

# from ..configuration import config

__all__ = [
    "Manager",
    "FluxCalibrationManager",
    "BoresightManager",
    "BoresightManagerByDate"
]

module_logger = logging.getLogger(__name__)


def timeit(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t0 = time.time()
        res = func(*args, **kwargs)
        module_logger.debug("{}: Took {:.2f} seconds".format(
            func.__name__, time.time() - t0
        ))
        return res
    return wrapper


class Manager(Pyro5Server):

    cache_file_name = "cache.json"
    cls_lookup = {
        "BoresightFileAnalyzer":
            BoresightFileAnalyzer,
        "FluxCalibrationFileAnalyzer":
            FluxCalibrationFileAnalyzer
    }

    def __init__(self, file_paths, **kwargs):
        super(Manager, self).__init__(obj=self)
        self.logger = module_logger.getChild(self.__class__.__name__)
        file_paths = Manager.squash_file_paths(file_paths)
        # self.logger.debug("__init__: file_paths: {}".format(file_paths))
        self.logger.debug("__init__: len(file_paths): {}".format(
            len(file_paths)))
        self.file_paths = file_paths
        self._analyzer_objects = self.create_analyzer_objects(self.file_paths,
                                                              **kwargs)

    @Pyro5.api.expose
    def iter(self):
        for obj in self:
            yield obj

    @Pyro5.api.expose
    @property
    def analyzer_objects(self):
        return self._analyzer_objects

    @timeit
    def _dump_analyzer_objects_to_cache(self, analyzer_objects,
                                        cache_file_path=None):
        if cache_file_path is None:
            cache_file_path = os.path.join(data_dir,
                                           self.cache_file_name)
        cache_dict = {}
        for file_name in analyzer_objects:
            obj = analyzer_objects[file_name]
            if hasattr(obj, "to_dict"):
                cache_dict[file_name] = obj.to_dict()
            else:
                cache_dict[file_name] = obj
        with open(cache_file_path, "w") as f:
            json.dump(cache_dict, f)

    @timeit
    def _load_cached_analyzer_objects(self):
        cache_file_path = os.path.join(tams_config.boresight_data_dir,
                                       self.cache_file_name)
        if not os.path.exists(cache_file_path):
            return {}
        with open(cache_file_path, "r") as f:
            try:
                cache_dict = json.load(f)
            except ValueError as err:
                module_logger.error(
                    ("_load_cached_analyzer_objects: "
                     "Couldn't load JSON object: {}".format(err)))
                return {}
        analyzer_objects = {}
        for file_name in cache_dict:
            obj_dict = cache_dict[file_name]
            if "__class__" not in obj_dict:
                obj = obj_dict.copy()
            else:
                cls_name = obj_dict["__class__"]
                cls_name = cls_name.split(".")[-1]
                cls = self.cls_lookup[cls_name]
                obj = cls.from_dict("", obj_dict)
            analyzer_objects[file_name] = obj
        return analyzer_objects

    def create_analyzer_objects(self, file_paths,
                                reload=True,
                                multiprocessed=False,
                                dump_cache=True):
        raise NotImplementedError("Subclasses should implement this method")

    def add_analyzer_object(self, file_path):

        objs = self.create_analyzer_objects(
            [file_path],
            reload=False,
            multiprocessed=False,
            dump_cache=False
        )

        self._analyzer_objects.update(objs)

    def generate_report(self, delimiter="\n"):
        self.logger.debug("generate_report: called")
        report = [
            obj.generate_report(header=False) for obj in self._analyzer_objects
        ]
        return delimiter.join(report)

    def save_markdown_report(self, **kwargs):
        """
        Takes as save directory the directory where the first file is located.
        """
        save_dir = os.path.dirname(self.file_paths[0])
        timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%Hh%Mm%Ss")
        markdown_file_name = "report_{}.md".format(timestamp)
        markdown_file_path = os.path.join(save_dir, markdown_file_name)
        report_str = self.generate_report()
        with open(markdown_file_path, "w") as md_file:
            md_file.write(report_str)

    def plot(self, **kwargs):
        self.logger.debug("plot: called")
        for obj in self:
            if hasattr(obj, "plot"):
                obj.plot(**kwargs)

    def _raise_analyzer_objects_not_initialized(self):
        raise RuntimeError(("analyzer_objects is not yet defined."))

    def __iter__(self):
        if self._analyzer_objects is None:
            self._raise_analyzer_objects_not_initialized()
        for file_name in self._analyzer_objects:
            yield self._analyzer_objects[file_name]

    def __getitem__(self, item):
        if self._analyzer_objects is None:
            self._raise_analyzer_objects_not_initialized()
        return self._analyzer_objects[item]

    def __contains__(self, item):
        if item in self._analyzer_objects:
            return True
        return False

    def __len__(self):
        if self._analyzer_objects is None:
            self._raise_analyzer_objects_not_initialized()
        return len(self._analyzer_objects)

    def keys(self, *args, **kwargs):
        if self._analyzer_objects is None:
            self._raise_analyzer_objects_not_initialized()
        return self._analyzer_objects.keys(*args, **kwargs)

    def __eq__(self, manager_obj):
        """
        This will attempt to determine if two manager objects are refering to
        the same files.
        """
        if len(self.file_paths) == len(self.manager_obj.file_paths):
            sorted_file_paths = sorted(self.file_paths)
            sorted_file_paths_other = sorted(manager_obj.file_paths)
            for p1, p2 in zip(sorted_file_paths, sorted_file_paths_other):
                if p1 != p2:
                    break
        return False

    def __ne__(self, manager_obj):
        return not self.__eq__(manager_obj)

    @staticmethod
    def squash_file_paths(file_paths):
        def _squash_file_paths(file_paths, res):
            if isinstance(file_paths, list):
                for file_path in file_paths:
                    res = _squash_file_paths(file_path, res)
                return res
            if os.path.isdir(file_paths):
                paths = [os.path.join(file_paths, f)
                         for f in os.listdir(file_paths)
                         if f.endswith(".hdf5")]
                res.extend(paths)
                return res
            if os.path.isfile(file_paths):
                if file_paths.endswith(".hdf5"):
                    res.append(file_paths)
                return res
        res = _squash_file_paths(file_paths, [])
        return res

    @staticmethod
    def detect_file_paths(base_dir, year=None, doy=None):
        """
        Get file paths for specified year(s) and doy(s). If neither year
        nor doy are provided, then find all file paths for every year and
        every doy.
        """
        def _detect_years(base_dir):
            for yr in os.listdir(base_dir):
                try:
                    datetime.datetime.strptime(yr, "%Y")
                    yield yr
                except ValueError as err:
                    continue

        def _detect_doys(base_dir, year):
            for doy in os.listdir(os.path.join(base_dir, year)):
                try:
                    datetime.datetime.strptime(doy, "%j")
                    yield doy
                except ValueError as err:
                    continue

        def _detect_hdf5(base_dir, year, doy):
            for file_name in os.listdir(os.path.join(base_dir, year, doy)):
                if file_name.endswith(".hdf5"):
                    yield (
                        file_name,
                        os.path.join(base_dir, year, doy, file_name)
                    )
                else:
                    continue

        if year is None:
            year = list(_detect_years(base_dir))
        if hasattr(year, "format"):
            year = [year]

        if doy is None:
            doy = {yr: list(_detect_doys(base_dir, yr)) for yr in year}
        elif hasattr(doy, "__getitem__"):
            if hasattr(doy, "format"):
                doy = [doy]
            if not hasattr(doy, "keys"):
                doy = {yr: doy for yr in year}

        file_names = {
            yr: {
                d: list(_detect_hdf5(base_dir, yr, d)) for d in doy[yr]
            } for yr in year
        }

        flattened = [
            i[-1] for yr in file_names
            for d in file_names[yr]
            for i in file_names[yr][d]
        ]

        return file_names, flattened


class FluxCalibrationManager(Manager):

    cache_file_name = "cache.flux_calibration.json"

    def __init__(self, file_paths, **kwargs):
        super(self.__class__, self).__init__(file_paths, **kwargs)

    def create_analyzer_objects(self, file_paths):
        flux_objects = []
        for f in file_paths:
            obj = FluxCalibrationAnalyzer(f)
            obj.load_data()
            obj.data.calculate_tsys_factors()
            flux_objects.append(obj)
        self._analyzer_objects = flux_objects
        return flux_objects

    @staticmethod
    def detect_file_paths(**kwargs):
        return Manager.detect_file_paths(
            tams_config.flux_calibration_data_dir, **kwargs)


class BoresightManager(Manager):

    cache_file_name = "cache.boresight.json"

    def __init__(self, file_paths, **kwargs):
        super(BoresightManager, self).__init__(file_paths, **kwargs)

    @timeit
    def create_analyzer_objects(self, file_paths,
                                multiprocessed=False,
                                reload=False,
                                dump_cache=True):

        def create_objects_multiprocessed(file_paths):
            boresight_objects = {}

            def process_data(file_paths, q):
                for i, f in enumerate(file_paths):
                    obj, _, file_name = self._detect_boresight_type(f)
                    q.put((file_name, obj))

            n_processors = multiprocessing.cpu_count()
            module_logger.debug(
                "create_analyzer_objects: Using {} processors".format(
                    n_processors))
            n_files = len(file_paths)
            n_files_per_processor = n_files / n_processors
            module_logger.debug(
                "create_analyzer_objects: Each processor taking care of {} files".format(
                    n_files_per_processor))
            manager = multiprocessing.Manager()
            q = manager.Queue()
            processes = []
            for i in range(n_processors):
                begin, end = i*n_files_per_processor, (i+1)*n_files_per_processor
                module_logger.debug(
                    "create_analyzer_objects: processor {} taking care of {} to {}".format(
                        i, begin, end))
                if i+1 == n_processors:
                    sub_file_paths = file_paths[begin:]
                else:
                    sub_file_paths = file_paths[begin:end]
                p = multiprocessing.Process(
                    target=process_data,
                    args=(sub_file_paths, q)
                )
                p.start()
                processes.append(p)

            for p in processes:
                p.join()

            while not q.empty():
                file_name, obj = q.get()
                boresight_objects[file_name] = obj
            return boresight_objects

        def create_objects_non_multiprocessed(file_paths):
            boresight_objects = {}
            for f in file_paths:
                obj, _, file_name = self._detect_boresight_type(f)
                boresight_objects[file_name] = obj

            return boresight_objects

        boresight_objects = {}
        remaining_boresight_files = []
        if not reload:
            boresight_objects.update(self._load_cached_analyzer_objects())
            existing_file_paths = {boresight_objects[file_name]["meta_data"]["file_path"]
                                   for file_name in boresight_objects}
            remaining_boresight_files = [
                file_path for file_path in file_paths
                if file_path not in existing_file_paths
            ]
        else:
            remaining_boresight_files = file_paths
        if len(remaining_boresight_files) < 200:
            multiprocessed = False

        if multiprocessed:
            boresight_objects.update(
                create_objects_multiprocessed(remaining_boresight_files))
        else:
            boresight_objects.update(
                create_objects_non_multiprocessed(remaining_boresight_files))
        if dump_cache:
            self._dump_analyzer_objects_to_cache(boresight_objects)
        file_paths = set(file_paths)
        ret_dict = {}
        for file_name in boresight_objects:
            obj = boresight_objects[file_name]
            if (not isinstance(obj, dict) and obj["meta_data"]["file_path"] in file_paths):
                ret_dict[file_name] = obj
            # else:
            #     self.logger.debug("couldn't put {} in ret_dict".format(obj))
            #     self.logger.debug(not isinstance(obj, dict))
            #     self.logger.debug(obj["meta_data"]["file_path"] in file_paths)
            #     self.logger.debug(obj["meta_data"]["file_path"])
            #     self.logger.debug(file_paths)
        return ret_dict

    def _detect_boresight_type(self, file_path):
        # self.logger.debug("_detect_boresight_type: {}".format(file_path))
        file_name = os.path.basename(file_path)
        cls = BoresightFileAnalyzer
        boresight_type = "stepping"
        if "scanning" in file_name:
            boresight_type = "scanning"
        default_obj = {"file_path": file_path}
        if cls is None:
            return default_obj, None, file_name
        else:
            try:
                obj = cls(file_path=file_path, boresight_type=boresight_type)
                obj.load()
                obj.meta_data["file_path"] = file_path
                obj.calculate_offsets(channels=["0",
                                                "1",
                                                "2",
                                                "3",
                                                ["0", "-", "2"],
                                                ["1", "-", "3"]])
                self.logger.debug("_detect_boresight_type: obj: {}".format(obj))
            except (IOError, KeyError, IndexError) as err:
                obj = default_obj
            self.logger.debug("_detect_boresight_type: obj: {}, cls: {}, file_name: {}".format(
                obj, cls, file_name
            ))
            return obj, cls, file_name

    def save_json_report(self, **kwargs):
        """
        """
        save_dir = tams_config.product_dir
        timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%Hh%Mm%Ss")
        json_file_name = "report_{}.json".format(timestamp)
        json_file_path = os.path.join(save_dir, json_file_name)
        self._dump_analyzer_objects_to_cache(cache_file_path=json_file_path)

    def get_paths(self, *paths, **kwargs):
        """
        Iterating through each boresight object, get a specific path in
        either meta_data or fits.
        """
        filter_fn = kwargs.get("filter_fn", None)

        def default_filter_fn(obj):
            return True

        def default_transport(x):
            return x

        if filter_fn is None:
            filter_fn = default_filter_fn
        key = kwargs.get("key", None)
        if key is None:
            objs = self._analyzer_objects
        else:
            objs = sorted(self._analyzer_objects, key=key)

        res = []
        for obj_name in objs:
            obj = objs[obj_name]
            if not filter_fn(obj):
                continue
            obj_res = []
            for i in range(len(paths)):
                path = paths[i]
                if isinstance(path, dict):
                    transport = path["transport"]
                    path = path["path"]
                else:
                    transport = default_transport
                traverse = obj[path[0]]
                for p in path[1:]:
                    if traverse is not None:
                        if p in traverse:
                            traverse = traverse[p]
                        elif p.upper() in traverse:
                            traverse = traverse[p.upper()]
                        else:
                            traverse = None
                obj_res.append(transport(traverse))
            if len(obj_res) == len(paths):
                res.append(obj_res)
        return res

    @staticmethod
    def detect_file_paths(**kwargs):
        return Manager.detect_file_paths(
            tams_config.boresight_data_dir, **kwargs)


class BoresightManagerByDate(BoresightManager):
    """
    Class that allows us to simple grab boresight by DOY and year.
    """
    def __init__(self, year, doys, **kwargs):
        self.year = year
        if not isinstance(doys, list):
            doys = [doys]
        doys = ["{:03d}".format(int(doy)) for doy in doys]
        self.doys = []
        file_paths = []
        for doy in doys:
            doy_path = os.path.join(tams_config.boresight_data_dir, year, doy)
            if os.path.exists(doy_path):
                self.doys.append(doy)
                file_paths.append(doy_path)
        self.products_dir = os.path.join(
            tams_config.boresight_data_dir, year, "products")
        module_logger.debug(
            "BoresightManagerByDate.__init__: products_dir: {}".format(
                self.products_dir))
        super(BoresightManagerByDate, self).__init__(file_paths, **kwargs)
