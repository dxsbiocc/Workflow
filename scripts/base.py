# -*- encoding: utf-8 -*-
# ============================================================
# File        : base.py
# Time        : 2022/08/17 17:15:11
# Author      : dengxsh 
# Version     : 1.0
# Contact     : 920466915@qq.com
# Copyright   : Copyright (c) 2022, dengxsh
# License     : MIT
# Description : The role of the current file 
# ============================================================


from abc import ABCMeta, abstractmethod
import logging
try:
    import colorlog
except:
    COLOR = False
else:
    COLOR = True
    
class Filter(logging.Filter):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        # color map
        self.color_map = {
            'ERROR'   : "\033[31mERROR\033[0m",
            'INFO'    : "\033[37mINFO\033[0m",
            'DEBUG'   : "\033[1mDEBUG\033[0m",
            'WARN'    : "\033[33mWARN\033[0m",
            'WARNING' : "\033[33mWARNING\033[0m",
            'CRITICAL': "\033[35mCRITICAL\033[0m",
        }

    def filter(self, record: logging.LogRecord) -> bool:
        record.levelname = self.color_map.get(record.levelname)
        return True
    
def get_logger(
    name: str, 
    level: int    = logging.INFO,
    fmt: str      = "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    fmt_date: str = "%Y-%m-%d %H:%M:%S %p",
    stream = None,
    filename = None
    ) -> logging.Logger:
    # color handler
    assert not (filename and stream), "Cannot set 'stream' and 'filename' parameter in the same time!"

    if filename:
        handler = logging.FileHandler(filename=filename, mode='w')
    else:
        handler = logging.StreamHandler(stream=stream)
    
    if COLOR:
        cfmt = '%(log_color)s' + fmt
        fmter = colorlog.ColoredFormatter(
            cfmt, fmt_date, 
            log_colors={
                'DEBUG':    'cyan',
                'INFO':     'green',
                'WARNING':  'yellow',
                'ERROR':    'red',
                'CRITICAL': 'purple'
            }
        )
        handler.setFormatter(fmter)
    else:
        fmter = logging.Formatter(fmt, fmt_date)
        handler.setLevel(level)
        handler.setFormatter(fmter)
        # filter
        my_filter = Filter()
        handler.addFilter(my_filter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    
    return logger

class WrapperBase(metaclass=ABCMeta):

    def __init__(self, snakemake) -> None:
        self.snakemake = snakemake
        # common parameters
        self.extra = snakemake.params.get("extra", "")
        self.tmp = snakemake.params.get("tmp_dir", "/tmp/")
        # log
        self.log = snakemake.log_fmt_shell(stdout=True, stderr=True)
        self.parser()
        self.run()
    
    @abstractmethod
    def parser(self):
        """Extract all informations from input、output and parameterss
        """
        pass

    @abstractmethod
    def run(self):
        """Run shell command
        """
        pass