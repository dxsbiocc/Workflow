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

    if filename and isinstance(stream, str):
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
    
def close_logger(logger):
    handlers = logger.handlers[:]
    for handler in handlers:
        logger.removeHandler(handler)
        handler.close()
    print("\033[35mlogger closed!\033[0m")