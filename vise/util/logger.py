# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import logging
import sys


def get_logger(name: str,
               level=logging.DEBUG,
               stream=sys.stdout,
               datetime: str = '%Y/%m/%d %H:%M:%S',
               logger_type: str = "simple",
               log_filename: str = None,
               file_handler_level: str = "WARNING"):

    if logger_type == "simple":
        log_format = "%(levelname)7s: %(message)s"
        formatter = logging.Formatter(log_format, datefmt=datetime)
    elif logger_type == "abundant":
        log_format = \
            "%(asctime)18s %(levelname)7s %(name)25s\n --> %(message)s"
        formatter = logging.Formatter(log_format, datefmt=datetime)
    else:
        raise ValueError

    logger = logging.getLogger(name)
    logger.setLevel(level)

    if log_filename:
        file_handler = logging.FileHandler(log_filename)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(file_handler_level)
        logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler(stream=stream)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger
