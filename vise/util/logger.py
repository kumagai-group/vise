# -*- coding: utf-8 -*-
import logging
import sys

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def get_logger(name: str,
               level=logging.DEBUG,
               stream=sys.stdout,
               log_format: str = None,
               log_filename: str = None):
#    log_filename: str = "vise.log"):

    log_format = log_format or \
                 '%(asctime)25s %(levelname)7s %(name)25s \n --> %(message)s'

    logger = logging.getLogger(name)
    logger.setLevel(level)

    formatter = logging.Formatter(log_format)

    if log_filename:
        fh = logging.FileHandler(log_filename)
        fh.setFormatter(formatter)
        fh.setLevel(logging.WARNING)
        logger.addHandler(fh)

    sh = logging.StreamHandler(stream=stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    return logger
