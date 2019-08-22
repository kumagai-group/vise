# -*- coding: utf-8 -*-
import logging
import sys

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


def get_logger(name, level=logging.DEBUG,
               log_format='%(asctime)s %(levelname)s %(name)s %(message)s',
               stream=sys.stdout):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter(log_format)
    sh = logging.StreamHandler(stream=stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger
