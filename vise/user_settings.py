# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import re
from pathlib import Path
from typing import List

import yaml
from vise.util.logger import get_logger
logger = get_logger(__name__)


class UserSettings:
    def __init__(self, yaml_filename: str):
        self._cwd = Path.cwd()
        self.yaml_filename = yaml_filename
        self.yaml_files_from_root_dir = self._make_yaml_file_list()

    def _make_yaml_file_list(self) -> List[Path]:
        result = []

        dirname = self._cwd
        while True:
            filenames = [self.yaml_filename, "." + self.yaml_filename]
            file_paths = [dirname / filename for filename in filenames]
            for file_path in file_paths:
                if file_path.exists():
                    result.append(file_path)

            if dirname == Path("/"):
                break
            else:
                dirname = dirname.parent

        return list(reversed(result))

    @property
    def user_settings(self) -> dict:
        result = {}

        for file_path in self.yaml_files_from_root_dir:
            logger.info(f"Setting file: {file_path} is parsed...")
            try:
                with open(str(file_path), "r") as fin:
                    settings = yaml.load(fin, Loader=yaml.SafeLoader)
                    for k, v in settings.items():
                        if k in result:
                            logger.info(
                                f"key {k} was overridden to {v} by {file_path}")
                            if isinstance(v, dict):
                                result[k].update(v)
                                continue

                        if self.is_path(v):
                            v = file_path.parent / v
                        result[k] = v
            except AttributeError:
                pass

        logger.info(f"-- Settings from {self.yaml_filename}:")
        logger.info(", ".join(f"{k}: {v} " for k, v in result.items()))
        return result

    @staticmethod
    def is_path(key):
        return isinstance(key, str) and re.match(r'\S*/\S*', key)
