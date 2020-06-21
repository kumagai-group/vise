# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import re
from copy import deepcopy
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
                    result.update(self._add_absolute_path(settings, file_path))
            except AttributeError:
                pass

        return result

    def _add_absolute_path(self, settings, file_path):
        result = deepcopy(settings)
        for key, value in settings.items():
            if self.is_path(value):
                result[key] = file_path.parent / value

        return result

    @staticmethod
    def is_path(key):
        return isinstance(key, str) and re.match(r'\S*/\S*', key)
