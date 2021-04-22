# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

from monty.serialization import loadfn


class ToJsonFileMixIn(ABC):
    def to_json_file(self, filename: Optional[str] = None) -> None:
        filename = filename or self._json_filename
        Path(filename).write_text(self.to_json())

    @abstractmethod
    def to_json(self):
        pass

    @property
    def _json_filename(self):
        """ ClassForThis -> class_for_this.json
        https://stackoverflow.com/questions/7322028/how-to-replace-uppercase-with-underscore
        """
        class_name = self.__class__.__name__
        return re.sub("(?<!^)(?=[A-Z])", "_", class_name).lower() + ".json"


class ToYamlFileMixIn(ABC):

    def to_yaml_file(self, filename: Optional[str] = None) -> None:
        filename = filename or self._yaml_filename()
        Path(filename).write_text(self.to_yaml())

    @abstractmethod
    def to_yaml(self):
        pass

    @classmethod
    def from_yaml(cls, filename: str = None):
        d = loadfn(filename or cls._yaml_filename())
        return cls(**d)

    @classmethod
    def _yaml_filename(cls):
        """ ClassForThis -> class_for_this.json
        https://stackoverflow.com/questions/7322028/how-to-replace-uppercase-with-underscore
        """
        class_name = cls.__name__
        return re.sub("(?<!^)(?=[A-Z])", "_", class_name).lower() + ".yaml"


