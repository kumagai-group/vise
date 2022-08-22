# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import csv
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Dict, List, Any

import yaml
from monty.serialization import loadfn


class ToFileMixIn(ABC):
    @property
    def _filename(self):
        """ ClassForThis -> class_for_this
        https://stackoverflow.com/questions/7322028/how-to-replace-uppercase-with-underscore
        """
        class_name = self.__class__.__name__
        return re.sub("(?<!^)(?=[A-Z])", "_", class_name).lower()


class ToJsonFileMixIn(ToFileMixIn, ABC):
    def to_json_file(self, filename: Optional[str] = None) -> None:
        filename = filename or self._json_filename
        Path(filename).write_text(self.to_json())

    @abstractmethod
    def to_json(self):
        pass

    @property
    def _json_filename(self):
        return self._filename + ".json"


class ToCsvFileMixIn(ToFileMixIn, ABC):

    def to_csv_file(self,
                    filename: Optional[str] = None,
                    ) -> None:
        filename = Path(filename or self._csv_filename)
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(self.csv_column_names)
            writer.writerows(self.csv_data)

    @property
    @abstractmethod
    def csv_column_names(self) -> List[Any]:
        pass

    @property
    @abstractmethod
    def csv_data(self) -> List[List[Any]]:
        pass

    @property
    def _csv_filename(self) -> str:
        return self._filename + ".csv"


class ToYamlFileMixIn(ToFileMixIn, ABC):

    def to_yaml_file(self, filename: Optional[str] = None) -> None:
        filename = filename or self._yaml_filename
        Path(filename).write_text(self.to_yaml())

    def to_yaml(self):
        return yaml.dump(self.as_dict())

    @abstractmethod
    def as_dict(self):
        pass

    @classmethod
    def from_yaml(cls, filename: str = None):
        d = loadfn(filename or cls._yaml_filename)
        if hasattr(cls, "from_dict"):
            return cls.from_dict(d)
        return cls(**d)

    @property
    def _yaml_filename(self):
        """ ClassForThis -> class_for_this.json
        https://stackoverflow.com/questions/7322028/how-to-replace-uppercase-with-underscore
        """
        return self._filename + ".yaml"


