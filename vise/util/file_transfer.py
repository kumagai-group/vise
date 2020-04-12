# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os
import shutil
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
import re
from typing import List, Dict

from monty.io import zopen

from vise.util.logger import get_logger

logger = get_logger(__name__)


class AFileTransfer(ABC):
    def __init__(self, abs_file_path: Path):
        self.abs_in_file = abs_file_path

    @property
    def file_name(self):
        return Path(self.abs_in_file).name

    @abstractmethod
    def transfer(self, abs_output_dir: Path):
        pass

    def to(self, abs_output_dir: Path) -> Path:
        return Path(abs_output_dir.absolute()) / self.file_name


class AFileMove(AFileTransfer):
    def transfer(self, abs_output_dir: Path):
        shutil.move(self.abs_in_file, self.to(abs_output_dir))


class AFileCopy(AFileTransfer):
    def transfer(self, abs_output_dir: Path):
        with zopen(self.abs_in_file, "rb") as fin, \
                zopen(self.to(abs_output_dir), "wb") as fout:
            shutil.copyfileobj(fin, fout)


class AFileLink(AFileTransfer):
    def transfer(self, abs_output_dir: Path):
        os.symlink(self.abs_in_file, self.to(abs_output_dir))


@dataclass
class FileTransfers:
    file_transfer_list: List[AFileTransfer]

    @classmethod
    def from_dict(cls, d: Dict[str, str], path: Path) -> "FileTransfers":
        file_transfers = []
        for filename, transfer_type in d.items():
            f = path.absolute() / filename
            initial = transfer_type[0].lower()
            if not f.is_file():
                logger.warning(f"{f} does not exist.")
            elif f.stat().st_size == 0:
                logger.warning(f"{f} is empty.")
            else:
                if initial == "m":
                    transfer_cls = AFileMove
                elif initial == "c":
                    transfer_cls = AFileCopy
                elif initial == "l":
                    transfer_cls = AFileLink
                else:
                    logger.warning(
                        f"{transfer_type} option for {filename} is invalid.")
                    continue
                file_transfers.append(transfer_cls(f))

        return cls(file_transfers)

    def delete_file_transfers_w_keywords(self, keywords: List[str]):
        pattern = re.compile("|".join(keywords))
        print(pattern)
        for a_file_transfer in list(self.file_transfer_list):
            if pattern.search(a_file_transfer.file_name):
                self.file_transfer_list.remove(a_file_transfer)

    def transfer(self, output_dir: Path) -> None:
        for transfer in self.file_transfer_list:
            transfer.transfer(output_dir.absolute())
