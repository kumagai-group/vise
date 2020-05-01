# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
import shutil
from abc import ABC, abstractmethod
from pathlib import Path
import re
from typing import List, Dict

from monty.io import zopen

from vise.util.logger import get_logger

logger = get_logger(__name__)


class FileTransfer(ABC):
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


class FileMove(FileTransfer):
    def transfer(self, abs_output_dir: Path):
        shutil.move(self.abs_in_file, self.to(abs_output_dir))


class FileCopy(FileTransfer):
    def transfer(self, abs_output_dir: Path):
        with zopen(self.abs_in_file, "rb") as fin, \
                zopen(self.to(abs_output_dir), "wb") as fout:
            shutil.copyfileobj(fin, fout)


class FileLink(FileTransfer):
    def transfer(self, abs_output_dir: Path):
        os.symlink(self.abs_in_file, self.to(abs_output_dir))


def transfer_instance(initial, filename):
    if initial == "m":
        result = FileMove
    elif initial == "c":
        result = FileCopy
    elif initial == "l":
        result = FileLink
    else:
        raise ViseFileTransferError
    return result(filename)


class FileTransfers:
    def __init__(self, transfer_manner_by_file: Dict[str, str], path: Path):

        file_transfers = []
        for filename, transfer_manner in transfer_manner_by_file.items():
            f = path.absolute() / filename
            initial = transfer_manner[0].lower()

            if not f.is_file():
                logger.warning(f"{f} does not exist.")
            elif f.stat().st_size == 0:
                logger.warning(f"{f} is empty.")
            else:
                file_transfers.append(transfer_instance(initial, f))

        self.file_transfers = file_transfers

    def delete_file_transfers_w_keywords(self, keywords: List[str]):
        pattern = re.compile("|".join(keywords))
        for a_file_transfer in list(self.file_transfers):
            if pattern.search(a_file_transfer.file_name):
                self.file_transfers.remove(a_file_transfer)

    def transfer(self, output_dir: Path) -> None:
        for transfer in self.file_transfers:
            transfer.transfer(output_dir.absolute())


class ViseFileTransferError(Exception):
    pass
