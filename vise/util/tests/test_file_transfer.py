# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import os
from pathlib import Path
import tempfile

import pytest

from vise.util.file_transfer import (
    AFileTransfer, AFileMove, AFileCopy, AFileLink, FileTransfers, A)


cwd = Path.cwd().absolute()


@pytest.fixture()
def file_transfer_instance():
    class Test(AFileTransfer):
        def transfer(self, abs_output_dir):
            pass

    return Test(abs_file_path=cwd / "a")


def test_file_transfer_abs_file_path(file_transfer_instance):
    assert file_transfer_instance.abs_in_file == cwd / "a"


def test_file_transfer_file_name(file_transfer_instance):
    assert file_transfer_instance.file_name == "a"


def test_file_transfer_to(file_transfer_instance):
    assert file_transfer_instance.to(Path("b")) == cwd / "b" / "a"


def test_file_move():
    print_string = "test"
    filename = "a"
    with tempfile.TemporaryDirectory() as tmp_from:
        file_move = AFileMove(abs_file_path=Path(tmp_from) / filename)
        os.chdir(tmp_from)

        with open(filename, "w") as f1:
            print(print_string, end="", file=f1)

        with tempfile.TemporaryDirectory() as tmp_to:
            file_move.transfer(Path(tmp_to))
            os.chdir(Path(tmp_to))

            with open(Path(tmp_to) / filename, 'r') as f2:
                assert f2.read() == print_string


def test_file_copy():
    print_string = "test"
    filename = "a"
    with tempfile.TemporaryDirectory() as tmp_from:
        file_copy = AFileCopy(abs_file_path=Path(tmp_from) / filename)
        os.chdir(tmp_from)

        with open(filename, "w") as f1:
            print(print_string, end="", file=f1)

        with tempfile.TemporaryDirectory() as tmp_to:
            file_copy.transfer(Path(tmp_to))
            os.chdir(Path(tmp_to))

            with open(Path(tmp_to) / filename, 'r') as f_to:
                with open(Path(tmp_from) / filename, 'r') as f_from:
                    assert f_to.read() == f_from.read()


def test_file_link():
    print_string = "test"
    filename = "a"
    with tempfile.TemporaryDirectory() as tmp_from:
        file_link = AFileLink(abs_file_path=Path(tmp_from) / filename)
        os.chdir(tmp_from)

        with open(filename, "w") as f1:
            print(print_string, end="", file=f1)

        with tempfile.TemporaryDirectory() as tmp_to:
            file_link.transfer(Path(tmp_to))

            assert (Path(tmp_to) / filename).is_symlink()

            os.chdir(Path(tmp_to))

            with open(Path(tmp_to) / filename, 'r') as f_to:
                with open(Path(tmp_from) / filename, 'r') as f_from:
                    assert f_to.read() == f_from.read()


def test_transfer_files():
    print_string = "test"
    with tempfile.TemporaryDirectory() as tmp_from:
        os.chdir(tmp_from)
        for i in ["a", "b", "c", "d"]:
            with open(i, "w") as f1:
                print(print_string, end="", file=f1)
        a = FileTransfers.from_dict({"a": "m",
                                     "b": "c",
                                     "c": "l",
                                     "d": "m"}, path=Path(tmp_from))

        with tempfile.TemporaryDirectory() as tmp_to:
            a.transfer(Path(tmp_to))
            for i in ["b", "c"]:
                with open(Path(tmp_to) / i, 'r') as f_to:
                    with open(Path(tmp_from) / i, 'r') as f_from:
                        assert f_to.read() == f_from.read()
            for i in ["a", "d"]:
                with open(Path(tmp_to) / i, 'r') as f2:
                    assert f2.read() == print_string


def test_transfer_files_logger_non_exist(mocker):
    mock = mocker.patch("vise.util.file_transfer.logger.warning")
    non_exist_filename = "a"

    with tempfile.TemporaryDirectory() as tmp_from:
        os.chdir(tmp_from)
        FileTransfers.from_dict({non_exist_filename: "m"}, path=Path(tmp_from))
        mock.assert_called_once_with(
            f"{Path(tmp_from) / non_exist_filename} does not exist.")


def test_transfer_files_logger_empty(mocker):
    mock = mocker.patch("vise.util.file_transfer.logger.warning")
    empty_filename = "b"

    with tempfile.TemporaryDirectory() as tmp_from:
        os.chdir(tmp_from)
        with open(empty_filename, "w") as f1:
            print("", end="", file=f1)
        FileTransfers.from_dict({empty_filename: "m"}, path=Path(tmp_from))
        mock.assert_called_once_with(
            f"{Path(tmp_from) / empty_filename} is empty.")


def test_transfer_files_logger_empty(mocker):
    mock = mocker.patch("vise.util.file_transfer.logger.warning")
    filename = "c"
    print_string = "test"

    with tempfile.TemporaryDirectory() as tmp_from:
        os.chdir(tmp_from)
        with open(filename, "w") as f2:
            print(print_string, end="", file=f2)
        FileTransfers.from_dict({filename: "x"}, path=Path(tmp_from))
        mock.assert_called_once_with(
            f"x option for {filename} is invalid.")

