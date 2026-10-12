"""Reject auxiliary distributions and unexpected upload artifacts before publishing."""
from pathlib import Path
import sys
import zipfile

import pytest

from scripts.check_wheel import check_wheel, main


def test_rejects_an_auxiliary_distribution(tmp_path):
    wheel = tmp_path / 'parascape_core-0.42.0-py3-none-any.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('parascape_core-0.42.0.dist-info/METADATA', 'Name: parascape-core\nVersion: 0.42.0\n')
    with pytest.raises(ValueError, match='Unexpected distribution'):
        check_wheel(wheel)


def test_rejects_missing_bundled_data(tmp_path):
    wheel = tmp_path / 'parascape-0.42.0-py3-none-any.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('parascape-0.42.0.dist-info/METADATA', 'Name: parascape\nVersion: 0.42.0\n')
    with pytest.raises(ValueError, match='Incomplete wheel'):
        check_wheel(wheel)


@pytest.mark.parametrize('extra', ['parascape_gpu-0.42.0.tar.gz', 'parascape-0.41.0.tar.gz', 'notes.txt'])
def test_rejects_extra_upload_artifacts(tmp_path, monkeypatch, extra):
    for name in ('parascape-0.42.0-py3-none-any.whl', 'parascape-0.42.0.tar.gz', extra):
        (tmp_path / name).touch()
    monkeypatch.setattr(sys, 'argv', ['check_wheel.py', str(tmp_path)])
    with pytest.raises(ValueError, match='only the matching'):
        main()
