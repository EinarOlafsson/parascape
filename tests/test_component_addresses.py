"""The address recorded for a newly built component part: Hugging Face first, GitHub second."""
import importlib.util
from pathlib import Path

import pytest

from parascape import data_components as C

spec = importlib.util.spec_from_file_location('component_builder_addresses', Path(__file__).resolve().parents[1] / 'scripts/build_data_components.py')
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)

IDENTITY = 'a' * 64
KIND = next(iter(C.KINDS))


def test_new_part_records_hugging_face_address():
    first, _ = B.part_addresses(KIND, 'part-000.tar', 'v0.0.13.0')
    assert first == f'{C.COMPONENT_HOST}/{KIND}/part-000.tar'
    assert first.startswith('https://huggingface.co/')


def test_fallback_order_is_hugging_face_then_github():
    first, second = B.part_addresses(KIND, 'part-000.tar', 'v0.0.9.0')
    assert 'huggingface.co' in first
    assert second == f'{C.GITHUB_RELEASES}/v0.0.9.0/part-000.tar'


def test_both_addresses_pass_the_component_check():
    for url in B.part_addresses(KIND, 'part-000.tar', 'v0.0.9.0'):
        assert C._part_url(url) == url


def test_unknown_host_and_plain_http_are_refused():
    for bad in ('http://huggingface.co/datasets/einarolafsson/parascape-data/resolve/main/components/x/y',
                f'https://example.com/{KIND}/part-000.tar'):
        with pytest.raises(ValueError):
            C._part_url(bad)


def test_recorded_address_is_tried_first_and_validates():
    url = B.part_addresses(KIND, 'p-000.tar', 'v0.0.13.0')[0]
    component = dict(directory=f'{KIND}/{IDENTITY}', files={'a.bin': dict(size=10, sha256='c' * 64)},
                     archive_size=10, archive_sha256='d' * 64,
                     parts=[dict(url=url, size=10, sha256='b' * 64)])
    C.validate(component)
    assert C.part_sources(component, component['parts'][0])[0] == url
