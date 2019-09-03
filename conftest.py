import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--suite", default="all", choices=["blast", "nonblast", "all"], help="test suite to run"
    )

def pytest_collection_modifyitems(config, items):
    suite = config.getoption("--suite")
    if suite == "all":
        # Don't skip any tests
        return
    skip = pytest.mark.skip(reason="skipping non-selected suite")
    for item in items:
        if suite != "blast" and "blast" in item.keywords:
            item.add_marker(skip)
        if suite == "blast" and "blast" not in item.keywords:
            item.add_marker(skip)
