import pytest

import runner


def test_find_tool_present():
    # python3 is guaranteed available in the test environment
    assert runner.find_tool("python3") is not None


def test_find_tool_missing():
    assert runner.find_tool("definitely-not-a-real-binary-xyz") is None


def test_find_tool_prefers_first():
    path = runner.find_tool("python3", "sh")
    assert path is not None and path.endswith("python3")


def test_default_seed_is_42_without_env(monkeypatch):
    monkeypatch.delenv("VS_SEED", raising=False)
    assert runner.default_seed() == 42


def test_default_seed_reads_vs_seed(monkeypatch):
    monkeypatch.setenv("VS_SEED", "7")
    assert runner.default_seed() == 7


def test_default_seed_rejects_non_integer(monkeypatch):
    monkeypatch.setenv("VS_SEED", "abc")
    with pytest.raises(SystemExit):
        runner.default_seed()
