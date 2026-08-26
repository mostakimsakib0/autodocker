import runner


def test_find_tool_present():
    # python3 is guaranteed available in the test environment
    assert runner.find_tool("python3") is not None


def test_find_tool_missing():
    assert runner.find_tool("definitely-not-a-real-binary-xyz") is None


def test_find_tool_prefers_first():
    path = runner.find_tool("python3", "sh")
    assert path is not None and path.endswith("python3")
