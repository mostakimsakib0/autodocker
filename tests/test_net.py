import pytest

import vspipeline.net as net


class FakeResp:
    def __init__(self, content, status_code=200):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        pass


class FakeRequests:
    def __init__(self, content):
        self.content = content
        self.timeouts = []

    def get(self, url, timeout=None):
        self.timeouts.append(timeout)
        return FakeResp(self.content)


def test_http_requests_uses_named_timeout(monkeypatch):
    fake = FakeRequests(b"data")
    monkeypatch.setattr(net, "requests", fake)
    assert net._http_get_bytes("http://x") == b"data"
    assert fake.timeouts == [net.HTTP_TIMEOUT]


def test_http_urllib_fallback(monkeypatch):
    monkeypatch.setattr(net, "requests", None)
    import io

    class Ctx:
        def __enter__(self):
            return io.BytesIO(b"u")

        def __exit__(self, *a):
            return False

    class FakeUrlopen:
        def urlopen(self, url, timeout=None):
            return Ctx()

    class FakeUrllib:
        request = FakeUrlopen()

    monkeypatch.setattr(net, "urllib", FakeUrllib)
    assert net._http_get_bytes("http://x") == b"u"


def test_programming_errors_are_not_retried(monkeypatch):
    calls = []

    def get_once(url, timeout=None):
        calls.append(url)
        raise RuntimeError("bug")

    monkeypatch.setattr(net, "_http_get_once", get_once)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    with pytest.raises(RuntimeError):
        net._http_get_bytes("http://x")
    assert len(calls) == 1


def test_http_404_is_not_retried(monkeypatch):
    calls = []

    class Resp:
        status_code = 404

    class NotFound(OSError):
        response = Resp()

    def get_once(url, timeout=None):
        calls.append(url)
        raise NotFound("404 Not Found")

    monkeypatch.setattr(net, "_http_get_once", get_once)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    with pytest.raises(NotFound):
        net._http_get_bytes("http://x")
    assert len(calls) == 1


def test_connection_errors_are_retried(monkeypatch):
    calls = []

    def get_once(url, timeout=None):
        calls.append(url)
        raise ConnectionRefusedError("refused")

    monkeypatch.setattr(net, "_http_get_once", get_once)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    with pytest.raises(ConnectionRefusedError):
        net._http_get_bytes("http://x", retries=3)
    assert len(calls) == 3


def test_429_is_retried_then_raised(monkeypatch):
    calls = []

    def get_once(url, timeout=None):
        calls.append(url)
        raise net._RetryableHTTPError("HTTP 429")

    monkeypatch.setattr(net, "_http_get_once", get_once)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    with pytest.raises(net._RetryableHTTPError):
        net._http_get_bytes("http://x", retries=2)
    assert len(calls) == 2


def test_is_network_failure():
    class Resp:
        status_code = 404

    class NotFound(OSError):
        response = Resp()

    assert net.is_network_failure(net._RetryableHTTPError("429"))
    assert net.is_network_failure(ConnectionResetError("reset"))
    assert not net.is_network_failure(NotFound("404"))
    assert not net.is_network_failure(ValueError("bad json"))
