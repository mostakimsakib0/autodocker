import vspipeline.net as net


class FakeResp:
    def __init__(self, content):
        self.content = content
        self.status_code = 200

    def raise_for_status(self):
        pass


class FakeRequests:
    def __init__(self, content):
        self.content = content

    def get(self, url, timeout=None):
        return FakeResp(self.content)


def test_http_requests(monkeypatch):
    monkeypatch.setattr(net, "requests", FakeRequests(b"data"))
    assert net._http_get_bytes("http://x") == b"data"


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


def test_http_no_lib(monkeypatch):
    monkeypatch.setattr(net, "requests", None)
    monkeypatch.setattr(net, "urllib", None)
    import pytest
    with pytest.raises(RuntimeError):
        net._http_get_bytes("http://x")


def test_http_404_is_not_retried(monkeypatch):
    import pytest
    calls = []

    class Resp:
        status_code = 404

    class NotFound(OSError):
        response = Resp()

    def get_once(url, timeout=30):
        calls.append(url)
        raise NotFound("404 Not Found")

    monkeypatch.setattr(net, "_http_get_once", get_once)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    with pytest.raises(NotFound):
        net._http_get_bytes("http://x")
    assert len(calls) == 1


def test_connection_errors_are_retried(monkeypatch):
    import pytest
    calls = []

    def get_once(url, timeout=30):
        calls.append(url)
        raise ConnectionRefusedError("refused")

    monkeypatch.setattr(net, "_http_get_once", get_once)
    monkeypatch.setattr(net.time, "sleep", lambda s: None)
    with pytest.raises(ConnectionRefusedError):
        net._http_get_bytes("http://x", retries=3)
    assert len(calls) == 3
