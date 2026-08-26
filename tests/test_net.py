import vspipeline.net as net


class FakeResp:
    def __init__(self, content):
        self.content = content

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
