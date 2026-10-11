"""Code conventions that keep runs reproducible across machines."""
import ast
import glob
import os

ROOT = os.path.join(os.path.dirname(__file__), "..")
SOURCES = sorted(glob.glob(os.path.join(ROOT, "autodocker", "**", "*.py"), recursive=True)
                 + glob.glob(os.path.join(ROOT, "scripts", "*.py")))


def _is_binary_mode(call):
    mode = call.args[1] if len(call.args) > 1 else next(
        (k.value for k in call.keywords if k.arg == "mode"), None)
    return isinstance(mode, ast.Constant) and "b" in str(mode.value)


def test_every_text_open_names_its_encoding():
    """open() without encoding= uses the locale's encoding, so the same
    input could be read differently on another system; every text-mode
    open() in the pipeline and scripts states UTF-8 explicitly."""
    missing = []
    for path in SOURCES:
        with open(path, encoding="utf-8") as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "open" and not _is_binary_mode(node)
                    and not any(k.arg == "encoding" for k in node.keywords)):
                missing.append(f"{os.path.relpath(path, ROOT)}:{node.lineno}")
    assert not missing, "open() without encoding=: " + ", ".join(missing)
