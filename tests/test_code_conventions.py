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


def test_readme_drug_list_size_matches_code():
    """The README states how many drugs the built-in list holds."""
    import re
    import runner  # noqa: F401  (sets up the package import path)
    from vspipeline.library import LibraryManager
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as fh:
        stated = re.search(r"fixed list of (\d+) well-known", fh.read())
    assert stated and int(stated.group(1)) == len(LibraryManager.FDA_DRUG_NAMES)


def test_report_ngl_version_matches_the_image_build():
    """report.py downloads/verifies the same NGL version the image installs."""
    import json
    import runner  # noqa: F401
    from vspipeline import report
    with open(os.path.join(ROOT, "tools", "ngl", "package.json"), encoding="utf-8") as fh:
        pinned = json.load(fh)["dependencies"]["ngl"]
    assert pinned == report.NGL_VERSION
