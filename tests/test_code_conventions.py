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
    from vspipeline.library import LibraryManager
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as fh:
        stated = re.search(r"fixed list of (\d+) well-known", fh.read())
    assert stated and int(stated.group(1)) == len(LibraryManager.FDA_DRUG_NAMES)


def test_report_ngl_version_matches_the_image_build():
    """report.py downloads/verifies the same NGL version the image installs."""
    import json
    from vspipeline import report
    with open(os.path.join(ROOT, "tools", "ngl", "package.json"), encoding="utf-8") as fh:
        pinned = json.load(fh)["dependencies"]["ngl"]
    assert pinned == report.NGL_VERSION


def test_readme_option_defaults_match_the_cli():
    """Each numeric default in README's option table equals the CLI default."""
    import re
    import runner
    from vspipeline import defaults as builtin
    defaults = {}
    for action in runner.build_parser()._actions:
        for opt in action.option_strings:
            defaults[opt] = action.default
    # --timeout follows $VS_COMMAND_TIMEOUT when set; the README states the
    # built-in value.
    defaults["--timeout"] = builtin.COMMAND_TIMEOUT
    with open(os.path.join(ROOT, "README.md"), encoding="utf-8") as fh:
        rows = [ln for ln in fh if ln.startswith("| `-")]
    checked = 0
    for row in rows:
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        stated = [float(n) for n in re.findall(r"-?\d+(?:\.\d+)?", cells[-1])]
        for opt in re.findall(r"`(--?[\w-]+)", cells[0]):
            value = defaults.get(opt)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                assert float(value) in stated, f"README default of {opt}: {cells[-1]} != {value}"
                checked += 1
    assert checked >= 10  # the table was parsed


def test_benchmark_conformers_use_the_pipeline_recipe():
    """The re-docking start conformers come from vspipeline.conformers."""
    with open(os.path.join(ROOT, "scripts", "_benchlib.py"), encoding="utf-8") as fh:
        src = fh.read()
    assert "conformers.embed_and_relax(" in src
    assert "AllChem.ETKDGv3()" not in src


def test_redock_sample_examples_match_the_shared_constants():
    """fetch_casf_rcsb.py and benchmark_redock.py must select the same
    complexes; their usage examples state the shared sample size and seed."""
    import re
    import sys
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import _benchlib
    for name in ("fetch_casf_rcsb.py", "benchmark_redock.py"):
        with open(os.path.join(ROOT, "scripts", name), encoding="utf-8") as fh:
            src = fh.read()
        assert "default=_benchlib.BENCH_SAMPLE_SIZE" in src, name
        for size, seed in re.findall(r"--sample (\d+)\s+--sample-seed (\d+)", src):
            assert int(size) == _benchlib.BENCH_SAMPLE_SIZE, name
            assert int(seed) == _benchlib.BENCH_SEED, name
