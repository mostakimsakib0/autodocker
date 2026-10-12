---
name: Tester report
about: Report from following docs/user_testing.md on your own machine
labels: user-test
---

**Machine**
- Operating system and version:
- CPU cores and RAM:
- Container engine and version (`docker --version` / `podman --version`), rootless or not:
- Your user ID (`id -u`, Linux only):

**How you got the image**: docker pull / built from source / could not get it
**Time to get the image** (minutes):
**Time for the first run** (minutes):

**Results**
- Did the first run finish? yes / no
- Needed `chmod -R a+rwX works`? yes / no
- Box on the inhibitor in `results_report.html`? yes / no / could not tell
- `ranking.csv` hashes of the two runs equal? yes / no
- `verify_run.py` exit status (one run / two runs):

**Ratings** (1 easy – 5 hard): getting the image: , running and reading results:

**Where did you get stuck or have to guess?**

**Error messages** (paste exactly)

**Suggestions**
