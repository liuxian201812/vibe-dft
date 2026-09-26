# Smoke test

The smoke test checks the public, local workflow without connecting to a
cluster or running a VASP executable.

It verifies:

1. POSCAR validation;
2. relax input generation;
3. band input generation;
4. independent submit-script generation.

## Run the privacy scan

```bash
bash tools/check_no_site_details.sh
```

The scanner uses generic rules rather than a repository-specific blacklist. It
checks tracked and non-ignored UTF-8 text for contact addresses, network
addresses, user-home or site mount paths, credential formats, private-key
material, embedded URL credentials, and literal scheduler account or partition
values. Descriptive angle-bracket placeholders remain allowed.

For safer CI logs, failures identify only the file and rule category; the
matched value is not printed. The scanner checks the current working tree. It
does not remove sensitive values from earlier Git commits or from external
forks, caches, workflow logs, or clones.

## Run the workflow smoke test

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
PYTHON=.venv/bin/python bash smoke_tests/run_smoke_test.sh
```

The test creates synthetic inputs under `smoke_tests/tmp/`, exports a simulated
VASP command, and uses a fake POTCAR library. It does not consume cluster time.

A successful run ends with:

```text
Smoke test passed.
```

GitHub Actions runs the privacy scan and workflow smoke test on pushes and pull
requests.

## Clean up

```bash
rm -rf smoke_tests/tmp
```

`smoke_tests/tmp/` is ignored by Git and must not be used to store permanent or
sensitive runtime configuration.
