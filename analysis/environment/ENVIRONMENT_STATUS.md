# Environment status

`requirements.in` records the direct Python dependencies. `requirements.lock`
is the resolved public Python 3.14 lock used for package smoke testing, while
`requirements-current-host.lock` and `RECORDED_RUNTIME_METADATA.json` preserve
the historical execution environment. The final Article 002 workflow is
Python-based; the retained R package files are provenance for legacy staged
utilities and are not required by `run_pipeline.py`.
