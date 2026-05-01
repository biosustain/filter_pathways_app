# Filter pathways app

## Development environment

Install package so that new code is picked up in a restared python interpreter:

```bash
pip install -e ".[dev]"
```

Then you want to run locally the commands which are check in the CI/CD pipeline 
(using GitHub Actions). You can type

```bash
# run unittests
pytest tests
# format code and sort imports
black .
isort .
# lint code and check for obvious errors
ruff check src
```

## Basic usage

> ToDo.

```python
import filter_pathways_app
print (filter_pathways_app.__version__)
```
