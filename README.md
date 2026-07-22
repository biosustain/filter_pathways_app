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

```python
import filter_pathways_app
print (filter_pathways_app.__version__)
```

You can also resolve KEGG KO terms to their common descriptions using acore

```python
from acore.io.kegg import fetch_kegg_ko_descriptions
df = fetch_kegg_ko_descriptions(["ko:K03007", "K02143"])
print(df[["ko_term", "common_description"]])
```
