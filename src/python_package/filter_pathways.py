"""Core logic for querying UniProt via acore and filtering annotations by keywords.

Typical usage::

    from python_package.filter_pathways import query_uniprot, filter_annotations

    df = query_uniprot(["P05067", "P12345"], fields="accession,go_p,go_c")
    filtered = filter_annotations(df, keywords=["apoptosis", "mitochondr"])
"""

from __future__ import annotations

import re as _re

import pandas as pd

from acore.io.uniprot import fetch_annotations, process_annotations

re_IGNORECASE = _re.IGNORECASE


def query_uniprot(
    uniprot_ids: list[str] | pd.Index,
    fields: str = "accession,go_p,go_c,go_f",
) -> pd.DataFrame:
    """Fetch annotations from UniProt for the given protein IDs and return a
    tidy long-format DataFrame.

    Parameters
    ----------
    uniprot_ids : list[str] | pd.Index
        UniProt accession IDs (e.g. ``["P05067", "P12345"]``).
    fields : str, optional
        Comma-separated UniProt return-field identifiers.
        ``accession`` is always included by the API.
        Default includes GO biological process, cellular component and
        molecular function.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns:
        ``identifier``, ``source``, ``annotation``.

    Examples
    --------
    >>> # Integration test – skipped in offline environments
    >>> query_uniprot(["P05067"]) # doctest: +SKIP
    """
    raw_df = fetch_annotations(uniprot_ids, fields=fields)
    long_df = process_annotations(raw_df, fields=fields)
    return long_df


def filter_annotations(
    annotations: pd.DataFrame,
    keywords: list[str],
    case_sensitive: bool = False,
) -> pd.DataFrame:
    """Filter a long-format annotations DataFrame to rows whose ``annotation``
    column contains **any** of the given keywords.

    Parameters
    ----------
    annotations : pd.DataFrame
        Long-format DataFrame with at least an ``annotation`` column (as
        returned by :func:`query_uniprot`).
    keywords : list[str]
        Keywords to search for.  A row is kept when at least one keyword
        appears in the annotation string.  Empty strings are ignored.
    case_sensitive : bool, optional
        Whether the search is case-sensitive.  Default ``False``.

    Returns
    -------
    pd.DataFrame
        Filtered copy of *annotations* where at least one keyword matched.
        Returns an empty DataFrame with the same columns when no keywords
        are provided or none match.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     "identifier": ["P1", "P1", "P2"],
    ...     "source": ["go_p", "go_p", "go_p"],
    ...     "annotation": [
    ...         "cell apoptosis [GO:0042981]",
    ...         "cell division",
    ...         "mitochondria",
    ...     ],
    ... })
    >>> result = filter_annotations(df, keywords=["apoptosis", "mitochondr"])
    >>> list(result["annotation"])
    ['cell apoptosis [GO:0042981]', 'mitochondria']
    >>> filter_annotations(df, keywords=[]).shape[0]
    0
    """
    # Drop empty / whitespace-only keywords
    clean_keywords = [kw.strip() for kw in keywords if kw.strip()]
    if not clean_keywords:
        return annotations.iloc[0:0].copy()

    flags = 0 if case_sensitive else re_IGNORECASE
    pattern = "|".join(
        _escape(kw) for kw in clean_keywords
    )
    mask = annotations["annotation"].str.contains(
        pattern, flags=flags, na=False
    )
    return annotations[mask].copy()


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _escape(keyword: str) -> str:
    return _re.escape(keyword)


def export_to_csv(df: pd.DataFrame) -> bytes:
    """Serialise a DataFrame to CSV bytes (UTF-8 with BOM for Excel compat).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to serialise.

    Returns
    -------
    bytes
        UTF-8-with-BOM encoded CSV bytes.

    Examples
    --------
    >>> import pandas as pd
    >>> data = export_to_csv(pd.DataFrame({"a": [1, 2]}))
    >>> data.startswith(b'\\xef\\xbb\\xbf')
    True
    """
    return df.to_csv(index=False).encode("utf-8-sig")
