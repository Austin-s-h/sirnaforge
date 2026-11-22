"""Isolated miRNA download worker, executed as a script to bound the entire request.

Keep this module standard-library-only: importing the data package before applying
an acquisition timeout would load unrelated annotation and alignment dependencies.
"""

import gzip
import html
import json
import sys
import urllib.request


def main() -> int:
    """Write a downloaded mature reference to stdout, or a diagnostic to stderr."""
    try:
        options = json.load(sys.stdin)
        request = urllib.request.Request(
            options["url"],
            headers={
                "User-Agent": "sirnaforge/1.0 (+https://github.com/austin-s-h/sirnaforge)",
                "Accept": "text/plain,application/octet-stream",
            },
        )
        # URLs come from the same configured ReferenceSource boundary as the base downloader.
        with urllib.request.urlopen(request, timeout=options["timeout"]) as response:  # nosec B310
            data = response.read()
        if options["compressed"] and options["url"].endswith(".gz"):
            data = gzip.decompress(data)
        try:
            content = data.decode("utf-8")
        except UnicodeDecodeError:
            content = data.decode("latin-1")
        content = html.unescape(content).replace("<br>", "\n").replace("<BR>", "\n")
        if not content.strip():
            raise ValueError("empty response")
        sys.stdout.write(content)
        return 0
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
