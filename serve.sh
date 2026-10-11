#!/bin/sh
# Serve this folder on http://localhost:8000 (the viewers fetch data files, so file:// does not work).
cd "$(dirname "$0")" && exec python3 -m http.server "${1:-8000}"
