"""Unit tests. Offline: no network, no Claude calls, no AMRFinderPlus runs.

Run from the project root:
    python -m unittest discover -s code/tests -t code
"""
import os

# Never start the background feed thread or touch the real database in tests.
os.environ.setdefault("OHP_DISABLE_FEED_WORKER", "1")
