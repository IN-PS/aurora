"""Permite executar ``python -m aurora``."""

import sys

from .interface.cli import main

sys.exit(main())
