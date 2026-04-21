from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "output" / "tqec_import"
GENERATED_ROOT = OUTPUT_ROOT / "generated"
INJECTED_ROOT = OUTPUT_ROOT / "injected"
LER_ROOT = OUTPUT_ROOT / "ler"

