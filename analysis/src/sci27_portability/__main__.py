from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import expand_environment_tokens


def main() -> int:
    parser = argparse.ArgumentParser(description="Render a SCI27 JSON configuration template.")
    parser.add_argument("template")
    parser.add_argument("output")
    args = parser.parse_args()
    source = Path(args.template)
    target = Path(args.output)
    payload = expand_environment_tokens(json.loads(source.read_text(encoding="utf-8")))
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
