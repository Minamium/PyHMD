#!/usr/bin/env python3
"""Check local Markdown targets and English/Japanese navigation; no network I/O."""
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit


def local_targets(page, root):
    """Read inline Markdown links outside fenced code blocks."""
    fence = None
    for number, line in enumerate(page.read_text(encoding="utf-8").splitlines(), 1):
        marker = re.match(r"^\s*(`{3,}|~{3,})(.*)$", line)
        if marker:
            token, suffix = marker.groups()
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence) and not suffix.strip():
                fence = None
            continue
        if fence is not None:
            continue
        for match in re.finditer(r"\[[^\]]*\]\(([^)]+)\)", line):
            destination = match.group(1).strip()
            if not destination:
                continue
            if destination.startswith("<"):
                destination = destination[1:destination.index(">")]
            else:
                destination = destination.split()[0]
            url = urlsplit(destination)
            if url.scheme or url.netloc or not url.path:
                continue
            path = unquote(url.path)
            target = root / path.lstrip("/") if path.startswith("/") else page.parent / path
            yield number, destination, target.resolve()


def check(root):
    root = Path(root).resolve()
    pages = sorted({*root.glob("README*.md"), root / "AGENTS.md",
                    *(root / "docs").rglob("*.md"), *(root / "tests").glob("*.md"),
                    *(root / ".agents/skills").rglob("SKILL.md")})
    errors, targets = [], {}
    for page in pages:
        if not page.is_file():
            errors.append(f"missing documentation: {page.relative_to(root)}")
            continue
        targets[page] = list(local_targets(page, root))
        for number, label, target in targets[page]:
            if not target.exists():
                errors.append(f"{page.relative_to(root)}:{number}: missing target {label}")
    pairs = [(root / "README.md", root / "README.ja.md")]
    pairs += [(page, root / "docs/ja" / page.name) for page in sorted((root / "docs").glob("*.md"))]
    for english, japanese in pairs:
        for page, counterpart in ((english, japanese), (japanese, english)):
            if not page.is_file():
                errors.append(f"missing language page: {page.relative_to(root)}")
            elif not any(target == counterpart.resolve() for _, _, target in targets[page]):
                errors.append(f"{page.relative_to(root)}: missing link to {counterpart.relative_to(root)}")
    for page in (root / "docs/ja").glob("*.md"):
        if not (root / "docs" / page.name).is_file():
            errors.append(f"missing English counterpart: {page.relative_to(root)}")
    return errors, len(pages), len(pairs), sum(len(items) for items in targets.values())


def main():
    errors, pages, pairs, links = check(Path(__file__).resolve().parents[1])
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"PASS {pages} Markdown files, {pairs} language pairs, {links} local links")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
