#!/usr/bin/env python3
"""Render a Marp-style markdown deck to PDF via markdown -> HTML -> Chromium.

This bypasses Marp CLI (which requires Node) by:
  1. Parsing the YAML frontmatter to extract inline `style:` CSS.
  2. Splitting body on `---` slide separators.
  3. Detecting per-slide HTML comments like `<!-- _class: lead -->`.
  4. Converting each slide's markdown to HTML (including raw HTML pass-through).
  5. Wrapping in a single document with one @page-sized section per slide.
  6. Printing to PDF via Playwright/Chromium.

Usage:
  python tools/render_slides.py slides/slides.md
  # produces slides/slides.pdf

Dependencies:
  pip install markdown playwright
  playwright install chromium
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import markdown
from playwright.sync_api import sync_playwright


if len(sys.argv) < 2:
    sys.exit("usage: render_slides.py <path/to/deck.md>")
SRC = Path(sys.argv[1]).resolve()
OUT_PDF = SRC.with_suffix(".pdf")
OUT_HTML = Path("/tmp") / (SRC.stem + ".html")

# 16:9 slide page (Marp default 1280x720 -> use that)
PAGE_W = "1280px"
PAGE_H = "720px"


def parse_frontmatter(text: str):
    """Extract YAML frontmatter (between leading --- pair) and body.

    Strips any leading <!-- ... --> comment block before checking for ---.
    """
    # Skip a leading HTML comment block, if any
    text = re.sub(r"\A\s*<!--.*?-->\s*", "", text, flags=re.S)
    if not text.startswith("---"):
        return {}, text
    end = text.find("\n---", 3)
    if end == -1:
        return {}, text
    fm = text[3:end].strip("\n")
    body = text[end + 4 :].lstrip("\n")
    # Naive YAML extraction of `style: |` block + simple top-level keys
    style_match = re.search(r"^style:\s*\|\s*\n(.+?)(?=^\S|\Z)", fm, re.M | re.S)
    style_css = ""
    if style_match:
        raw = style_match.group(1)
        # Dedent by smallest leading-space count of non-empty lines
        lines = [ln for ln in raw.splitlines() if ln.strip()]
        indent = min((len(ln) - len(ln.lstrip(" "))) for ln in lines) if lines else 0
        style_css = "\n".join(ln[indent:] for ln in raw.splitlines())
    return {"style": style_css}, body


def split_slides(body: str):
    """Split on lines that are exactly '---' (slide separator)."""
    parts = re.split(r"^---\s*$", body, flags=re.M)
    return [p.strip("\n") for p in parts if p.strip()]


def extract_class(slide_md: str):
    """Pull `<!-- _class: foo -->` directive (Marp per-slide class)."""
    m = re.search(r"<!--\s*_class:\s*([\w\- ]+)\s*-->", slide_md)
    if not m:
        return "", slide_md
    cls = m.group(1).strip()
    cleaned = re.sub(r"<!--\s*_class:\s*[\w\- ]+\s*-->", "", slide_md, count=1).strip("\n")
    return cls, cleaned


def strip_html_comments(s: str) -> str:
    return re.sub(r"<!--.*?-->", "", s, flags=re.S)


_MARP_IMG_RE = re.compile(r"!\[([whWH]):(\d+)\]\(([^)]+)\)")


def _marp_img_sub(m: re.Match) -> str:
    dim = m.group(1).lower()
    val = m.group(2)
    src = m.group(3)
    style = f"width:{val}px;" if dim == "w" else f"height:{val}px;"
    return f'<img src="{src}" style="{style}">'


def md_to_html(md_text: str) -> str:
    # Convert Marp-style image-width directives to plain <img> tags first
    md_text = _MARP_IMG_RE.sub(_marp_img_sub, md_text)
    md = markdown.Markdown(
        extensions=["tables", "fenced_code", "md_in_html", "attr_list"],
    )
    return md.convert(md_text)


def render_html(frontmatter_css: str, slides: list[tuple[str, str]]) -> str:
    base_css = f"""
@page {{ size: {PAGE_W} {PAGE_H}; margin: 0; }}
html, body {{ margin: 0; padding: 0; background: #fff; }}
section {{
    width: {PAGE_W};
    height: {PAGE_H};
    box-sizing: border-box;
    page-break-after: always;
    break-after: page;
    overflow: hidden;
    position: relative;
}}
section:last-child {{ page-break-after: auto; }}
img {{ max-width: 100%; height: auto; }}
"""
    base_href = SRC.parent.as_uri() + "/"
    parts = [
        "<!doctype html>",
        f'<html><head><meta charset="utf-8"><base href="{base_href}"><style>',
        base_css,
        frontmatter_css,
        "</style></head><body>",
    ]
    for cls, html in slides:
        cls_attr = f' class="{cls}"' if cls else ""
        parts.append(f"<section{cls_attr}>{html}</section>")
    parts.append("</body></html>")
    return "\n".join(parts)


def main() -> int:
    text = SRC.read_text()
    fm, body = parse_frontmatter(text)

    slides_md = split_slides(body)
    rendered = []
    for s in slides_md:
        s = strip_html_comments(s) if not re.search(r"_class:", s) else s
        cls, cleaned = extract_class(s)
        cleaned = strip_html_comments(cleaned)
        # markdown-in-html: wrap raw <div> blocks so inner markdown still parses
        cleaned = re.sub(r"(<div[^>]*>)", r"\1\n\n", cleaned)
        cleaned = re.sub(r"(</div>)", r"\n\n\1", cleaned)
        html = md_to_html(cleaned)
        rendered.append((cls, html))

    full_html = render_html(fm.get("style", ""), rendered)
    OUT_HTML.write_text(full_html)
    print(f"[html] wrote {OUT_HTML} ({len(rendered)} slides)")

    # Render to PDF via Chromium
    with sync_playwright() as p:
        browser = p.chromium.launch()
        ctx = browser.new_context(viewport={"width": 1280, "height": 720})
        page = ctx.new_page()
        page.goto(OUT_HTML.as_uri(), wait_until="networkidle")
        page.emulate_media(media="print")
        page.pdf(
            path=str(OUT_PDF),
            width=PAGE_W,
            height=PAGE_H,
            margin={"top": "0", "right": "0", "bottom": "0", "left": "0"},
            print_background=True,
            prefer_css_page_size=True,
        )
        browser.close()

    print(f"[pdf]  wrote {OUT_PDF} ({OUT_PDF.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
