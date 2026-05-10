#!/usr/bin/env python3
"""
onepot_query.py - Playwright skeleton to submit SMILES queries to onepot.ai
=============================================================================

WHAT THIS DOES
--------------
For each (id, SMILES) row in an input CSV, drive a real Chromium browser
(via Playwright) to https://www.onepot.ai/, paste the SMILES into the
search interface, choose a search depth (quick / balanced / comprehensive),
submit, wait for the asynchronous results to render, scrape a few summary
fields out of the DOM, and append a row to an output CSV. Designed so a
crash or Ctrl-C does not lose progress (per-row flush + fsync, --resume).

WHY A REAL BROWSER (not requests/httpx)
---------------------------------------
onepot.ai is a Next.js SPA whose JS bundle performs end-to-end encryption
(RSA-OAEP + AES-256-GCM) on the SMILES *in the user's browser* before
talking to the backend. There is no documented HTTP API, no /api route,
no /docs, no /robots.txt, no /terms. Reverse-engineering the crypto is
out of scope, so we drive the real client.

SETUP (run on the user's laptop, not on a headless HPC node)
------------------------------------------------------------
    # in the project conda env (e.g. <conda_env_prefix>)
    pip install playwright
    playwright install chromium

EXAMPLE INVOCATIONS
-------------------
    # Quick smoke test on the bundled 4-5 SMILES test input, watch it run:
    python scripts/team/onepot_query.py \
        --input  scripts/team/onepot_test_input.csv \
        --output data/onepot_results_test.csv \
        --depth  quick \
        --rate-limit-s 12

    # Full batch on the candidate pool (resume-safe; can be re-run):
    python scripts/team/onepot_query.py \
        --input  data/full_pool_input.csv \
        --output data/onepot_results.csv \
        --depth  balanced

    # Background / headless once selectors are confirmed:
    python scripts/team/onepot_query.py --headless --depth balanced

ETHICAL NOTE - PLEASE READ BEFORE RUNNING AT SCALE
--------------------------------------------------
onepot.ai went out of its way to encrypt user SMILES end-to-end so that
even *their own* server cannot see the molecules a user is querying. That
is a strong product signal that they consider user molecules sensitive
and that they care about respectful use of their service. There is no
public terms-of-service page, but absence of a ToS is not consent to
heavy scraping; if anything it suggests the product is pre-launch.
Before running a 60+ SMILES batch, you should:
  (1) Email the onepot.ai team, describe the chordoma/TBXT hackathon use
      case, and ask whether they offer a CSV bulk-evaluation endpoint
      or an API key for academic / non-commercial work.
  (2) If you proceed without hearing back, keep --rate-limit-s >= 10s,
      run during off-peak hours, and stop immediately if the site shows
      any rate-limit / abuse warning.
  (3) Do not redistribute scraped route data outside the hackathon team.

USAGE WITHIN THE TBXT PIPELINE
------------------------------
Output CSV columns: id, smiles, status, n_routes, top_route_score,
                    raw_html_snippet, timestamp, error
Downstream consumers (e.g. scripts/team/onepot_filter.py,
scripts/team/filter_onepot_candidates.py) should treat status=='ok' rows
with n_routes >= 1 as "synthetically reachable via onepot's 7-reaction
toolkit" and use top_route_score for tie-breaking.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

# Playwright is imported lazily inside main() so that --help works even
# when Playwright is not yet installed in the user's env.


# ---------------------------------------------------------------------------
# SELECTORS - TODO(user): confirm these against the live site before running
# ---------------------------------------------------------------------------
# How to confirm:
#   1. Open https://www.onepot.ai/ in Chrome (logged out is fine).
#   2. Open DevTools (F12) -> Elements panel.
#   3. Click each UI element below, right-click the highlighted node in
#      the Elements panel, choose "Copy" -> "Copy selector", and paste
#      the result here, replacing the placeholder string.
#   4. Prefer stable attributes (id, data-testid, aria-label) over deeply
#      nested CSS paths, which break on the next Next.js rebuild.

# Verified against live https://www.onepot.ai/ on 2026-05-09.
# The text input that holds the SMILES string.
SMILES_INPUT_SELECTOR = 'input[placeholder*="SMILES" i]'

# Submit button. Site has no Quick/Balanced/Comprehensive selector on the
# homepage - depth defaults to "Balanced" (returns up to 50 results).
SUBMIT_BUTTON_SELECTOR = 'button:has-text("Start Synthesis")'

# Depth selector - NOT exposed on the homepage; site default is "balanced"
# (50 results). The --depth flag is currently a no-op; left in place for
# forward-compat if the site adds an explicit selector later.
DEPTH_BUTTON_SELECTOR_TEMPLATE = ''  # not used

# After clicking submit the URL stays at "/", but the page DOM updates
# in-place. We detect results via the body text "molecules" / "Results found".
RESULTS_READY_TEXT_RE = r"\d+\s+molecules|Results found"

# Empty-state text (when nothing matches at the chosen depth).
# Use a word boundary so "0 molecules" doesn't match the "0" in "50 molecules".
NO_RESULTS_TEXT_RE = r"\bNo results\b|\b0\s+molecules\b"

# Validation-error text the site shows when a SMILES fails its sanity check.
INVALID_SMILES_TEXT_RE = r"Please enter a valid SMILES"


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------
OUTPUT_FIELDS = [
    "id",
    "smiles",
    "status",            # ok | no_routes | error | skipped_placeholder
    "n_routes",
    "top_route_score",
    "raw_html_snippet",  # first 500 chars of results container, for debugging
    "timestamp",
    "error",
]


def read_input(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or "id" not in reader.fieldnames or "smiles" not in reader.fieldnames:
            raise SystemExit(
                f"Input {path} must have 'id' and 'smiles' columns; got {reader.fieldnames}"
            )
        for r in reader:
            rid = (r.get("id") or "").strip()
            smi = (r.get("smiles") or "").strip()
            if not rid:
                continue
            rows.append({"id": rid, "smiles": smi})
    return rows


def load_done_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    done: set[str] = set()
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rid = (r.get("id") or "").strip()
            status = (r.get("status") or "").strip()
            # Only treat terminal-success rows as "done"; allow retry of errors.
            if rid and status in {"ok", "no_routes", "skipped_placeholder"}:
                done.add(rid)
    return done


def open_output(path: Path) -> tuple:
    """Open the output CSV for append, writing a header if new."""
    path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not path.exists() or path.stat().st_size == 0
    fh = path.open("a", newline="")
    writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS)
    if new_file:
        writer.writeheader()
        fh.flush()
        os.fsync(fh.fileno())
    return fh, writer


def write_row(fh, writer, row: dict) -> None:
    # Defensive: enforce schema so ad-hoc fields don't sneak in.
    safe = {k: row.get(k, "") for k in OUTPUT_FIELDS}
    writer.writerow(safe)
    fh.flush()
    os.fsync(fh.fileno())


# ---------------------------------------------------------------------------
# Playwright interactions
# ---------------------------------------------------------------------------
def now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def is_placeholder_smiles(smi: str) -> bool:
    if not smi:
        return True
    s = smi.strip().upper()
    return s.startswith("REPLACE_") or s in {"TODO", "PLACEHOLDER", "FIXME"}


def submit_one(page, smiles: str, depth: str, max_wait_s: int) -> dict:
    """
    Drive one SMILES through the UI. Returns a partial output-row dict
    (status, n_routes, top_route_score, raw_html_snippet, error).
    Caller fills in id/smiles/timestamp.
    """
    import re

    # 1. Navigate fresh each time so we don't carry state between queries.
    page.goto("https://www.onepot.ai/", wait_until="domcontentloaded")

    # 2. Wait for the SPA to be interactive - the SMILES input must exist.
    page.wait_for_selector(SMILES_INPUT_SELECTOR, timeout=20_000)

    # 3. React-friendly fill sequence (verified). Direct .fill() leaves the
    #    input visually filled but React's controlled-state stays empty if an
    #    autocomplete dropdown is open. Sequence: focus -> Escape any popup
    #    -> refocus -> fill -> Tab to blur and trigger debounced validator.
    inp = page.locator(SMILES_INPUT_SELECTOR)
    inp.click()
    page.wait_for_timeout(800)
    page.keyboard.press("Escape")
    page.wait_for_timeout(300)
    inp.click()
    page.wait_for_timeout(300)
    inp.fill(smiles)
    page.keyboard.press("Tab")
    page.wait_for_timeout(2500)  # debounce + async validator

    # 4. Check the validator didn't reject the SMILES.
    body_text = page.evaluate("() => document.body.innerText")
    if re.search(INVALID_SMILES_TEXT_RE, body_text, re.I):
        return {
            "status": "error",
            "n_routes": 0,
            "top_route_score": "",
            "raw_html_snippet": body_text[:500].replace("\n", " "),
            "error": "site rejected SMILES as invalid",
        }

    # 5. Submit the form.
    page.locator(SUBMIT_BUTTON_SELECTOR).click(timeout=10_000)

    # 6. Wait for results - page renders in-place (no URL change).
    deadline = time.time() + max_wait_s
    found = False
    no_results = False
    while time.time() < deadline:
        try:
            txt = page.evaluate("() => document.body.innerText")
        except Exception:
            time.sleep(0.5)
            continue
        if re.search(NO_RESULTS_TEXT_RE, txt, re.I):
            no_results = True
            break
        if re.search(RESULTS_READY_TEXT_RE, txt, re.I):
            found = True
            break
        time.sleep(1.0)

    raw_snippet = ""
    try:
        raw_snippet = (page.evaluate("() => document.body.innerText") or "")[:500].replace("\n", " ")
    except Exception:
        pass

    if no_results:
        return {
            "status": "no_routes",
            "n_routes": 0,
            "top_route_score": "",
            "raw_html_snippet": raw_snippet,
            "error": "",
        }

    if not found:
        return {
            "status": "error",
            "n_routes": 0,
            "top_route_score": "",
            "raw_html_snippet": raw_snippet,
            "error": f"timeout waiting {max_wait_s}s for results",
        }

    # 7. Parse the molecule count from body text. Site reports e.g. "50 molecules".
    n_routes = 0
    m = re.search(r"(\d+)\s+molecules", raw_snippet, re.I)
    if m:
        n_routes = int(m.group(1))
    # Best-effort grab of the first similarity % shown (proxy for top hit quality).
    top_score = ""
    m2 = re.search(r"(\d{1,3}%)\s*similarity", raw_snippet, re.I) \
         or re.search(r"(\d{1,3})%", raw_snippet)
    if m2:
        top_score = m2.group(1)

    return {
        "status": "ok",
        "n_routes": n_routes,
        "top_route_score": top_score,
        "raw_html_snippet": raw_snippet,
        "error": "",
    }


def save_debug_screenshot(page, debug_dir: Path, rid: str) -> str:
    debug_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = debug_dir / f"{rid}_{ts}.png"
    try:
        page.screenshot(path=str(out), full_page=True)
        return str(out)
    except Exception as exc:
        return f"<screenshot failed: {exc}>"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]  # .../TBXT
    p = argparse.ArgumentParser(
        description="Submit SMILES to onepot.ai via Playwright and write a CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input",
        type=Path,
        default=repo_root / "scripts" / "team" / "onepot_test_input.csv",
        help="Input CSV with columns: id, smiles",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=repo_root / "data" / "onepot_results.csv",
        help="Output CSV (appended to; safe to re-run with --resume).",
    )
    p.add_argument(
        "--depth",
        choices=["quick", "balanced", "comprehensive"],
        default="balanced",
        help="onepot search depth (quick=10, balanced=50, comprehensive=200 max items).",
    )
    p.add_argument(
        "--headless",
        action="store_true",
        default=False,
        help="Run Chromium headless (default: visible so you can watch).",
    )
    p.add_argument(
        "--rate-limit-s",
        type=float,
        default=10.0,
        help="Sleep between submissions (be polite - site uses E2E crypto).",
    )
    p.add_argument(
        "--max-wait-s",
        type=int,
        default=90,
        help="Max seconds to wait for results to render per query.",
    )
    p.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip ids already present (with terminal status) in --output.",
    )
    p.add_argument(
        "--debug-dir",
        type=Path,
        default=repo_root / "data" / "onepot_debug",
        help="Where to dump failure screenshots.",
    )
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)

    if not args.input.exists():
        print(f"[onepot_query] ERROR: input not found: {args.input}", file=sys.stderr)
        return 2

    rows = read_input(args.input)
    done = load_done_ids(args.output) if args.resume else set()
    todo = [r for r in rows if r["id"] not in done]

    print(
        f"[onepot_query] input={args.input}  total={len(rows)}  "
        f"already_done={len(done)}  to_run={len(todo)}  depth={args.depth}  "
        f"headless={args.headless}  rate_limit_s={args.rate_limit_s}",
        flush=True,
    )
    if not todo:
        print("[onepot_query] nothing to do.")
        return 0

    # Lazy import so --help works without playwright installed.
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print(
            "[onepot_query] ERROR: playwright not installed.\n"
            "  pip install playwright\n"
            "  playwright install chromium",
            file=sys.stderr,
        )
        return 3

    # On Ubuntu 26.04 (e.g. WSL2 Resolute Raccoon) Playwright doesn't ship
    # an official chromium build. The platform-override pulls the
    # ubuntu24.04-x64 binary as a fallback. Harmless on supported distros.
    os.environ.setdefault("PLAYWRIGHT_HOST_PLATFORM_OVERRIDE", "ubuntu24.04-x64")

    # The site does crypto + WebGL work that hangs the headless-shell variant.
    # Force the full chromium binary if available; safe to omit otherwise.
    full_chrome = Path.home() / ".cache/ms-playwright/chromium-1217/chrome-linux64/chrome"
    launch_kwargs = {"headless": args.headless, "args": ["--no-sandbox"]}
    if full_chrome.exists():
        launch_kwargs["executable_path"] = str(full_chrome)

    out_fh, writer = open_output(args.output)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(**launch_kwargs)
            context = browser.new_context(
                viewport={"width": 1400, "height": 900},
                # Use a standard desktop Chrome UA so the SPA renders normally.
                user_agent=(
                    "Mozilla/5.0 (X11; Linux x86_64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/131.0.0.0 Safari/537.36"
                ),
            )
            page = context.new_page()

            for i, row in enumerate(todo, 1):
                rid = row["id"]
                smi = row["smiles"]
                base = {
                    "id": rid,
                    "smiles": smi,
                    "status": "",
                    "n_routes": 0,
                    "top_route_score": "",
                    "raw_html_snippet": "",
                    "timestamp": now_iso(),
                    "error": "",
                }

                if is_placeholder_smiles(smi):
                    base.update(
                        status="skipped_placeholder",
                        error="smiles field is a placeholder; fill it in",
                    )
                    write_row(out_fh, writer, base)
                    print(f"[{i}/{len(todo)}] {rid} SKIPPED (placeholder smiles)")
                    continue

                print(f"[{i}/{len(todo)}] {rid} -> submitting ({smi[:60]}...)", flush=True)
                try:
                    result = submit_one(page, smi, args.depth, args.max_wait_s)
                    base.update(result)
                    if result["status"] == "error":
                        shot = save_debug_screenshot(page, args.debug_dir, rid)
                        base["error"] = f"{result['error']}; screenshot={shot}"
                except Exception as exc:
                    shot = save_debug_screenshot(page, args.debug_dir, rid)
                    base.update(
                        status="error",
                        error=f"{type(exc).__name__}: {exc}; screenshot={shot}",
                    )
                    traceback.print_exc()

                base["timestamp"] = now_iso()
                write_row(out_fh, writer, base)
                print(
                    f"    status={base['status']} n_routes={base['n_routes']} "
                    f"top_score={base['top_route_score']!r}",
                    flush=True,
                )

                # Polite pacing - don't hammer the SPA.
                if i < len(todo):
                    time.sleep(args.rate_limit_s)

            context.close()
            browser.close()
    finally:
        out_fh.close()

    print(f"[onepot_query] done. wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
