#!/usr/bin/env python3
"""
NotebookLM integration for single-cell RNA-seq courses.

Uploads course content to Google NotebookLM so you can:
  - Generate audio overviews (podcasts) to listen while commuting
  - Ask questions and get answers from your course materials
  - Generate quizzes and flashcards for self-testing

Requirements:
  pip install "notebooklm-py[browser]"
  playwright install chromium

First-time setup:
  notebooklm login   # Opens browser for Google auth

Usage:
  python scripts/notebooklm_upload.py prerequisites
  python scripts/notebooklm_upload.py 1.single_cell_processing_course
  python scripts/notebooklm_upload.py all
  python scripts/notebooklm_upload.py all --generate audio
"""

import argparse
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "notebooklm_courses.yaml"


def load_config():
    if not yaml:
        print("Install PyYAML: pip install pyyaml")
        sys.exit(1)
    if not CONFIG_PATH.exists():
        print(f"Config not found: {CONFIG_PATH}")
        sys.exit(1)
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def resolve_sources(course_id: str, config: dict) -> list[Path]:
    """Resolve all source file paths for a course."""
    courses = config.get("courses", {})
    if course_id not in courses:
        print(f"Unknown course: {course_id}")
        print(f"Available: {', '.join(courses.keys())}")
        sys.exit(1)

    course = courses[course_id]
    paths = []

    for rel in course.get("sources", []):
        p = REPO_ROOT / rel
        if p.exists():
            paths.append(p)
        else:
            print(f"  (skip missing: {rel})")

    for pattern in course.get("sources_glob", []):
        base = REPO_ROOT / pattern.split("*")[0].rstrip("/")
        if not base.exists():
            continue
        for p in base.rglob("*.md"):
            if p.is_file() and p not in paths:
                paths.append(p)

    return sorted(paths)


def run_notebooklm(args: list[str]) -> bool:
    """Run notebooklm CLI command."""
    try:
        result = subprocess.run(
            ["notebooklm"] + args,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(result.stderr or result.stdout)
            return False
        if result.stdout:
            print(result.stdout)
        return True
    except FileNotFoundError:
        print(
            "notebooklm CLI not found. Install with:\n"
            "  pip install 'notebooklm-py[browser]'\n"
            "  playwright install chromium\n"
            "  notebooklm login  # first-time auth"
        )
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Upload course content to NotebookLM for audio/quiz generation"
    )
    parser.add_argument(
        "course",
        choices=["prerequisites", "1.single_cell_processing_course", "2.statistics_de_course",
                 "3.clustering_annotation_course", "4.data_integration_course",
                 "5.trajectory_analysis_course", "6.cell_cell_communication_course", "all"],
        help="Course to upload (or 'all')",
    )
    parser.add_argument(
        "--generate",
        choices=["audio", "quiz", "flashcards", "slide-deck"],
        help="Generate artifact after upload (e.g. audio for podcast)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only list files, don't upload",
    )
    args = parser.parse_args()

    config = load_config()
    courses_config = config.get("courses", {})

    if args.course == "all":
        course_ids = [k for k in courses_config if k != "prerequisites"]
    else:
        course_ids = [args.course]

    for course_id in course_ids:
        paths = resolve_sources(course_id, config)
        if not paths:
            print(f"No sources found for {course_id}")
            continue

        name = courses_config[course_id].get("name", course_id)
        print(f"\n=== {name} ({len(paths)} files) ===")
        for p in paths:
            print(f"  {p.relative_to(REPO_ROOT)}")

        if args.dry_run:
            continue

        nb_name = f"scRNA-seq: {name}"
        print(f"\nCreating notebook: {nb_name}")
        if not run_notebooklm(["create", nb_name]):
            continue

        # notebooklm use selects the most recent notebook
        run_notebooklm(["use"])

        for p in paths:
            print(f"Adding: {p.name}")
            run_notebooklm(["source", "add", str(p)])

        if args.generate:
            print(f"\nGenerating {args.generate}...")
            if args.generate == "audio":
                run_notebooklm(["generate", "audio", "--wait"])
            elif args.generate == "quiz":
                run_notebooklm(["generate", "quiz", "--wait"])
            elif args.generate == "flashcards":
                run_notebooklm(["generate", "flashcards", "--wait"])
            elif args.generate == "slide-deck":
                run_notebooklm(["generate", "slide-deck", "--wait"])

    print("\nDone. Open https://notebooklm.google.com to chat, generate more content, or download.")


if __name__ == "__main__":
    main()
