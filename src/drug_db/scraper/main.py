import subprocess
import sys
import time
from pathlib import Path

SCRAPER_DIR = Path(__file__).parent.resolve()

SCRAPERS = [
    "withdrawn.py",
    "drugbank.py",
    "ames_muta.py",
    "resp_tox.py",
    "chembl.py",
    "lincs.py",
]


def run_pipeline():
    print("--- Starting DrugDB Scraper Pipeline ---")
    print(f"Found {len(SCRAPERS)} scrapers to execute in sequence.\n")

    for scraper_file in SCRAPERS:
        script_path = SCRAPER_DIR / scraper_file

        if not script_path.exists():
            print(f"ERROR: Could not find {scraper_file} at {script_path}")
            continue

        print(f"\n{'=' * 30}")
        print(f"RUNNING: {scraper_file}")
        print(f"{'=' * 30}")

        start_time = time.time()

        try:
            _ = subprocess.run(["uv", "run", str(script_path)], check=True, text=True)

            elapsed = time.time() - start_time
            print(f"\n[{scraper_file}] completed in {elapsed:.2f} seconds.")

        except subprocess.CalledProcessError as e:
            print(f"\nERROR: {scraper_file} crashed with exit code {e.returncode}.")
            print("Aborting...")
            sys.exit(1)

        except KeyboardInterrupt:
            sys.exit(0)

    print("\nALL SCRAPERS COMPLETED SUCCESSFULLY")


if __name__ == "__main__":
    run_pipeline()
