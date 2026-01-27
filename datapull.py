from pathlib import Path

from huggingface_hub import snapshot_download

CUR_DIR = Path(__file__).parent
ROOT_DIR = CUR_DIR.parent.parent.parent
DATASTORE_DIR = ROOT_DIR / "datastore"

REPO_ID = "ethch19/drugdb"


def download_repo(repo_id: str, local_dir: Path):
    local_dir.mkdir(parents=True, exist_ok=True)

    try:
        snapshot_download(
            repo_id=repo_id, repo_type="dataset", local_dir=local_dir, etag_timeout=30
        )
    except Exception as e:
        print(f"ERROR: cannot download dataset {e}")


if __name__ == "__main__":
    download_repo(REPO_ID, DATASTORE_DIR)
    print(f"Dataset download complete, in {DATASTORE_DIR}")
