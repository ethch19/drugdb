# DrugDB

## Project Structure
- **Code** stored here on [Github](https://github.com/ethch19/drugdb)
- **Processed data** stored on [Huggingface](https://huggingface.co/datasets/ethch19/drugdb)
- **Raw data** not hosted

## Prerequisites
- Install [uv](https://docs.astral.sh/uv/getting-started/installation/) (or `pip`)
- Python 3.14+

## Installation
Run these commands
```sh
git clone https://github.com/ethch19/drugdb.git # or ssh alternatively
cd drugdb
uv sync
uv run datapull.py # pulls processed data from Huggingface
```

## How to use
To run a specific scraper (stored in `drugs.db` under `~/datastore` directory):
```sh
uv run src/drug_db/scraper/withdrawn.py
```

To run model inferences:
```sh
uv run src/drug_db/inferences/protox.py
```

Or: use programmatically by importing modules directly in your own scripts

## Database Schema (`drugs.db`)
### `drugs` table
| Column        | Type           | Description  |
| ------------- |:-------------:| -----:|
| id      | INT(PK) | Unique ID |
| inchi_key      | TEXT | 27-char hashed InChI |
| generic_name      | TEXT      |  INN name  |
| synonyms | JSON      |  List of drug names   |
| smiles | TEXT      |  Canonical/Isomeric SMILES   |
| inchi | TEXT      |  InChI string   |
| chem_formula | TEXT      |  Molecular formula   |
| mol_weight | FLOAT      |  Molecular weight in Daltons   |
| drugbank_id | TEXT      |  N/A   |
| chembl_id | TEXT      |  N/A   |


### `sources` table
| Column        | Type           | Description  |
| ------------- |:-------------:| -----:|
| id      | INT(PK) | Unique ID |
| name      | TEXT | Name of source |
| version      | TEXT     |  Source version  |
| date_accessed | TEXT      |  ISO8601 date  |
| url | TEXT      |  For online sources  |


### `source_records` table
| Column        | Type           | Description  |
| ------------- |:-------------:| -----:|
| id      | INT(PK) | Unique ID |
| drug_inchi_key      | TEXT(FK) | References `drugs` table |
| source_id      | TEXT(FK)      |  References `sources` table  |
| data_json | JSON      |  Filtered data stored as JSON  |