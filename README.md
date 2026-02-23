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
To run entire scraping pipeline (stored in `drugs.db` under `~/datastore` directory):
```sh
uv run src/drug_db/scraper/main.py
```

To run a specific scraper 
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

## Sources Metadata (`source_records*/data_json`)
### Drugbank ([link](https://go.drugbank.com/releases/latest))
| Key | Type |  Description |
| --- | --- | --- |
| unii | str | Unique Ingredient Identifier (UNII), searchable from `unii.py` or `UNII_Records*.csv` |
| description | str | Non-uniform paragraph about the drug |
| toxicity | str | Non-uniform paragraph about toxicity |
| moa | str | Non-uniform paragraph about drug's mechanism of action (MoA) |
| synonyms | list[str] | List of strings containing other names of the drug. Stored as a square-bracketed comma separated string |
| groups | list[str] | List containing the following: approved, illicit, experimental, withdrawn, nutraceutical, investigational, vet_approved | 
| categories | list[dict[str, any]] | List of different categories that the drug falls under. Pretty random strings. Stored as a dict, key=category name, value=mesh ID for that category |
| atc_codes | list[str] | List of 7-character alphanumerical Anatomical Therapeutic Chemical (ATC) classification codes from the WHO |

### Withdrawn 2.0 ([link](https://bioinformatics.charite.de/withdrawn_3/index.php))
| Key | Type | Description |
| --- | --- | --- |
| withdrawn_id | str | Unique ID specific to Withdrawn 2.0 |
| atc_code | str | List of 7-character alphanumerical Anatomical Therapeutic Chemical (ATC) classification codes from the WHO |
| pubchem_cid | str | PubChem Compounds ID, searchable through `pubchem/retrieval.py` |
| ctd_id | str | ID of compound from Comparative Toxicogenomics Database (CTD) |
| inchikey_jchem | str | InChIKey computed/sourced from JChem |
| ld50 | float | Lethal dose at which 50% of population died (`mg/kg`) - INFERENCED FROM PROTOX-II |
| protox_toxclass | int | Globally Harmonized System(GHS, 2005 edition) toxicity class in rodents, ranging I to VI - INFERENCED FROM PROTOX-II |
| first_withdrawn_year| int | First year of drug being withdrawn from any country |
| last_withdrawn_year | int | Last/Most recent? year of drug being withdrawn from any country |
| first_approval_year | int | First year of drug being approved for medical use in any country |
| withdrawal_countries | list[str] | List of ISO-alpha-3 country codes where drug has been withdrawn. Stored as a square-bracketed comma separated string |
| toxicity_types | list[str] | Endpoint toxicity (e.g hepatotoxicity). Not always available |

### Resp Tox ([link](https://www.mdpi.com/1999-4923/14/4/832))
> TableS1 from Supplementary Table's xlsx file

| Key | Type | Description |
| --- | --- | --- |
| respiratory_toxicity | bool | Binary classification of resp toxicity |

### Ames Mutagenicity V2 ([link](https://doc.ml.tu-berlin.de/toxbenchmark/))
| Key | Type | Description |
| --- | --- | --- |
| ames_mutagenic | bool | Binary classification of Ames Mutagenicity |

### LINCS L1000 Phase II 2020 ([link](https://clue.io/data/CMap2020#LINCS2020))
> Level 5 data used, reduced to landmark assays yielding 978 dimensions
- Check [documentation](https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit?usp=sharing) for details about how z-scores are calculated, and the definition of levels.

| Key | Type | Description |
| --- | --- | --- |
| sig_id | str | Signature ID, unique to LINCS |
| lincs_zscores_array | list[float] | Fixed 978-long vector storing ALL L1000 assays as z-scores |

### ChEMBL DB ([link](https://chembl.gitbook.io/chembl-interface-documentation/downloads))
> Only assay data extracted
Filtered for assays performed on:
- Minimium 100 molecules
- `organism` = 'Homo sapiens'
- `assay-type` = 'Binding' or 'Functional'
- `confidence-level` >= 8

| Key | Type | Description |
| --- | --- | --- |
| chembl_pchembl_array | list[float] | Fixed 1283-long vector storing PChemBl assay values associated with a particular molecule/drug |

