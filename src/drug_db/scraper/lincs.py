from typing import Iterator

import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
from sqlalchemy.orm import Session

from drug_db.models.schemas import DrugSchema, RecordSchema, ScraperPayload
from drug_db.scraper.adapters import BaseAdapter
from drug_db.scraper.base import BaseScraper


class LincsAdapter(BaseAdapter):
    def standardise(self, raw_data: dict) -> ScraperPayload:
        drug = DrugSchema(
            inchi_key=raw_data.get("inchi_key"),
            smiles=raw_data.get("canonical_smiles"),
            generic_name=raw_data.get("pert_iname"),
        )

        record = RecordSchema(
            data_json={
                "sig_id": raw_data.get("sig_id"),
                "lincs_zscores_array": raw_data.get("zscores"),
            }
        )
        return ScraperPayload(drug=drug, record=record)


class LincsScraper(BaseScraper):
    def __init__(
        self,
        session: Session,
        gctx_path: str,
        sig_info_path: str,
        comp_info_path: str,
        gene_info_path: str,
        chunk_size: int = 5000,
    ):
        (
            self.gctx_path,
            self.sig_info_path,
            self.comp_info_path,
            self.gene_info_path,
            self.chunk_size,
        ) = (
            gctx_path,
            sig_info_path,
            comp_info_path,
            gene_info_path,
            chunk_size,
        )
        super().__init__(
            session, source_name="LINCS_L1000_PhaseII", version="CMAP_2020"
        )

    def _get_adapter(self) -> BaseAdapter:
        return LincsAdapter()

    def _extract(self) -> Iterator[dict]:
        sig_info = pd.read_csv(self.sig_info_path, sep="\t", low_memory=False)
        comp_info = pd.read_csv(self.comp_info_path, sep="\t", low_memory=False)
        gene_info = pd.read_csv(self.gene_info_path, sep="\t", dtype=str)

        # CMAP 2020 Beta different column names with CMAP 2017
        if "feature_space" in gene_info.columns:
            landmark_ids = gene_info[gene_info["feature_space"] == "landmark"][
                "gene_id"
            ].tolist()
        elif "pr_is_lm" in gene_info.columns:
            landmark_ids = gene_info[gene_info["pr_is_lm"] == "1"][
                "pr_gene_id"
            ].tolist()
        else:
            raise KeyError(
                f"Could not find landmark column. Available columns: {gene_info.columns.tolist()}"
            )

        merged_info = (
            pd.merge(
                sig_info[sig_info["pert_type"] == "trt_cp"],
                comp_info,
                on="pert_id",
                how="left",
            )
            .dropna(subset=["inchi_key"])
            .query("inchi_key != '-666' and inchi_key != 'restricted'")
        )

        all_col_meta = parse(str(self.gctx_path), col_meta_only=True)
        all_col_ids = all_col_meta.index.tolist()
        id_to_idx = {sig_id: i for i, sig_id in enumerate(all_col_ids)}

        unique_valid_sig_ids = [
            sid for sid in merged_info["sig_id"].unique() if sid in id_to_idx
        ]
        print(f"Processing {len(unique_valid_sig_ids)} unique compound signatures...")

        for i in range(0, len(unique_valid_sig_ids), self.chunk_size):
            batch_ids = unique_valid_sig_ids[i : i + self.chunk_size]
            batch_indices = sorted(list(set(id_to_idx[sid] for sid in batch_ids)))

            print(f"Vector-parsing LINCS batch {i // self.chunk_size + 1}...")

            gct_obj = parse(str(self.gctx_path), rid=landmark_ids, cidx=batch_indices)
            matrix = gct_obj.data_df

            batch_meta = merged_info[merged_info["sig_id"].isin(matrix.columns)]

            for sig_id, meta_group in batch_meta.groupby("sig_id"):
                zscores = matrix[sig_id].tolist()
                meta = meta_group.iloc[0]  # ignore duplicates/aliases

                cmap_name = meta.get("cmap_name")
                pert_iname = meta.get("pert_iname")

                if pd.notna(cmap_name):
                    generic_name = cmap_name
                elif pd.notna(pert_iname):
                    generic_name = pert_iname
                else:
                    generic_name = None

                smiles = meta.get("canonical_smiles")
                if pd.isna(smiles) or smiles == ["-666", "restricted"]:
                    smiles = None

                yield {
                    "sig_id": sig_id,
                    "pert_iname": generic_name,
                    "inchi_key": meta["inchi_key"],
                    "canonical_smiles": smiles,
                    "zscores": zscores,
                }


if __name__ == "__main__":
    import sys
    from pathlib import Path

    from drug_db.database.manager import db

    BASE_DIR = Path(__file__).parent.parent.parent.resolve()
    DATA_DIR = BASE_DIR / ".." / "data"
    gctx_path = DATA_DIR / "level5_beta_trt_cp_n720216x12328.gctx"
    compoundinfo_path = DATA_DIR / "compoundinfo_beta.txt"
    siginfo_path = DATA_DIR / "siginfo_beta.txt"
    geneinfo_path = DATA_DIR / "geneinfo_beta.txt"

    db.create_tables()
    session = db.get_session()
    try:
        with session:
            scraper = LincsScraper(
                session,
                gctx_path.resolve(),
                siginfo_path.resolve(),
                compoundinfo_path.resolve(),
                geneinfo_path.resolve(),
            )
            scraper.run()

    except KeyboardInterrupt:
        sys.exit(0)
    finally:
        session.close()
