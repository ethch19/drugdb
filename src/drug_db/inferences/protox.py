"""
TO DO:
- Graphs donwload (dist curves, toxicity radar, network chart)
"""

import re
from dataclasses import dataclass, field

import requests
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords


@dataclass
class ProToxConfig:
    ld50: bool = False
    tox_class: bool = False
    avg_sim: bool = False
    pred_acc: bool = False
    mol_props: bool = False
    sim_compounds: bool = False
    tox_targets: bool = False
    models: list[str] = field(default_factory=list)


class ProToxBuilder:
    def __init__(self):
        self.config = ProToxConfig()

    def with_ld50(self):
        self.config.ld50 = True
        return self

    def with_tox_class(self):
        self.config.tox_class = True
        return self

    def with_avg_sim(self):
        self.config.avg_sim = True
        return self

    def with_pred_acc(self):
        self.config.pred_acc = True
        return self

    def with_mol_props(self):
        self.config.mol_props = True
        return self

    def with_tox_targets(self):
        self.config.tox_targets = True
        return self

    def with_sim_compounds(self):
        self.config.sim_compounds = True
        return self

    def with_models(self, models: list[str]):
        self.config.models.extend(models)
        return self

    def with_all_models(self):
        self.config.models.extend(
            [
                "dili",
                "neuro",
                "nephro",
                "respi",
                "cardio",
                "carcino",
                "immuno",
                "mutagen",
                "cyto",
                "bbb",
                "eco",
                "clinical",
                "nutri",
                "nr_ahr",
                "nr_ar",
                "nr_ar_lbd",
                "nr_aromatase",
                "nr_er",
                "nr_er_lbd",
                "nr_ppar_gamma",
                "sr_are",
                "sr_hse",
                "sr_mmp",
                "sr_p53",
                "sr_atad5",
                "mie_thr_alpha",
                "mie_thr_beta",
                "mie_ttr",
                "mie_ryr",
                "mie_gabar",
                "mie_nmdar",
                "mie_ampar",
                "mie_kar",
                "mie_ache",
                "mie_car",
                "mie_pxr",
                "mie_nadhox",
                "mie_vgsc",
                "mie_nis",
                "CYP1A2",
                "CYP2C19",
                "CYP2C9",
                "CYP2D6",
                "CYP3A4",
                "CYP2E1",
            ]
        )
        return self

    def build(self) -> ProToxConfig:
        return self.config


class ProToxModel:
    BASE_URL = "https://tox.charite.de/protox3"

    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/142.0.0.0 Safari/537.36"
            }
        )

    @staticmethod
    def get_mol_block(smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        Compute2DCoords(mol)
        raw_block = Chem.MolToMolBlock(mol, forceV3000=False)
        lines = raw_block.split("\n")
        chemdoodle_header = [  # spoofer
            "Molecule from ChemDoodle Web Components",
            "",
            "http://www.ichemlabs.com",
        ]
        new_lines = chemdoodle_header + lines[3:]
        return "\n".join(new_lines)

    def execute(self, smiles: str, config: ProToxConfig):
        mol_block = ProToxModel.get_mol_block(smiles)
        if not mol_block:
            print("Unable to get mol_block")
            return None

        meta_payload = {
            "smilesString": mol_block,
            "defaultName": "Placeholder",
            "smiles": "Placeholder",
            "pubchem_name": "Placeholder",
        }

        results = {}

        # index.php?site=compound_search_similarity
        if any(
            [
                config.ld50,
                config.tox_class,
                config.mol_props,
                config.tox_targets,
                config.pred_acc,
                config.avg_sim,
            ]
        ):
            try:
                response_1 = self.session.post(
                    f"{self.BASE_URL}/index.php?site=compound_search_similarity",
                    data=meta_payload,
                    timeout=30,
                )
                response_1.raise_for_status()

                match = re.search(r"var server_id='(\d+)';", response_1.text)
                if match:
                    server_id = match.group(1)
                else:
                    print("ERROR: Could not find server_id in response")
                    return None

                soup = BeautifulSoup(response_1.text, "lxml")

                if config.ld50:
                    ld50_text = soup.find(string=re.compile("Predicted LD50"))
                    results["ld50"] = (
                        ld50_text.split(":")[1].strip() if ld50_text else None
                    )

                if config.tox_class:
                    tox_class_text = soup.find(
                        string=re.compile("Predicted Toxicity Class")
                    )  # GHS toxicity classification
                    results["tox_class"] = (
                        tox_class_text.split(":")[1].strip() if tox_class_text else None
                    )

                if config.avg_sim:
                    avg_sim_text = soup.find(string=re.compile("Average similarity"))
                    results["avg_sim"] = (
                        avg_sim_text.split(":")[1].strip() if avg_sim_text else None
                    )

                if config.pred_acc:
                    pred_acc_text = soup.find(
                        string=re.compile("Prediction accuracy")
                    )  # depends on the similarity of the input compound to compounds with known LD50 values as well as the hit rates obtained in a cross-validation study.
                    results["pred_acc"] = (
                        pred_acc_text.split(":")[1].strip() if pred_acc_text else None
                    )

                tables = soup.find_all(
                    "table", {"id": "table-out"}
                )  # compound props, similar compounds
                print(len(tables))
                if config.mol_props:
                    mol_properties = {}
                    if tables[0]:
                        rows = tables[0].find_all("tr")
                        for row in rows:
                            cols = row.find_all("td")
                            if len(cols) == 2:
                                key = cols[0].get_text(strip=True)
                                val = cols[1].get_text(strip=True)
                                mol_properties[key] = val
                    results["mol_props"] = mol_properties

                if config.sim_compounds:
                    sim_compounds = []
                    if len(tables) > 2:
                        for i in range(1, len(tables)):
                            rows = tables[i].find_all("tr")
                            comp_prop = {}
                            for j in range(2, len(rows)):
                                cols = rows[j].find_all("td")
                                key = cols[0].get_text(strip=True)
                                val = cols[1].get_text(strip=True)
                                comp_prop[key] = val
                            sim_compounds.append(comp_prop)
                    results["sim_compounds"] = sim_compounds

                if config.tox_targets:
                    targets = []
                    target_div = soup.find("div", {"id": "print_targets"})
                    tox_tables = target_div.find_all("table")
                    if tox_tables[1]:
                        rows = tox_tables[1].find_all("tr")
                        target_keys = []
                        for index, row in enumerate(rows):
                            if index == 0:
                                cols = row.find_all("th")
                                cols.pop(0)
                                for col in cols:
                                    target_keys.append(col.get_text(strip=True))
                                continue
                            cols = row.find_all("td")
                            if len(cols) == 4:
                                cols.pop(0)
                                target_vals = []
                                for col in cols:
                                    target_vals.append(col.get_text(strip=True))
                                targets.append(dict(zip(target_keys, target_vals)))
                    results["tox_targets"] = targets

            except Exception as e:
                print(f"ERROR: Cannot retrieve metadata: {e}")

        if config.models:
            try:
                response_2 = self.session.post(
                    f"{self.BASE_URL}/src/run_models.php",
                    data={
                        "models": " ".join(config.models),
                        "sdfile": "empty",
                        "mol": mol_block,
                        "id": server_id,
                    },
                    headers={
                        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/142.0.0.0 Safari/537.36",
                        "Accept": "*/*",
                        "Accept-Language": "en-GB,en-US;q=0.9,en;q=0.8",
                        "Origin": "https://tox.charite.de",
                        "Referer": "https://tox.charite.de/protox3/index.php?site=compound_search_similarity",
                        "X-Requested-With": "XMLHttpRequest",
                        "Sec-Fetch-Dest": "empty",
                        "Sec-Fetch-Mode": "cors",
                        "Sec-Fetch-Site": "same-origin",
                    },
                    timeout=45,
                )

                try:
                    results["models"] = response_2.json()
                except Exception:
                    print(
                        f"Failed to extract JSON from model predictions.\nStatus Code: {response_2.status_code}"
                    )
                    print("--- RAW RESPONSE DUMP START ---")
                    print(repr(response_2.text))
                    print("--- RAW RESPONSE DUMP END ---")

            except Exception as e:
                print(f"Error: Cannot retrieve model predictions: {e}")

        return results


if __name__ == "__main__":
    model = ProToxModel()

    request_full = (
        ProToxBuilder()
        .with_ld50()
        .with_avg_sim()
        .with_mol_props()
        .with_pred_acc()
        .with_tox_class()
        .with_tox_targets()
        .with_sim_compounds()
        # .with_all_models()
        .with_models(["dili", "CYP1A2"])
        .build()
    )

    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    data = model.execute(aspirin_smiles, request_full)

    import json

    print(json.dumps(data, indent=2))


"""
Organ Toxicity

dili Hepatotoxicity
neuro Neurotoxicity
nephro Nephrotoxicity
respi Respiratory toxicity 
cardio Cardiotoxicity


Toxicity end points

carcino Carcinogenicity
immuno Immunotoxicity
mutagen Mutagenicity
cyto Cytotoxicity
bbb BBB-barrier
eco Ecotoxicity 
clinical Clinical toxicity
nutri Nutritional toxicity


Tox21 Nuclear receptor signalling pathways

nr_ahr Aryl hydrocarbon Receptor (AhR)
nr_ar Androgen Receptor (AR)
nr_ar_lbd Androgen Receptor Ligand Binding Domain (AR-LBD)
nr_aromatase Aromatase
nr_er Estrogen Receptor Alpha (ER)
nr_er_lbd Estrogen Receptor Ligand Binding Domain (ER-LBD)
nr_ppar_gamma Peroxisome Proliferator Activated Receptor Gamma (PPAR-Gamma)

Tox21 Stress response pathways

sr_are Nuclear factor (erythroid-derived 2)-like 2/antioxidant responsive element (nrf2/ARE)
sr_hse Heat shock factor response element (HSE)
sr_mmp Mitochondrial Membrane Potential (MMP)
sr_p53 Phosphoprotein (Tumor Suppressor) p53
sr_atad5 ATPase family AAA domain containing protein 5 (ATAD5)

Molecular Initiating Events

mie_thr_alpha Thyroid hormone receptor alpha (THRa)
mie_thr_beta Thyroid hormone receptor beta (THRβ)
mie_ttr Transtyretrin (TTR)
mie_ryr Ryanodine receptor (RYR)
mie_gabar GABA receptor (GABAR)
mie_nmdar Glutamate N-methyl-D-aspartate receptor (NMDAR)
mie_ampar alpha-amino-3-hydroxy-5-methyl-4-isoxazolepropionate receptor (AMPAR)
mie_kar Kainate receptor (KAR)
mie_ache Achetylcholinesterase (AChE)
mie_car Constitutive androstane receptor (CAR)
mie_pxr Pregnane X receptor (PXR)
mie_nadhox NADH-quinone oxidoreductase (NADHOX)
mie_vgsc Voltage gated sodium channel (VGSC)
mie_nis Na+/I- symporter (NIS)


Metabolism

CYP1A2 Cytochrome CYP1A2
CYP2C19 Cytochrome CYP2C19
CYP2C9 Cytochrome CYP2C9
CYP2D6 Cytochrome CYP2D6
CYP3A4 Cytochrome CYP3A4
CYP2E1 Cytochrome CYP2E1
"""
