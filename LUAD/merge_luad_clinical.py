from pathlib import Path

import pandas as pd


BASE_DIR = Path("/mnt/disk1/qiuzerui/downloads/LUAD")
GSE123904_PATH = BASE_DIR / "GSE123904" / "GSE123904_LUAD_Patient_Clinical_Data.csv"
GSE131907_PATH = BASE_DIR / "GSE131907" / "GSE131907_Lung_Cancer_Feature_Summary.csv"
OUTPUT_PATH = Path("/home/zerui/code/files/LUAD_clinical_merged_from_image.csv")


IMAGE_MAPPING = [
    ("P1", "LUNG_T09"),
    ("P2", "LUNG_T08"),
    ("P3", "LX676"),
    ("P4", "LX661"),
    ("P5", "LUNG_T25"),
    ("P6", "LX682"),
    ("P7", "LUNG_T06"),
    ("P8", "LUNG_T34"),
    ("P9", "LUNG_T31"),
    ("P10", "LX653"),
    ("P11", "LUNG_T30"),
    ("P12", "LUNG_T19"),
    ("P13", "LX675"),
    ("P14", "LX684"),
    ("P15", "LUNG_T20"),
    ("P16", "LUNG_T18"),
    ("P17", "LUNG_T28"),
    ("M1", "NS_06"),
    ("M2", "NS_16"),
    ("M3", "NS_02"),
    ("M4", "NS_19"),
    ("M5", "NS_17"),
    ("M6", "NS_04"),
    ("M7", "NS_13"),
    ("M8", "NS_03"),
    ("M9", "NS_12"),
    ("M10", "NS_07"),
    ("C1", "LX701"),
    ("C2", "LX255B"),
    ("C3", "LX681"),
]


def normalize_lx_id(value):
    value = str(value).strip()
    if value.startswith("MSK-"):
        value = value.removeprefix("MSK-")
    return value


def first_present(row, candidates):
    for col in candidates:
        if col in row and pd.notna(row[col]):
            return row[col]
    return pd.NA


def prefix_row(row, prefix, columns):
    return {f"{prefix}__{col}": row[col] if col in row else pd.NA for col in columns}


def choose_gse123904_row(df, study_id, source_id):
    hits = df[df["normalized_patient"] == source_id]
    if hits.empty:
        return None

    if study_id.startswith("P"):
        primary = hits[hits["Tissue Site"].astype(str).str.lower().eq("primary")]
        if not primary.empty:
            return primary.iloc[0]

    return hits.iloc[0]


def main():
    gse123904 = pd.read_csv(GSE123904_PATH, encoding="utf-8-sig")
    gse131907 = pd.read_csv(GSE131907_PATH, encoding="utf-8-sig")

    gse123904["normalized_patient"] = gse123904["Patient"].map(normalize_lx_id)

    gse123904_cols = [col for col in gse123904.columns if col != "normalized_patient"]
    gse131907_cols = list(gse131907.columns)

    rows = []
    for study_id, source_id in IMAGE_MAPPING:
        source_dataset = pd.NA
        matched = None

        if source_id.startswith(("LUNG_", "NS_")):
            hits = gse131907[gse131907["Samples"].astype(str).str.strip() == source_id]
            if not hits.empty:
                source_dataset = "GSE131907"
                matched = hits.iloc[0]
        else:
            matched = choose_gse123904_row(gse123904, study_id, source_id)
            if matched is not None:
                source_dataset = "GSE123904"

        out = {
            "Patients ID in this study": study_id,
            "Patients ID in GSE131907 / GSE123904": source_id,
            "Source_Dataset": source_dataset,
            "Dataset_Patient_ID": pd.NA,
            "Dataset_Sample_ID": pd.NA,
            "Tissue": pd.NA,
            "Smoking": pd.NA,
            "Diagnosis_or_Histology": pd.NA,
            "Stage": pd.NA,
            "EGFR_or_Known_Oncogenic_Mutations": pd.NA,
        }

        out.update({f"GSE131907__{col}": pd.NA for col in gse131907_cols})
        out.update({f"GSE123904__{col}": pd.NA for col in gse123904_cols})

        if matched is None:
            rows.append(out)
            continue

        if source_dataset == "GSE131907":
            out["Dataset_Patient_ID"] = matched["Patient id"]
            out["Dataset_Sample_ID"] = matched["Samples"]
            out["Tissue"] = matched["Tissue origins"]
            out["Smoking"] = matched["Smoking"]
            out["Diagnosis_or_Histology"] = matched["Histology"]
            out["Stage"] = matched["Stages"]
            out["EGFR_or_Known_Oncogenic_Mutations"] = matched["EGFR"]
            out.update(prefix_row(matched, "GSE131907", gse131907_cols))
        elif source_dataset == "GSE123904":
            out["Dataset_Patient_ID"] = matched["Patient"]
            out["Dataset_Sample_ID"] = source_id
            out["Tissue"] = matched["Tissue Site"]
            out["Smoking"] = matched["Smoking History"]
            out["Diagnosis_or_Histology"] = matched["Diagnosis"]
            out["Stage"] = matched["Stage"]
            out["EGFR_or_Known_Oncogenic_Mutations"] = matched[
                "Known Oncogenic Mutations"
            ]
            out.update(prefix_row(matched, "GSE123904", gse123904_cols))

        rows.append(out)

    merged = pd.DataFrame(rows)
    merged = merged.fillna("NA")

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUTPUT_PATH, index=False)
    print(f"Wrote {merged.shape[0]} rows x {merged.shape[1]} columns to {OUTPUT_PATH}")

    unmatched = merged[merged["Source_Dataset"].eq("NA")]
    if not unmatched.empty:
        print("Unmatched rows:")
        print(unmatched[["Patients ID in this study", "Patients ID in GSE131907 / GSE123904"]])


if __name__ == "__main__":
    main()
