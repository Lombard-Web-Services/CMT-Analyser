#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cmt_vcf_analyzer.py
~~~~~~~~~~~~~~~~~~~

Analyse un fichier de variantes MyHeritage (CSV) contre le VCF ClinVar
à la recherche de variants pathogènes liés à la maladie de Charcot‑Marie‑Tooth (CMT).

Version 3.1 - Avec débogage pour les gènes critiques
"""

import re
import sys
import gzip
import pandas as pd
from collections import defaultdict

# ----------------------------------------------------------------------
# -------------------------- CONFIGURATION -----------------------------
# ----------------------------------------------------------------------
# Sexe de la personne testée (modifier selon le cas)
SEX = "male"  # "male" ou "female"

# Activer le débogage pour les gènes sensibles
DEBUG_GENES = ["BRCA1", "BRCA2", "MLH1", "APC", "GJB1", "SH3TC2"]

GENES_CMT = {
    "PMP22", "MPZ", "MFNR", "GJB1", "EGR2", "LITAF", "HSPB1",
    "HSPB8", "DNM2", "YARS", "GDAP1", "NEFL", "FIG4", "SBF2",
    "SH3TC2", "GDAP1", "LRSAM1", "TRPV4", "MED25", "FGD4"
}

CLNSIG_SCORE = {
    "pathogenic": 5,
    "likely_pathogenic": 4,
    "pathogenic/likely_pathogenic": 4,
    "conflicting_interpretations": 2,
    "uncertain_significance": 1,
    "benign": 0,
    "likely_benign": 0,
}

GENE_INHERITANCE = defaultdict(
    lambda: "UNKNOWN",
    {
        "PMP22": "AD", "MPZ": "AD", "MFNR": "AD", "GJB1": "XLR",
        "EGR2": "AD", "LITAF": "AD", "HSPB1": "AD", "HSPB8": "AD",
        "DNM2": "AD", "YARS": "AD", "GDAP1": "AR", "NEFL": "AD",
        "FIG4": "AD", "SBF2": "AR", "SH3TC2": "AR", "LRSAM1": "AR",
        "TRPV4": "AD", "MED25": "AR", "FGD4": "AR",
    }
)

# ----------------------------------------------------------------------
# -------------------------- UTILITAIRES -------------------------------
# ----------------------------------------------------------------------
def norm_chrom(chrom: str) -> str:
    """Normalise les représentations de chromosome."""
    if not isinstance(chrom, str):
        chrom = str(chrom)
    chrom = chrom.strip().upper()
    if chrom.startswith("CHR"):
        chrom = chrom[3:]
    if chrom == "23":
        return "X"
    if chrom == "24":
        return "Y"
    if chrom in ("MT", "M"):
        return "MT"
    return chrom


def complement_base(b: str) -> str:
    """Renvoie le complément d'une base."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return comp.get(b.upper(), b)


def complement_sequence(seq: str) -> str:
    """Renvoie le complément d'une séquence."""
    return ''.join(complement_base(b) for b in reversed(seq))


def parse_genotype(gt: str):
    """
    Parse un génotype MyHeritage en allèles normalisés.
    Pour un homme sur le chromosome X, on ne duplique PAS l'allèle.
    """
    if not isinstance(gt, str):
        return []
    
    gt = gt.upper().strip()
    
    # Cas spéciaux : insertions/délétions
    if gt in ["II", "DD", "--", ".."]:
        return [gt, gt]
    
    # Cas avec séparateurs
    if "/" in gt:
        alleles = gt.split("/")
    elif "|" in gt:
        alleles = gt.split("|")
    elif ":" in gt:
        alleles = gt.split(":")
    else:
        # Pas de séparateur -> un seul allèle (hémizygote) ou deux identiques
        if len(gt) == 1:
            return [gt]  # Un seul allèle (hémizygote)
        elif len(gt) == 2:
            return [gt[0], gt[1]]
        else:
            # Indel
            if gt in ["INS", "DEL", "INS/DEL"]:
                return [gt, gt]
            half = len(gt) // 2
            return [gt[:half], gt[half:]]
    
    # Nettoyer les allèles
    alleles = [a.strip() for a in alleles if a.strip()]
    
    return alleles


def normalize_allele(allele: str) -> str:
    """Normalise un allèle pour la comparaison."""
    if not allele:
        return ""
    return allele.upper()


def compare_alleles(user_allele: str, ref: str, alt: str) -> str:
    """
    Compare un allèle de l'utilisateur avec REF/ALT.
    Retourne : "REF", "ALT", "ALT_COMP", "OTHER"
    """
    user_norm = normalize_allele(user_allele)
    ref_norm = normalize_allele(ref)
    alt_norm = normalize_allele(alt)
    
    # Comparaison exacte
    if user_norm == ref_norm:
        return "REF"
    if user_norm == alt_norm:
        return "ALT"
    
    # Comparaison avec complément (brin opposé)
    alt_comp = complement_sequence(alt_norm)
    if user_norm == alt_comp:
        return "ALT_COMP"
    
    # Gestion des indels
    if alt_norm in ["INS", "DEL", "I", "D"]:
        if user_norm in ["INS", "DEL", "I", "D"]:
            return "ALT"
    
    return "OTHER"


def get_zygosity(user_gt: str, ref: str, alt: str, chrom: str):
    """
    Détermine la zygosité réelle en comparant avec REF et ALT.
    Gère correctement les hémizygotes (chromosome X chez l'homme).
    """
    alleles = parse_genotype(user_gt)
    if not alleles:
        return ("UNKNOWN", 0)
    
    # Détection hémizygote : un seul allèle sur le chromosome X
    is_hemizygous = (len(alleles) == 1 and chrom == "X" and SEX == "male")
    
    alt_count = 0
    allele_status = []
    
    for allele in alleles:
        status = compare_alleles(allele, ref, alt)
        allele_status.append(status)
        if status in ["ALT", "ALT_COMP"]:
            alt_count += 1
        elif status == "REF":
            alt_count += 0
    
    # Cas hémizygote
    if is_hemizygous:
        if alt_count == 1:
            return ("HEM_ALT", 1)
        elif alt_count == 0:
            return ("HEM_REF", 0)
    
    # Cas diploïde normal
    if alt_count == 2:
        return ("ALT_ALT", 2)
    elif alt_count == 1:
        return ("HET", 1)
    elif alt_count == 0:
        return ("REF_REF", 0)
    else:
        return ("UNKNOWN", 0)


def carries_alt(user_gt: str, ref: str, alt: str, chrom: str) -> bool:
    """Vérifie si l'utilisateur porte l'allèle alternatif."""
    zygosity_type, alt_count = get_zygosity(user_gt, ref, alt, chrom)
    return alt_count > 0


def parse_info_field(info: str, tag: str):
    """Extrait la valeur d'un champ INFO du VCF."""
    pattern = rf"{tag}=([^;]+)"
    m = re.search(pattern, info)
    if not m:
        return []
    raw = m.group(1)
    return [v.strip() for v in raw.split(",") if v.strip()]


def get_clnsig_score(clnsig_list):
    """Retourne le score le plus sévère pour CLNSIG."""
    best_score = 0
    best_label = None
    for token in clnsig_list:
        token_low = token.lower().replace("/", "_")
        for label, score in CLNSIG_SCORE.items():
            if label in token_low and score > best_score:
                best_score = score
                best_label = label
    return best_label, best_score


def get_review_status_score(review_list):
    """Retourne un score de 0 à 5 pour CLNREVSTAT."""
    if not review_list:
        return 0
    score = 0
    for token in review_list:
        t = token.lower().replace("_", " ")
        if "practice guideline" in t:
            score = max(score, 5)
        elif "expert panel" in t:
            score = max(score, 4)
        elif "multiple submitters" in t:
            score = max(score, 3)
        elif "single submitter" in t:
            score = max(score, 2)
        else:
            score = max(score, 1)
    return score


def open_vcf(path):
    """Ouvre un fichier VCF compressé ou non."""
    try:
        return gzip.open(path, "rt")
    except OSError:
        return open(path, "rt")


# ----------------------------------------------------------------------
# -------------------------- CHARGEMENT MYHERITAGE --------------------
# ----------------------------------------------------------------------
def load_myheritage(csv_path: str):
    """Charge le fichier MyHeritage CSV."""
    try:
        df = pd.read_csv(
            csv_path,
            dtype=str,
            sep=",",
            comment="#",
            engine="c",
            quoting=1,
            keep_default_na=False
        )
    except Exception as e:
        sys.exit(f"❌ Erreur lors de la lecture du fichier CSV : {e}")

    df.columns = [c.strip().upper() for c in df.columns]

    required = ["RSID", "CHROMOSOME", "POSITION", "RESULT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"❌ Colonnes manquantes : {', '.join(missing)}")

    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].str.strip('"')

    df["CHROMOSOME"] = df["CHROMOSOME"].apply(norm_chrom)
    df["POSITION"] = pd.to_numeric(df["POSITION"], errors="coerce").astype("Int64")
    df["KEY_POS"] = df["CHROMOSOME"] + ":" + df["POSITION"].astype(str)

    df = df.drop_duplicates(subset=["KEY_POS"], keep="first")

    dict_by_pos = df.set_index("KEY_POS")[["RSID", "RESULT", "CHROMOSOME", "POSITION"]].to_dict("index")

    dict_by_rsid = {}
    for _, row in df.iterrows():
        rsid = row["RSID"]
        if pd.notna(rsid) and rsid != "" and rsid != ".":
            dict_by_rsid.setdefault(rsid, []).append(row.to_dict())

    print(f"   ⚠️ {len(dict_by_pos)} positions uniques")
    return df, dict_by_pos, dict_by_rsid


# ----------------------------------------------------------------------
# -------------------------- ANALYSE VCF -------------------------------
# ----------------------------------------------------------------------
def analyze_vcf(vcf_path, dict_by_pos, dict_by_rsid):
    """Analyse le VCF ClinVar."""
    results = []
    debug_count = 0

    with open_vcf(vcf_path) as f:
        # Lire l'en-tête
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                idx_chrom = cols.index("#CHROM")
                idx_pos = cols.index("POS")
                idx_id = cols.index("ID")
                idx_ref = cols.index("REF")
                idx_alt = cols.index("ALT")
                idx_info = cols.index("INFO")
                break

        # Traiter les variants
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            chrom_raw = parts[idx_chrom]
            pos_raw = parts[idx_pos]
            rsid_raw = parts[idx_id]
            ref = parts[idx_ref]
            alt_list = parts[idx_alt].split(",")
            info_field = parts[idx_info]

            chrom = norm_chrom(chrom_raw)
            try:
                pos = int(pos_raw)
            except ValueError:
                continue

            # Recherche dans les index
            key_pos = f"{chrom}:{pos}"
            mh_entry = dict_by_pos.get(key_pos)

            if not mh_entry and rsid_raw not in (".", ""):
                mh_entry_list = dict_by_rsid.get(rsid_raw)
                if mh_entry_list:
                    mh_entry = mh_entry_list[0]

            if not mh_entry:
                continue

            # ---- Extraction des champs INFO ----
            clnsig_list = parse_info_field(info_field, "CLNSIG")
            review_list = parse_info_field(info_field, "CLNREVSTAT")
            clndn_list = parse_info_field(info_field, "CLNDN")
            gene_info = parse_info_field(info_field, "GENEINFO")

            # ---- Détermination du gène ----
            gene = "UNKNOWN"
            if gene_info and isinstance(gene_info, list) and len(gene_info) > 0:
                gene_str = gene_info[0]
                if "|" in gene_str:
                    first = gene_str.split("|")[0]
                else:
                    first = gene_str
                if ":" in first:
                    gene = first.split(":")[0]
                else:
                    gene = first

            is_cmt_gene = gene in GENES_CMT

            # ---- DEBUG : Afficher les gènes critiques ----
            if gene in DEBUG_GENES:
                debug_count += 1
                vos_gt = mh_entry["RESULT"]
                parsed_gt = parse_genotype(vos_gt)
                zygosity_type, alt_count = get_zygosity(vos_gt, ref, alt_list[0], chrom) if alt_list else ("UNKNOWN", 0)
                
                print(f"\n🔍 DEBUG #{debug_count} - Gène: {gene}")
                print(f"   Chromosome: {chrom}, Position: {pos}")
                print(f"   RSID: {rsid_raw}")
                print(f"   REF: {ref}, ALT: {','.join(alt_list)}")
                print(f"   Génotype MyHeritage: {vos_gt}")
                print(f"   Allèles parsés: {parsed_gt}")
                print(f"   Zygosité calculée: {zygosity_type} (alt_count={alt_count})")
                print(f"   Sexe: {SEX}, Hémizygote? {chrom == 'X' and SEX == 'male'}")
                
                # Analyse détaillée des allèles
                if parsed_gt:
                    for i, allele in enumerate(parsed_gt):
                        status = compare_alleles(allele, ref, alt_list[0]) if alt_list else "N/A"
                        print(f"   Allèle {i+1}: {allele} -> {status}")
                
                print(f"   {'✅' if alt_count > 0 else '❌'} Allèle alternatif présent: {alt_count > 0}")

            # ---- Scores ----
            clnsig_label, clnsig_score = get_clnsig_score(clnsig_list)
            review_score = get_review_status_score(review_list)

            # ---- Zygosité RÉELLE avec gestion du chromosome X ----
            vos_gt = mh_entry["RESULT"]

            best_zygosity = "UNKNOWN"
            best_alt_count = 0
            best_alt_match = False
            best_alt = None

            for alt in alt_list:
                zygosity_type, alt_count = get_zygosity(vos_gt, ref, alt, chrom)
                if alt_count > best_alt_count:
                    best_alt_count = alt_count
                    best_zygosity = zygosity_type
                    best_alt_match = alt_count > 0
                    best_alt = alt

            # ---- Scores de zygosité adaptés ----
            if best_zygosity == "ALT_ALT":
                zygo_score = 2
            elif best_zygosity == "HEM_ALT":
                zygo_score = 2  # Pour un homme, hémizygote alt = aussi grave que homozygote
            elif best_zygosity == "HET":
                zygo_score = 1
            else:
                zygo_score = 0

            # ---- Inheritance match adapté ----
            inher_score = 0
            if is_cmt_gene and best_zygosity != "UNKNOWN":
                exp = GENE_INHERITANCE[gene]
                if exp != "UNKNOWN":
                    if exp in ("AR", "MT"):
                        # Récessif autosomique : besoin de ALT_ALT
                        if best_zygosity == "ALT_ALT":
                            inher_score = 2
                    elif exp == "XLR":
                        # Récessif lié à l'X : ALT_ALT ou HEM_ALT
                        if best_zygosity in ("ALT_ALT", "HEM_ALT"):
                            inher_score = 2
                    else:  # AD, XLD
                        # Dominant : HET, ALT_ALT, ou HEM_ALT
                        if best_zygosity in ("HET", "ALT_ALT", "HEM_ALT"):
                            inher_score = 2

            # ---- Score total ----
            total_score = (
                clnsig_score * 2 +
                review_score * 1.5 +
                zygo_score * 1 +
                inher_score * 1.5 +
                (2 if best_alt_match else 0)
            )

            # ---- Construction du résultat ----
            result = {
                "CHROM": chrom,
                "POS": pos,
                "RSID": rsid_raw if rsid_raw not in (".", "") else None,
                "REF": ref,
                "ALT": ",".join(alt_list),
                "GENE": gene,
                "IS_CMT_GENE": is_cmt_gene,
                "CLNSIG": ";".join(clnsig_list) if clnsig_list else None,
                "CLNSIG_LABEL": clnsig_label,
                "CLNSIG_SCORE": clnsig_score,
                "CLNREVSTAT": ";".join(review_list) if review_list else None,
                "REVIEW_SCORE": review_score,
                "CLNDN": ";".join(clndn_list) if clndn_list else None,
                "MYHERITAGE_RESULT": vos_gt,
                "ZYGOSITY_TYPE": best_zygosity,
                "ZYGOSITY_SCORE": zygo_score,
                "ALT_COUNT": best_alt_count,
                "INHERITANCE_SCORE": inher_score,
                "ALT_PRESENT": best_alt_match,
                "TOTAL_SCORE": round(total_score, 2),
                "SEX": SEX,
            }
            results.append(result)

    print(f"\n📊 DEBUG: {debug_count} variants dans les gènes critiques analysés")
    return results


# ----------------------------------------------------------------------
# -------------------------- MAIN ---------------------------------------
# ----------------------------------------------------------------------
def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python cmt_vcf_analyzer.py <myheritage.csv> <clinvar.vcf.gz>")

    myheritage_file = sys.argv[1]
    vcf_file = sys.argv[2]

    print(f"🧬 Sexe configuré : {SEX}")
    print(f"🔍 Debug pour les gènes : {', '.join(DEBUG_GENES)}")
    print("🔄 Chargement du fichier MyHeritage...")
    df_mh, dict_by_pos, dict_by_rsid = load_myheritage(myheritage_file)
    print(f"   → {len(df_mh)} variants chargés.")

    print("🔄 Analyse du VCF ClinVar...")
    results = analyze_vcf(vcf_file, dict_by_pos, dict_by_rsid)
    print(f"   → {len(results)} variants trouvés dans le VCF et présents dans MyHeritage.")

    if not results:
        print("⚠️ Aucun variant correspondant trouvé.")
        return

    # Filtrer les vrais ALT
    results_filtered = [r for r in results if r["ALT_PRESENT"]]
    print(f"   → {len(results_filtered)} variants avec allèle alternatif présent.")

    if not results_filtered:
        print("⚠️ Aucun variant avec allèle alternatif trouvé.")
        return

    # Filtrer les gènes CMT en priorité
    results_cmt = [r for r in results_filtered if r["IS_CMT_GENE"]]
    print(f"   → {len(results_cmt)} variants dans les gènes CMT.")

    # Tri par score
    results_sorted = sorted(results_filtered, key=lambda x: x["TOTAL_SCORE"], reverse=True)

    # Affichage - priorité aux gènes CMT
    print("\n📊 Top 10 variantes (par score total) :")
    print("{:<6} {:<8} {:<10} {:<6} {:<12} {:<10} {:<15} {:<8}".format(
        "CHROM", "POS", "RSID", "GENE", "CLNSIG", "SCORE", "ZYGOSITE", "CMT?"))

    count = 0
    # Afficher d'abord les gènes CMT
    for r in results_sorted:
        if count >= 10:
            break
        if r["IS_CMT_GENE"]:
            print("{:<6} {:<8} {:<10} {:<6} {:<12} {:<10} {:<15} {:<8}".format(
                r["CHROM"],
                r["POS"],
                r["RSID"] if r["RSID"] else ".",
                r["GENE"],
                r["CLNSIG_LABEL"] if r["CLNSIG_LABEL"] else ".",
                r["TOTAL_SCORE"],
                r["ZYGOSITY_TYPE"],
                "✅ CMT"
            ))
            count += 1
    
    # Compléter avec les non-CMT si besoin
    if count < 10:
        for r in results_sorted:
            if count >= 10:
                break
            if not r["IS_CMT_GENE"]:
                print("{:<6} {:<8} {:<10} {:<6} {:<12} {:<10} {:<15} {:<8}".format(
                    r["CHROM"],
                    r["POS"],
                    r["RSID"] if r["RSID"] else ".",
                    r["GENE"],
                    r["CLNSIG_LABEL"] if r["CLNSIG_LABEL"] else ".",
                    r["TOTAL_SCORE"],
                    r["ZYGOSITY_TYPE"],
                    ""
                ))
                count += 1

    # Affichage spécifique des gènes CMT
    if results_cmt:
        print("\n📊 Variantes dans les gènes CMT :")
        print("{:<6} {:<8} {:<10} {:<6} {:<12} {:<10} {:<15}".format(
            "CHROM", "POS", "RSID", "GENE", "CLNSIG", "SCORE", "ZYGOSITE"))
        for r in sorted(results_cmt, key=lambda x: x["TOTAL_SCORE"], reverse=True)[:20]:
            print("{:<6} {:<8} {:<10} {:<6} {:<12} {:<10} {:<15}".format(
                r["CHROM"],
                r["POS"],
                r["RSID"] if r["RSID"] else ".",
                r["GENE"],
                r["CLNSIG_LABEL"] if r["CLNSIG_LABEL"] else ".",
                r["TOTAL_SCORE"],
                r["ZYGOSITY_TYPE"]
            ))

    # Export
    out_csv = "resultats_cmt_vrais_positifs.csv"
    out_df = pd.DataFrame(results_sorted)
    out_df.to_csv(out_csv, index=False, sep="\t")
    print(f"\n✅ Résultats détaillés enregistrés dans : {out_csv}")
    
    # Export spécifique CMT
    if results_cmt:
        out_csv_cmt = "resultats_cmt_uniquement.csv"
        out_df_cmt = pd.DataFrame(results_cmt)
        out_df_cmt.to_csv(out_csv_cmt, index=False, sep="\t")
        print(f"✅ Résultats CMT uniquement enregistrés dans : {out_csv_cmt}")


if __name__ == "__main__":
    main()
