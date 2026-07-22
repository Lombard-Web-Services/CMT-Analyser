# CMT-Analyse

Outil d'analyse genetique locale pour le depistage de variants pathogenes associes a la maladie de Charcot-Marie-Tooth (CMT) a partir de donnees de genotypage MyHeritage et du VCF ClinVar.

---

## Table des matieres

1. [Description](#description)
2. [Prerequis](#prerequis)
3. [Installation](#installation)
4. [Utilisation](#utilisation)
5. [Architecture](#architecture)
6. [Fichiers](#fichiers)
7. [Methodologie](#methodologie)
8. [Scores et ponderation](#scores-et-ponderation)
9. [Limites](#limites)
10. [Licence](#licence)

---

## Description

Ce projet fournit deux scripts d'analyse complementaires permettant de croiser un fichier de variants genetiques issu d'un test MyHeritage avec la base ClinVar (NCBI) afin d'identifier d'eventuels variants pathogenes impliques dans la maladie de Charcot-Marie-Tooth (CMT).

L'outil est entierement local, ne necessite aucun envoi de donnees genetiques sur un serveur distant et s'execute en ligne de commande.

---

## Prerequis

- Python 3.8 ou superieur
- pandas
- Fichier de donnees MyHeritage au format CSV
- Fichier ClinVar au format VCF.GZ (GRCh37 recommande)

### Telechargement de ClinVar

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
```

---

## Installation

```bash
git clone https://github.com/Lombard-Web-Services/CMT-Analyse.git
cd CMT-Analyse
pip install pandas
```

---

## Utilisation

### Script v1 : Analyse simplifiee

```bash
python analyse_cmt_vcf.py
```

Le script attend les fichiers `MH.csv` et `clinvar.vcf.gz` dans le repertoire courant.

### Script v2 : Analyse avancee avec scoring

```bash
python cmt_vcf_analyzer.py <myheritage.csv> <clinvar.vcf.gz>
```

Exemple :

```bash
python cmt_vcf_analyzer.py donnees_MH.csv clinvar.vcf.gz
```

---

## Fichiers

| Fichier | Description |
|---------|-------------|
| `analyse_cmt_vcf.py` | Script principal v1, analyse par position et RSID avec verification de brin |
| `cmt_vcf_analyzer.py` | Script v2, analyse avancee avec systeme de scoring multi-criteres |
| `MH.csv` | Fichier de donnees MyHeritage (fourni par l'utilisateur) |
| `clinvar.vcf.gz` | Base ClinVar compressee (a telecharger separement) |
| `resultats_cmt_vrais_positifs.csv` | Export des variants confirmes (genere) |
| `resultats_cmt_uniquement.csv` | Export filtre sur les genes CMT uniquement (genere) |

---

## Architecture

```
CMT-Analyse/
|-- test.py          # Script v1
|-- test2.py         # Script v2
|-- MH.csv                      # Donnees MyHeritage (input)
|-- clinvar.vcf.gz              # Base ClinVar (input)
|-- resultats_cmt_vrais_positifs.csv   # Resultats (output)
|-- resultats_cmt_uniquement.csv         # Resultats filtres (output)
|-- README.md
```

---

## Methodologie

### 1. Indexation des donnees MyHeritage

Les donnees sont indexees de maniere bidirectionnelle :
- Par position chromosomique (`CHR:POS`)
- Par identifiant RSID

Cette double indexation permet de maximiser le taux de correspondance meme en cas de decalage de position entre les deux bases.

### 2. Filtrage des genes CMT

Les scripts ciblent les genes suivants, associes a la CMT :

| Gene | Chromosome | Transmission |
|------|------------|--------------|
| PMP22 | 17 | Autosomique dominante |
| MPZ | 1 | Autosomique dominante |
| MFN2 | 1 | Autosomique dominante |
| GJB1 | X | Liee a l'X |
| EGR2 | 10 | Autosomique dominante |
| LITAF | 16 | Autosomique dominante |
| HSPB1 | 7 | Autosomique dominante |
| HSPB8 | 12 | Autosomique dominante |
| DNM2 | 19 | Autosomique dominante |
| YARS | 1 | Autosomique dominante |
| GDAP1 | 8 | Autosomique recessive |
| NEFL | 8 | Autosomique dominante |
| FIG4 | 6 | Autosomique dominante |
| SBF2 | 11 | Autosomique recessive |
| SH3TC2 | 5 | Autosomique recessive |
| LRSAM1 | 9 | Autosomique recessive |
| TRPV4 | 12 | Autosomique dominante |
| MED25 | 19 | Autosomique recessive |
| FGD4 | 12 | Autosomique recessive |

### 3. Classification ClinVar

Les variants sont filtres selon leur classification ClinVar :
- Pathogenic
- Likely pathogenic
- Pathogenic / Likely pathogenic

Les classifications Benign, Likely benign, Uncertain significance et Conflicting interpretations sont exclues du rapport final.

### 4. Verification des alleles

Le systeme compare les alleles de l'utilisateur avec la reference (REF) et l'allele alternatif (ALT) du VCF :
- Comparaison directe
- Comparaison par complement (inversion de brin)
- Gestion des indels (insertions / deletions)
- Gestion des hemizygotes (chromosome X chez l'homme)

### 5. Determination de la zygosite

| Type | Description |
|------|-------------|
| REF_REF | Homozygote reference |
| HET | Heterozygote alternatif |
| ALT_ALT | Homozygote alternatif |
| HEM_ALT | Hemizygote alternatif (X, male) |
| HEM_REF | Hemizygote reference (X, male) |
| UNKNOWN | Non determinable |

---

## Scores et ponderation

Le script v2 attribue un score total base sur plusieurs criteres ponderes :

| Critere | Poids | Description |
|---------|-------|-------------|
| CLNSIG | x2 | Classification ClinVar (pathogenic = 5, likely pathogenic = 4) |
| CLNREVSTAT | x1.5 | Niveau de relecture (practice guideline = 5, expert panel = 4, multiple submitters = 3) |
| Zygosite | x1 | Homozygote / hemizygote = 2, heterozygote = 1 |
| Heredite | x1.5 | Concordance entre zygosite et mode de transmission attendu |
| Presence ALT | +2 | Allele alternatif effectivement detecte |

Le score total permet de hierarchiser les resultats par niveau de confiance.

---

## Limites

1. **Duplications / CNV** : La duplication du gene PMP22 (responsable d'environ 70% des cas de CMT1A) n'est pas detectable par ce type d'analyse ponctuelle.

2. **Variants rares** : Seuls les variants repertories dans ClinVar sont analyses. Les variants de novo ou tres rares peuvent echapper a la detection.

3. **Precision des puces** : Les donnees de genotypage MyHeritage peuvent comporter des erreurs ou des non-concordances de position.

4. **Build genomique** : Les coordonnees doivent correspondre au build GRCh37 (hg19). Un decalage de build (GRCh38) peut induire des faux negatifs.

5. **Inversion de brin** : Bien que prise en compte, la verification du brin correct n'est pas garantie a 100% pour tous les variants.

6. **Non diagnostic** : Cet outil est a vocation informative et de recherche. Il ne constitue en aucun cas un diagnostic medical.

---

## Avertissement medical

Les resultats produits par cet outil sont strictement informatifs. Ils ne remplacent en aucun cas une consultation medicale, un avis de geneticien ou un test genetique clinique valide. En cas de detection d'un variant pathogene, il est imperatif de consulter un professionnel de sante qualifie pour confirmation et conseil genetique.

---

## Licence

MIT License

Copyright (c) 2025 Lombard Web Services

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
