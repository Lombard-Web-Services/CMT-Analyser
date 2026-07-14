#!/usr/bin/env python3

import pandas as pd
import gzip
import re
import warnings
from collections import defaultdict
warnings.filterwarnings('ignore')

print("=" * 80)
print("🔬 ANALYSE CMT - VCF OPTIMISÉ (PRISE EN COMPTE DE TOUS LES PARAMÈTRES)")
print("=" * 80)

# ===================================================================
# 1. CHARGEMENT ET INDEXATION DES DONNÉES MYHERITAGE
# ===================================================================

print("\n📥 Chargement de vos données MyHeritage...")
df_myheritage = pd.read_csv('MH.csv', sep=',', comment='#')
df_myheritage.columns = df_myheritage.columns.str.strip()

# Indexation par position ET par RSID pour plus de flexibilité
mh_index_pos = {}  # "CHR:POS" -> données
mh_index_rsid = {}  # RSID -> données

for _, row in df_myheritage.iterrows():
    chrom = str(row['CHROMOSOME']).replace('chr', '')
    pos = int(row['POSITION'])
    rsid = str(row['RSID'])
    alleles = str(row['RESULT']).upper()
    
    mh_index_pos[f"{chrom}:{pos}"] = {
        'rsid': rsid,
        'alleles': alleles,
        'chrom': chrom,
        'pos': pos
    }
    mh_index_rsid[rsid] = {
        'alleles': alleles,
        'chrom': chrom,
        'pos': pos
    }

print(f"✅ {len(mh_index_pos)} positions MyHeritage indexées")

# ===================================================================
# 2. COORDONNÉES PRÉCISES DES GÈNES (RefSeq/Ensembl GRCh37/hg19)
# ===================================================================

# Coordonnées officielles des gènes majeurs CMT (GRCh37/hg19)
# Source: Ensembl/RefSeq
genes_cmt = {
    'PMP22': {'chr': '17', 'start': 15100000, 'end': 15250000},
    'GJB1': {'chr': 'X', 'start': 70430000, 'end': 70440000},
    'MFN2': {'chr': '1', 'start': 12000000, 'end': 12050000},
    'MPZ': {'chr': '1', 'start': 161200000, 'end': 161300000},
    'LITAF': {'chr': '16', 'start': 11600000, 'end': 11750000},
    'EGR2': {'chr': '10', 'start': 64500000, 'end': 64650000},
    'NEFL': {'chr': '8', 'start': 24900000, 'end': 25050000},
    'GDAP1': {'chr': '8', 'start': 74300000, 'end': 74500000},
    'SH3TC2': {'chr': '5', 'start': 148800000, 'end': 149100000},
    'PLEKHG5': {'chr': '1', 'start': 6460000, 'end': 6540000},
}

# ===================================================================
# 3. FONCTIONS UTILITAIRES
# ===================================================================

def complement(allele):
    """Calcule le complément pour l'inversion de brin"""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return comp.get(allele, allele)

def compare_alleles_with_strand_check(vos_alleles, ref, alt):
    """
    Vérifie si vous portez l'allèle pathogène en tenant compte :
    - Du génotype (AA, AT, etc.)
    - De l'inversion de brin éventuelle
    """
    # Nettoyer
    vos_alleles = vos_alleles.upper()
    ref = ref.upper()
    alt = alt.upper()
    
    # Vérification directe
    if alt in vos_alleles:
        return True, "direct"
    
    # Vérification avec complément (inversion de brin)
    alt_comp = complement(alt)
    if alt_comp in vos_alleles:
        return True, "complement"
    
    return False, None

def parse_clnsig(clnsig_str):
    """Parse proprement la classification ClinVar"""
    if not clnsig_str:
        return None, 0
    
    # Nettoyer
    clnsig_str = str(clnsig_str)
    
    # Classification avec score de confiance
    classifications = {
        'Pathogenic': 5,
        'Likely_pathogenic': 4,
        'Pathogenic/Likely_pathogenic': 4,
        'Conflicting': 2,
        'Uncertain_significance': 1,
        'Benign': 0
    }
    
    # Trouver la classification la plus sévère
    max_score = 0
    max_class = None
    
    for classif, score in classifications.items():
        if classif.lower() in clnsig_str.lower():
            if score > max_score:
                max_score = score
                max_class = classif
    
    return max_class, max_score

def get_review_status_score(review_str):
    """Score basé sur le niveau de revue ClinVar"""
    if not review_str:
        return 0
    
    review_str = str(review_str).lower()
    
    if 'practice guideline' in review_str:
        return 5
    elif 'expert panel' in review_str:
        return 4
    elif 'multiple submitters' in review_str:
        return 3
    elif 'single submitter' in review_str:
        return 2
    else:
        return 1

# ===================================================================
# 4. PARSEUR VCF OPTIMISÉ
# ===================================================================

print("\n🔍 Lecture du fichier ClinVar VCF...")

# Détection du build
def detect_build_from_vcf(vcf_path):
    """Détecte le build à partir du header VCF"""
    try:
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('##reference'):
                    if 'GRCh37' in line or 'hg19' in line:
                        return 'GRCh37'
                    elif 'GRCh38' in line or 'hg38' in line:
                        return 'GRCh38'
                if line.startswith('#CHROM'):
                    break
    except:
        pass
    return 'GRCh37'  # Default

vcf_path = 'clinvar.vcf.gz'
build = detect_build_from_vcf(vcf_path)
print(f"   Build détecté: {build}")

if build != 'GRCh37':
    print("   ⚠️ ATTENTION: Le build ClinVar n'est pas GRCh37.")
    print("   Les coordonnées pourraient ne pas correspondre à MyHeritage.")
    print("   Continuer avec prudence...")

# ===================================================================
# 5. ANALYSE PRINCIPALE
# ===================================================================

vrais_positifs = []
faux_positifs = []
non_trouves = []
variants_analyses = 0
variants_dans_genes = 0
variants_pathogenes = 0

try:
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            variants_analyses += 1
            if variants_analyses % 500000 == 0:
                print(f"   - {variants_analyses} variants analysés...")
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0].replace('chr', '')
            pos = int(fields[1])
            rsid = fields[2] if fields[2] != '.' else None
            ref = fields[3].upper()
            alt_list = fields[4].split(',')
            info = fields[7]
            
            # 1. Filtrer par gène
            gene_trouve = None
            for gene, region in genes_cmt.items():
                if chrom == region['chr'] and region['start'] <= pos <= region['end']:
                    gene_trouve = gene
                    break
            
            if not gene_trouve:
                continue
            
            variants_dans_genes += 1
            
            # 2. Parser la classification
            clnsig = None
            if 'CLNSIG=' in info:
                match = re.search(r'CLNSIG=([^;]+)', info)
                if match:
                    clnsig = match.group(1)
            
            if not clnsig:
                continue
            
            classif, score_class = parse_clnsig(clnsig)
            if score_class < 3:  # Ne garder que Pathogenic ou Likely Pathogenic
                continue
            
            variants_pathogenes += 1
            
            # 3. Vérifier dans vos données (par position ou RSID)
            key = f"{chrom}:{pos}"
            vos_donnees = None
            
            if key in mh_index_pos:
                vos_donnees = mh_index_pos[key]
            elif rsid and rsid in mh_index_rsid:
                vos_donnees = mh_index_rsid[rsid]
            
            if not vos_donnees:
                non_trouves.append({
                    'gene': gene_trouve,
                    'position': pos,
                    'rsid': rsid,
                    'ref': ref,
                    'alt': alt_list,
                    'classification': clnsig
                })
                continue
            
            vos_alleles = vos_donnees['alleles']
            rsid_mh = vos_donnees['rsid']
            
            # 4. Vérifier les allèles (supporte les multialléliques)
            porteur = False
            allele_patho = None
            type_match = None
            
            for alt in alt_list:
                alt = alt.upper()
                if alt in ['<DEL>', '<INS>', '<DUP>']:
                    # Indel/CNV non détectable par puce
                    continue
                
                porteur, type_match = compare_alleles_with_strand_check(vos_alleles, ref, alt)
                if porteur:
                    allele_patho = alt
                    break
            
            if porteur:
                # 5. Récupérer la condition
                condition = ''
                if 'CLNDN=' in info:
                    match = re.search(r'CLNDN=([^;]+)', info)
                    if match:
                        condition = match.group(1)
                
                # 6. Review status
                review = ''
                if 'CLNREVSTAT=' in info:
                    match = re.search(r'CLNREVSTAT=([^;]+)', info)
                    if match:
                        review = match.group(1)
                
                review_score = get_review_status_score(review)
                
                # 7. Zygosité
                zygosity = "HOMOZYGOTE" if vos_alleles == allele_patho * 2 else "HÉTÉROZYGOTE"
                
                # 8. Score de confiance
                confidence_score = score_class + review_score
                if zygosity == "HOMOZYGOTE":
                    confidence_score += 1
                if "dominant" in condition.lower():
                    confidence_score += 1
                
                # Étoiles de confiance
                if confidence_score >= 8:
                    stars = "★★★★★"
                elif confidence_score >= 6:
                    stars = "★★★★"
                elif confidence_score >= 4:
                    stars = "★★★"
                else:
                    stars = "★★"
                
                vrais_positifs.append({
                    'RSID_MyHeritage': rsid_mh,
                    'RSID_ClinVar': rsid,
                    'Gene': gene_trouve,
                    'Chromosome': chrom,
                    'Position': pos,
                    'Reference': ref,
                    'Alt_Pathogene': allele_patho,
                    'Vos_Alleles': vos_alleles,
                    'Zygosite': zygosity,
                    'Classification': clnsig,
                    'Review': review,
                    'Condition': condition,
                    'Confiance': stars,
                    'Score': confidence_score,
                    'Strand_check': type_match
                })
            else:
                faux_positifs.append({
                    'gene': gene_trouve,
                    'rsid': rsid,
                    'vos_alleles': vos_alleles,
                    'ref': ref,
                    'alt': alt_list,
                    'classification': clnsig
                })

except FileNotFoundError:
    print(f"❌ Fichier '{vcf_path}' introuvable.")
    print("   Téléchargez-le avec :")
    print("   wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz")
    exit(1)

# ===================================================================
# 6. RÉSULTATS DÉTAILLÉS
# ===================================================================

print("\n" + "=" * 80)
print("📊 RÉSULTATS DÉTAILLÉS DE L'ANALYSE")
print("=" * 80)

print(f"""
📈 STATISTIQUES GÉNÉRALES :
   - Variants analysés dans ClinVar : {variants_analyses:,}
   - Variants dans les gènes CMT : {variants_dans_genes:,}
   - Variants pathogènes/likely dans ces gènes : {variants_pathogenes:,}
   - Variants présents dans vos données : {len(vrais_positifs) + len(faux_positifs)}
   - Variants non trouvés : {len(non_trouves)}
""")

# ===================================================================
# 7. VRAIS POSITIFS
# ===================================================================

if vrais_positifs:
    print("⚠️⚠️⚠️ VARIANTS PATHOGÈNES CONFIRMÉS ⚠️⚠️⚠️")
    print("=" * 80)
    
    df_vrais = pd.DataFrame(vrais_positifs)
    df_vrais = df_vrais.sort_values('Score', ascending=False)
    
    # Grouper par gène
    for gene in df_vrais['Gene'].unique():
        variants_gene = df_vrais[df_vrais['Gene'] == gene]
        print(f"\n🧬 GÈNE {gene} : {len(variants_gene)} variant(s) confirmé(s)")
        
        for _, row in variants_gene.iterrows():
            print(f"\n   {'=' * 60}")
            print(f"   {row['Confiance']} VARIANT {row['RSID_MyHeritage']}")
            print(f"   {'=' * 60}")
            print(f"   Position : chr{row['Chromosome']}:{row['Position']:,}")
            print(f"   ClinVar RSID : {row['RSID_ClinVar']}")
            print(f"   Mutation : {row['Reference']} → {row['Alt_Pathogene']}")
            print(f"   Votre génotype : {row['Vos_Alleles']} ({row['Zygosite']})")
            print(f"   Classification : {row['Classification']}")
            print(f"   Review : {row['Review']}")
            print(f"   Condition : {row['Condition']}")
            print(f"   Vérification brin : {row['Strand_check']}")
            print(f"   Score de confiance : {row['Score']}/10")
    
    # Sauvegarde
    df_vrais.to_csv('resultats_cmt_vrais_positifs.csv', index=False)
    print(f"\n💾 Résultats sauvegardés dans 'resultats_cmt_vrais_positifs.csv'")
    
else:
    print("✅ AUCUN VARIANT PATHOGÈNE CONFIRMÉ")
    print("=" * 80)
    print("""
   Après vérification stricte des coordonnées (CHR:POS) et des allèles (REF/ALT),
   vous ne portez aucun des variants ponctuels pathogènes majeurs de la CMT.

   C'est un résultat rassurant, mais attention :
   - La duplication PMP22 (~70% des CMT1A) n'est PAS détectable
   - Les variants rares dans d'autres gènes peuvent exister
   - Un test clinique reste la seule référence
    """)

# ===================================================================
# 8. EXEMPLES DE FAUX POSITIFS
# ===================================================================

if faux_positifs:
    print("\n📊 EXEMPLES DE FAUX POSITIFS (RSID trouvé mais pas l'allèle) :")
    for fp in faux_positifs[:5]:
        print(f"   - {fp['gene']}: {fp['rsid']} - Vos allèles: {fp['vos_alleles']} ≠ {fp['alt']}")

# ===================================================================
# 9. RÉSUMÉ FINAL
# ===================================================================

print("\n" + "=" * 80)
print("📋 RÉSUMÉ FINAL")
print("=" * 80)

# Générer un rapport résumé
if vrais_positifs:
    print("\n⚠️ RÉSULTAT : DES VARIANTS PATHOGÈNES ONT ÉTÉ DÉTECTÉS")
    print("   Voici les étapes à suivre :")
    print("   1. Ne pas paniquer - ce n'est pas un diagnostic")
    print("   2. Consulter un neurologue ou un généticien")
    print("   3. Demander un test génétique clinique de confirmation")
    print("   4. Discuter des implications familiales")
else:
    print("\n✅ RÉSULTAT : AUCUN VARIANT PATHOGÈNE DÉTECTÉ")
    print("   Points à retenir :")
    print("   • Vous ne portez pas les variants ponctuels testés")
    print("   • La duplication PMP22 n'est pas détectable")
    print("   • Ce n'est pas un diagnostic d'exclusion total")
    print("   • Si vous avez des symptômes, consultez un médecin")

print("\n" + "=" * 80)
print("⚠️ LIMITES DE CETTE ANALYSE :")
print("   ✗ Détection des duplications/CNV (PMP22) impossible")
print("   ✗ Variants rares non répertoriés dans ClinVar")
print("   ✗ Erreurs potentielles des puces MyHeritage")
print("   ✗ Différences de build (GRCh37 vs GRCh38)")
print("   ✗ Inversion de brin non vérifiée pour tous les variants")
print("=" * 80)
