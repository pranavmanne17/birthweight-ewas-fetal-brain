"""
Functional classification of the 729 birthweight-associated genes.

Manne P. Maternal Micronutrient Deficiency, DNA Methylation, and Fetal
Neurodevelopment in Low-Resource Settings: A Secondary Analysis of
Growth-Associated EWAS Signatures. (Methods 2.5; Supplementary Table S1)

For each gene, Gene Ontology Biological Process (GO BP) terms are retrieved
from MyGene.info. Each term is lower-cased and searched for the keywords of
four categories. A gene matching several categories goes to the highest-
priority one (growth > metabolic > immune > neurodevelopment); a gene with
no GO BP terms or no keyword match goes to "(v) Other / unclassified".
Genes that match several categories or none are flagged for manual review.

Usage:
    pip install -r requirements.txt
    python classify_genes.py                      # reads data/genes_729.txt
    python classify_genes.py path/to/genes.txt    # any one-symbol-per-line list

Output: results/go_gene_classification_keyword.csv
Note: GO annotations change over time. The analysis in the paper used
annotations retrieved on 21 June 2026, so a rerun today can differ slightly.
"""
import sys
from pathlib import Path

import mygene
import pandas as pd

# Gene symbols that spreadsheet software converts to dates, mapped back to
# their current HGNC symbols. Three such symbols were in the original list.
DATE_FIXES = {
    "2018-03-03 00:00:00": "MARCHF3", "3-Mar": "MARCHF3", "MARCH3": "MARCHF3",
    "2018-03-09 00:00:00": "MARCHF9", "9-Mar": "MARCHF9", "MARCH9": "MARCHF9",
    "2018-09-11 00:00:00": "SEPTIN11", "11-Sep": "SEPTIN11", "SEPT11": "SEPTIN11",
}

# Category keywords in priority order (first match wins). Supplementary Table S1.
GO_KEYWORDS = [
    ("(i) Growth factor & hormone signalling", [
        "growth factor", "insulin-like growth", "insulin receptor signaling",
        "response to insulin", "hormone", "somatotropin", "growth hormone",
        "regulation of growth", "developmental growth", "skeletal system",
        "epidermal growth factor", "fibroblast growth factor", "steroid hormone",
        "thyroid hormone", "androgen", "estrogen"]),
    ("(ii) Metabolic regulation (incl. one-carbon)", [
        "one-carbon", "folic acid", "folate", "methionine", "homocysteine",
        "s-adenosylmethionine", "tetrahydrofolate", "dna methylation", "methylation",
        "transmethylation", "lipid metabolic", "fatty acid", "cholesterol",
        "glucose metabolic", "carbohydrate metabolic", "metabolic process",
        "oxidation-reduction", "cellular respiration", "gluconeogenesis"]),
    ("(iii) Immune & inflammatory (incl. JAK-STAT)", [
        "immune", "inflammat", "cytokine", "interferon", "interleukin", "jak-stat",
        "innate immune", "adaptive immune", "tumor necrosis", "chemokine",
        "leukocyte", "t cell", "b cell", "nf-kappab", "toll-like receptor",
        "antigen", "complement activation"]),
    ("(iv) Neurodevelopment & synaptic function", [
        "neuron", "axon", "synap", "dendrit", "nervous system development",
        "neurogenesis", "neural", "glial", "gliogenesis", "myelin",
        "brain development", "cognition", "learning", "memory",
        "neurotransmitter", "neuromuscular", "forebrain", "cerebral cortex"]),
]
OTHER = "(v) Other / unclassified"


def load_genes(path):
    genes = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    genes = [DATE_FIXES.get(g, g) for g in genes]
    return list(dict.fromkeys(genes))  # de-duplicate, keep order


def classify(terms):
    """Return (assigned category, all matching categories joined by ';')."""
    tl = [t.lower() for t in terms if t]
    matches = [cat for cat, kws in GO_KEYWORDS if any(kw in t for kw in kws for t in tl)]
    return (matches[0] if matches else OTHER), ";".join(matches)


def go_terms(hit):
    go = hit.get("go")
    if not isinstance(go, dict):
        return []
    bp = go.get("BP", [])
    if isinstance(bp, dict):
        bp = [bp]
    return [d.get("term", "") for d in bp]


def main():
    gene_file = sys.argv[1] if len(sys.argv) > 1 else "data/genes_729.txt"
    genes = load_genes(gene_file)
    print(f"Looking up GO BP terms for {len(genes)} genes ...")

    res = mygene.MyGeneInfo().querymany(
        genes, scopes="symbol", fields="go.BP", species="human", returnall=True)

    rows = []
    for hit in res["out"]:
        terms = go_terms(hit)
        cat, allm = classify(terms)
        rows.append({
            "gene": hit.get("query"),
            "n_GO_BP_terms": len(terms),
            "keyword_category": cat,
            "all_matching_categories": allm,
            "needs_review": cat == OTHER or len(allm.split(";")) > 1,
        })
    df = pd.DataFrame(rows).drop_duplicates("gene")

    summary = (df["keyword_category"].value_counts()
               .rename_axis("category").reset_index(name="n_genes"))
    summary["pct"] = (100 * summary["n_genes"] / len(df)).round(1)
    print("\nKeyword-only category counts (before manual review):")
    print(summary.to_string(index=False))
    print(f"\nGenes flagged for manual review: {int(df['needs_review'].sum())}")
    print(f"Genes with no MyGene.info match: {len(res.get('missing', []))}")

    Path("results").mkdir(exist_ok=True)
    out = Path("results/go_gene_classification_keyword.csv")
    df.to_csv(out, index=False)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
