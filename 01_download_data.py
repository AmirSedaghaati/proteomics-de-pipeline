import os
import urllib.request

os.makedirs("data", exist_ok=True)

# Source: Merrihew et al. 2023, Scientific Data, "A peptide-centric quantitative
# proteomics dataset for the phenotypic assessment of Alzheimer's disease"
# ProteomeXchange: PXD034525 | DOI: 10.6069/wefm-vv52 | License: CC BY 4.0
# Hosted publicly on Panorama Public (University of Washington) - no login required.
# Region used: Hippocampus (matches the brain region used in the candidate's own
# Scientific African transcriptomics paper, at the protein level instead of mRNA).

PROTEIN_ABUNDANCE_URL = (
    "https://panoramaweb.org/Panorama%20Public/2022/MacCoss%20-%20Human%20AD%20"
    "Clean%20Diagnosis%20DIA%20Data/Hippocampus/filecontent-sendFile.view?"
    "fileName=Level%203B/Hipp-protGrp.batchadj.csv"
)

# Sample-to-diagnosis metadata, bundled with the paper's own companion R package
# (uw-maccosslab/ADBrainCleanDiagDIA on GitHub), fetched from its raw GitHub content.
METADATA_URL = (
    "https://raw.githubusercontent.com/uw-maccosslab/ADBrainCleanDiagDIA/main/"
    "inst/extdata/HuAD-Clean-Hipp-B1-B3-meta.csv"
)

print("Downloading Hippocampus protein-group abundance table (Level 3B, batch-adjusted)...")
urllib.request.urlretrieve(PROTEIN_ABUNDANCE_URL, "data/Hipp-protGrp.batchadj.csv")
print("Saved to data/Hipp-protGrp.batchadj.csv")

print("Downloading sample metadata (diagnosis group per sample)...")
urllib.request.urlretrieve(METADATA_URL, "data/Hipp-meta.csv")
print("Saved to data/Hipp-meta.csv")

print("\nDone. Both files are real, published data (CC BY 4.0) - not mock/fabricated.")
