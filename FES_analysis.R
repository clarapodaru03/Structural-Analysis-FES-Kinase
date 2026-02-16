# Structural Characterization of Fes/Fps Tyrosine Kinase
# Main goal: Extracting structural features from AphaFold2 models

# 1. Load Bio3D library
if (!require(bio3d)) install.packages("bio3d")
library(bio3d)

# 2. Path Protein Structure
# Path to the PDB file
pdb_file <- "AF-P07332-F1-model_v6.pdb"
pdb <- read.pdb(pdb_file)

# 3. Sequence Extraction and Quality Assessment
# Export sequence for multiple sequence alignment (MSA)
seq <- pdbseq(pdb)
write.fasta(seqs = seq, ids = "FES_human", file = "FES_sequence_native.fasta")

# Calculate pLDDT (confidence score stored in B factor column)
cat("Average pLDDT (AI Confidence):", mean(pdb$atom$b), "\n")
# Average pLDDT (AI Confidence): 89.62793ç

# 4. Structural Descriptors (HPC-ready metrics)

# Calculate Radius of Gyration (protein compactness)
# 1. Extract XYZ coordinates
coords <- matrix(pdb$xyz, ncol=3, byrow=TRUE)

# 2. Calculate the Geometric Center (Center of Mass assuming equal weights)
center <- colMeans(coords)

# 3. Calculate squared distances from each atom to the center
sq_dist <- rowSums(sweep(coords, 2, center)^2)

# 4. Rg is the square root of the mean of squared distances
rg <- sqrt(mean(sq_dist))

cat("Radius of Gyration (Manual Calc):", round(rg, 2), "Å\n")
#Radius of Gyration (Manual Calc): 51.08 Å

# 5. Residue Flexibility Analysis

# We select only C-Alpha atoms to ensure a 1:1 mapping between data points and residues
ca_indices <- atom.select(pdb, "calpha")
# In AlphaFold models, the B-factor column contains the pLDDT (Confidence score)
plddt_scores <- pdb$atom$b[ca_indices$atom]

# Saving the professional plot
png("FES_Structural_Confidence.png", width=1000, height=700, res=120)

# plot.bio3d is used to visualize the pLDDT profile across the residue index
plot.bio3d(plddt_scores, 
           resno = pdb$atom$resno[ca_indices$atom], 
           typ = "l", 
           lwd = 2, 
           col = "darkblue",
           ylab = "pLDDT Score (AI Confidence)", 
           xlab = "Residue Number", 
           main = "FES Structural Confidence Profile (AlphaFold2)")

# Adding a threshold line for high-confidence regions (pLDDT > 70)
abline(h=70, col="darkred", lty=2, lwd=1.5)
legend("bottomright", legend=c("pLDDT", "Confidence Threshold (70)"), 
       col=c("darkblue", "darkred"), lty=c(1,2), bty="n")
dev.off()

# 6. Export Processed Data
# Save C-alpha only model for simplified molecular dynamics or docking setup
ca_pdb <- trim.pdb(pdb, ca_indices)
write.pdb(ca_pdb, file="FES_CAlpha_Model.pdb")

