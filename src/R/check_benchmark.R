# --- Step 1: Load aaSim function ---
aaSim <- readRDS("./dat/aaSim.4.1.Rds")

# --- Step 2: Define the standard genetic code ---
std_code <- c(
  "TTT"="F","TTC"="F","TTA"="L","TTG"="L",
  "TCT"="S","TCC"="S","TCA"="S","TCG"="S",
  "TAT"="Y","TAC"="Y","TAA"="*","TAG"="*",
  "TGT"="C","TGC"="C","TGA"="*","TGG"="W",
  "CTT"="L","CTC"="L","CTA"="L","CTG"="L",
  "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
  "CAT"="H","CAC"="H","CAA"="Q","CAG"="Q",
  "CGT"="R","CGC"="R","CGA"="R","CGG"="R",
  "ATT"="I","ATC"="I","ATA"="I","ATG"="M",
  "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
  "AAT"="N","AAC"="N","AAA"="K","AAG"="K",
  "AGT"="S","AGC"="S","AGA"="R","AGG"="R",
  "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
  "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
  "GAT"="D","GAC"="D","GAA"="E","GAG"="E",
  "GGT"="G","GGC"="G","GGA"="G","GGG"="G"
)

# --- Step 3: Define function to get 1-step neighbors ---
neighbors_1mut <- function(codon) {
  bases <- c("A","C","G","T")
  cs <- strsplit(codon, "")[[1]]
  out <- character(0)
  for (i in 1:3) {
    for (b in setdiff(bases, cs[i])) {
      tmp <- cs
      tmp[i] <- b
      out <- c(out, paste0(tmp, collapse=""))
    }
  }
  out
}

# --- Step 4: Compute benchmark score ---
codons <- names(std_code)
score <- 0

for (c in codons) {
  aa1 <- std_code[[c]]
  for (n in neighbors_1mut(c)) {
    aa2 <- std_code[[n]]
    score <- score + aaSim(aa1, aa2)
  }
}

cat("Benchmark score =", score, "\n")
