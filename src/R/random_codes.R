# ---- CSB195: random_codes.R ----
# Generates random genetic codes, scores them, and compares to Standard.
# Outputs: ./out/random_scores.csv and ./out/random_scores_hist.png

# --- CONFIG ---
set.seed(195)                 # reproducible
N <- 10000                    # number of random codes to simulate
model <- "cover_only" # "degeneracy_preserve" or "cover_only"

# --- Load aaSim() ---
aaSim <- readRDS("./dat/aaSim.4.1.Rds")

# --- Standard Genetic Code (SGC) ---
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
stopifnot(length(std_code) == 64)

# --- Helpers ---
neighbors_1mut <- function(codon) {
  bases <- c("A","C","G","T")
  cs <- strsplit(codon, "", fixed = TRUE)[[1]]
  out <- character(0)
  for (i in 1:3) {
    for (b in setdiff(bases, cs[i])) {
      tmp <- cs
      tmp[i] <- b
      out <- c(out, paste0(tmp, collapse = ""))
    }
  }
  out
}

score_code <- function(code_map) {
  codons <- names(code_map)
  s <- 0
  for (c in codons) {
    aa1 <- code_map[[c]]
    for (n in neighbors_1mut(c)) {
      aa2 <- code_map[[n]]
      s <- s + aaSim(aa1, aa2)
    }
  }
  s
}

# --- Random code generators (two models) ---

# (A) Degeneracy-preserving model:
# Shuffle which amino acid each codon encodes while preserving
# the exact number of codons per amino acid (and number of stops).
gen_code_deg_preserve <- function(sgc) {
  codons <- names(sgc)
  counts <- table(unname(sgc))     # how many codons per AA / stop
  # Build a label vector with repeated AA symbols according to counts
  labels <- unlist(mapply(rep, names(counts), as.integer(counts), SIMPLIFY = FALSE), use.names = FALSE)
  # Shuffle labels and assign to codons
  shuffled <- sample(labels, length(labels), replace = FALSE)
  structure(shuffled, names = codons)
}

# (B) Coverage-only model:
# Only enforce: all 20 amino acids appear at least once, and >=1 stop codon.
# Remaining codons are filled uniformly at random from AA set + stop.
gen_code_cover_only <- function() {
  codons <- names(std_code)
  aa20 <- sort(unique(unname(std_code[std_code != "*"])))  # 20 amino acids
  # Ensure one of each AA and one stop
  base_labels <- c(aa20, "*")  # length 21
  # Fill remaining 43 codons by sampling from 21 symbols
  extra <- sample(c(aa20, "*"), size = 64 - length(base_labels), replace = TRUE)
  labels <- sample(c(base_labels, extra), size = 64, replace = FALSE)
  structure(labels, names = codons)
}

gen_random_code <- function(model = c("degeneracy_preserve","cover_only")) {
  model <- match.arg(model)
  if (model == "degeneracy_preserve") return(gen_code_deg_preserve(std_code))
  else return(gen_code_cover_only())
}

# --- Compute standard score (control) ---
std_score <- score_code(std_code)
cat("Standard code score:", std_score, "\n")  # should be 9856.116

# --- Simulate ---
scores <- numeric(N)
for (i in seq_len(N)) {
  rc <- gen_random_code(model)
  scores[i] <- score_code(rc)
  if (i %% 1000 == 0) cat("Simulated:", i, "of", N, "\n")
}

# --- Summaries ---
cat("\nModel:", model, "\n")
cat("N:", N, "\n")
cat("Mean:", mean(scores), " SD:", sd(scores), "\n")
qs <- quantile(scores, c(0.01, 0.05, 0.5, 0.95, 0.99))
print(qs)
prop_better <- mean(scores >= std_score)  # how often random >= standard
cat("Proportion random >= standard:", prop_better, "\n")

# --- Save outputs ---
# --- Save outputs ---
dir.create("./out", showWarnings = FALSE)

# include model name in filenames
write.csv(data.frame(score = scores),
          paste0("./out/random_scores_", model, ".csv"), row.names = FALSE)

png(paste0("./out/random_scores_hist_", model, ".png"), width = 1200, height = 900, res = 150)
hist(scores, breaks = 60, main = paste("Random code scores (", model, ")", sep = ""),
     xlab = "Score (sum of aaSim over all 1-nt neighbors)")
abline(v = std_score, lwd = 3)   # mark the standard code
dev.off()


cat("\nSaved:\n - ./out/random_scores.csv\n - ./out/random_scores_hist.png\n")
