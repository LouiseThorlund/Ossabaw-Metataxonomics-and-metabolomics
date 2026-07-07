library(readxl)
library(stringr)
library(readxl)
library(stringr)

DP3 <- as.data.frame(DP3)

# ----------------------------
# Functions
# ----------------------------
get_metabolites <- function(sample) {
  vec <- DP3[sample, , drop = TRUE]
  names(vec)[vec > 0 & !is.na(vec)]
}

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

# ----------------------------
# Setup
# ----------------------------
samples <- rownames(DP3)
pigs <- unique(str_extract(samples, "G\\d+"))

results <- data.frame()

# ----------------------------
# Loop per pig
# ----------------------------
for (pig in pigs) {
  
  pig_samples <- samples[str_detect(samples, pig)]
  
  H1 <- pig_samples[str_detect(pig_samples, "HR1")]
  H2 <- pig_samples[str_detect(pig_samples, "HR2")]
  L1 <- pig_samples[str_detect(pig_samples, "LR1")]
  L2 <- pig_samples[str_detect(pig_samples, "LR2")]
  
  # skip incomplete pigs
  if (length(H1) == 0 || length(H2) == 0 ||
      length(L1) == 0 || length(L2) == 0) next
  
  # metabolite sets
  H1_met <- get_metabolites(HR1)
  H2_met <- get_metabolites(HR2)
  L1_met <- get_metabolites(LR1)
  L2_met <- get_metabolites(LR2)
  
  # store results
  results <- rbind(
    results,
    data.frame(
      Pig = pig,
      
      HR1 = HR1,
      HR2 = HR2,
      LR1 = LR1,
      LR2 = LR2,
      
      H_shared = length(intersect(H1_met, H2_met)),
      H_total  = length(union(H1_met, H2_met)),
      H_jaccard = jaccard(H1_met, H2_met),
      
      L_shared = length(intersect(L1_met, L2_met)),
      L_total  = length(union(L1_met, L2_met)),
      L_jaccard = jaccard(L1_met, L2_met)
    )
  )
}

# ----------------------------
# Output
# ----------------------------
results


#RSD
DP3 <- as.data.frame(DP3)

samples <- rownames(DP3)
pigs <- unique(stringr::str_extract(samples, "G\\d+"))

get_value <- function(sample) {
  sum(DP3[sample, , drop = TRUE] > 0, na.rm = TRUE)
}

rsd_fun <- function(x, y) {
  m <- mean(c(x, y))
  sd(c(x, y)) / m * 100
}

results <- data.frame()

for (pig in pigs) {
  
  pig_samples <- samples[stringr::str_detect(samples, pig)]
  
  HR1 <- pig_samples[stringr::str_detect(pig_samples, "HR1")]
  HR2 <- pig_samples[stringr::str_detect(pig_samples, "HR2")]
  LR1 <- pig_samples[stringr::str_detect(pig_samples, "LR1")]
  LR2 <- pig_samples[stringr::str_detect(pig_samples, "LR2")]
  
  if (length(HR1)==0 || length(HR2)==0 || length(LR1)==0 || length(LR2)==0) next
  
  results <- rbind(
    results,
    data.frame(
      Pig = pig,
      H_RSD = rsd_fun(get_value(HR1), get_value(HR2)),
      L_RSD = rsd_fun(get_value(LR1), get_value(LR2))
    )
  )
}

results


results <- data.frame()

for (pig in pigs) {
  
  pig_samples <- samples[stringr::str_detect(samples, pig)]
  
  HR1 <- pig_samples[stringr::str_detect(pig_samples, "HR1")]
  HR2 <- pig_samples[stringr::str_detect(pig_samples, "HR2")]
  LR1 <- pig_samples[stringr::str_detect(pig_samples, "LR1")]
  LR2 <- pig_samples[stringr::str_detect(pig_samples, "LR2")]
  
  h_rsd <- NA
  l_rsd <- NA
  
  if (length(HR1) > 0 & length(HR2) > 0) {
    h_rsd <- rsd_fun(get_value(HR1), get_value(HR2))
  }
  
  if (length(LR1) > 0 & length(LR2) > 0) {
    l_rsd <- rsd_fun(get_value(LR1), get_value(LR2))
  }
  
  results <- rbind(
    results,
    data.frame(
      Pig = pig,
      H_RSD = h_rsd,
      L_RSD = l_rsd
    )
  )
}

results
