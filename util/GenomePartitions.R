
library(package = "data.table")
library(package = "HDInterval")

ReadGenomePartitionsSamples = function(Files, nSamp) {
  # Read genome partitions samples
  # Files - character, files with genome partitions samples
  # nSamp - numeric, number of samples to read
  nPart = length(Files)
  if (nPart < 2) {
    stop("ERROR: need at least two genome partitions!")
  }
  nInd = length(scan(file = Files[1], what = numeric(), nlines = 1, quiet = TRUE))
  Ret = array(dim = c(nInd, nPart, nSamp))
  for (Part in 1:nPart) {
    # Part = 1
    cat("Reading partition: ", Part, "\n", sep = "")

    # Via base::scan()
    # To understand this code note that:
    #   - scan reads line by line, in this case one sample of nInd values
    #   - objects in R are filled by inner dimension first, in this case a scanned
    #     line of nInd values goes into Ind column
    # Ret[, Part, ] = scan(file = Files[Part], what = numeric(), nlines = nSamp, n = nSamp * nInd, quiet = TRUE)

    # Via data.table::fread() to cut down memory usage
    Tmp = fread(file = Files[Part], nrows = nSamp, header = FALSE)
    for (Ind in 1:nInd) {
      # Ind = 1
      Ret[Ind, Part, ] = Tmp[[Ind]]
    }
    rm(Tmp); gc(verbose = FALSE)

  }
  Ret
}

CalculateCovForMultiVariableSamples = function(x) {
  # Calculate covariance matrix samples from multiple variable samples
  # x - array, with dimmensions (n, nVar, nSamp)
  Tmp = dim(x)
  nVar = Tmp[2]
  nSamp = Tmp[3]
  Ret = array(dim = c(nVar, nVar, nSamp))
  for (Samp in 1:nSamp) {
    # Samp = 1
    Ret[, , Samp] = cov(x[, , Samp])
  }
  Ret
}

Cov2CorSamples = function(x) {
  # Convert covariance matrix samples to correlation matrix samples
  # x - array, with dimmensions (nVar, nVar, nSamp)
  Tmp = dim(x)
  Ret = array(dim = Tmp)
  nSamp = Tmp[3]
  for (Samp in 1:nSamp) {
    # Samp = 1
    Ret[, , Samp] = cov2cor(x[, , Samp])
  }
  Ret
}

PartVarVsTotalVarSamples = function(x) {
  # Express partitioned variance (as a matrix) relative to total variance on samples
  # x - array, of dimension (nVar, nVar, nSamp)
  Tmp = dim(x)
  Ret = array(dim = Tmp)
  nSamp = Tmp[3]
  for (Samp in 1:nSamp) {
    # Samp = 1
    Ret[, , Samp] = x[, , Samp] / sum(x[, , Samp])
  }
  Ret
}

SummarizeMatrixSamples = function(x) {
  # Summarize covariances samples
  # x - array, with dimmensions (nVar, nVar, nSamp)
  Ret = vector(mode = "list", length = 4)
  names(Ret) = c("Mean", "Sd", "HpdL", "HpdU")
  nVar = dim(x)[1]
  Ret$Mean = Ret$Sd = Ret$HpdL = Ret$HpdU = array(dim = c(nVar, nVar))
  for (Var1 in 1:nVar) {
    # Var1 = 1
    for (Var2 in Var1:nVar) {
      # Var2 = Var1
      Ret$Mean[Var1, Var2] = mean(x[Var1, Var2, ])
      Ret$Sd[Var1, Var2]   = sd(x[Var1, Var2, ])
      Tmp = hdi(x[Var1, Var2, ])
      Ret$HpdL[Var1, Var2] = Tmp["lower"]
      Ret$HpdU[Var1, Var2] = Tmp["upper"]
    }
  }
  LowTri = lower.tri(Ret$Mean)
  UppTri = upper.tri(Ret$Mean)
  Ret$Mean[LowTri] = Ret$Mean[UppTri]
  Ret$Sd[LowTri]   = Ret$Sd[UppTri]
  Ret$HpdL[LowTri] = Ret$HpdL[UppTri]
  Ret$HpdU[LowTri] = Ret$HpdU[UppTri]
  Ret
}

# setwd("~/Gitbox/AlphaSuite/AlphaAnalyse/example/")
GenoPartSamp = ReadGenomePartitionsSamples(Files = c("GenomePartition1Samples.txt", "GenomePartition2Samples.txt", "GenomePartitionRemainderSamples.txt"),
                                           nSamp = 900)
CovGenoPartSamp = CalculateCovForMultiVariableSamples(x = GenoPartSamp)
CorGenoPartSamp = Cov2CorSamples(x = CovGenoPartSamp)
CovPropGenoPartSamp = PartVarVsTotalVarSamples(x = CovGenoPartSamp)

(CovGenoPart = SummarizeMatrixSamples(x = CovGenoPartSamp))
(CorGenoPart = SummarizeMatrixSamples(x = CorGenoPartSamp))
(CoVPropGenoPart = SummarizeMatrixSamples(x = CovPropGenoPartSamp))
