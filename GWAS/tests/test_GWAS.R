options(import.path="/homes/hannah/projects/GWAS")
ldfilter <- modules::import("LDfilter")

context('Read sample data')
sig <- as.matrix(read.table("../sample_data/sigChr2.txt", header=FALSE,
                  stringsAsFactors=FALSE))
ld <- read.table("../sample_data/sigChr2_ld.txt", header=FALSE,
                  stringsAsFactors=FALSE)
colnames(ld) <- c("SNP","CHR","BP", "NTAG", "LEFT", "RIGHT", "KBSPAN", "TAGS")

context('Test LDfilter functions')
testthat::test_that('Output of findLD contains all input SNPs', {
                        snp2ld <- lapply(sig, ldfilter$findLD,  LD=ld)
                        snp2ld <- do.call(rbind, snp2ld)
                        testthat::expect_true(all(snp2ld$SNP == sig))
                  })
testthat::test_that('findLD fails with column error', {
                        testthat::expect_error(ldfilter$findLD(sig[1], LD=ld,
                                                               snp_name="rsID"),
                                               )
                  })
