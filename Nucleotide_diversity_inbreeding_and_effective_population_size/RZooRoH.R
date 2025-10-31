# Script: RZooRoH analysis
# Description: This script is launched by 5_launcher_RZooRoH_analyses.sh

library(RZooRoH)

args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]        # Population or pool name
base_name <- args[2]  # Base dataset name (e.g., "masked_HWE_10covperpool")

# Load genotype likelihood data
genofile_PL <- paste0("GL/", pop, "_", base_name, "_GL2PL.txt")
raw_data <- read.table(genofile_PL, header = FALSE)


# Load allele frequency data
freq_file <- paste0("freq_alleliq/frequence_", pop, "_", base_name, ".txt")
freq_data <- read.table(freq_file, header = FALSE)

# Convert physical to genetic coordinates
# Assumes recombination rate = 4 cM/Mb (i.e., multiply by factor of 4)
raw_data[, 2] <- raw_data[, 2] * 4

# Write to a temporary file for RZooRoH input
temp_file <- tempfile(fileext = ".txt")
write.table(raw_data, file = temp_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Prepare data for RZooRoH
ROH_data <- zoodata(
  genofile = temp_file,
  zformat = "gl",        
  chrcol = 1,           
  poscol = 2,            
  supcol = 2,            
  allelefreq = freq_data$V1
)

# Define RZooRoH model
model <- zoomodel(K = 14, base_rate = 2, layers = FALSE)

# Run RZooRoH for the current pool
cat(paste0("Running RZooRoH for population: ", pop, "\n"))
result <- zoorun(model, ROH_data)

# Save results
output_file <- paste0("results/zoorun_result_", pop, ".RData")
save(result, file = output_file)
cat(paste0("Results saved to: ", output_file, "\n"))

# Clean up temporary files
unlink(temp_file)

cat("RZooRoH analysis completed successfully.\n")
