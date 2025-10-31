# Script for convert GL in PL for RZooRoH

# Function conversion GL into PL
convert_gl_to_pl <- function(gl) {
  if (gl == 0) {
    return(0)  # Si GL est 0, PL doit aussi être 0
  } else {
    return(round(-10 * gl / log(10), 0))  # Conversion GL à PL
  }
}

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
  
# Read GL data
input_file <- paste0("GL/",file,"_GL.txt")
df <- read.table(input_file)
  
# Convert GL into PL 
df_pl <- df
df_pl[, 3:ncol(df)] <- lapply(df[, 3:ncol(df)], function(x) sapply(x, convert_gl_to_pl))
  
# Write results
output_file <- paste0("GL/",file,"_GL2PL.txt")
write.table(df_pl, output_file, row.names = FALSE, col.names = FALSE, quote = F)
