#-------------------------------------------------------------------------------
# The following code is an example of how dplyr can be used to interface with
# the SQLite database.
# Please see the schema_db folder for description of column headers
#-------------------------------------------------------------------------------
stopifnot(require("dplyr"))  # requisite library
db <- src_sqlite("cage.sqlite", create = FALSE)
# View data base description and tables available
db
# src:  sqlite 3.8.6 [cage.sqlite]
# tbls: associations_1e6, eqtl_cojo, gemma_top_snps,
#  gene_transcription_start_stop, heritability, probe_genome_location,
#  westra_cis, westra_trans
#-------------------------------------------------------------------------------
# tbl() grabs the table from the database without loading it into R
# this is fast and cheap!
tbl.assoc <- tbl(db, "associations_1e6", n = Inf)  # get a table
head(tbl.assoc)
str(tbl.assoc)
# List of 2
 # $ src:List of 2
  # ..$ con :Formal class 'SQLiteConnection' [package "RSQLite"] with 5 slots
  # .. .. ..@ Id                 :<externalptr> 
  # .. .. ..@ dbname             : chr "cage.sqlite"
  # .. .. ..@ loadable.extensions: logi TRUE
  # .. .. ..@ flags              : int 6
  # .. .. ..@ vfs                : chr ""
  # ..$ path: chr "cage.sqlite"
  # ..- attr(*, "class")= chr [1:3] "src_sqlite" "src_sql" "src"
 # $ ops:List of 3
  # ..$ src :List of 2
  # .. ..$ con :Formal class 'SQLiteConnection' [package "RSQLite"] with 5 slots
  # .. .. .. ..@ Id                 :<externalptr> 
  # .. .. .. ..@ dbname             : chr "cage.sqlite"
  # .. .. .. ..@ loadable.extensions: logi TRUE
  # .. .. .. ..@ flags              : int 6
  # .. .. .. ..@ vfs                : chr ""
  # .. ..$ path: chr "cage.sqlite"
  # .. ..- attr(*, "class")= chr [1:3] "src_sqlite" "src_sql" "src"
  # ..$ x   :Classes 'ident', 'sql', 'character'  chr "associations_1e6"
  # ..$ vars: chr [1:13] "SNP" "CHR" "BP" "GEN_POS" ...
  # ..- attr(*, "class")= chr [1:3] "op_base_remote" "op_base" "op"
 # - attr(*, "class")= chr [1:4] "tbl_sqlite" "tbl_sql" "tbl_lazy" "tbl"
subset <- select(tbl.assoc, c(SNP, PROBE, BETA, SE))  # select() some columns
res    <- dplyr::filter(subset, BETA <= -0.5)  # filter() using logical predicates
# collect() brings data back from the database and puts it into an R data.frame
res <- collect(res)
# To collect all you need to use n = Inf. You can collect everything into
# a data frame all in one go if you prefer
res <- as.data.frame(collect(tbl(db, "associations_1e6"), n = Inf))
# The pipe operator in dplyr is very powerful for subsetting the data 
# For example get the sentinel SNP for each probe
sent <- res %>% group_by(PROBE) %>% arrange(P_BOLT_LMM) %>% slice(1)
# ------------------------------------------------------------------------------
# The illuminaHumanv4.db package from Bioconductor is an excellent resource for 
# annotating the results
# Install if necessary
# ------------------------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("illuminaHumanv4.db")
library(illuminaHumanv4.db)
library(DBI)
anno.f <- illuminaHumanv4fullReannotation()
# Pull in the heritability table
h2      <- as.data.frame(collect(tbl(db, "heritability")))
# ------------------------------------------------------------------------------
# Null results are stored as -9 so that they are numeric in the sqlite data base
# Set all -9s to NAs 
# ------------------------------------------------------------------------------
h2[(h2==-9)] <- NA
head(h2)
# ------------------------------------------------------------------------------
# COJO eQTL 
# ------------------------------------------------------------------------------
eqtl <- as.data.frame(collect(tbl(db, "eqtl_cojo"), n = Inf))
head(eqtl)
# -------------------------------------------------------------------------------
# Take it down to only good probes
# -------------------------------------------------------------------------------
cojo.db <- dplyr::filter(eqtl, PROBE_QUALITY != "Bad" & PROBE_QUALITY != "NA")
gen.loc  <- as.data.frame(collect(tbl(db, "probe_genome_location")))
# Inner join results with genomic location
cojo.gl <- inner_join(cojo.db, gen.loc , by = c("PROBE"))
dim(cojo.gl) # 14995  35 -- 17608 - 14995  -- 2613 
length(unique(cojo.gl$PROBE))
# -------------------------------------------------------------------------------
# Process the probes with genomic location
# -------------------------------------------------------------------------------
gen.loc  <- as.data.frame(collect(tbl(db, "probe_genome_location")))
# Inner join results with genomic location
cojo.gl <- inner_join(cojo.db, gen.loc , by = c("PROBE"))
dim(cojo.gl) # 14995    35
# Paste a chr in front of the eQTL results in order to match by chromosome first
cojo.gl$CHR_CAGE <- paste0("chr", cojo.gl$CHR.x)
head(cojo.gl$SymbolReannotated)
length(unique(cojo.gl$PROBE)) # cojo.gl 9967
# -------------------------------------------------------------------------------
# Take down to chromosomes and put haploype chromosomes as being on
# the chromosome they are mapped to. Unplaced contigs stay as Un 
# -------------------------------------------------------------------------------
cojo.gl$CHR_SUB <-  sub("_\\S*", "", cojo.gl$CHR.y)
# Which values are on the same chromosome cis classification
cojo.gl.cis     <- dplyr::filter(cojo.gl, CHR_SUB == CHR_CAGE)
dim(cojo.gl.cis) #  11204    37
length(unique(cojo.gl.cis$PROBE)) # 7857
cojo.gl.trans <- dplyr::filter(cojo.gl, CHR_SUB != CHR_CAGE)
dim(cojo.gl.trans) # 3791   37
length(unique(cojo.gl.trans$PROBE)) # 3338





















