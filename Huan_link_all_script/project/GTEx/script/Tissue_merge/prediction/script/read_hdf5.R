library(rhdf5)


interactions  <- rhdf5::h5read("ENCFF497EDU.h5", "interactions")
bin_p <- rhdf5::h5read("ENCFF497EDU.h5", "bin_positions")

aaa <-t(bin_p)
chr_bin_range <- rhdf5::h5read("ENCFF497EDU.h5", "chr_bin_range")
chrs <- rhdf5::h5read("ENCFF497EDU.h5", "chrs")




library(rhdf5)

file <-


H5Fget_name("ENCFF497EDU.h5")
H5Rget_name("ENCFF497EDU.h5")

H5Fopen("ENCFF497EDU.h5")

library(hdf5r)
library(Seurat) 
Read10X_h5("ENCFF497EDU.h5")

h5file <- system.file("testfiles", "ENCFF497EDU.h5", package = "rhdf5")
fid <- H5Fopen(h5file)
H5Fget_name(fid)
