require(asdreader)
require(dplyr)
require(tidyr)

# import every spectra file 
setwd("C:/nico")
list_asd = list.files("./FIELDSPEC4_data2",
                      pattern = "\\.asd$",
                      full.names = TRUE)

spec = get_spectra(list_asd, type = "reflectance")

# remove spectralon reading (blank reference)
spec = spec[2:nrow(spec ), ]

# tidy up the spectra names and convert it to data.frame
rownames(spec) = sprintf("spec%03d", 301:(300 + nrow(spec)))
spec = as.data.frame(spec)

# isolate references
soil_spec = spec[61, ]
PP_spec = spec[62, ] 
PVC_spec = spec[63, ]
PET_spec = spec[64, ]
PE_spec = spec[65, ]

# set contaminated soil data
spec = spec[1:60, ]

# load the sample design spreadsheet
design = read.csv2("./Dissertação/SENSNEXUS_data/sample_design2.csv")

# merge datasets
design$spc = spec
design = design[, -1]

# save compiled data
saveRDS(design, "./Dissertação/SENSNEXUS_data/preprocessed_data/datsoil2.rds")

# save references
ref_spec2 = list(
  soil = soil_spec,
  PP   = PP_spec,
  PVC  = PVC_spec,
  PET  = PET_spec,
  PE   = PE_spec
)

saveRDS(ref_spec2, "./Dissertação/SENSNEXUS_data/preprocessed_data/ref_spec2.rds")

 