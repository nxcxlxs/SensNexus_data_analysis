require(asdreader)


# import every spectra file
list_asd = list.files("../raw_spectra/FIELDSPEC4_data2",
                      pattern = "\\.asd$",
                      full.names = TRUE)

spec = get_spectra(list_asd, type = "reflectance")

# remove spectralon reading (blank reference)
spec = spec[2:nrow(spec ), ]

# isolate contaminated soil data (see `raw_spectra/notes [04.08.2025].txt`)
spec = spec[1:60, ]

# tidy up the spectra names and convert it to data.frame
rownames(spec) = sprintf("spec%03d", 289:(288 + nrow(spec)))
spec = as.data.frame(spec)

# load the sample design spreadsheet
design2 = read.csv2("../raw_spectra//colored_sample_design.csv")

# merge datasets
design2$spc = spec
colored = design2

# save compiled data
saveRDS(colored, "../raw_spectra/raw_colored.rds")
