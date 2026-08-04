require(dplyr)


# define experimental factors
polymers = c("ALL")
sizes_mm = c("4.80–1.00", "1.00–0.60", "0.60–0.053")
masses_mg = c("3750", "750", "150", "30", "6")
replicates = 1:4

# create full factorial design
design2 = expand.grid(
  POLYMER = polymers,
  SIZE_INTERVALS_mm = sizes_mm,
  MASS_mg = masses_mg,
  REPLICATE = replicates,
  stringsAsFactors = FALSE
)

# map size codes
size_code_map = c("4.80–1.00" = "A", "1.00–0.60" = "B", "0.60–0.053" = "C")

design2 = design2 |>
  mutate(
    SIZE_CODE = size_code_map[as.character(SIZE_INTERVALS_mm)],
    SAMPLE_CODE = paste0("MIX_", SIZE_CODE, MASS_mg, "_", REPLICATE)
  )

# define sort order for sample_ID assignment
design2 = design2 |> 
  arrange(SIZE_CODE, desc(as.numeric(as.character(MASS_mg))), REPLICATE)

design2$SAMPLE_ID = 289:348

# reorder columns
design2 = design2 |>
  select(SAMPLE_ID, SAMPLE_CODE, POLYMER, SIZE_INTERVALS_mm,
         SIZE_CODE, MASS_mg, REPLICATE)

# save spreadsheet
write.csv2(design2, "../raw_spectra/colored_sample_design.csv", row.names = F)
