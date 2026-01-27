require(dplyr)
require(writexl)

# define experimental factors
polymers = c("PP", "PVC", "PET", "PE")
sizes_mm = c("4.80–1.00", "1.00–0.60", "0.60–0.053")
masses_mg = c("3750", "750", "150", "30", "6", "1.2")
replicates = 1:4

# create full factorial design
design = expand.grid(
  POLYMER = polymers,
  SIZE_INTERVALS_mm = sizes_mm,
  MASS_mg = masses_mg,
  REPLICATE = replicates,
  stringsAsFactors = FALSE
)

# identify prepared combinations (PP and PVC, two size intervals, all masses)
prepared_combos = expand.grid(
  POLYMER = c("PP", "PVC"),
  SIZE_INTERVALS_mm = sizes_mm[1:2],  # "4.80–1.00 mm" and "1.00–0.60 mm"
  MASS_mg = masses_mg,
  stringsAsFactors = FALSE
) |> 
  mutate(PREPARED = TRUE)

# mark prepared combinations in design
design = design |>
  left_join(prepared_combos, by = c("POLYMER",
                                    "SIZE_INTERVALS_mm",
                                    "MASS_mg")) |>
  mutate(PREPARED = if_else(is.na(PREPARED), FALSE, TRUE))

# define sort order for Sample_ID assignment
design = design |>
  mutate(
    POLYMER = factor(POLYMER, levels = polymers),
    SIZE_INTERVALS_mm = factor(SIZE_INTERVALS_mm, levels = sizes_mm),
    MASS_mg = factor(MASS_mg, levels = masses_mg)
  )

# assign sample_IDs: prepared samples first, ordered by polymer > size > mass > replicate
prepared_samples = design |> 
  filter(PREPARED) |>
  arrange(POLYMER, SIZE_INTERVALS_mm, MASS_mg, REPLICATE) |>
  mutate(SAMPLE_ID = row_number())

unprepared_samples = design |>
  filter(!PREPARED) |>
  arrange(POLYMER, SIZE_INTERVALS_mm, MASS_mg, REPLICATE) |>
  mutate(SAMPLE_ID = max(prepared_samples$SAMPLE_ID) + row_number())

# recombine full dataset
design = bind_rows(prepared_samples, unprepared_samples) |>
  arrange(SAMPLE_ID)

# map size codes
size_code_map = c("4.80–1.00" = "A", "1.00–0.60" = "B", "0.60–0.053" = "C")
design = design |>
  mutate(
    SIZE_CODE = size_code_map[as.character(SIZE_INTERVALS_mm)],
    SAMPLE_CODE = paste0(POLYMER, "_", SIZE_CODE, MASS_mg, "_", REPLICATE)
  )

# reorder columns
design = design |>
  select(SAMPLE_ID, SAMPLE_CODE, POLYMER, SIZE_INTERVALS_mm,
         SIZE_CODE, MASS_mg, REPLICATE, PREPARED)


# save spreadsheet
setwd("C:/Users/nicolas/Documents/Dissertação/SENSNEXUS_data")
write.csv2(design, "sample_design.csv")