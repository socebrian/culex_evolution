## Sonia cebrian camison
## proyecto gen families
#21 octubre 2025 / modified 29/12/25 for round 2 of bitacora

# ============================================================
# Script: Build curation summary tables for GR candidates
# - Part 1 (efor BITACORA round1): merge BITACORA gene_table + TOPCONS
# - Part 2 (new for Round2): keep ONLY genes that are new in R2 vs R1,
#   then merge those new genes with TOPCONS (R2)
#
# Outputs (TSV + XLSX), with spp + species included in:
#   - filename
#   - columns inside the output table
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(openxlsx)

# =========================
# 0) Directories (WORKING FROM GR/)
# =========================
# TOPCONS directories
dir_topcons <- "TOPCONS/1st_round"   # if you don't have R1 TOPCONS, keep but it will be skipped
dir_topcons_r2 <- "TOPCONS/2nd_round"

# BITACORA gene_table directories
dir_classif <- "bitacora-primeros-filtros/1st_round"
dir_bitacora_r2 <- "bitacora-primeros-filtros/2nd_round"

# Output directory
dir_out <- "curado/1st_round"
dir_out_r2 <- "curado/2nd_round"

# =========================
# 0.5) File naming conventions (EDIT IF NEEDED)
# =========================
# Given spp/base_name like "1pip", we will look for:
#   R1 gene_table: <spp>_GR_gene_table.tsv
#   R2 gene_table: <spp>_GR2_gene_table.tsv
#
# And TOPCONS tables ending with:  TOPCONS.table.tsv. GR or GR2 depending on the round


#========================
# 1) FIRST ROUND OF BITACORA
#========================
# Listar archivos TOPCONS
top_files <- list.files(dir_topcons, pattern = "\\.table\\.tsv$", full.names = TRUE)
class_files <- list.files(dir_classif, pattern = "_gene_table\\.tsv$", full.names = TRUE)

# Procesar cada par de archivos
for (top_path in top_files) {
  
  # Obtener el identificador (ej: 1pip)
  base_name <- sub("GR_TOPCONS\\.table\\.tsv$", "", basename(top_path))
  
  # Buscar el classif correspondiente
  class_path <- file.path(dir_classif, paste0(base_name, "GR_gene_table.tsv"))
  
  if (!file.exists(class_path)) {
    message("⚠️ No se encontró el archivo de classif para ", base_name)
    next
  }
  
  # Leer tablas
  topcons <- read.table(top_path, sep = "\t", header = TRUE)
  classif <- read.table(class_path, sep = "\t", header = TRUE)
  
  # Crear columna 80%
  classif$`80%` <- ifelse(classif$Length > 296,
                          "Complete",
                          classif$Status.Complete.Partial.Pseudogene)
  
  # Merge
  tabla <- merge(classif, topcons, by.x = "Gene.name", by.y = "GeneID", all.x = TRUE)
  
  # Extraer especie (parte del nombre después del número, ej: 1pip -> pipiens)
  species <- gsub("^[0-9]+", "", base_name)
  tabla$species <- species
  
  # Exportar resultado
  # 1. TSV (texto plano, compatible con R y Excel)
  out_file <- file.path(dir_out, paste0(base_name, "_GR_curado.tsv"))
  write.table(tabla, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 2. Excel real (.xlsx)
  out_file2 <- file.path(dir_out, paste0(base_name, "_GR_curado.xlsx"))
  write.xlsx(tabla, out_file2, rowNames = FALSE)
  
  message("✅ Procesado: ", base_name)
}


# =========================
# 2) SECOND ROUND OF BITACORA (NEW-ONLY vs ROUND1)
# =========================
# Goal:
#   - For each TOPCONS R2 table:
#       (a) read BITACORA R2 gene_table
#       (b) read BITACORA R1 gene_table (same spp) to define "already known"
#       (c) keep only genes that are in R2 but NOT in R1
#       (d) merge those NEW genes with TOPCONS R2, like in Round1
#
# Outputs per spp:
#   curado/2nd_round/<base_name>_OR_curado_NEWONLY.tsv
#   curado/2nd_round/<base_name>_OR_curado_NEWONLY.xlsx

top_files_r2 <- list.files(dir_topcons_r2, pattern = "\\.table\\.tsv$", full.names = TRUE)

for (top_path_r2 in top_files_r2) {
  #infer basename (spp)
    base_name_r2 <- sub("GR2\\_TOPCONS\\.table\\.tsv$", "", basename(top_path_r2))
  
  # --- locate BITACORA R2 gene_table ---
  bitacora_r2_path <- file.path(dir_bitacora_r2, paste0(base_name_r2, "GR2_gene_table.tsv"))
  if (!file.exists(bitacora_r2_path)) {
    # fallback: some pipelines keep the same suffix as R1
    bitacora_r2_path2 <- file.path(dir_bitacora_r2, paste0(base_name_r2, "GR_gene_table.tsv"))
    if (file.exists(bitacora_r2_path2)) bitacora_r2_path <- bitacora_r2_path2
  }
  
  if (!file.exists(bitacora_r2_path)) {
    message("⚠️ No se encontró BITACORA R2 gene_table para ", base_name_r2, " | esperado: ", bitacora_r2_path)
    next
  }
  
  # --- locate BITACORA R1 gene_table (to subtract known genes) ---
  bitacora_r1_path <- file.path(dir_classif, paste0(base_name_r2, "GR_gene_table.tsv"))
  if (!file.exists(bitacora_r1_path)) {
    message("⚠️ No se encontró BITACORA R1 gene_table para ", base_name_r2, " (no puedo calcular NEW-only). Saltando.")
    next
  }
  
  # --- read tables ---
  topcons_r2 <- read.table(top_path_r2, sep = "\t", header = TRUE)
  classif_r2 <- read.table(bitacora_r2_path, sep = "\t", header = TRUE)
  classif_r1 <- read.table(bitacora_r1_path, sep = "\t", header = TRUE)
  
  # --- add 80% flag in R2 (same rule as Round1) ---
  classif_r2$`80%` <- ifelse(classif_r2$Length > 296,
                             "Complete",
                             classif_r2$Status.Complete.Partial.Pseudogene)
  
  # --- define NEW genes: Gene.name in R2 but not in R1 ---
  # Normalize IDs lightly (trim + remove \r) to prevent false mismatches
  r1_ids <- classif_r1 %>%
    mutate(.id = str_replace_all(str_trim(as.character(Gene.name)), "\r", "")) %>%
    filter(!is.na(.id), .id != "") %>%
    distinct(.id)
  
  classif_r2_new <- classif_r2 %>%
    mutate(.id = str_replace_all(str_trim(as.character(Gene.name)), "\r", "")) %>%
    anti_join(r1_ids, by = ".id") %>%     # keep only R2 \ R1
    select(-.id)
  
  # --- merge NEW-only BITACORA R2 with TOPCONS R2 ---
  # Same structure as Round1 merge, but on the filtered R2 table
  tabla_r2_new <- merge(classif_r2_new, topcons_r2, by.x = "Gene.name", by.y = "GeneID", all.x = TRUE)
  
  # --- add species column ---
  species <- gsub("^[0-9]+", "", base_name_r2)
  tabla_r2_new$species <- species
  
  # --- export outputs ---
  #out_file_r2  <- file.path(dir_out_r2, paste0(base_name_r2, "_GR_curado_NEWONLY.tsv"))
  #out_file_r22 <- file.path(dir_out_r2, paste0(base_name_r2, "_GR_curado_NEWONLY.xlsx"))
  
  #write.table(tabla_r2_new, out_file_r2, sep = "\t", quote = FALSE, row.names = FALSE)
  #write.xlsx(tabla_r2_new, out_file_r22, rowNames = FALSE)
  
  message("✅ Procesado R2 NEW-only: ", base_name_r2,
          " | NEW genes: ", nrow(classif_r2_new))
}


# =========================
# 2B) NEWONLY "GR*" genes: map to Round1 by Scaffold.location
#     Goal:
#       - Keep ONLY genes in NEWONLY whose Gene.name starts with "GR"
#       - For each of those, bring the corresponding Round1 info (if any)
#         by matching Scaffold.location (same locus).
#     Output:
#       - Excel with:
#           * MAPPED_FULL: all columns from R2 NEWONLY (incl TOPCONS R2) + R1 columns
#           * OVERVIEW: compact mapping table (R2 name vs R1 name)
# =========================

# Column name used as coordinate key (exists already in your tables)
coord_col <- "Scaffold.location"

# 1) Keep only NEWONLY genes whose name starts with "GR"
tabla_r2_new_GR <- tabla_r2_new %>%
  dplyr::filter(grepl("^GR", Gene.name))

if (nrow(tabla_r2_new_GR) == 0) {
  message("ℹ️  No hay genes NEWONLY que empiecen por 'GR' para ", base_name_r2,
          " (se omite el mapeo por Scaffold.location).")
} else {
  
  # 2) Read TOPCONS Round1 for this spp (optional but recommended)
  topcons_r1_path <- file.path(dir_topcons, paste0(base_name_r2, "GR_TOPCONS.table.tsv"))
  
  if (!file.exists(topcons_r1_path)) {
    message("⚠️ No se encontró TOPCONS R1 para ", base_name_r2,
            " | esperado: ", topcons_r1_path,
            " (se hará join solo con gen_table R1).")
    topcons_r1 <- NULL
  } else {
    topcons_r1 <- read.table(topcons_r1_path, sep = "\t", header = TRUE)
  }
  
  # 3) Build FULL Round1 table (gen_table R1 + TOPCONS R1 if available)
  #    Note: classif_r1 is already loaded earlier in your loop.
  if (!is.null(topcons_r1)) {
    tabla_r1_all <- merge(classif_r1, topcons_r1, by.x = "Gene.name", by.y = "GeneID", all.x = TRUE)
  } else {
    tabla_r1_all <- classif_r1
  }
  
  # 4) Check coordinate column exists in both tables
  if (!(coord_col %in% names(tabla_r2_new_GR)) || !(coord_col %in% names(tabla_r1_all))) {
    message("⚠️ Falta la columna ", coord_col, " en R1 o R2 para ", base_name_r2,
            ". Revisa el header.")
  } else {
    
    # Ensure consistent type (avoid factor mismatch)
    tabla_r2_new_GR[[coord_col]] <- as.character(tabla_r2_new_GR[[coord_col]])
    tabla_r1_all[[coord_col]]    <- as.character(tabla_r1_all[[coord_col]])
    
    # Add explicit gene-name columns for clarity in the final table
    r2_side <- tabla_r2_new_GR %>%
      mutate(Gene.name.R2 = Gene.name)
    
    r1_side <- tabla_r1_all %>%
      mutate(Gene.name.R1 = Gene.name)
    
    # 5) LEFT JOIN by Scaffold.location:
    #    Keep all NEWONLY_GR (R2) rows; add R1 columns where coord matches.
    mapped_GR <- dplyr::left_join(
      r2_side,
      r1_side,
      by = coord_col,
      suffix = c(".R2", ".R1")
    ) %>%
      mutate(
        spp = base_name_r2,
        species = gsub("^[0-9]+", "", base_name_r2),
        has_match_in_R1 = !is.na(Gene.name.R1),
        name_changed = ifelse(!is.na(Gene.name.R1) & Gene.name.R1 != Gene.name.R2, TRUE, FALSE)
      )
    
    # 6) Compact overview (easy to inspect in Excel)
    overview <- mapped_GR %>%
      select(
        spp, species,
        !!coord_col,
        Gene.name.R2, Gene.name.R1,
        has_match_in_R1, name_changed
      ) %>%
      arrange(desc(name_changed), desc(has_match_in_R1), .data[[coord_col]])
    
    # 7) Export to the same directory as other Round2 outputs
    species_tag <- gsub("^[0-9]+", "", base_name_r2)
    
    out_map_xlsx <- file.path(dir_out_r2, paste0(base_name_r2, "_", species_tag, "_GR_NEWONLY_mapped_to_R1.xlsx"))
    out_map_tsv  <- file.path(dir_out_r2, paste0(base_name_r2, "_", species_tag, "_GR_NEWONLY_mapped_to_R1.tsv"))
    
    wb <- createWorkbook()
    addWorksheet(wb, "MAPPED_FULL")
    writeData(wb, "MAPPED_FULL", mapped_GR)
    
    addWorksheet(wb, "OVERVIEW")
    writeData(wb, "OVERVIEW", overview)
    
    saveWorkbook(wb, out_map_xlsx, overwrite = TRUE)
    write.table(overview, out_map_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    message("✅ Mapeo NEWONLY GR* -> R1 por ", coord_col, ": ", base_name_r2,
            " | GR* NEWONLY: ", nrow(tabla_r2_new_GR),
            " | matches in R1: ", sum(overview$has_match_in_R1, na.rm = TRUE),
            " | name_changed: ", sum(overview$name_changed, na.rm = TRUE))
  }
}
