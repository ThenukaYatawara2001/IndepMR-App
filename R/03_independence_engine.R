# R/03_independence_engine.R
# ============================================================
# IndepMR — Independence engine (ancestry / extracted terms)
#
# ✅ RULE (UPDATED):
#   - If either profile contains main group "Other",
#     then the pair is NOT independent.
#   - "Other" cannot be independent from any group (even "Other").
# ============================================================

# ----------------------------
# 1) Split helpers
# ----------------------------
split_extracted_term <- function(term) {
  t <- tolower(trimws(as.character(term)))
  if (nchar(t) == 0) return(character(0))
  
  t <- stringr::str_replace_all(t, "\\s+and\\s+", ",")
  t <- stringr::str_replace_all(t, "\\s+or\\s+", ",")
  t <- stringr::str_replace_all(t, "/", ",")
  t <- stringr::str_replace_all(
    t,
    "\\bafrican american\\s+afro caribbean\\b",
    "african american, afro caribbean"
  )
  
  parts <- trimws(unlist(strsplit(t, ",")))
  parts <- parts[nchar(parts) > 0]
  parts
}

TERM_MAP <- c(
  "other european"="european",
  "european"="european",
  "european ancestry"="european",
  "white"="european",
  "other white"="european",
  "whole exome sequenced european"="european",
  "whole genome sequenced european"="european",
  "north western european"="northern european",
  "bristish"="british",
  "white british"="british",
  "metropolitan french"="french",
  "old order amish"="amish",
  "british central"="british",
  "native american ancestries"="native american",
  "han chinese"="chinese",
  "chinese han"="chinese",
  "southern han chinese"="chinese",
  "indian asian"="indian",
  "greater middle eastern"="middle eastern",
  "hispanic american"="hispanic",
  "latino american"="latino",
  "latin american"="latino",
  "pima indian"="indian",
  "central"="central asian",
  "native"="native american",
  "black"="other",
  "admixed"="other",
  "unknown"="other",
  "other"="other",
  "other admixed"="other",
  "admixed american"="other",
  "other admixed ancestry"="other",
  "mixed"="other",
  "nr"="other",
  "admixed african"="other"
)

normalize_term <- function(term) {
  t <- tolower(trimws(as.character(term)))
  t <- stringr::str_replace_all(t, "\\s+", " ")
  t <- stringr::str_replace_all(t, "^(and|or)\\s+", "")
  t <- stringr::str_replace_all(t, "[\\.;]+$", "")
  if (t %in% names(TERM_MAP)) return(unname(TERM_MAP[t]))
  t
}

split_dual_identity <- function(term) {
  t <- tolower(trimws(as.character(term)))
  t <- stringr::str_replace_all(t, "\\s+", " ")
  
  dual_map <- list(
    "european american"=c("european","american"),
    "african american"=c("african","american"),
    "east asian american"=c("east asian","american"),
    "japanese american"=c("japanese","american"),
    "african british"=c("african","british"),
    "british african"=c("african","british"),
    "asian british"=c("asian","british"),
    "american indian"=c("american","indian")
  )
  if (t %in% names(dual_map)) return(dual_map[[t]])
  
  m <- stringr::str_match(t, "^(.+)\\s+american$")
  if (!is.na(m[1,2])) {
    left <- trimws(m[1,2])
    if (nchar(left) > 0 && !(left %in% c("non","other"))) return(c(left, "american"))
  }
  c(t)
}

# ----------------------------
# 2) Group dictionaries
# ----------------------------
UMBRELLA_MAIN <- list(
  "European"=c("european","caucasian","non hispanic white"),
  "African"=c("african","sub saharan african"),
  "Asian"=c("asian"),
  "Middle Eastern"=c("middle eastern"),
  "American"=c("american"),
  "Pacific Islander"=c("pacific islander"),
  "Other"=c("other")  # label exists, but now "Other" blocks independence
)

SUBGROUP_MEMBERS <- list(
  "east asian"=c("chinese","japanese","korean","taiwanese","tibetan"),
  "south asian"=c("indian","pakistani","bangladeshi"),
  "south east asian"=c("malay","filipino"),
  "central asian"=character(0),
  "scandinavian"=c("norwegian","swedish","danish"),
  "iberian"=c("spanish","basque"),
  "british"=c("scottish","orcadian","shetland isles origin","black british"),
  "dutch"=c("erasmus rucphen"),
  "italian"=c("cilento"),
  "greek"=c("mylopotamos"),
  "west african"=c("yoruban"),
  "hispanic"=c("latino")
)
SUBGROUP_LABELS <- names(SUBGROUP_MEMBERS)
MEMBER_TO_SUBGROUP <- unlist(lapply(names(SUBGROUP_MEMBERS), function(sg) {
  setNames(rep(sg, length(SUBGROUP_MEMBERS[[sg]])), SUBGROUP_MEMBERS[[sg]])
}), use.names = TRUE)

EUROPEAN_TERMS <- tolower(c(
  "european","british","icelandic","sardinian","finnish","german","shetland isles origin",
  "caucasian","greek","french","amish","scandinavian","orcadian","norwegian","scottish","dutch",
  "erasmus rucphen","non finnish european","basque","swedish","iberian","spanish","irish",
  "black british","italian","danish","estonian","non hispanic white","barbadian","cilento",
  "korculan","mylopotamos","pomak"
))
AFRICAN_TERMS <- tolower(c("african","sub saharan african","west african","yoruban","ugandan","malagasy","afro caribbean"))
ASIAN_TERMS   <- tolower(c("asian","east asian","south asian","central asian","south east asian",
                           "chinese","japanese","korean","taiwanese","tibetan",
                           "indian","pakistani","bangladeshi","malay","filipino"))
MENA_TERMS    <- tolower(c("middle eastern","maghrebian","arab","turkish","iranian","lebanese","moroccan","qatari"))
AMERICAN_TERMS<- tolower(c("american","hispanic","latino","brazilian","peruvian","surinamese","native american","kalinago","alaskan native"))
PACIFIC_TERMS <- tolower(c("pacific islander","native hawaiian","samoan"))
OTHER_TERMS   <- tolower(c("other","unknown","nr","admixed","mixed","black"))

# ----------------------------
# 3) Main group + type mapping
# ----------------------------
get_main_group <- function(term) {
  t <- tolower(trimws(as.character(term)))
  
  if (t %in% EUROPEAN_TERMS) return("European")
  if (t %in% AFRICAN_TERMS)  return("African")
  if (t %in% ASIAN_TERMS)    return("Asian")
  if (t %in% MENA_TERMS)     return("Middle Eastern")
  if (t %in% AMERICAN_TERMS) return("American")
  if (t %in% PACIFIC_TERMS)  return("Pacific Islander")
  if (t %in% OTHER_TERMS)    return("Other")
  
  if (t %in% tolower(UMBRELLA_MAIN$European)) return("European")
  if (t %in% tolower(UMBRELLA_MAIN$African))  return("African")
  if (t %in% tolower(UMBRELLA_MAIN$Asian))    return("Asian")
  if (t %in% tolower(UMBRELLA_MAIN$`Middle Eastern`)) return("Middle Eastern")
  if (t %in% tolower(UMBRELLA_MAIN$American)) return("American")
  if (t %in% tolower(UMBRELLA_MAIN$`Pacific Islander`)) return("Pacific Islander")
  if (t %in% tolower(UMBRELLA_MAIN$Other)) return("Other")
  
  "Other"
}

type_map <- function(term) {
  t <- tolower(trimws(term))
  mg <- get_main_group(t)
  if (t %in% tolower(UMBRELLA_MAIN[[mg]])) return(list(type="UMBRELLA", val=mg))
  if (t %in% SUBGROUP_LABELS) return(list(type="SUBGROUP", val=t))
  if (t %in% names(MEMBER_TO_SUBGROUP)) return(list(type="MEMBER", val=unname(MEMBER_TO_SUBGROUP[t])))
  list(type="OTHER", val=NA_character_)
}

row_profile <- function(term_set) {
  term_set <- unique(tolower(trimws(term_set)))
  term_set <- term_set[nchar(term_set) > 0]
  
  prof <- list(
    main_groups       = character(0),
    umbrellas_by_main = character(0),
    subgroup_labels   = character(0),
    member_subgroups  = character(0),
    exact_terms       = term_set
  )
  
  for (term in term_set) {
    mg <- get_main_group(term)
    prof$main_groups <- unique(c(prof$main_groups, mg))
    
    tt <- type_map(term)
    if (tt$type == "UMBRELLA") prof$umbrellas_by_main <- unique(c(prof$umbrellas_by_main, tt$val))
    if (tt$type == "SUBGROUP") prof$subgroup_labels   <- unique(c(prof$subgroup_labels, tt$val))
    if (tt$type == "MEMBER" && !is.na(tt$val)) prof$member_subgroups <- unique(c(prof$member_subgroups, tt$val))
  }
  
  prof
}

# ----------------------------
# 4) Independence logic  (UPDATED)
# ----------------------------
pair_is_independent <- function(p1, p2) {
  if (is.null(p1) || is.null(p2)) return(NA)
  
  # ✅ NEW RULE:
  # If either profile contains "Other" main group -> NOT independent
  if ("Other" %in% p1$main_groups || "Other" %in% p2$main_groups) {
    return(FALSE)
  }
  
  mg1 <- unique(p1$main_groups)
  mg2 <- unique(p2$main_groups)
  
  common_main <- intersect(mg1, mg2)
  
  # If no main-group overlap -> independent
  if (length(common_main) == 0) return(TRUE)
  
  for (mg in common_main) {
    # umbrella blocks independence
    if (mg %in% p1$umbrellas_by_main || mg %in% p2$umbrellas_by_main) return(FALSE)
    
    # subgroup/member conflicts
    if (length(intersect(p1$subgroup_labels, p2$member_subgroups)) > 0) return(FALSE)
    if (length(intersect(p2$subgroup_labels, p1$member_subgroups)) > 0) return(FALSE)
    if (length(intersect(p1$subgroup_labels, p2$subgroup_labels)) > 0) return(FALSE)
  }
  
  TRUE
}

parse_extracted_terms_string <- function(x) {
  if (is.na(x) || trimws(x) == "") return(character(0))
  raw <- tolower(as.character(x))
  chunks <- trimws(unlist(strsplit(raw, "[,:;|]")))
  chunks <- chunks[nchar(chunks) > 0]
  
  terms <- character(0)
  for (seg in chunks) {
    for (ft in split_extracted_term(seg)) {
      tnorm <- normalize_term(ft)
      for (st in split_dual_identity(tnorm)) {
        st <- normalize_term(st)
        st <- trimws(st)
        if (nchar(st) < 2) next
        terms <- c(terms, st)
      }
    }
  }
  unique(tolower(trimws(terms)))
}

# ----------------------------
# 5) Cache + public functions
# ----------------------------
PROFILE_CACHE <- new.env(parent = emptyenv())

get_profile_for_gcst <- function(gcst_id, merged_df) {
  gcst_id <- stringr::str_trim(as.character(gcst_id))
  if (is.na(gcst_id) || gcst_id == "") return(NULL)
  
  key <- paste0("P__", gcst_id)
  if (exists(key, envir = PROFILE_CACHE, inherits = FALSE)) {
    return(get(key, envir = PROFILE_CACHE, inherits = FALSE))
  }
  
  row <- merged_df %>% dplyr::filter(`STUDY ACCESSION` == gcst_id) %>% dplyr::slice(1)
  if (nrow(row) == 0) {
    assign(key, NULL, envir = PROFILE_CACHE)
    return(NULL)
  }
  
  terms <- parse_extracted_terms_string(row$EXTRACTED_TERMS[[1]])
  if (length(terms) == 0) {
    assign(key, NULL, envir = PROFILE_CACHE)
    return(NULL)
  }
  
  prof <- row_profile(terms)
  assign(key, prof, envir = PROFILE_CACHE)
  prof
}

is_independent_from_all <- function(candidate_id, reference_ids, merged_df) {
  candidate_prof <- get_profile_for_gcst(candidate_id, merged_df)
  if (is.null(candidate_prof)) return(NA)
  
  reference_ids <- unique(stringr::str_trim(as.character(reference_ids)))
  reference_ids <- reference_ids[reference_ids != "" & !is.na(reference_ids)]
  if (length(reference_ids) == 0) return(NA)
  
  for (rid in reference_ids) {
    rprof <- get_profile_for_gcst(rid, merged_df)
    if (is.null(rprof)) return(NA)
    
    ok <- pair_is_independent(candidate_prof, rprof)
    if (is.na(ok)) return(NA)
    if (!ok) return(FALSE)
  }
  
  TRUE
}

# ============================================================
# Add cohort overlap column
# ============================================================

add_overlap_column <- function(df, exposure_row) {
  
  if (is.null(exposure_row) || nrow(exposure_row) == 0) {
    df$Overlap <- "UNKNOWN"
    return(df)
  }
  
  exp_cohort <- exposure_row$COHORT[1]
  
  df$Overlap <- vapply(df$COHORT, function(out_cohort) {
    detect_cohort_overlap(exp_cohort, out_cohort)
  }, character(1))
  
  df
}