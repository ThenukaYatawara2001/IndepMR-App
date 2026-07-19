# R/03_independence_engine.R
# ============================================================
# IndepMR — Ancestry Independence Engine
#
# DESIGN PRINCIPLE (simplified for reliability and clarity)
# ─────────────────────────────────────────────────────────────
# Ancestry independence is determined at the MAIN GROUP level only.
#
# Five main ancestry groups are recognised:
#   1. European
#   2. African
#   3. Asian
#   4. Middle Eastern
#   5. American
#   (+ Pacific Islander and Other as secondary categories)
#
# RULE: Two studies are ancestry-independent if and only if their
#       EXTRACTED_TERMS map to DIFFERENT main groups.
#       Same main group → same ancestry (not independent).
#       Different main group → different ancestry (independent).
#
# WITHIN-ANCESTRY INDEPENDENCE (second checkbox)
# ─────────────────────────────────────────────────────────────
# Within the same main group, finer subgroup-level independence
# is assessed using the known subgroup hierarchy:
#
#   European subgroups  : Scandinavian, British, Iberian, Dutch,
#                         Italian, Greek, ...
#   Asian subgroups     : East Asian, South Asian, South East Asian,
#                         Central Asian
#   African subgroups   : West African
#   American subgroups  : Hispanic
#
# Studies from DIFFERENT known subgroups within the same main
# group are treated as sub-group independent (distinct LD
# structure). Studies from the SAME subgroup, or where only a
# broad continental label is available, are not separable at
# the subgroup level.
#
# "Other" / unknown ancestry → independence cannot be determined.
# ============================================================

# ─────────────────────────────────────────────────────────────
# 1.  Alias / normalisation map
# ─────────────────────────────────────────────────────────────

TERM_MAP <- c(
  # European aliases
  "other european"                  = "european",
  "european ancestry"               = "european",
  "white"                           = "european",
  "other white"                     = "european",
  "whole exome sequenced european"  = "european",
  "whole genome sequenced european" = "european",
  "north western european"          = "northern european",
  "bristish"                        = "british",
  "white british"                   = "british",
  "metropolitan french"             = "french",
  "old order amish"                 = "amish",
  "british central"                 = "british",
  "non finnish european"            = "european",
  # Asian aliases
  "han chinese"                     = "chinese",
  "chinese han"                     = "chinese",
  "southern han chinese"            = "chinese",
  "indian asian"                    = "indian",
  "central"                         = "central asian",
  "native american ancestries"      = "native american",
  "pima indian"                     = "indian",
  # Middle Eastern
  "greater middle eastern"          = "middle eastern",
  # American
  "hispanic american"               = "hispanic",
  "latino american"                 = "latino",
  "latin american"                  = "latino",
  # Unknown / Other
  "black"                           = "other",
  "admixed"                         = "other",
  "unknown"                         = "other",
  "other admixed"                   = "other",
  "admixed american"                = "other",
  "other admixed ancestry"          = "other",
  "mixed"                           = "other",
  "nr"                              = "other",
  "admixed african"                 = "other",
  "native"                          = "native american"
)

normalize_term <- function(term) {
  t <- tolower(trimws(as.character(term)))
  t <- stringr::str_replace_all(t, "\\s+", " ")
  t <- stringr::str_replace_all(t, "^(and|or)\\s+", "")
  t <- stringr::str_replace_all(t, "[\\.;]+$", "")
  if (t %in% names(TERM_MAP)) return(unname(TERM_MAP[t]))
  t
}

# ─────────────────────────────────────────────────────────────
# 2.  Main group lookup tables
# ─────────────────────────────────────────────────────────────

EUROPEAN_TERMS <- tolower(c(
  "european", "british", "icelandic", "sardinian", "finnish", "german",
  "shetland isles origin", "caucasian", "greek", "french", "amish",
  "scandinavian", "orcadian", "norwegian", "scottish", "dutch",
  "erasmus rucphen", "basque", "swedish", "iberian", "spanish", "irish",
  "black british", "italian", "danish", "estonian", "non hispanic white",
  "barbadian", "cilento", "korculan", "mylopotamos", "pomak",
  "northern european", "non finnish european"
))

AFRICAN_TERMS <- tolower(c(
  "african", "sub saharan african", "west african", "yoruban",
  "ugandan", "malagasy", "afro caribbean"
))

ASIAN_TERMS <- tolower(c(
  "asian", "east asian", "south asian", "central asian", "south east asian",
  "chinese", "japanese", "korean", "taiwanese", "tibetan",
  "indian", "pakistani", "bangladeshi", "malay", "filipino"
))

MIDDLE_EASTERN_TERMS <- tolower(c(
  "middle eastern", "maghrebian", "arab", "turkish", "iranian",
  "lebanese", "moroccan", "qatari"
))

AMERICAN_TERMS <- tolower(c(
  "american", "hispanic", "latino", "brazilian", "peruvian",
  "surinamese", "native american", "kalinago", "alaskan native"
))

PACIFIC_TERMS <- tolower(c(
  "pacific islander", "native hawaiian", "samoan"
))

OTHER_TERMS <- tolower(c(
  "other", "unknown", "nr", "admixed", "mixed", "black"
))

# ─────────────────────────────────────────────────────────────
# 3.  Map a single normalised term → main group
# ─────────────────────────────────────────────────────────────

get_main_group <- function(term) {
  t <- tolower(trimws(as.character(term)))
  if (t %in% EUROPEAN_TERMS)       return("European")
  if (t %in% AFRICAN_TERMS)        return("African")
  if (t %in% ASIAN_TERMS)          return("Asian")
  if (t %in% MIDDLE_EASTERN_TERMS) return("Middle Eastern")
  if (t %in% AMERICAN_TERMS)       return("American")
  if (t %in% PACIFIC_TERMS)        return("Pacific Islander")
  if (t %in% OTHER_TERMS)          return("Other")
  "Other"   # fallback for unrecognised terms
}

# ─────────────────────────────────────────────────────────────
# 4.  Parse EXTRACTED_TERMS string → set of main groups
# ─────────────────────────────────────────────────────────────

# Split a raw EXTRACTED_TERMS string into individual terms
.split_terms_string <- function(x) {
  if (is.na(x) || trimws(x) == "") return(character(0))
  raw    <- tolower(as.character(x))
  # Split on common delimiters
  chunks <- trimws(unlist(strsplit(raw, "[,:;|/]")))
  chunks <- chunks[nchar(chunks) > 0]
  # Further split on " and " / " or "
  chunks <- trimws(unlist(strsplit(chunks, "\\s+and\\s+|\\s+or\\s+")))
  chunks <- chunks[nchar(chunks) > 1]
  unique(chunks)
}

# Parse a raw EXTRACTED_TERMS value → character vector of main groups
get_main_groups_for_terms <- function(terms_string) {
  chunks <- .split_terms_string(terms_string)
  if (length(chunks) == 0) return(character(0))
  groups <- vapply(chunks, function(ch) get_main_group(normalize_term(ch)), character(1))
  unique(groups[groups != "Other"])   # exclude "Other" from the useful groups
}

# ─────────────────────────────────────────────────────────────
# 5.  Subgroup hierarchy (for within-ancestry independence)
# ─────────────────────────────────────────────────────────────

SUBGROUP_MEMBERS <- list(
  # European
  "scandinavian"    = c("norwegian", "swedish", "danish"),
  "iberian"         = c("spanish", "basque"),
  "british"         = c("scottish", "orcadian", "shetland isles origin", "black british"),
  "dutch"           = c("erasmus rucphen"),
  "italian"         = c("cilento"),
  "greek"           = c("mylopotamos"),
  # Asian
  "east asian"      = c("chinese", "japanese", "korean", "taiwanese", "tibetan"),
  "south asian"     = c("indian", "pakistani", "bangladeshi"),
  "south east asian"= c("malay", "filipino"),
  "central asian"   = character(0),
  # African
  "west african"    = c("yoruban"),
  # American
  "hispanic"        = c("latino")
)

SUBGROUP_LABELS    <- names(SUBGROUP_MEMBERS)
MEMBER_TO_SUBGROUP <- unlist(
  lapply(names(SUBGROUP_MEMBERS), function(sg) {
    setNames(rep(sg, length(SUBGROUP_MEMBERS[[sg]])), SUBGROUP_MEMBERS[[sg]])
  }), use.names = TRUE
)

# Broad continental umbrella terms (cannot distinguish sub-populations)
UMBRELLA_TERMS <- tolower(c(
  "european", "african", "asian", "american",
  "middle eastern", "pacific islander",
  "caucasian", "non hispanic white", "sub saharan african"
))

# ─────────────────────────────────────────────────────────────
# 6.  Profile: structured summary of a study's ancestry
# ─────────────────────────────────────────────────────────────
#
# A study's EXTRACTED_TERMS may contain MULTIPLE ancestry labels,
# e.g. "East Asian, European" or "British, Irish".
# Each label is normalised, assigned a main group, and classified
# as either an umbrella term, a known subgroup, or a member term.
#
# KEY DESIGN — per-main-group umbrella tracking:
# ─────────────────────────────────────────────────────────────
# $umbrella_groups  — set of main groups for which THIS study uses
#                     only a broad/umbrella label (e.g. "European"
#                     with no finer British/Scandinavian term).
#
# This is critical for multi-term studies.  Example:
#   EXTRACTED_TERMS = "East Asian, European"
#   → main_groups = c("Asian", "European")
#   → subgroups   = c("east asian")          (from the Asian side)
#   → umbrella_groups = c("European")        (the European term is broad)
#
# When checking subgroup independence against a British exposure:
#   • Shared main group = "European"
#   • Outcome has "European" in umbrella_groups → UMBRELLA_ONLY
#   (The East Asian subgroup label is irrelevant here because it
#    belongs to a different main group than the shared one.)
#
# Returns a list:
#   $main_groups      — unique main groups (excl. "Other")
#   $has_other        — TRUE if any term maps to "Other"
#   $is_umbrella_only — TRUE if ALL main groups use only umbrella terms
#   $umbrella_groups  — main groups where only umbrella labels are present
#   $subgroups        — known subgroup labels present in the study
#   $member_subgroups — subgroups inferred from member terms
#   $subgroup_by_main — named list: main_group → subgroups resolved for it
#   $all_terms        — all normalised terms
build_ancestry_profile <- function(terms_string) {
  chunks <- .split_terms_string(terms_string)
  if (length(chunks) == 0) {
    return(list(
      main_groups      = character(0),
      has_other        = FALSE,
      is_umbrella_only = TRUE,
      umbrella_groups  = character(0),
      subgroups        = character(0),
      member_subgroups = character(0),
      subgroup_by_main = list(),
      all_terms        = character(0)
    ))
  }

  norm_terms <- unique(vapply(chunks, function(ch) normalize_term(ch), character(1)))
  main_groups_all <- unique(vapply(norm_terms, get_main_group, character(1)))
  has_other       <- "Other" %in% main_groups_all
  useful_groups   <- main_groups_all[main_groups_all != "Other"]

  # For each term, determine its subgroup or member status
  subgroups       <- norm_terms[norm_terms %in% SUBGROUP_LABELS]
  member_terms    <- norm_terms[norm_terms %in% names(MEMBER_TO_SUBGROUP)]
  member_subgroups <- unique(unname(MEMBER_TO_SUBGROUP[member_terms]))

  # ── Per-main-group umbrella / subgroup tracking ──────────────────────────
  # For each useful main group, collect which terms belong to it and
  # classify whether those terms are umbrella-only or have subgroup resolution.
  umbrella_groups  <- character(0)
  subgroup_by_main <- list()

  for (mg in useful_groups) {
    # Terms belonging to this main group
    mg_terms <- norm_terms[vapply(norm_terms, function(t) get_main_group(t) == mg,
                                  logical(1))]

    # Subgroup labels explicitly for this main group
    mg_subs <- mg_terms[mg_terms %in% SUBGROUP_LABELS]
    # Member terms → their parent subgroup name
    mg_mem_raw  <- mg_terms[mg_terms %in% names(MEMBER_TO_SUBGROUP)]
    mg_mem_subs <- unique(unname(MEMBER_TO_SUBGROUP[mg_mem_raw]))
    # Specific named populations that are neither a subgroup label nor a member
    # (e.g. "Finnish", "Sardinian", "Estonian") — treat each as its own
    # unique sub-population.  This allows them to be recognised as
    # subgroup-independent from other named sub-populations (e.g. British).
    mg_specific <- mg_terms[
      !mg_terms %in% UMBRELLA_TERMS &
      !mg_terms %in% SUBGROUP_LABELS &
      !mg_terms %in% names(MEMBER_TO_SUBGROUP)
    ]
    mg_all_subs <- unique(c(mg_subs, mg_mem_subs, mg_specific))

    subgroup_by_main[[mg]] <- mg_all_subs

    # Is this main group "contaminated" by an umbrella term?
    # If ANY term contributing to this main group is a broad umbrella label
    # (e.g. the study says "Italian, European"), then the main group CANNOT be
    # treated as subgroup-specific.  The broad "European" label means this
    # study may include participants from ANY European sub-group, even if
    # "Italian" is also listed.  Therefore, if any mg_term is an umbrella
    # term, we classify this main group as umbrella.
    mg_has_umbrella <- any(mg_terms %in% UMBRELLA_TERMS)
    if (mg_has_umbrella) umbrella_groups <- c(umbrella_groups, mg)
  }

  # Overall is_umbrella_only: TRUE when EVERY main group has at least one
  # umbrella term (i.e. no main group has ONLY specific subgroup terms).
  # Used to hide the second checkbox when the exposure cannot be resolved
  # to a specific sub-population in any of its main groups.
  is_umbrella_only <- length(useful_groups) > 0 &&
                      length(umbrella_groups) == length(useful_groups)

  list(
    main_groups      = useful_groups,
    has_other        = has_other,
    is_umbrella_only = is_umbrella_only,
    umbrella_groups  = umbrella_groups,
    subgroups        = subgroups,
    member_subgroups = member_subgroups,
    subgroup_by_main = subgroup_by_main,
    all_terms        = norm_terms
  )
}

# ─────────────────────────────────────────────────────────────
# 7.  Profile cache
# ─────────────────────────────────────────────────────────────

PROFILE_CACHE <- new.env(parent = emptyenv())

get_profile_for_gcst <- function(gcst_id, merged_df) {
  gcst_id <- stringr::str_trim(as.character(gcst_id))
  if (is.na(gcst_id) || gcst_id == "") return(NULL)

  key <- paste0("P__", gcst_id)
  if (exists(key, envir = PROFILE_CACHE, inherits = FALSE))
    return(get(key, envir = PROFILE_CACHE, inherits = FALSE))

  row <- merged_df %>%
    dplyr::filter(`STUDY ACCESSION` == gcst_id) %>%
    dplyr::slice(1)
  if (nrow(row) == 0) {
    assign(key, NULL, envir = PROFILE_CACHE)
    return(NULL)
  }

  terms_str <- as.character(row$EXTRACTED_TERMS[[1]])
  prof      <- build_ancestry_profile(terms_str)
  assign(key, prof, envir = PROFILE_CACHE)
  prof
}

# ─────────────────────────────────────────────────────────────
# 8.  Core independence test  (MAIN GROUP LEVEL)
# ─────────────────────────────────────────────────────────────

# Returns:
#   TRUE   — profiles belong to DIFFERENT main groups (ancestry-independent)
#   FALSE  — profiles share at least one main group (same ancestry)
#   NA     — cannot determine (missing profile, Other ancestry, or no groups found)
pair_is_independent <- function(p1, p2) {
  if (is.null(p1) || is.null(p2))               return(NA)
  if (p1$has_other || p2$has_other)             return(NA)
  if (length(p1$main_groups) == 0 ||
      length(p2$main_groups) == 0)              return(NA)
  length(intersect(p1$main_groups, p2$main_groups)) == 0
}

# Convenience: test one candidate against multiple reference IDs
# Returns TRUE / FALSE / NA
is_independent_from_all <- function(candidate_id, reference_ids, merged_df) {
  prof_c <- get_profile_for_gcst(candidate_id, merged_df)
  if (is.null(prof_c)) return(NA)

  reference_ids <- unique(stringr::str_trim(as.character(reference_ids)))
  reference_ids <- reference_ids[reference_ids != "" & !is.na(reference_ids)]
  if (length(reference_ids) == 0) return(NA)

  for (rid in reference_ids) {
    prof_r <- get_profile_for_gcst(rid, merged_df)
    if (is.null(prof_r)) return(NA)
    ok <- pair_is_independent(prof_c, prof_r)
    if (is.na(ok)) return(NA)
    if (!ok)       return(FALSE)
  }
  TRUE
}

# ─────────────────────────────────────────────────────────────
# 9.  Within-ancestry sub-group independence
# ─────────────────────────────────────────────────────────────
#
# Classify a candidate study vs. a reference study at the
# subgroup level, for each main group they share.
#
# RULES (applied per shared main group)
# ─────────────────────────────────────────────────────────────
# Rule A — Other / unknown ancestry → OTHER_PRESENT
#
# Rule B — Reference (exposure) is umbrella-only in a shared
#           main group:
#           e.g. exposure terms include plain "European" with no
#           specific sub-group.  We cannot find any subgroup-
#           independent outcome because the exposure's sub-group
#           within that main group is unknown.
#           → UMBRELLA_ONLY
#
# Rule C — Candidate (outcome) is umbrella-only in a shared
#           main group:
#           e.g. outcome = "East Asian, European" and exposure =
#           "British" (European sub-group).  The shared main group
#           is "European".  The outcome's European term is the broad
#           umbrella "European" which INCLUDES British.  Even though
#           the outcome also has "East Asian" (a different main group
#           with a specific sub-group label), that does NOT make it
#           subgroup-independent with respect to the shared European
#           main group.
#           → UMBRELLA_ONLY
#
# Rule D — Both have specific sub-group resolution in the shared
#           main group:
#           Compare the sub-groups for that main group only.
#           Same sub-group → SUBGROUP_OVERLAPPING
#           Different sub-groups → SUBGROUP_INDEPENDENT
#
# When multiple main groups are shared, the MOST RESTRICTIVE
# result across all shared groups is returned.
#
# Returns one of:
#   "SUBGROUP_INDEPENDENT"  — all shared main groups show different
#                             known sub-groups
#   "SUBGROUP_OVERLAPPING"  — at least one shared main group has the
#                             same sub-group
#   "UMBRELLA_ONLY"         — at least one shared main group uses a
#                             broad label in either the candidate or
#                             the reference (cannot resolve)
#   "OTHER_PRESENT"         — ancestry is broad / unknown
#   "UNKNOWN"               — profile missing or no shared main group
get_within_ancestry_status <- function(p_candidate, p_reference) {
  if (is.null(p_candidate) || is.null(p_reference)) return("UNKNOWN")

  # Rule A
  if (p_candidate$has_other || p_reference$has_other) return("OTHER_PRESENT")

  shared_main <- intersect(p_candidate$main_groups, p_reference$main_groups)
  if (length(shared_main) == 0) return("UNKNOWN")  # caller should check this

  # Priority order for final result (most restrictive wins across shared groups)
  priority <- c(
    "OTHER_PRESENT"       = 5,
    "UNKNOWN"             = 4,
    "SUBGROUP_OVERLAPPING"= 3,
    "UMBRELLA_ONLY"       = 2,
    "SUBGROUP_INDEPENDENT"= 0
  )

  worst <- "SUBGROUP_INDEPENDENT"

  for (mg in shared_main) {
    # ── Rule B: exposure umbrella in this shared main group ──────────────────
    if (mg %in% p_reference$umbrella_groups) {
      st <- "UMBRELLA_ONLY"
      if (priority[st] > priority[worst]) worst <- st
      next
    }

    # ── Rule C: candidate umbrella in this shared main group ─────────────────
    # A broad continental label (e.g. "European") contains ALL sub-groups,
    # so it can never be subgroup-independent from any specific sub-group.
    if (mg %in% p_candidate$umbrella_groups) {
      st <- "UMBRELLA_ONLY"
      if (priority[st] > priority[worst]) worst <- st
      next
    }

    # ── Rule D: both have sub-group resolution for this main group ───────────
    # Use per-main-group sub-group lists (not the global ones).
    c_subs <- p_candidate$subgroup_by_main[[mg]]
    r_subs <- p_reference$subgroup_by_main[[mg]]

    # If neither side has sub-group info for this main group
    # (they have non-umbrella terms that are not recognised sub-groups),
    # treat conservatively as umbrella-level
    if (length(c_subs) == 0 || length(r_subs) == 0) {
      st <- "UMBRELLA_ONLY"
      if (priority[st] > priority[worst]) worst <- st
      next
    }

    if (length(intersect(c_subs, r_subs)) > 0) {
      st <- "SUBGROUP_OVERLAPPING"
    } else {
      st <- "SUBGROUP_INDEPENDENT"
    }
    if (priority[st] > priority[worst]) worst <- st
  }

  worst
}

# Aggregate over multiple reference IDs (most restrictive wins)
# Priority: OTHER_PRESENT > UNKNOWN > SUBGROUP_OVERLAPPING > UMBRELLA_ONLY > SUBGROUP_INDEPENDENT
get_within_ancestry_status_from_all <- function(candidate_id, reference_ids, merged_df) {
  prof_c <- get_profile_for_gcst(candidate_id, merged_df)
  if (is.null(prof_c)) return("UNKNOWN")

  reference_ids <- unique(stringr::str_trim(as.character(reference_ids)))
  reference_ids <- reference_ids[reference_ids != "" & !is.na(reference_ids)]
  if (length(reference_ids) == 0) return("UNKNOWN")

  priority <- c(
    "OTHER_PRESENT"       = 5,
    "UNKNOWN"             = 4,
    "SUBGROUP_OVERLAPPING"= 3,
    "UMBRELLA_ONLY"       = 2,
    "SUBGROUP_INDEPENDENT"= 0
  )

  worst <- "SUBGROUP_INDEPENDENT"
  for (rid in reference_ids) {
    prof_r <- get_profile_for_gcst(rid, merged_df)
    st     <- get_within_ancestry_status(prof_c, prof_r)
    if (priority[st] > priority[worst]) worst <- st
  }
  worst
}

# ─────────────────────────────────────────────────────────────
# 10.  UI helper: describe exposure ancestry for messages
# ─────────────────────────────────────────────────────────────

# Returns a list describing the ancestry of the given exposure GCST IDs.
# Used by the server to decide which messages to show and whether to
# display the second (sub-group independence) checkbox.
#
# Fields returned:
#   $main_groups       — unique main-group labels (e.g. "European")
#   $subgroups         — known subgroup labels present in the exposure
#                        (e.g. "british", "scandinavian")
#   $has_other         — TRUE if any term maps to "Other" / unknown
#   $is_umbrella_only  — TRUE if the exposure uses only broad continental
#                        labels with no specific subgroup information.
#                        In this case the second checkbox is hidden because
#                        we cannot find subgroup-independent outcomes for an
#                        exposure whose subgroup is not specified.
#
# IMPORTANT: even when is_umbrella_only = FALSE (i.e. the exposure has a
# known subgroup such as "British"), some candidate outcomes may still be
# classified as UMBRELLA_ONLY because their own ancestry label is broad
# (e.g. "European").  Such outcomes are shown as "Broad only" in the
# Sub-group column and are excluded when the second checkbox is ON.
describe_exposure_ancestry <- function(gcst_ids, merged_df) {
  gcst_ids <- gcst_ids[!is.na(gcst_ids) & gcst_ids != ""]
  if (length(gcst_ids) == 0) return(NULL)

  all_groups    <- character(0)
  all_subgroups <- character(0)
  has_other     <- FALSE
  is_umb        <- TRUE   # assume umbrella until proven otherwise

  for (gid in gcst_ids) {
    prof <- get_profile_for_gcst(gid, merged_df)
    if (is.null(prof)) next
    if (prof$has_other) has_other <- TRUE
    all_groups    <- unique(c(all_groups,    prof$main_groups))
    all_subgroups <- unique(c(all_subgroups, prof$subgroups, prof$member_subgroups))
    if (!prof$is_umbrella_only) is_umb <- FALSE
  }

  list(
    main_groups      = all_groups,
    subgroups        = all_subgroups,     # specific subgroups present in exposure
    has_other        = has_other,
    is_umbrella_only = is_umb && length(all_subgroups) == 0
  )
}

# ─────────────────────────────────────────────────────────────
# 11. Legacy aliases kept for compatibility
# ─────────────────────────────────────────────────────────────

# Older function names used in earlier versions of the server;
# map them to the new implementations so no server code breaks.

get_subgroup_independence_from_all <- function(candidate_id, reference_ids, merged_df) {
  get_within_ancestry_status_from_all(candidate_id, reference_ids, merged_df)
}

get_subgroup_independence_status <- function(p_candidate, p_reference) {
  get_within_ancestry_status(p_candidate, p_reference)
}
