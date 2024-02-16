#Testing first CL ontologies
renv::restore()

library(rols)
library(data.table)

#loading IHEC df
main_csv <- 'IHEC_metadata_harmonization.v1.2.extended.csv'
original_dt <- fread(main_csv)
separator <- '::'


###Ontologies
cl <- Ontology("cl")
efo <- Ontology("efo")
uberon <- Ontology("uberon")
ncit <- Ontology("ncit")

# helper function for different url for EFO
my_Term <-  function(object, id) {
  ## See https://github.com/EBISPOT/ols4/issues/621
  url <- paste0(
    "https://www.ebi.ac.uk/ols4/api/ontologies/",
    object,
    "/terms/http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252F")
  url <- paste0(url, sub(":", "_", id))
  httr2::request(url) |>
    httr2::req_perform() |>
    httr2::resp_body_json() |>
    rols:::termFromJson()
}

hierarchicalAncestors <- function(object) {
  stopifnot(inherits(object, "Term"))
  if (object@is_root)
    return(NULL)
  url <- termLinks(object)[["hierarchicalAncestors"]]
  ans <- lapply(rols:::ols_requests(url, "terms"),
                rols:::termFromJson)
  names(ans) <- sapply(ans, termId)
  rols:::.Terms(x = ans)
}

# this adds one column per entry in ls to the main_dt containing higher_level annotations with overall at most l terms (for each l in ls)
add_higher_order <- function(main_dt, rows_of_interest=rep(TRUE, times=main_dt[, .N]), ls = c(15, 30)){
  
  ls <- sort(ls, decreasing = TRUE)
  
  ancestor_occurences <- as.data.table(main_dt[rows_of_interest, sort(table(unlist(ancestors)))])
  setnames(ancestor_occurences, c('V1', 'N'), c('term', 'counts'))
  unique_ancestor2occurences <- ancestor_occurences[, counts]
  names(unique_ancestor2occurences) <- ancestor_occurences[, term]
  
  main_dt[rows_of_interest, ancestors:=sapply(ancestors, function(a) list(a[order(unique_ancestor2occurences[a])]))]
  
  last_k <- 1L
  for (l in ls) {
    for (k in seq.int(last_k, sum(rows_of_interest, na.rm = TRUE))) {
      new_col <-
        main_dt[rows_of_interest, sapply(ancestors, function(ancestor_vector){
          if(is.null(ancestor_vector)) return(NA)
          ancestor_vector[tryCatch(
            max(which(unique_ancestor2occurences[ancestor_vector] <= k)),
            warning = function(w)
              ifelse(w$message == 'no non-missing arguments to max; returning -Inf', 1, w)
          )]
        }
        )]
      if (l >= uniqueN(new_col)) {
        if (sum(rows_of_interest, na.rm = TRUE) != length(new_col))
          stop('new column has different length than number of rows')
        main_dt[rows_of_interest, (paste('order', l, sep='_')):=new_col]
        last_k <- k
        break
      }
    }
  }

  main_dt
}

ncit_terms <- Terms(ncit, pagesize = 10000)
ncit2ncim <- lapply(ncit_terms, function(x) {
    unique(c(unlist(x@annotation$UMLS_CUI), x@annotation$NCI_META_CUI))
  })
rm(ncit_terms)
gc()
ncit2ncim_atomic <- setNames(unlist(ncit2ncim, use.names=F), rep(names(ncit2ncim), lengths(ncit2ncim)))

for (ncim_col in c('harmonized_sample_disease_ontology_curie', 'harmonized_donor_health_status_ontology_curie')) {
  new_col <- paste('automated', ncim_col, 'ncit', sep = '_')
  new_column_content <- sapply(original_dt[, get(ncim_col)], function(ncim){
    ncits <- unlist(lapply(tstrsplit(ncim, separator, fixed=TRUE), function(single_curie) {
      cui <- unlist(tstrsplit(single_curie, ':', fixed=TRUE, keep = 2))
      result <- names(ncit2ncim_atomic)[cui == ncit2ncim_atomic]
      if (length(result) == 0)
        result <- NULL
      result
    }))
    paste(ncits[!is.null(ncits)], collapse = separator)
  })
  # print(original_dt[get(new_col) != get(paste(sub('^harmonized_', '', ncim_col), 'ncit', sep = '_'))])
  original_dt[, (new_col) := NULL]
  original_dt[, (new_col) := new_column_content]
}

# now go trough the curie cols ----
cols_to_drop <- names(original_dt)[startsWith(names(original_dt), prefix = 'automated') & grepl('ontology', names(original_dt), fixed = TRUE) & !endsWith(names(original_dt), '_ncit')]
# remove the columns in cols_to_drop from original_dt
metadata_dt <- copy(original_dt[, !cols_to_drop, with = FALSE])

ncit_cols <- names(original_dt)[endsWith(names(original_dt), '_ncit')]
for (curie_col in c(ncit_cols, 'harmonized_sample_ontology_curie')){
  print(curie_col)
  spread_curies <- original_dt[, .(single_curie=unlist(tstrsplit(get(curie_col), separator, fixed=TRUE))), by=mget(curie_col)]
  # get the terms from the ontology curie
  if (curie_col == 'harmonized_sample_ontology_curie') {
    spread_curies[, automated_harmonized_sample_ontology:=tstrsplit(single_curie, ':', fixed=TRUE, keep = 1)]
    spread_curies[automated_harmonized_sample_ontology=='CL', term:=sapply(single_curie, Term, object=cl)]
    spread_curies[automated_harmonized_sample_ontology=='EFO', term:=sapply(single_curie, my_Term, object=olsNamespace(efo))]
    spread_curies[automated_harmonized_sample_ontology=='UBERON', term:=sapply(single_curie, Term, object=uberon)]
    spread_curies[, automated_harmonized_sample_ontology:=as.factor(automated_harmonized_sample_ontology)]
  } else {
    spread_curies[, term:=sapply(single_curie, Term, object=ncit)]
  }
  # get the term_name and the ancestors
  spread_curies[, term_name:=sapply(term, function(t) unname(termLabel(t)))]
  spread_curies[, ancestors:=lapply(term, function(t) {
    ancestors <- termLabel(hierarchicalAncestors(t))
    ancestors <- ancestors[startsWith(names(ancestors), toupper(termOntology(t)))]
    c(termLabel(t), ancestors)
  })]
  # aggregate ancestors and add higher level on the unique terms
  if (curie_col == 'harmonized_sample_ontology_curie') {
    aggregated_ancestors <- spread_curies[, .(term_name=paste(sort(term_name), collapse = separator), ancestors=list(unique(unlist(Reduce(intersect, ancestors))))),
                                          by=.(automated_harmonized_sample_ontology, harmonized_sample_ontology_curie)]
    metadata_dt[aggregated_ancestors, on=.(harmonized_sample_ontology_curie), (c('automated_harmonized_sample_ontology', 'automated_harmonized_sample_ontology_term')):=mget(c('automated_harmonized_sample_ontology', 'term_name'))]
    metadata_dt[, automated_harmonized_sample_ontology:=as.factor(automated_harmonized_sample_ontology)]
    for (this_ontology in aggregated_ancestors[, levels(automated_harmonized_sample_ontology)]) {
      add_higher_order(aggregated_ancestors, rows_of_interest = aggregated_ancestors[, automated_harmonized_sample_ontology == this_ontology])
    }
    print(aggregated_ancestors[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30)), by=automated_harmonized_sample_ontology])
  } else {
    aggregated_ancestors <- spread_curies[, .(term_name=paste(sort(term_name), collapse = separator), ancestors=list(unique(unlist(Reduce(intersect, ancestors))))),
                                          by=mget(curie_col)]
    add_higher_order(aggregated_ancestors)
    print(aggregated_ancestors[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30))])
  }
  # now merge the new cols with higher level annotation to the main dt
  cols_to_add <- names(aggregated_ancestors)[startsWith(names(aggregated_ancestors), 'order')]
  curie_col_term <- sub('curie', 'term', curie_col, fixed = TRUE)
  metadata_dt[aggregated_ancestors, on=curie_col, (c(paste(curie_col_term, cols_to_add, 'unique', sep = '_'), 'ancestors')):=mget(c(cols_to_add, 'ancestors'))]

  # now add the higher level while considering all entries, also check for inconsistencies
  if (curie_col == 'harmonized_sample_ontology_curie') {
    # check for incosistencies START
    inconsistent_dt <- metadata_dt[, .(EpiRR, project, harmonized_biomaterial_type, harmonized_sample_ontology_curie, automated_harmonized_sample_ontology, harmonized_cell_type, harmonized_tissue_type, harmonized_cell_line, automated_harmonized_sample_ontology_term)]
    inconsistent_dt[harmonized_biomaterial_type %in% c('primary cell', 'primary cell culture'), sample_ontology_term_manual:=harmonized_cell_type]
    inconsistent_dt[harmonized_biomaterial_type == 'primary tissue', sample_ontology_term_manual:=harmonized_tissue_type]
    inconsistent_dt[harmonized_biomaterial_type == 'cell line', sample_ontology_term_manual:=harmonized_cell_line]
    if (inconsistent_dt[, any(automated_harmonized_sample_ontology_term != sample_ontology_term_manual, na.rm = TRUE)]) {
      print('Some automatically extracted sample ontology terms are not in line with the manual free text:')
      print(inconsistent_dt[sample_ontology_term_manual != automated_harmonized_sample_ontology_term, .(n_EpiRR=uniqueN(EpiRR), n_project=uniqueN(project)), by=.(harmonized_biomaterial_type, harmonized_sample_ontology_curie, automated_harmonized_sample_ontology, harmonized_cell_type, harmonized_tissue_type, harmonized_cell_line, automated_harmonized_sample_ontology_term, sample_ontology_term_manual)])
    }
    # check for incosistencies END
    # check if manual anno in ancestors START
    non_matching_manual <- metadata_dt[!mapply(function(term, a) term %in% a, harmonized_sample_ontology_intermediate, ancestors)]
    if(nrow(non_matching_manual) > 0){
      print('Some harmonized_sample_ontology_intermediate are not in the list of ancestors:')
      print(non_matching_manual[, .(n_EpiRR = uniqueN(EpiRR), n_project = uniqueN(project)), by=.(harmonized_biomaterial_type, harmonized_sample_ontology_curie, harmonized_cell_type, harmonized_tissue_type, harmonized_cell_line, harmonized_sample_ontology_intermediate)])
    }
    # check if manual anno in ancestors END
    for (this_ontology in metadata_dt[, levels(automated_harmonized_sample_ontology)]) {
      add_higher_order(metadata_dt, rows_of_interest = metadata_dt[, automated_harmonized_sample_ontology == this_ontology])
    }
    print(metadata_dt[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30)), by=automated_harmonized_sample_ontology])
  } else {
    if (curie_col == 'automated_harmonized_sample_disease_ontology_curie_ncit') {
      # check for incosistencies START
      inconsistent_dt <- metadata_dt[, .(
        EpiRR,
        project,
        harmonized_sample_disease_ontology_curie,
        automated_harmonized_sample_disease_ontology_curie_ncit,
        harmonized_sample_disease,
        automated_harmonized_sample_ontology_term = unlist(lapply(ancestors, function(a) ifelse(length(a) == 0, '', a[1]))) # this doesn't work for lists
      )]
      inconsistent_dt[grep(separator,
                           automated_harmonized_sample_disease_ontology_curie_ncit,
                           fixed = TRUE),
                      automated_harmonized_sample_ontology_term :=
                        sapply(strsplit(
                             automated_harmonized_sample_disease_ontology_curie_ncit,
                             separator,
                             fixed = TRUE
                           ), function(curies)
                             paste(sort(sapply(curies, function(curie)
                               termLabel(Term(ncit, curie)))), collapse = separator))]
      if (inconsistent_dt[, any(harmonized_sample_disease != automated_harmonized_sample_ontology_term, na.rm = TRUE)]) {
        print('Some automatically extracted sample disease terms are not in line with the manual free text:')
        print(inconsistent_dt[harmonized_sample_disease != automated_harmonized_sample_ontology_term, .(n_EpiRR = uniqueN(EpiRR), n_project =
                                                                uniqueN(project)), by = .(
                                                                  harmonized_sample_disease_ontology_curie,
                                                                  automated_harmonized_sample_disease_ontology_curie_ncit,
                                                                  harmonized_sample_disease,
                                                                  automated_harmonized_sample_ontology_term
                                                                )])
      }
      # check for incosistencies END
      # check if manual anno in ancestors START
      non_matching_manual <- metadata_dt[!mapply(function(term, a) term %in% a, harmonized_sample_disease_intermediate, ancestors)]
      if(nrow(non_matching_manual)>0){
        print('Some harmonized_sample_disease_intermediate are not in the list of ancestors:')
        print(non_matching_manual[, .(n_EpiRR = uniqueN(EpiRR), n_project =
                                        uniqueN(project)), by=.(harmonized_sample_disease_ontology_curie, automated_harmonized_sample_disease_ontology_curie_ncit, harmonized_sample_disease, harmonized_sample_disease_intermediate, harmonized_sample_disease_high)])
      }
      # check if manual anno in ancestors END
    } else if (curie_col == 'automated_harmonized_donor_health_status_ontology_curie_ncit') {
      # check for incosistencies START
      inconsistent_dt <- metadata_dt[, .(
        EpiRR,
        project,
        harmonized_donor_health_status_ontology_curie,
        automated_harmonized_donor_health_status_ontology_curie_ncit,
        harmonized_donor_health_status,
        automated_donor_health_status_ontology_term = unlist(lapply(ancestors, function(a) ifelse(length(a) == 0, '', a[1])))
      )]
      inconsistent_dt[grep(separator,
                           automated_harmonized_donor_health_status_ontology_curie_ncit,
                           fixed = TRUE),
                      automated_donor_health_status_ontology_term :=
                        sapply(strsplit(
                          automated_harmonized_donor_health_status_ontology_curie_ncit,
                          separator,
                          fixed = TRUE
                        ), function(curies)
                          paste(sort(sapply(curies, function(curie)
                            termLabel(Term(ncit, curie)))), collapse = separator))]
      if (inconsistent_dt[, any(harmonized_donor_health_status != automated_donor_health_status_ontology_term, na.rm = TRUE)]) {
        print('Some automatically extracted donor disease ontology terms are not in line with the manual free text:')
        print(inconsistent_dt[harmonized_donor_health_status != automated_donor_health_status_ontology_term, .(n_EpiRR = uniqueN(EpiRR), n_project =
                                                                             uniqueN(project)), by = .(
                                                                               harmonized_donor_health_status_ontology_curie,
                                                                               automated_harmonized_donor_health_status_ontology_curie_ncit,
                                                                               harmonized_donor_health_status,
                                                                               automated_donor_health_status_ontology_term
                                                                             )])
      }
      # check for incosistencies END
    }
    add_higher_order(metadata_dt)
    print(metadata_dt[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30))])
  }
  # remove ancestors and rename cols
  metadata_dt[, ancestors:=NULL]
  cols_to_rename <- names(metadata_dt)[startsWith(names(metadata_dt), 'order')]
  setnames(metadata_dt, cols_to_rename, paste(curie_col_term, cols_to_rename, sep = '_'))
}


names(metadata_dt) <- gsub('(ncit_)?order_15', 'high_order', names(metadata_dt))
names(metadata_dt) <- gsub('(ncit_)?order_30', 'intermediate_order', names(metadata_dt))

higher_annos <- names(metadata_dt)[(endsWith(names(metadata_dt), 'order')|endsWith(names(metadata_dt), 'order_unique')) & !startsWith(names(metadata_dt), 'automated_')]
setnames(metadata_dt, higher_annos, paste('automated', higher_annos, sep='_'))

fwrite(metadata_dt, file = main_csv)
