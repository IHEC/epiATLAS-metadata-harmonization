#Testing first CL ontologies
renv::restore()

library(rols)
library(data.table)
library(ggplot2)

#loading IHEC df
original_dt <- fread('../IHEC_metadata_harmonization.v0.9.ncit.csv')
separator <- '::'

###Ontologies
#cl ontology
cl <- Ontology("cl")
efo <- Ontology("efo")
uberon <- Ontology("uberon")
ncit <- Ontology("ncit")

# this adds one column per entry in ls to the main_dt containing higher_level annotations with overall at most l terms (for each l in ls)
add_higher_order <- function(main_dt, rows_of_interest=rep(TRUE, times=main_dt[, .N]), ls = c(15, 30)){
  
  ls <- sort(ls, decreasing = TRUE)
  
    # df_ihec[aggregated_ancestors, on=.(sample_ontology_curie), sample_ontology_ancestors:=ancestors]
    # ancestor_occurences <- as.data.table(df_ihec[, sort(table(unlist(sample_ontology_ancestors)))])
    # setnames(ancestor_occurences, c('V1', 'N'), c('term', 'counts'))
    # ancestor2occurences <- ancestor_occurences[, counts]
    # names(ancestor2occurences) <- ancestor_occurences[, term]
    # cutoff <- 100
    # ggplot(ancestor_occurences[counts >= cutoff], aes(y = reorder(term, counts), x = counts)) + geom_bar(stat = 'identity') + labs(y = 'term', title=paste('ancestor terms with >=', cutoff, 'occurences in all sample_ontology_entries'))
    # ggplot(ancestor_occurences, aes(x = counts)) + geom_histogram() + labs(title=paste('occurence distribution of terms in all sample_ontology_entries'))
    
    
    
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

# now go trough the curie cols ----
metadata_dt <- copy(original_dt)

ncit_cols <- names(original_dt)[endsWith(names(original_dt), '_ncit')]
for (curie_col in c('sample_ontology_curie', ncit_cols)){
  message(curie_col)
  spread_curies <- original_dt[, .(single_curie=unlist(tstrsplit(get(curie_col), separator, fixed=TRUE))), by=mget(curie_col)]
  if (curie_col == 'sample_ontology_curie') {
    spread_curies[, sample_ontology:=tstrsplit(single_curie, ':', fixed=TRUE, keep = 1)]
    spread_curies[sample_ontology=='CL', term:=sapply(single_curie, term, object=cl)]
    spread_curies[sample_ontology=='EFO', term:=sapply(single_curie, term, object=efo)]
    spread_curies[sample_ontology=='UBERON', term:=sapply(single_curie, term, object=uberon)]
    spread_curies[, sample_ontology:=as.factor(sample_ontology)]
  } else {
    spread_curies[, term:=sapply(single_curie, term, object=ncit)]
  }
  spread_curies[, term_name:=sapply(term, function(t) unname(termLabel(t)))]
  spread_curies[, ancestors:=lapply(term, function(t) c(unname(termLabel(t)), unname(termLabel(unique(ancestors(t))))))]
  if (curie_col == 'sample_ontology_curie') {
    aggregated_ancestors <- spread_curies[, .(term_name=paste(sort(term_name), collapse = separator), ancestors=list(unique(unlist(Reduce(intersect, ancestors))))), 
                                          by=.(sample_ontology, sample_ontology_curie)]
    metadata_dt[aggregated_ancestors, on=.(sample_ontology_curie), (c('sample_ontology', 'sample_ontology_term')):=mget(c('sample_ontology', 'term_name'))]
    for (this_ontology in aggregated_ancestors[, levels(sample_ontology)]) {
      add_higher_order(aggregated_ancestors, rows_of_interest = aggregated_ancestors[, sample_ontology == this_ontology])
    }
    print(aggregated_ancestors[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30)), by=sample_ontology])
  } else {
    aggregated_ancestors <- spread_curies[, .(term_name=paste(sort(term_name), collapse = separator), ancestors=list(unique(unlist(Reduce(intersect, ancestors))))), 
                                          by=mget(curie_col)]
    add_higher_order(aggregated_ancestors)
    print(aggregated_ancestors[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30))])
  }
  cols_to_add <- names(aggregated_ancestors)[startsWith(names(aggregated_ancestors), 'order')]
  cure_col_term <- sub('curie', 'term', curie_col, fixed = TRUE)
  metadata_dt[aggregated_ancestors, on=curie_col, (c(paste(cure_col_term, cols_to_add, 'unique', sep = '_'), 'ancestors')):=mget(c(cols_to_add, 'ancestors'))]
  if (curie_col == 'sample_ontology_curie') {
    for (this_ontology in metadata_dt[, levels(sample_ontology)]) {
      add_higher_order(metadata_dt, rows_of_interest = metadata_dt[, sample_ontology == this_ontology])
    }
    print(metadata_dt[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30)), by=sample_ontology])
    # add the manual annotations
    jeffrey_anno <- fread('../IHEC_metadata_harmonization.v0.8.JeffreyHyacinthe.csv', header = TRUE, select = c('EpiRR', 'category'))
    jonathan_anno <- fread('../IHEC_metadata_harmonization.v0.9.JonathanSteif.csv', header = TRUE, select = c('EpiRR', 'aggregate_group'))
    stopifnot(all(sapply(list(jonathan_anno[, EpiRR], jeffrey_anno[, EpiRR]), FUN = identical, metadata_dt[, EpiRR])))
    metadata_dt[jeffrey_anno, on=.(EpiRR), sample_ontology_term_high_order_JeffreyHyacinthe:=category]
    metadata_dt[jonathan_anno, on=.(EpiRR), sample_ontology_term_high_order_JonathanSteif:=aggregate_group]
    manual_cutoff <- fread('../manual_higher_order_annotation.tsv')
    metadata_dt[manual_cutoff, on=.(sample_ontology_term), sample_ontology_term_high_order_manual:=prelim_classification_v0.9]
  } else {
    add_higher_order(metadata_dt)
    print(metadata_dt[, .(unique_terms_higher_order=uniqueN(order_15), unique_terms_intermediate_order=uniqueN(order_30))])
  }
  metadata_dt[, ancestors:=NULL]
  cols_to_rename <- names(metadata_dt)[startsWith(names(metadata_dt), 'order')]
  setnames(metadata_dt, cols_to_rename, paste(cure_col_term, cols_to_rename, sep = '_'))
}


# metadata_dt[aggregated_ancestors, on=.(sample_ontology_curie), sample_ontology_ancestors:=ancestors]


# metadata_dt[, sapply(.SD, function(x) paste(sort(unique(x)), collapse = separator), simplify=FALSE), by=eval(names(metadata_dt)[-grep('order', names(metadata_dt))])]
metadata_melt <- melt(metadata_dt, id.vars = c('EpiRR', 'sample_ontology'), measure.vars = names(metadata_dt)[grep('order', names(metadata_dt))], variable.name = 'annotation_col', value.name = 'ontology_term', variable.factor = FALSE)
metadata_melt[, unique_terms_only:=endsWith(annotation_col, '_unique')]
metadata_melt[grepl('order_15', annotation_col), annotation:='high order']
metadata_melt[grepl('order_30', annotation_col), annotation:='intermediate order']
metadata_melt[, 'ontology column':=gsub("_order.*$","", annotation_col)]

for (ontology_column in metadata_melt[ , unique(`ontology column`)]) {
  if (ontology_column == 'sample_ontology_term') {

      ggplot(metadata_melt[ontology_column == `ontology column` &
                             !is.na(ontology_term)], aes(
                               y = reorder(ontology_term, ontology_term, length),
                               fill = unique_terms_only
                             )) + geom_bar(position = 'dodge', color = 'black') + facet_wrap(
                               . ~ annotation + sample_ontology,
                               scales = 'free',
                               nrow = 2
                             ) + labs(title= ontology_column)

  } else {

      ggplot(metadata_melt[ontology_column == `ontology column` &
                             !is.na(ontology_term)], aes(
                               y = reorder(ontology_term, ontology_term, length)
                             )) + geom_bar(position = 'dodge', color = 'black') + facet_wrap(
                               . ~ annotation + unique_terms_only,
                               scales = 'free',
                               nrow = 2
                             ) + labs(title= ontology_column)

  }
}

names(metadata_dt) <- gsub('(ncit_)?order_15', 'high_order', names(metadata_dt))
names(metadata_dt) <- gsub('(ncit_)?order_30', 'intermediate_order', names(metadata_dt))

fwrite(metadata_dt, file = '../IHEC_metadata_harmonization.v0.9.extended.csv')
