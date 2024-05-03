checkMetaData <- function(meta, counts){
  # If the row names of the meta data and column names of the count data do not match perfectly (exact same names in the exact same order), the downstream analysis will have errors.

  # Check if the ow names of the meta data and column names of the count data are identical
  if(identical(row.names(meta),colnames(counts))){
    message("Data matches!") # if they are identical, we print "data matches" and dont do anything else
  } else{ # otherwise...

    # if there are sample names in one set that are not in the other, this set of if statements will notify you and print any mismatches.
    if(any(!row.names(meta)%in%colnames(counts)) |
       any(!colnames(counts)%in%row.names(meta))){

      if(any(!row.names(meta)%in%colnames(counts))){
        warning("There are sample names in the meta data without a match in the count data:")
        warning(row.names(meta)[!row.names(meta)%in%colnames(counts)])
      }
      if(any(!colnames(counts)%in%row.names(meta))){
        warning("There are sample names in the count data without a match in the meta data!")
        warning(colnames(counts)[!colnames(counts)%in%row.names(meta)])
      }

    } else { # otherwise, the two sets are not identical because they are in different orders. This statement will reorder the meta data to match the order of the count data.
      message("Reordering meta data...")
      meta[order(match(rownames(meta), colnames(counts))), , drop = FALSE] -> meta
    }

  }
  return(meta) # this returns the meta data. If it was reordered, it returns the reordered object. Otherwise it is the same as the input.
}
