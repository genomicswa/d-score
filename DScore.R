# Generate the aggr-compatible Categories csv file after cell de-hashing.

# Try to load Seurat, required for all D-Score functions.
if ( !require(Seurat) ) {
  print( "ERROR: D-Score requires the Seurat package to be installed. Please install Seurat first then load D-Score." )
  return(0)
}

###########################################################################

# A one-size-fits-all function to (i) return a SeuratObject, and (ii) write csv file of hashtag assignments.
dscore<-function( 
  path=NA,
  seurat.obj=NA,
  seurat.fbassay="HTO",
  features=c("SBO01", "SBO02", "SBO03", "SBO04", "SBO05"),
  algorithm="dscore",
  cat.filename=NA,
  cat.name="Cell Assignment",
  datasetnumber=1 ) {

  algorithms.all<-c("seurat-default","dscore")

  # D-Score values for minreads and "next largest" dscore threshold
  d.minreads<-50
  d.threshold<-3

  # Helper function: Vector-based d-score calc
  aux_calcD<-function(x) {
    total<-sqrt(sum(x^2)+0.00001)
    x^2/total
  }
  # Helper function: find max D from a given row
  aux_getMaxSBO<-function(x,rows) { rows[which( x==max(x) )][1] }

  # Mode 1: Load data from CR output
  tenx.hashtag<-NA
  if ( !is.na(path) ) {

    print("Loading 10X data from path")
    
    tenx.data <- Read10X(path)
    tenx.umis <- tenx.data$'Gene Expression'
    tenx.htos <- tenx.data$'Antibody Capture'

    # Check features are in assay - only use the intersect of 
    assay.features<-rownames(tenx.htos)
    tenx.htos <-  tenx.htos[intersect(assay.features, features),]

    # Subset RNA and HTO counts by joint cell barcodes
    joint.bcs <- intersect(colnames(tenx.umis), colnames(tenx.htos))
    tenx.umis <- tenx.umis[, joint.bcs]
    tenx.htos <- as.matrix(tenx.htos[, joint.bcs])

    # Setup Seurat object
    tenx.hashtag <- CreateSeuratObject(counts = tenx.umis)
    # Add HTO data as a new assay independent from RNA
    tenx.hashtag[[seurat.fbassay]] <- CreateAssayObject(counts = tenx.htos)
  }

  # Mode 2: Given a Seurat Object, calculate (and add) the D-Score
  else if ( class(seurat.obj) == "Seurat" ) {
    # Simply assign the given object to tenx.hashtag
    tenx.hashtag<-seurat.obj

    # Check if the seurat.fbassay assay exist?
    if ( !is.element(seurat.fbassay, names(tenx.hashtag)) ) {
      stop("Cannot find assay called '",seurat.fbassay,"'. Please check the Seurat Object contains feature barcoding information (use 'CreateAssayObject' to create assay), and set the name of the FB assay in 'seurat.fbassay'.\nAvailable assays: ", paste( names(tenx.hashtag), collapse =", "))
    }

  } else {  # Throw error (need path or seurat obj)
    stop("Require cellranger output directory (path) or a Seurat Object (seurat.obj). Cannot continue, returning..")
  }


  ## Check overlap between features and Feature Barcoding assay
  # In assay but not in features list:
  check1<-setdiff(rownames(tenx.hashtag[[seurat.fbassay]][]), features)
  if (length(check1)>0) {
    print( paste0("WARNING: Following feature(s) found in assay but not in supplied features list - continuing but please check if they should be in 'features': ", paste(check1, collapse=", ") ))
  }
  # In features but not in assay:
  check2<-setdiff(features, rownames(tenx.hashtag[[seurat.fbassay]][]))
  if (length(check2)>0) {
    print( paste0("WARNING: Following feature(s) supplied but not found in assay - continuing but please remove from 'features': ", paste(check2, collapse=", ") ))
  }


  ## TODO: Check the UMI and HTO cell IDs match?
  

  # Run the required algorithm: D-Score (default) or Seurat's HTODemux

  if ( algorithm == algorithms.all[1] ) {
    # Seurat default algorithm
    # CLR normalise and get Seurat categories
    tenx.hashtag <- NormalizeData(tenx.hashtag, assay = seurat.fbassay, normalization.method = "CLR")
    tenx.hashtag <- HTODemux(tenx.hashtag, assay = seurat.fbassay, positive.quantile = 0.99)

    # Standardise the output and specify doublets - nothing to do here, just use hash.ID

  } else if ( algorithm == algorithms.all[2] ) {

    # Use the polarisation D-score algorithm
    raw.counts<-GetAssayData(object = tenx.hashtag, assay=seurat.fbassay, slot = "counts")
    Dscores<-t(apply( raw.counts, 2, aux_calcD ))

    ## Create an intermediate dataframe with all the necessary info.
    Ddf<-data.frame(
      maxD.SBO=apply( Dscores, 1, aux_getMaxSBO, features ),
      maxD=apply(  Dscores, 1, function(x) max(x)       ),
      nextD=apply( Dscores, 1, function(x) sort(x,T)[2] ),
      restD=apply( Dscores, 1, function(x) sum(sort(x,T)[-1]) ),
      nCount_HTO=tenx.hashtag[['nCount_HTO']]
    )
    # Calc ratio of largest and next largest dscores
    Ddf$Rn<-Ddf$maxD/(Ddf$nextD+0.00001)

    # Set filters: Min reads and the D-Score ratio.
    Ddf$filter1<- Ddf$nCount_HTO>d.minreads
    Ddf$filter2<- Ddf$Rn>=d.threshold

    # Assign the classification
    Ddf$hash.ID<-ifelse( Ddf$filter1 & Ddf$filter2, Ddf$maxD.SBO, "Unassigned" )
    tenx.hashtag[["hash.ID"]]<-Ddf$hash.ID
    tenx.hashtag[["dscore"]]<-Ddf$hash.ID
    rm(Ddf)

  } # else error


  # If a category filename has been given, create the output object and write to file.
  if (!is.na(cat.filename)) {
    dehash.df<-data.frame(
      Barcodes=rownames(tenx.hashtag[[]]),
      Assignment=tenx.hashtag$hash.ID
    )
    rownames(dehash.df)<-NULL
    colnames(dehash.df)<-c("Barcodes",cat.name)

    # Correct the barcodes with the right GEM ID suffix
    dehash.df$Barcodes<-gsub("-\\d$", paste0("-",datasetnumber), dehash.df$Barcodes, perl=T)

    # Write out to csv - append if dataset number >1
    write.table(dehash.df, file=cat.filename, sep=",", append=(datasetnumber>1), row.names=F, col.names=(datasetnumber==1))
  }

  # Some stats
  print("Summary of Assignments")
  print(table(tenx.hashtag[['hash.ID']]))
  if ( algorithm == algorithms.all[2] ) {
    print(paste0("D-Score assignment rate: ",round((nrow(tenx.hashtag[[]])-sum(tenx.hashtag[['hash.ID']]=="Unassigned"))/nrow(tenx.hashtag[[]]),4)*100,"%"))
  }

  # Return the seurat object
  return(tenx.hashtag)
}
