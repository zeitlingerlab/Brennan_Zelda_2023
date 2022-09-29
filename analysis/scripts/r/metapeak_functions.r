#Melanie Weilert
library(viridis)
library(DECIPHER)
library(testit)
#Functions:
#matrix_obtain [subset of matrix_obtain]
#multiplot
#matrix_format
############################################################################################################################

#Format and Develop Data from Peak Calls

#Purpose: Imports the bigiwg +/- strands, reads the peaks of interest, obtains a matrix [peaksxselectedbp] with selected parameters.
#Inputs:
#   pos_file: filepath string of the sequencing data aligned to the positive strand (bw)
#   neg_file: filepath string of the sequencing data aligned to the negative strand (bw)
#   region_path: filepath string of peak regions in the form of an .rds file (seqnames,start,end,width,strand,sig) 
#               [start and end must be same at 0bp location]
#   peak_count: range of desired peaks from this .rds file, it can only handle less than 5000 peaks
#               [must pre-sort peaks to desired order prior to this]
#   background: desired cutoff for background, defaults to NA

matrix_format<-function(pos_file, neg_file, region_path, peak_count, upstream=100, downstream=100, background=NA){
  
  #Check for peak_count limit.
  if (length(peak_count)>5010) stop("Peak_count cannot process numbers greater than 5000.")
  
  #import data
  pos <- import(pos_file)
  neg <- import(neg_file)
  
  #parse down peak regions
  region<-readRDS(region_path)
  region_parsed<-region[peak_count]
  
  fmatrix <- matrix_obtain(pos_file = pos, neg_file = neg, gr_object = region_parsed, background=background, upstream=upstream, downstream=downstream)
  
}

#####################################################################################################################################
# Function to Read and Format .BW into exo-matrix
# Starting with normalized BigWig files, we must format the data, taking the GRanges tss locations that are already loaded. 
# The options for this function are as follows:

# - pos_file: the sequencing data aligned to the positive strand
# - neg_file: the sequencing data aligned to the negative strand
# - gr_object: GRanges object that contains peaks of interest and proposed tss
# - upstream: how many base pairs upstream (left of tss) should be incorperated *default is 100
# - downstream: how many base pairs downstream (right of tss) should be incorperated *default is 100
# - background: predefined background cutoff which is eliminated from the pos and neg strand.

# The function returns a gRanges object which includesa list plus and minus strand matrices of every peak.

matrix_obtain <- function(pos_file, neg_file, gr_object, upstream=100, downstream=100, background=NA) {
  
  # if there is a background cutoff
  if (is.na(background)==F){
    # remove background from pos_file and neg_file
    pos_file$score <- ifelse(pos_file$score<background, 0, pos_file$score - background)
    neg_file$score <- ifelse(abs(neg_file$score)<background, 0, neg_file$score + background)
  }
  
  # convert to an Rle list with positive and negative factor levels
  RleList <- list(pos= coverage(pos_file, weight=pos_file$score), neg =  coverage(neg_file, weight=neg_file$score))
  
  peak_matrix <- exo_metapeak_matrix(regions.gr=gr_object, sample=RleList, upstream=upstream, downstream=downstream)
  
  # return data frames
  peak_matrix
}


#############################################################################################################################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#<http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/>

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#############################################################################################################################
#Normalize reads of ChIP-nexus samples: NORMALIZED OVER TOTAL READS in READS PER MILLION

generate_normalized_coverage <- function(gr, chip_type="chipnexus", total_reads=NA, path=NA){ #Generate normalized coverage: 
  message("processing ", path)
  if(is.na(total_reads)){total_reads<-length(gr)}
  if(chip_type=="chipnexus"){
    gr <- resize(gr, 1, "start")
    gr_p <- gr[strand(gr) == "+"| strand(gr) == "*"]
    gr_n <- gr[strand(gr) == "-"]
    cov_list <- list(pos = coverage(gr_p) / total_reads * 1000000, 
                     neg = coverage(gr_n)/ total_reads * 1000000 * (-1))
    if(!is.na(path)){saveRDS(cov_list, path)}
    return(cov_list)
  }
  else{	
    cov <- coverage(gr) / total_reads * 1000000
    return(cov)
  }
}

#############################################################################################################################
# Plot Multiple Metapeaks on a Single GGPLOT
#
# Inputs: A list with the layer format: [Layer 1] List of sample names [Layer 2] pos=positive metapeak vector, neg=negative metapeak vector
# For more detailed structure information, go to the tutorial at the Zeitlinger Lab Drive under Analysis/ChIP-nexus/overlay_multiple_metagene_plots
# Outputs: A ggplot that plots multiple metapeaks on the sample plot

plot_multiple_metapeaks<-function(sample_metapeak_list, palette=viridis(length(sample_metapeak_list))){
  
  #Convert to data frame
  df <- ldply(sample_metapeak_list, data.frame)
  df$.id<-factor(df$.id)
  
  x<-ggplot(df, aes(x=position, ymin=0))+
    geom_line(aes(y=pos, color=.id))+
    geom_line(aes(y=neg, color=.id))+
    geom_area(aes(y=pos, fill=.id), alpha=.2, position="identity")+
    geom_area(aes(y=neg, fill=.id), alpha=.2, position="identity")+
    geom_vline(xintercept = 0, color="gray50", linetype="dashed")+
    scale_fill_manual(values=palette, name="Samples")+
    scale_color_manual(values=palette, name="Samples")+
    xlab("BP Position")+
    ylab("Signal")+
    ggtitle("ChIP Signals between samples")+
    theme_classic()+
    theme(text=element_text(size=18, family="Times"), plot.title=element_text(hjust = .5))
  print(x)
  return(x)
  
} 

########################################################################################################

#Plot heatmap
normalize_standard_matrix<-function(matrix, removal_threshold, normalize_threshold){
  
  max.per.gene.pos <- apply(matrix, 1, function(x){quantile(x, normalize_threshold)})
  min.per.gene.pos <- apply(matrix, 1, function(x){quantile(x, removal_threshold)})

  matrix.p <- matrix
  matrix.p [matrix.p <= min.per.gene.pos] <- NA #Remove all values that are below 50th percentile
  matrix.p <- pmin(matrix.p / max.per.gene.pos,1) #Re-normalize values that remain.
  
  matrix.p
}

format_standard_matrix_for_heatmap<-function(matrix, downstream, upstream){
  mat<-as.data.frame(matrix)
  colnames(mat)<-as.numeric((-upstream+1):downstream)
  mat$region<-1:nrow(mat)
  df<-mat %>% data.table %>% melt.data.table(mat, id.vars = c("region"), variable.name="position", 
                                             value.name="signal", measure.vars = paste0(as.numeric((-upstream+1):downstream)))
  #Remove NA values
  df<-df[!is.na(df$signal),]
  df$position<-as.numeric(as.character(df$position))
  return(df)
}


#Associated functions
normalize_each_matrix<-function(matrix_list, removal_threshold, normalize_threshold){
  
  max.per.gene.pos <- apply(matrix_list$pos, 1, function(x){quantile(x, normalize_threshold)})
  min.per.gene.pos <- apply(matrix_list$pos, 1, function(x){quantile(x, removal_threshold)})
  max.per.gene.neg <- apply(matrix_list$neg, 1, function(x){quantile(x, normalize_threshold)})
  min.per.gene.neg <- apply(matrix_list$neg, 1, function(x){quantile(x, removal_threshold)})
  
  matrix.p <- matrix_list$pos
  matrix.p [matrix.p <= min.per.gene.pos] <- NA #Remove all values that are below 50th percentile
  matrix.p <- pmin(matrix.p / max.per.gene.pos,1) #Re-normalize values that remain.
  
  matrix.n <- matrix_list$neg
  matrix.n [matrix.n <= min.per.gene.neg] <- NA
  matrix.n <- pmin(matrix.n / max.per.gene.neg,1)
  list(pos=matrix.p, neg=matrix.n)
}

format_matrix_for_heatmap<-function(matrix, strand, downstream, upstream){
  mat<-as.data.frame(matrix)
  colnames(mat)<-as.numeric((-upstream+1):downstream)
  mat$region<-1:nrow(mat)
  mat$original_idx<-rownames(mat)
  df<-mat %>% data.table %>% melt.data.table(mat, id.vars = c("region", "original_idx"), variable.name="position", 
                                            value.name="signal", measure.vars = paste0(as.numeric((-upstream+1):downstream)))
  #Remove NA values
  df<-df[!is.na(df$signal),]
  df$strand<-as.character(strand)
  df$position<-as.numeric(as.character(df$position))
  return(df)
}

#Required inputs:
# regions.gr = GRanges of regions desired to plot across
# sample = list(pos="", neg="") with appropriate filepaths to the nexus bws

#Optional inputs
# order = option that can either be "sum" or "clustering" that will determine how the heatmap gets arranged. 
#         Any other value will return the order of the granges as input.
# output_file_name = option that will write the plot to a PDF if a filepath is provided
# reduce = option to reduce the GRanges to remove redundancies
# upstream/downstream = how far you want your heatmap to extend relative to the center of the provided GRanges
# alpha_value = current option to determine color balancing (in progress)
# removal_threshold = quantile value to remove background regions
# normalize_threshold = upper quantile value to keep background regions
# reverse = option to reverse strand (only useful if you are plotting reverse complement of a motif, say)

#Outputs: ggplot with heatmap of regions (xaxis=position, yaxis=regions, fill=quantile normalized chipnexus value)

plot_nexus_heatmap<-function(regions.gr, sample, title="Heatmap of ChIP-nexus Signals", order="sum", output_file_name=NA, reduce=F, 
                             upstream=50, downstream=50, 
                             alpha_value=.5, removal_threshold=.5, normalize_threshold=0.99, reverse=F, return_only_df = F){
  library(testit)
  testit::assert('Reducing the GRanges and applying a custom ordering is not compatible settings.', 
                 ifelse(reduce, ifelse((order=='sum' | order=='clustering'), TRUE, FALSE), TRUE))
  
  if(reverse){strand(regions.gr) <- ifelse(strand(regions.gr) == '+', '-', '+')}
  if(reduce){regions.gr<-GenomicRanges::reduce(regions.gr)}
  print(length(regions.gr))
  
  #Get signals across regions
  mat<-exo_metapeak_matrix(regions.gr = resize(regions.gr, 1, "start"), sample=sample, upstream=upstream, downstream=downstream)
  
  #Order region rows
  if(order=="clustering"){
    order_of_rows<-hclust(d = (dist(mat$pos+mat$neg))/2)$order
  }
  else if (order=="sum"){
    nrows<-length(mat$pos[,1])
    sum_by_rows_vec<-(apply(mat$pos, 1, sum) + abs(apply(mat$neg, 1, sum)))/2
    order_of_rows<-order(sum_by_rows_vec, decreasing=F)
    # #TODO: Fix row orders
    # print(data.frame(order=order_of_rows, sum_amt=summed_by_row[order_of_rows]))
  }
  else{order_of_rows<-1:length(regions.gr)}

  rownames(mat$pos)<-1:nrow(mat$pos)
  rownames(mat$neg)<-1:nrow(mat$neg)

  mat$pos<-mat$pos[order_of_rows, ]
  mat$neg<-mat$neg[order_of_rows, ]
  
  #Normalize matrix
  norm_mat_matrix<-normalize_each_matrix(matrix_list = mat, removal_threshold =  removal_threshold, normalize_threshold = normalize_threshold)
  
  #Format matrix
  mat_pos.df<-format_matrix_for_heatmap(matrix = norm_mat_matrix$pos, strand = "pos", downstream = downstream, upstream = upstream)
  mat_neg.df<-format_matrix_for_heatmap(matrix = norm_mat_matrix$neg, strand = "neg", downstream = downstream, upstream = upstream)
  mat_final_df<-rbind(mat_pos.df, mat_neg.df)
  
  if(return_only_df){
    return(mat_final_df)
  }
  else{
    #Plot matrix
    g<-ggplot()+
      geom_tile(data=mat_pos.df, aes(x=position, y=region, fill=signal))+
      geom_tile(data=mat_neg.df, aes(x=position, y=region, fill=-signal), alpha=alpha_value)+
      ggtitle(title, subtitle = length(regions.gr))+
      scale_x_continuous("Position (bp)", expand = expand_scale(mult = c(0, 0)))+
      scale_y_reverse("Regions", expand = expand_scale(mult = c(0, 0)))+
      scale_fill_gradientn(colors = c('#08306b', '#08519c' , '#2171b5', '#4292c6', '#6baed6', 
                                      '#9ecae1', '#c6dbef','#deebf7', '#f7fbff', "white", '#fff5f0',
                                      '#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'), name = "norm. signal") + 
      # scale_fill_gradient2(low="#08306b", mid="white", high="#67000d")+
      theme_classic()+
      theme(text=element_text(size=14), legend.position = "none", panel.grid = element_blank(), panel.border = element_blank())
    
    if(!is.na(output_file_name)){ggsave(filename = paste(getwd(), output_file_name, sep="/"), plot = g, width = 12, height=12)}
    return(g)
  }
}

#Required inputs:
# regions.gr = GRanges of regions desired to plot across
# sample filepath 

#Optional inputs
# order = option that can either be "sum" or "clustering" that will determine how the heatmap gets arranged. 
#         Any other value will return the order of the granges as input.
# output_file_name = option that will write the plot to a PDF if a filepath is provided
# reduce = option to reduce the GRanges to remove redundancies
# upstream/downstream = how far you want your heatmap to extend relative to the center of the provided GRanges
# removal_threshold = quantile value to remove background regions
# normalize_threshold = upper quantile value to keep background regions
# reverse = option to reverse strand (only useful if you are plotting reverse complement of a motif, say)

#Outputs: ggplot with heatmap of regions (xaxis=position, yaxis=regions, fill=quantile normalized chipnexus value)

plot_standard_heatmap<-function(regions.gr, sample, title="Heatmap of ChIP-seq Signals", order="sum", output_file_name=NA, reduce=F, 
                             upstream=50, downstream=50, 
                             removal_threshold=.5, normalize_threshold=0.99, reverse=F, return_only_df = F){
  library(testit)
  testit::assert('Reducing the GRanges and applying a custom ordering is not compatible settings.', 
                 ifelse(reduce, ifelse((order=='sum' | order=='clustering'), TRUE, FALSE), TRUE))
  
  if(reverse){strand(regions.gr) <- ifelse(strand(regions.gr) == '+', '-', '+')}
  if(reduce){regions.gr<-GenomicRanges::reduce(regions.gr)}
  print(length(regions.gr))
  
  #Get signals across regions
  mat<-standard_metapeak_matrix(regions.gr = resize(regions.gr, 1, "start"), sample=sample, upstream=upstream, downstream=downstream)
  
  #Order region rows
  if(order=="clustering"){
    order_of_rows<-hclust(d = dist(mat))$order
  }
  else if (order=="sum"){
    nrows<-length(mat[,1])
    sum_by_rows_vec<-apply(mat, 1, sum) 
    order_of_rows<-order(sum_by_rows_vec, decreasing=F)
    # #TODO: Fix row orders
    # print(data.frame(order=order_of_rows, sum_amt=summed_by_row[order_of_rows]))
  }
  else{order_of_rows<-1:length(regions.gr)}
  
  rownames(mat)<-1:nrow(mat)
  mat<-mat[order_of_rows, ]

  #Normalize matrix
  norm_mat_matrix<-normalize_standard_matrix(mat, removal_threshold =  removal_threshold, normalize_threshold = normalize_threshold)
  
  #Format matrix
  mat_final_df<-format_standard_matrix_for_heatmap(matrix = norm_mat_matrix, downstream = downstream, upstream = upstream)

  if(return_only_df){
    return(mat_final_df)
  }
  else{
    #Plot matrix
    g<-ggplot()+
      geom_tile(data=mat_final_df, aes(x=position, y=region, fill=signal))+
      ggtitle(title, subtitle = length(regions.gr))+
      scale_x_continuous("Position (bp)", expand = expand_scale(mult = c(0, 0)))+
      scale_y_reverse("Regions", expand = expand_scale(mult = c(0, 0)))+
      scale_fill_gradientn(colors = c("white", '#fff5f0', '#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'), name = "norm. signal") + 
      # scale_fill_gradient2(low="#08306b", mid="white", high="#67000d")+
      theme_classic()+
      theme(text=element_text(size=14), legend.position = "none", panel.grid = element_blank(), panel.border = element_blank())
    
    if(!is.na(output_file_name)){ggsave(filename = paste(getwd(), output_file_name, sep="/"), plot = g, width = 12, height=12)}
    return(g)
  }
}
########################################################################################################

#Function to plot sequence

# Required Inputs:
#   gr = regions you with to plot sequence across
#   title = title of sequence plot
#   genome = BSgenome Object to provide the sequences across the input ranges
#   
# Optional Inputs:
#   subset = if T, will subset the plot by the name of the granges
#   cluster = if T, will cluster the GRanges before plotting to group sequences by their similarity to each other using the DECPHER package
#   window = number of bases to plot across
#   show = if T, print the plot regardless of output
# Output: list that will provide the plotting data frame, the plot object, and the DNAStringSet object of the sequences returned.


plot_sequence<-function(gr, title, genome, subset=F, cluster=F, show=T, window=201, y_axis_name = "Regions", x_axis_name = "Motif distance (bp)", cores = 8){
  
  #if(window %% 2==0){warning("Window might need to be an odd number")}
  #if(width(gr)[1] %% 2 != 0){gr[which(strand(gr)=="-")]<-IRanges::shift(gr[which(strand(gr)=="-")], shift=1)}
  # else{
  #   gr[which(strand(gr)=="+")]<-IRanges::shift(gr[which(strand(gr)=="+")], shift=1)
  # }
  
  gr<-resize(gr, 1, "start")
  gr$row<-1:length(gr)
  resized_to_desired_window<-resize(gr, width=window, fix="center")
  sequences<-getSeq(genome, resized_to_desired_window) %>% as.character
  
  #reassign to a data frame
  sequence.df<-lapply(1:length(sequences), function(x){
    split_vec<-strsplit(sequences[x], split="") %>% unlist
    df<-data.frame(nt=split_vec, position=(1:length(split_vec))-(ceiling(window/2)+1), row=x, stringsAsFactors = F)
    return(df)
  }) %>% rbindlist
  gr.df<-as.data.frame(gr)
  
  #Cluster the sequences if needed
  if(cluster){
    library(DescTools)
    seqs<-getSeq(genome, resized_to_desired_window)
    dna_dist <- DistanceMatrix(seqs, type="dist") # returns an object of class 'dist'
    dna_clust<-hclust(d = dna_dist)
    dna_order<-dna_clust$order
    gr.df$plot_row<-dna_order
    
  }
  else{gr.df$plot_row<-gr.df$row}  
  
  #Merge data together for plotting
  sequence.df<-merge(sequence.df, gr.df, by.x="row", by.y="row")
  
  #Remove N-containing sequences
  N_idx<-which(sequence.df$nt=="N")
  if(length(N_idx)>0){sequence.df<-sequence.df[-N_idx,]}
  
  #Plot
  #sequence.df$plot_row<-sequence.df$plot_row %>% as.character %>% as.integer
  sequence.df$position<-sequence.df$position %>% as.character %>% as.integer
  
  g<-ggplot(sequence.df, aes(x=position, y=plot_row, fill=nt))+
    geom_tile()+
    scale_y_continuous(name = y_axis_name, expand = c(0, 0))+
    scale_x_continuous(name = x_axis_name, expand = c(0, 0))+
    scale_fill_manual(values = c("#36982F", "#402CFD", "#FFB530", "#FC3437"), name="Nucleotide")+
    ggtitle(title)+
    theme_classic()+
    theme(text=element_text(size=14), legend.position="bottom", panel.background = element_blank())
  
  if(subset){g<-g+facet_grid(subset_name ~ ., scales="free")}
  if(show){print(g)}
  return(list(plot=g, df=sequence.df, dna=getSeq(genome, resized_to_desired_window)))
}

########################################################################################################

#Plot metapeak

#Required inputs:
# regions.gr = GRanges of regions desired to plot across 
# sample = list(pos="", neg="") with appropriate filepaths to the nexus bws

#Optional inputs
# output_file_name = option that will write the plot to a PDF if a filepath is provided
# upstream/downstream = how far you want your heatmap to extend relative to the center of the provided GRanges

#Outputs: ggplot with heatmap of regions (xaxis=position, yaxis=regions, fill=quantile normalized chipnexus value)

plot_nexus_metapeak<-function(regions.gr, sample, title="Metapeak of ChIP-nexus Signals", output_file_name=NA, 
                              upstream=50, downstream=50, return_only_df = F){
  assert("GRanges needs to be resized to 1bp (either centered or positioned on start/end site.", width(regions.gr[1])==1)
  
  #Get signals across regions
  df<-exo_metapeak(gr = regions.gr, sample=sample, upstream=upstream, downstream=downstream)
  
  #Either return only DF or plot metapeak
  if(return_only_df){
    return(df)
  }
  else{
    #Plot matrix
    g<-ggplot(df, aes(tss_distance, reads))+
      geom_area(aes(group=sample, fill = strand))+
      ggtitle(title, subtitle = length(regions.gr))+
      scale_x_continuous(name = "Position (bp)")+
      scale_y_continuous(name = "Reads")+
      scale_fill_manual(values = c('#a50f15', '#08519c'), labels = c("+", "-"), name = "strand") + 
      theme_classic()+
      theme(text=element_text(size=14))
    if(!is.na(output_file_name)){ggsave(filename = paste(getwd(), output_file_name, sep="/"), plot = g, width = 12, height=12)}
    return(g)
  }
}



########################################################################################################
#Align plots kudos = https://stackoverflow.com/questions/16255579/how-can-i-make-consistent-width-plots-in-ggplot-with-legends
########################################################################################################

same.size.ggplot <- function(vector.string.graph, # a vector of strings which correspond to Robject ggplot graphs
                             reference.string.graph, # a string of a  Robject ggplot graphs where height and/or height will be taken for reference
                             width = T, # if you wanna adapat only the width
                             height = F # if you wanna adapat only the height
) {
  
  # example: same.size.ggplot(p0rep(c("a", "b"), thre), "a30") 
  
  
  which(vector.string.graph %in% reference.string.graph)
  
  newref <- ggplotGrob(get(reference.string.graph))
  ref.width <- newref$widths
  ref.height <- newref$heights
  
  assign(reference.string.graph, newref, env = parent.frame(1))
  
  for(i in seq_along(vector.string.graph)) {
    if(vector.string.graph[i] != reference.string.graph) {
      new <- ggplotGrob(get(vector.string.graph[i]))
      if( width ) {
        new$widths <- ref.width
      }
      if( height ) {
        new$heights <- ref.height
      }
      assign(vector.string.graph[i], new, env = parent.frame(1))
    }
  }
}



  