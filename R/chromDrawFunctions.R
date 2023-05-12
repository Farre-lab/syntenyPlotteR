#' Draw Evolution Highway Plots
#'
#' This function draws Evolution Highway style plots. It requieres as input the syntenic blocks following this
#' format: chr,start,end,targetChr,targetStart,targetEnd,orient,targetSpecies separated by tabs
#' It also requieres the output file name and the range of chromosomes of the reference species.
#' Example: draw.eh(input.csv,goat,"1:29")
#'
#' @title Evolution Highway style plot
#' @param infile Path to the syntenic blocks file
#' @param output file name
#' @param chrRange range of chromosome numbers in the reference "1:29"
#' @return A pdf file with the comparative drawings
#' @export
draw.eh<-function(infile,output,chrRange) {
  outfile<-chr<-start<-end<-tarChr<-tarSt<-tarEnd<-orient<-tar<-text_size2<-NULL
  data<-read.table(infile, header=FALSE)
  colnames(data) = c("chr","start","end","tarChr","tarSt","tarEnd","orient","tar")
  data$orient = factor(data$orient, levels=c("1","-1"))
  data$text_size2=80*((data$end-data$start)/100000)

  pdf(paste0(outfile,".pdf"),width=5.5, height =10, pointsize = 10)
  for (ID in c(chrRange)) {
    #ID=chrRange
    print(ID)
    subsetChr1<-subset(data,chr==ID, select=c(chr,start,end,tarChr,tarSt,tarEnd,orient,tar,text_size2))
    min=min(subsetChr1$start)
    max=max(subsetChr1$end)
    print(ggplot2::ggplot() +
            ggplot2::geom_rect(data=subsetChr1, mapping=ggplot2::aes(xmin=start, xmax=end, ymin=0, ymax=0.5, fill=orient, group=tar), color="white",
                      alpha = 1, size = 0.1 ) +
            ggplot2::geom_rect(data=subsetChr1, ggplot2::aes(xmin=min,xmax=max,ymin=0,ymax=0.5), size=0.3, color="black", fill="NA") +
            ggplot2::facet_grid(~ tar) +
            ggplot2::coord_flip() +
            ggplot2::scale_x_reverse() +
            ggplot2::scale_y_discrete(expand=c(0,0)) +
            ggplot2::theme(
                panel.spacing.y = ggplot2::unit(c(-0.5,-0.5), "lines"),
                panel.spacing.x = ggplot2::unit(0,"lines"),
                panel.background = ggplot2::element_blank(),
                strip.background = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                legend.position="none"
           ) +
            ggplot2::scale_fill_manual(values = c("1" = "lightblue", "-1" = "lightpink")) +
            ggplot2::geom_text(data=subsetChr1,ggplot2::aes(x=start+(end-start)/2,y=0+(0.5-0)/2,label=tarChr,size=text_size2))
    )


  }
  dev.off()
}


#' Draw Pairwise Synteny Plots
#'
#' This function draws pairwise synteny plots.
#' It requires:
#' 1. the output file name;
#' 2. the format which to save the image in e.g. PDF or PNG
#' 3. a file with all chromosomes, chromosome lengths, and species identifiers for all species in the synteny analysis in this format:
#' chromosome ID, chromosome length, species identifier
#' 4. files containing the syntenic blocks (one file per alignment, in order from most recent species alignment file to ancestor alignment file) following this format:
#' reference chromosome, reference start position, reference end position, target chromosome,
#' target start position, target end position, orient, reference species identifier, target species identifier
#'
#' Please separate files by tab and ensure any species identifiers used between length and alignment files are matching (same identifiers and caseing)
#'
#'
#' Example: draw.pairwise("outputName","pdf", "ChromosomeSizefile", "syntenyfile1", "syntenyfile2", ...)
#'
#' @title Pairwise synteny plot
#' @param output file name
#' @param output file format
#' @param ChromosomeSize file
#' @param synteny files (any number of alignment files can be entered)
#' @return A pdf file with comparative drawings
#' @export
#'
#'
draw.pairwise <- function(output,fileformat,sizefile,...){
  #The below function converts coordinates to linear genome and creates synteny polygon coordinates
  synteny.data.reframing <- function(data,tar.y,ref.y,compiled.size){
    synteny <- data.frame()
    for (i in c(1:nrow(data))){
      reference = data[i,"ref.species"]
      target = data[i,"tar.species"]
      tar_chr = data[i,"tarchr"]
      ref_chr = data[i,"refchr"]
      dir = data[i, "dir"]
      tar_sizes = compiled.size[compiled.size$species == target,]
      names(tar_sizes) <- c("tarchr","size","species","xstart","xend")
      ref_sizes = compiled.size[compiled.size$species == reference,]
      names(ref_sizes) <- c("refchr","size","species","xstart","xend")
      tar_add = tar_sizes[as.character(tar_sizes$tarchr)==as.character(tar_chr),]$xstart
      ref_add = ref_sizes[as.character(ref_sizes$refchr)==as.character(ref_chr),]$xstart
      tar_y = tar.y
      ref_y = ref.y
      tar_xstart = data[i,"tarstart"] + tar_add
      tar_xend = data[i,"tarend"] + tar_add
      ref_xstart = data[i,"refstart"] + ref_add
      ref_xend = data[i,"refend"] + ref_add

      inverted = grepl("-", dir, fixed = TRUE)
      if(inverted == TRUE){
        df = data.frame(x = c(tar_xstart, tar_xend, ref_xstart, ref_xend), y = c(tar_y, tar_y, ref_y, ref_y),
                        fill = ref_chr, group = paste0("s",i),ref = reference, tar = target)
      } else {
        df = data.frame(x = c(tar_xstart, ref_xstart, ref_xend,  tar_xend), y = c(tar_y, ref_y, ref_y, tar_y),
                        fill = ref_chr, group = paste0("s",i),ref = reference, tar = target)
      }
      synteny = rbind(synteny,df)
    }
    return(synteny)
  }
  xstart<-xend<-refchr<-tarchr<-x<-y<-group<-fill<-NULL
  sizes <-read.delim(sizefile, header=FALSE) #to be consistent with naming in EH
  names(sizes) <- c("chromosome","size","species")
  sizes$size <- as.numeric(gsub(",","",sizes$size))

  count = 0
  compiled.size <- data.frame()
  #This adds gap in between  chromosomes and convert to "linear" genome
  for(i in unique(sizes$species)){
    size.intermediate <- sizes[sizes$species == i,]
    for (x in c(1:nrow(size.intermediate))){
      #print(i)
      if (x == 1){
        total_start = 1
        total_end = size.intermediate[x, "size"]
      } else {
        total_start = total_end + 6000000
        total_end = total_start + size.intermediate[x, "size"]
      }
      size.intermediate[x,"xstart"] = total_start
      size.intermediate[x, "xend"] = total_end
    }
    compiled.size <- rbind(compiled.size,size.intermediate)
  }

  #This calculates a position for each genome on the y axis
  for(z in unique(compiled.size$species)){
    compiled.size$y[compiled.size$species == z] <- count
    count = count + 2 }

  #makes a list of all synteny files input to function
  list.of.files <-  list()
  for(i in list(...)){
    list.of.files[[i]] <- i
  }

  #for each file in the list of synteny files prepare the polygon coordinates
  listsynt <- list()
  for(i in 1:length(list.of.files)){
    num <- i
    file <- list.of.files[[num]]
    dataTMP <- read.delim(file, header=FALSE)
    data2 <-dataTMP[,c(4,5,6,1,2,3,7,8,9)]
    colnames(data2) = c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "dir", "ref.species","tar.species")
    data2$tarstart <- as.numeric(gsub(",","",data2$tarstart))
    data2$tarend <- as.numeric(gsub(",","",data2$tarend))
    data2$refstart <- as.numeric(gsub(",","",data2$refstart))
    data2$refend <- as.numeric(gsub(",","",data2$refend))
    reference <- data2[1,"ref.species"]
    target <- data2[1,"tar.species"]
    ref_y <- compiled.size[compiled.size$species == reference,"y"]
    tar_y <- compiled.size[compiled.size$species == target,"y"]
    ref_y <- ref_y[1]
    tar_y <- tar_y[1] + 0.1
    x <- synteny.data.reframing(data2,tar_y,ref_y,compiled.size)
    x$fill <- as.factor(x$fill)
    listsynt[[i]] <- x
  }

  #ensure chromosomes in chromosome column are factors
  compiled.size$chromosome <-as.factor(compiled.size$chromosome)

  #prepare plot
  p <- ggplot()

  #for each file input to function which has been prepared for the polygon plot, plot onto graph
  for(i in 1:length(listsynt)){
    data <- listsynt[[i]]
    reference <- data[1,"ref"]
    target <- data[1,"tar"]
    ref_sizes <- compiled.size[compiled.size$species == reference,]
    tar_sizes <- compiled.size[compiled.size$species == target,]
    p = p + ggplot2::geom_rect(data=ref_sizes, mapping=ggplot2::aes(xmin=xstart, xmax=xend, ymin=y, ymax=y+0.10, fill=chromosome),
                               color="black", alpha = 0.85, size = 0.2 ) +
      ggplot2::geom_text(data=ref_sizes,ggplot2::aes(x=(xstart+xend)/2,y=y+0.2,label=chromosome),size=2,angle=45) +
      ggplot2::geom_text(data=ref_sizes,mapping=ggplot2::aes(x=2,y=y, label=species),size=3,hjust = 1) +
      ggplot2::geom_rect(data=tar_sizes, mapping=ggplot2::aes(xmin=xstart, xmax=xend, ymin=y, ymax=y+0.10),fill="grey85",
                         color="black", alpha = 0.85, size = 0.2 ) +
      ggplot2::geom_text(data=tar_sizes,ggplot2::aes(x=(xstart+xend)/2,y=y+0.2,label=chromosome),size=2,angle=45) +
      ggplot2::geom_text(data=tar_sizes,mapping=ggplot2::aes(x=2,y=y, label=species),size=3,hjust = 1) +
      ggplot2::geom_polygon(data = data, alpha = .5, ggplot2::aes(x = x, y = y, group = group, fill = fill))
  }



  #edit graph to keep colours constant and 'tidy' axis
  p = p +  ggplot2::scale_fill_manual(values = c("1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
                                                 "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
                                                 "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
                                                 "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
                                                 "13"= "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
                                                 "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
                                                 "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
                                                 "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
                                                 "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
                                                 "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
                                                 "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
                                                 "LGE22" = "grey", "LGE64" = "gray64",
                                                 "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
                                                 "Gap" = "white", "LG2" = "black", "LG5" = "#CC99CC")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank(),
                   legend.position="none")

  #save plot as image
  ggsave(paste0(output,".",fileformat),p,device = fileformat,scale = 1.5)
}


#' Draw synteny ideograms in inferCARS style
#'
#' This function draws pairwise synteny plots in inferCARS style.
#' Inputs are tab separated files, one with the synteny blocks and two files with target and reference chromosome sizes.
#' Synteny blocks file should be as: targetChr targetStart targetEnd referenceChr referenceStart referenceEnd Orientation
#' Output will be a pdf file with the ideogram.
#'
#' Target is the species which chromosomes will be painted. Reference will be used for painting and diagonals.
#' Chromosomes will be in the same order as in the target sizes file.
#'
#' Example: draw.ideogram(synteny_file, target_chr_size, reference_chr_size)
#' @title Draw ideograms in inferCARs style
#' @param file_data Path to the syntenic blocks file
#' @param file_tarsize Path to the target chromosomes length file
#' @param file_refsize Path to the reference chromosomes length file
#' @return A pdf file with the ideogram
#' @export

draw.ideogram <- function(file_data, file_tarsize, file_refsize) {
  # To make the rectangles wider, change the height of the pdf
  # Refchr refstart ref end tarchr tarstart tarend - ref is ideogram target is used for painting and diagonals

  # Read input files
  size<-tarstart<-tarend<-refchr<-ystart<-yend<-NULL
  tar_sizes = read.delim(file_tarsize, header = FALSE)
  ref_sizes = read.delim(file_refsize, header = FALSE)
  data = read.delim(file_data, header = FALSE)

  colnames(tar_sizes) = c("tarchr", "size")
  colnames(ref_sizes) = c("refchr", "size")
  colnames(data) = c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "orien","notUsed")

  data$tarchr = factor(data$tarchr, levels = tar_sizes$tarchr)
  data$refchr = factor(data$refchr, levels = ref_sizes$refchr)

  for (i in c(1:nrow(data))){
    dir = data[i, "orien"]
    chr = data[i, "refchr"]
    full_len = ref_sizes[ref_sizes$refchr==chr,2]
    y1 = round(data[i, "refstart"]/full_len, digits = 4)
    y2 = round(data[i, "refend"]/full_len, digits = 4)

    inverted = grepl("-", dir, fixed = TRUE)
    if(inverted == TRUE){
      data[i,"ystart"] = y2
      data[i, "yend"] = y1
    } else{
      data[i,"ystart"] = y1
      data[i,"yend"] = y2
    }
  }

  pdf(paste0(file_data,".pdf"), width = 8.5, height = 10, pointsize = 5)

  print(ggplot2::ggplot(size = 0.2, font = 10, data = data) +
          ggplot2::geom_rect(data=tar_sizes, mapping=ggplot2::aes(xmin=1, xmax=size, ymin=-0.1, ymax=1.1),
                             fill="white", color="black", alpha = 0.85, size = 0.2 ) +
          ggplot2::geom_rect(data=data, mapping=ggplot2::aes(xmin=tarstart, xmax=tarend, ymin=-0.1, ymax=1.1, fill=refchr),
                             color="black", alpha = 0.85, size = 0.2 ) +
          ggplot2::geom_segment(data=data, mapping=ggplot2::aes(x=tarstart, y=ystart, xend=tarend, yend=yend), size = 0.2) +
          ggplot2::facet_grid(tarchr ~ .) +
          ggplot2::labs(fill = "Reference", x = "Chomosome length (Mb)", size = 10) +
          ggplot2::theme(
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 10),
            axis.ticks.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_line(size = 0.2),
            strip.text.y = ggplot2::element_text(angle = 0, face = "bold", size = 10),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            legend.title.align = 0.5) +
          ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
          ggplot2::scale_fill_manual(values = c("1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
                                                "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
                                                "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
                                                "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
                                                "13"= "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
                                                "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
                                                "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
                                                "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
                                                "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
                                                "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
                                                "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
                                                "LGE22" = "grey", "LGE64" = "gray64",
                                                "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
                                                "Gap" = "white")) +
          ggplot2::scale_x_continuous(breaks = c(0,2.5e+07,5e+07,7.5e+07,1e+08,1.25e+08,1.5e+08,1.75e+08,2e+08),
                                      labels = c("0","25","50","75","100","125","150","175","200")))

  dev.off()
}
