#' Draw Evolution Highway Plots
#'
#' This function draws Evolution Highway style plots. It requieres as input the syntenic blocks following this
#' format: referenceSpecies,chr,start,end,targetChr,targetStart,targetEnd,orient,targetSpecies separated by tabs
#' It also requieres the output file name and the range of chromosomes of the reference species.
#' Example: draw.eh(input.csv,goat,"1:29")
#'
#' @title Evolution Highway style plot
#' @param infile Path to the syntenic blocks file
#' @param output file name
#' @param chrRange range of chromosome numbers in the reference "1:29,X"
#' @return A pdf file with the comparative drawings
#' @export
draw.eh<-function(infile,outfile,chrRange) {
  data<-read.table(infile, header=TRUE)
  colnames(data)
  data$orient = factor(data$orient, levels=c("1","-1"))
  data$text_size2=80*((data$end-data$start)/100000)

  pdf(paste0(outfile,".pdf"),width=5.5, height =10, pointsize = 10)
  for (ID in c(chrRange,"X")) {
    #ID=chrRange
    print(ID)
    subsetChr1<-subset(data,chr==ID, select=c(ref,chr,start,end,tarChr,tarSt,tarEnd,orient,tar,text_size2))
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
#' This function draws pairwise synteny plots. It requieres as input the syntenic blocks following this
#' format: referenceSpecies,chr,start,end,targetChr,targetStart,targetEnd,orient,targetSpecies. Please separate bu tabs.
#' It also requieres the output file name and the range of chromosomes of the reference and target species.
#' Example: draw.pairwise(input.txt,outputName, refSizes, tarSizes, refName, tarName)
#'
#' @title Pairwise synteny plot
#' @param infile Path to the syntenic blocks file
#' @param output file name
#' @param tarSizes file with chromosome sizes for the target species
#' @param refSizes file with the chromosome sizes for the reference species
#' @param refName name of the reference species
#' @param tarName name of the target species
#' @return A pdf file with the comparative drawings
#' @export
draw.pairwise <- function(infile,output,refSizes,tarSizes,refName,tarName) {
  data<- read.delim(infile, header=FALSE)
  ref_sizes <-read.delim(refSizes, header=FALSE)
  tar_sizes <-read.delim(tarSizes, header=FALSE)
  colnames(data) = c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "dir")
  colnames(ref_sizes) = c("refchr", "size")
  colnames(tar_sizes) = c("tarchr", "size")

  #This adds gap in between reference chromosomes and convert to "linear" genome
  for (i in c(1:nrow(ref_sizes))){
    #print(i)
    if (i == 1){
      total_start = 1
      total_end = ref_sizes[i, "size"]
    } else {
      total_start = total_end + 6000000
      total_end = total_start + ref_sizes[i, "size"]
    }
    ref_sizes[i,"xstart"] = total_start
    ref_sizes[i, "xend"] = total_end
  }

  #This adds gap in between target chromosomes
  for (i in c(1:nrow(tar_sizes))){
    #print(i)
    if (i == 1){
      total_start = 1
      total_end = tar_sizes[i, "size"]
    } else {
      total_start = total_end + 6000000
      total_end = total_start + tar_sizes[i, "size"]
    }
    tar_sizes[i,"xstart"] = total_start
    tar_sizes[i, "xend"] = total_end
  }

  #This converts coordinates to linear genome and creates synteny polygon coordinates
  synteny = data.frame()
  for (i in c(1:nrow(data))){
    tar_chr = data[i,"tarchr"]
    ref_chr = as.factor(data[i,"refchr"])
    dir = data[i, "dir"]
    tar_add = tar_sizes[as.character(tar_sizes$tarchr)==as.character(tar_chr),]$xstart
    ref_add = ref_sizes[as.character(ref_sizes$refchr)==as.character(ref_chr),]$xstart
    tar_y = 0.10
    ref_y = 2
    tar_xstart = data[i,"tarstart"] + tar_add
    tar_xend = data[i,"tarend"] + tar_add
    ref_xstart = data[i,"refstart"] + ref_add
    ref_xend = data[i,"refend"] + ref_add

    inverted = grepl("-", dir, fixed = TRUE)
    if(inverted == TRUE){
      df = data.frame(x = c(tar_xstart, tar_xend, ref_xstart, ref_xend), y = c(tar_y, tar_y, ref_y, ref_y),
                      fill = ref_chr, group = paste0("s",i))
    } else {
      df = data.frame(x = c(tar_xstart, ref_xstart, ref_xend,  tar_xend), y = c(tar_y, ref_y, ref_y, tar_y),
                      fill = ref_chr, group = paste0("s",i))
    }
    synteny = rbind(synteny,df)
  }

  #making sure chr columns are factors
  tar_sizes$tarchr<-as.factor(tar_sizes$tarchr)
  ref_sizes$refchr<-as.factor(ref_sizes$refchr)


  pdf(paste0(output,".pdf"),width=20, height =5, pointsize = 10)
  #This prints plot
  ggplot2::ggplot(size = 0.2, font = 10, data = data) +
    ggplot2::geom_rect(data=ref_sizes, mapping=ggplot2::aes(xmin=xstart, xmax=xend, ymin=2, ymax=2.10, fill=refchr),
                       color="black", alpha = 0.85, size = 0.2 ) +
    ggplot2::geom_text(data=ref_sizes,ggplot2::aes(x=(xstart+xend)/2,y=2.15,label=refchr),size=2,angle=45) +
    ggplot2::geom_text(mapping=ggplot2::aes(x=2,y=2.3, label=refName),size=3,hjust = 1) +
    ggplot2::geom_rect(data=tar_sizes, mapping=ggplot2::aes(xmin=xstart, xmax=xend, ymin=0, ymax=0.10),fill="grey85",
                       color="black", alpha = 0.85, size = 0.2 ) +
    ggplot2::geom_text(data=tar_sizes,ggplot2::aes(x=(xstart+xend)/2,y=-0.05,label=tarchr),size=2,angle=45) +
    ggplot2::geom_text(mapping=ggplot2::aes(x=2,y=-0.20, label=tarName),size=3,hjust = 1) +
    ggplot2::geom_polygon(data = synteny, alpha = .5, ggplot2::aes(x = x, y = y, group = group, fill = fill)) +
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

  dev.off()
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
#' Example: draw_ideogram(synteny_file, target_chr_size, reference_chr_size)
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
  tar_sizes = read.delim(file_tarsize, header = FALSE)
  ref_sizes = read.delim(file_refsize, header = FALSE)
  data = read.delim(file_data, header = FALSE)

  colnames(tar_sizes) = c("tarchr", "size")
  colnames(ref_sizes) = c("refchr", "size")
  colnames(data) = c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "orien")

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

  pdf(paste0(file_data,".pdf"), width = 8.5, height = 11, pointsize = 10)

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
