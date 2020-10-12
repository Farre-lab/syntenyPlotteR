#' Draw Evolution Highway Plots
#'
#' This function draws Evolution Highway style plots. It requieres as input the syntenic blocks following this
#' format: referenceSpecies,chr,start,end,targetChr,targetStart,targetEnd,orient,targetSpecies
#' It also requieres the output file name and the range of chromosomes of the reference species.
#' Example: draw.eh(input.csv,goat,"1:29,X")
#'
#' @param infile Path to the syntenic blocks file
#' @param output file name
#' @param chrRange range of chromosome numbers in the reference "1:29,X"
#' @return A pdf file with the comparative drawings
#' @export
draw.eh<-function(infile,outfile,chrRange) {
  data<-read.csv(infile, header=TRUE)
  colnames(data)
  data$orient = factor(data$orient, levels=c("1","-1"))
  data$text_size2=80*((data$end-data$start)/100000)

  pdf(outfile,width=5.5, height =10, pointsize = 10)
  for (ID in c(chrRange)) {
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
