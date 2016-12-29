args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(ggplot2))
datain<-read.table(args[1],header=F)
colnames(datain)<-c('Chrom','Pos','Type','Freq')
datain$Type = factor(datain$Type,
                     levels(datain$Type)[c(1,3,4,2)])
hist_plot<-ggplot(datain,aes(x=Freq, fill=Type)) +
  theme_bw()+
  geom_histogram(binwidth = 0.5, alpha=0.4) +
  ggtitle(args[1]) +
  ylab("Count positions") +
  xlab("Allele Freq") +
  scale_x_continuous(limits=c(1,100))+
  scale_fill_manual(values=c("#920000", "#24FF24", "#6DB6FF", "#FFFF6D")) #Friendly to colour-blind users: http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
pdf(file=args[2])
suppressWarnings(print(hist_plot))
invisible(dev.off())
