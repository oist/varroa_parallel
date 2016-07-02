#---------------------- Morphometrics of Arz-Txs populations ------------------------
#------------------------------- arora@evolbio.mog.de -------------------------------
rm(list=ls())
library(ggplot2)
library(shapes)
library(asbio)
library(Hmisc)

# convert from 2d matrix of coordinates to 3d matrix
make3dMatrix = function(dat){

   k = 19 # number of landmarks
   m = 2  # number of dimensions
   ns = nrow(dat)	 # number of samples

   input = array(0, dim=c(k,m,ns))
   for(n in 1:ns){
      i = 1
      for(k in 1:19){
         for(m in 1:2){
            input[k,m,n] =  unlist(dat[n,], use.names=F)[i]
            i = i + 1
         }
      }
   }

   return(input)
}

#---------------------------------------------------------------------------------
# load africanisation data
admix = read.table("admixture_K3.txt",header=T)

file.create("stat.txt")
to.plot = c()
for(pop in c("arizona", "texas")){

   # load raw coordinates data in _format and _db formats
   raw = read.delim(paste(pop, "_coord_Jun_Format.txt", sep=""), header = F, sep = "\t", skip = 1)
   db = read.delim(paste(pop, "_coord_Jun_DB.txt", sep=""))
   db = db[,-ncol(db)]

   # load samples year information
   yr = read.delim(paste(pop, "_yrs.txt", sep=""), header = F)
#   if(pop == "texas"){
#      yr = yr[,-2]
#      colnames(yr) = c("V1","V2")
#   }

   # convert coordinate matrix into 3d one and do procrustes anaylsis
   mat = make3dMatrix(dat = raw)
   p = procGPA(x = mat, scale = T, distances = T)

   # extract centroid sizes
   centroid.size = data.frame(size = p$size)
   centroid.size$id = sub("-[1-2]\\..*","",db$image)
   centroid.size$id = sub("\\.[1-2]\\..*","",centroid.size$id)

   # for correlation with africanisation, extract the ids of samples
#   if(pop == "texas"){
#      admix.pop = admix[grep("-", admix$id),]
#      admix.pop$id = sub("-.*","",admix.pop$id)
#   } else{
#      admix.pop = admix[-grep("-", admix$id),]
#   }

   # merge centroid size and admixture, and check correlation between them

   centroid.admix = merge(x = centroid.size, y = admix, by = "id", all.x = T, all.y = F, sort = F)
   oldnames <- colnames(centroid.admix)
   centroid.admix <- aggregate(centroid.admix[,c(2:5)],by=list(centroid.admix$id),FUN=mean)
   colnames(centroid.admix) <- oldnames
   t = cor.test(centroid.admix$size,centroid.admix$POP2, method = "s", exact = F)
   print(paste(pop,t$estimate,t$p.value))

   # store centroid size thereby excluding drones
   frame = merge(x = centroid.size, y = yr, by.x = "id", by.y = "V1", all.x = T, sort = F)
   frame = frame[frame$size < 900,]
   frame$pop = capitalize(pop)
   to.plot = rbind(to.plot, frame)

   # calcualte pairwise wilcoxon rank sum test, paired, between populations
#   write.table(x = pop, file = "stat.txt", append = T, quote = F, row.names = F, col.names = F)
#   write.table(x = paste(c("year1", "year2", "wilcoxon.paired.pvalue"), collapse = "\t"), file = "stat.txt", append = T, quote = F, row.names = F, col.names = F)
#   yrs = sort(unique(frame$V2))
#   for(i in 1:(length(yrs)-1)){
#      wilcox.t = wilcox.test(frame[frame$V2 == yrs[i],"size"], frame[frame$V2 == yrs[i+1],"size"])
#      write.table(x = paste(c(as.character(c(yrs[i], yrs[i+1])), round(wilcox.t$p.value,4)), collapse = "\t"), file = "stat.txt", append = T, quote = F, row.names = F, col.names = F)
#   }

}

to.plot2 <- na.omit(aggregate(to.plot$size, by=list(to.plot$id),FUN=mean))
colnames(to.plot2) <- c("id", "size")
metadata <- unique(na.omit(to.plot[,c(1,3,4)]))
to.plot2 <- na.omit(merge(x = to.plot2, y = metadata, by = "id", all.x = T, all.y = T, sort = F))

with(subset(to.plot2, pop == "Arizona"), cor.test(size,V2,method="s"))
with(subset(to.plot2, pop == "Texas"), cor.test(size,V2,method="s"))

# plot centroid size, faceted by populations
p <- ggplot(data = to.plot2, aes(x = factor(V2), y = size)) + geom_boxplot() + xlab("Year") + ylab("Centroid size") + theme_bw()
p <- p + facet_wrap( ~ pop, ncol = 1, ) + theme(panel.grid.major = element_blank(), strip.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(p)

ggplot(na.omit(join(to.plot2,admix)),aes(as.factor(V2),POP2))+facet_grid(pop~.)+geom_boxplot()


#Figures S1
to.plot3 <- na.omit(join(to.plot2,admix))
to.plot3$pop <- factor(to.plot3$pop, levels = c("Texas", "Arizona"))
colnames(to.plot3)[3] <- "year"
colnames(to.plot3)[6] <- "African"
to.plot3$epoch <- "pre"
to.plot3$epoch[to.plot3$year>1995] <- "post"

grid.newpage()
grid.raster(brewer.pal(7,"PRGn"), 0.25, 0.5, 0.4, 1)
ggplot(to.plot3,aes(as.factor(year),size,color=African))+facet_grid(pop~.,scales = "free_y")+geom_jitter(width = .3)+ xlab("Year") + ylab("Centroid size") + theme_bw()+theme(panel.grid.major = element_blank(), strip.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),legend.position="top")+stat_smooth(color="blue",method="lm", aes(group=epoch)) + scale_colour_gradient()
ggsave("size over time.pdf", width = 5, height = 3,device=cairo_pdf())
ggsave("size over time.svg", width = 5, height = 3)


size.pop.lm <- lm(size~POP2+pop, data=subset(na.omit(join(to.plot2,admix)),V2>1996))

summary(size.pop.lm)
partial.R2(lm(size~pop, data=subset(na.omit(join(to.plot2,admix))),V2>1996), size.pop.lm)

summary(lm(size~POP2, data=subset(na.omit(join(to.plot2,admix)),V2==1995)))
ggplot(subset(na.omit(join(to.plot2,admix)),V2>=1994 & V2<=1995),aes(POP2,size))+geom_point()
