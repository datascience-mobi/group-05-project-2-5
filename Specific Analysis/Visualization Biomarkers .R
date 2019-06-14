####################################  visualization biomarkers ################################################

install.packages("rbokeh")
library(rbokeh)

figure() %>%
  ly_points(biomarkers_FC_values, color = "red", size = 20)

n <- nrow(biomarkers_FC_values)
ramp <- colorRampPalette(c("red", "blue"))(n)
figure() %>%
  ly_points(biomarkers_FC_values, color = ramp, size = seq_len(n))

################################################################################################################
################################################################################################################

library(ggplot2)

x=cbind(biomarkers_FC_values,biomarkers_FC_genes)
test=as.data.frame(x)

p=ggplot(test)+aes(biomarkers_FC_genes,FC_meanrow)+theme(axis.text.x=element_blank(),axis.text.y=element_blank())
p + geom_point(aes(size = FC_meanrow),show.legend = FALSE,col='blue')
# point size according to p.Value 

################################################################################################################

### more colors 
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
library(RColorBrewer)

pal <- wes_palette("Zissou1", 100, type = "continuous")
p + geom_point(aes(size = FC_meanrow),show.legend = FALSE)+scale_color_brewer(palette="Dark2")

# does not work 

################################################################################################################
# Install
install.packages("wesanderson")
# Load
library(wesanderson)
names(wes_palettes)
p + geom_point(aes(size = FC_meanrow),show.legend = FALSE)+scale_color_manual(values=wes_palette(n=3, name="Royal1"))

# does not work 
