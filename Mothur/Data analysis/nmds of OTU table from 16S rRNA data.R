
# nmds of OTU table from 16S rRNA data

# set wd and load libraries


library(ggplot2)
library(grid)
library(reshape2)
library(vegan)
library(RColorBrewer)

# ==================== read in data & do NMDS ==============================================# 

# 1 read in OTU table -> put samples as rows, variables as columns
# should only contain numeric data. if it has taxonomy in, subset to remove
df <- read.table("data.txt", sep = "\t", header=T, row.names=1)

mdf <- as.matrix(df) # make sure it's in correct form

# read in your metadata -> sample names must match those in OTU table
meta <- read.table("metadata.txt", header = TRUE, row.names = 1, sep='\t')



# create distance matrix. Change method to which ever is better for your data

df.dist <- vegdist(mdf, method="bray", binary=FALSE, diag=FALSE, upper=FALSE)


# now do nmds - for distance = put which ever method used in vegdist() step
# first time, try k = 2 (ie a 2D nmds). 
NMDS <- metaMDS(df.dist, distance = "euclidean", k = 2, trymax = 20, autotransform =TRUE,
                noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
                old.wa = FALSE)

# Copy and paste what's printed to console into notepad for ref.
# ie note down stress value. If it's larger than 0.2 (or even 0.1)
# consider re-running metaMDS with k = 3 to make a 3d nmds

# check the fit of your data 
stressplot(NMDS) #again note down R2 
                # check r2 value and scattered-ness of points 
              # from the line shouldnt be big

# then basic plot of nmds to check it worked
plot(NMDS3, type="t")

## now we will read data out of here so it can be plotted in
# ggplot to look prettier. 

data.scores <- as.data.frame(scores(NMDS)) 
data.scores$id<- rownames(data.scores) # create a column called id which is sample names
# now add info re treatment and sample day as two extra colums in data.scores

data.scores$treatment<- meta$Treatment # add treatment column
data.scores$time<- meta$timepoint # add timepoint column

# now run setFactorOrder.R script to order your factors so that 
# when you plot they are in correct order instead of alphabetical order

data.scores[["treatment"]] <- setFactorOrder(data.scores[["treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
data.scores[["time"]] <- setFactorOrder(data.scores[["time"]], c("T0", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13"))

# ==================== and now to plotting  ==============================================# 

P <- ggplot(data.scores, aes(x= NMDS1, y= NMDS2, colour=treatment))+
  geom_point(size=7, stroke =2, shape = 1)  + # size and line thickness of plotted points
  #scale_shape_manual(values = c(1, 2)) + if you want diff shapes uncomment add shapes=
  scale_colour_manual(values = colorRampPalette(brewer.pal(8, "Set1"))(20)) =
  #set your colour palette (colourbrewer) - change number in brackets
  # to number of samples at the moment its 20 samples
  theme(legend.key.size=unit(0.3,"cm")) + # change legend icon size
  theme_bw()  # remove ugly grey background
P

# add labels to centre of each point 
# leave this step out if you dont want labelled points 
P2 <- P +geom_text(aes(label=time),hjust=0.45, vjust=0.3, 
                   size=3, colour="black", fontface = "bold", show.legend=FALSE) 

# now sorting font size etc. 
Phase1 = P2 + 
  labs(shape="Time point", colour="Treatment", x = "Axis 1", y = "Axis 2") + 
  # update legend titles.
  guides(color = guide_legend(override.aes = list(size=7))) +
  # change legend icon size for Treatment
  theme(legend.key.size = unit(1.8,"line"),
        # change spacing of legend icons. 
        axis.text.x=element_text(size=14, colour = "black"), # x axis text
        axis.title.x=element_text(size=16, colour = "black"),
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y=element_text(size=16, colour="black"),
        legend.title = element_text(size=16, colour = "black"),
        legend.text = element_text(size=14, colour = "black"))

# ======================= make your own colout palette ================ # 

# use ggthemr to plot with my WP3 colour theme so 4 treatments
# labelled as in all my other plots. 

library(ggthemr)
# define you colours (based on names off ggplot colours)
WP3_colsA <- c("black", "chocolate4", "slateblue", "olivedrab")
# add background
WP3_cols <- c("#555555", WP3_colsA)

# define palette
WP3Cols <- define_palette(
  swatch = WP3_cols, # colours for plotting points and bars
  gradient = c(lower = WP3_cols[1L], upper = WP3_cols[2L]), #upper and lower colours for continuous colours
  background = "white" #defining a grey-ish background 
)

# set new col theme as current

ggthemr(WP3Cols)    

# and replot with new color scheme by printing name you assigned to plot. 

