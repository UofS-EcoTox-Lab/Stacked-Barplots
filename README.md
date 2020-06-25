### Stacked barplots

```
library(vegan) # Used for decostand

rm(list = ls())
```

#### 1. Read in the data

The data are available [here](https://www.dropbox.com/s/ns1wlpg7zkdkxca/ulrich.df.csv?dl=0).

`
ulrich.df <- read.csv("~/Dropbox/r code repository/Data/ulrich.df.csv", header = TRUE)
`

#### 2. Convert raw counts to relative abundances

``
ulrich.prop <- decostand(ulrich.df[,c(8:25)], "total", MARGIN = 2)*100 
colSums(ulrich.prop) # Check to see if the math was right (all should sum to 100)
``

#### 3. Calculate mean abundance for each Genus within each sample

``````
ph <- rowMeans(ulrich.prop[,c(1,3,5)])
pl <- rowMeans(ulrich.prop[,c(2,4,6)])
uh <- rowMeans(ulrich.prop[,c(7,9,11)])
ul <- rowMeans(ulrich.prop[,c(8,10,12)])
ch <- rowMeans(ulrich.prop[,c(13,15,17)])
cl <- rowMeans(ulrich.prop[,c(14,16,18)])
``````

#### 4. Create the df to use for plotting

Convert to a data frame sorted by Genus - this is a bit of roundabout way to do this (sum of one cell is just that cell value):

``````
ph1 <- as.data.frame(tapply(ph, INDEX = ulrich.df$Genus, sum))
pl1 <- as.data.frame(tapply(pl, INDEX = ulrich.df$Genus, sum))
uh1 <- as.data.frame(tapply(uh, INDEX = ulrich.df$Genus, sum))
ul1 <- as.data.frame(tapply(ul, INDEX = ulrich.df$Genus, sum))
ch1 <- as.data.frame(tapply(ch, INDEX = ulrich.df$Genus, sum))
cl1 <- as.data.frame(tapply(cl, INDEX = ulrich.df$Genus, sum))
``````

Combine the sorted relative abundances into a matrix:

`
ulrich.bp <- as.matrix(cbind(uh1, ul1, ch1, cl1, ph1, pl1)) # Order is unlabelled, carbon, phosphate
`

Add column names:

`
colnames(ulrich.bp) <- c("uh","ul", "ch","cl","ph","pl")
`

Assign the genera as row names:

`
rownames(ulrich.bp) <- levels(as.factor(rownames(ph1)))
`

Convert to a data frame

`
ulrich.bp1 <- as.data.frame(ulrich.bp)
`

Calculate the means of each genera among all samples and use this to sort the df in order of increasing abundance:

`
ulrich.bp1$mean <- rowMeans(ulrich.bp)
`

Sort by mean abundance

`
ulrich.bp1 <- ulrich.bp1[order(ulrich.bp1$mean),]
`

Subset to the 20 most abundant genera

`
ulrich.bp2 <- ulrich.bp1[-c(1:78),]
`

Here are the remainders to use to calculate the total abundance of other, non-top-20 abundant genera:

`
ulrich.bp3 <- ulrich.bp1[c(1:78),]
`

Calculate total abundance of the "others":

`
ulrich.others <- as.data.frame(t(colSums(ulrich.bp3)))
`

Add row name and then add to the top-20 data frame:

`
rownames(ulrich.others) <- "Others"
`

Adding "others" to top-20 abundance data frame:

`
ulrich.bp4 <- rbind(ulrich.bp2, ulrich.others)
`

Convert to matrix and remove the mean column added earlier:

`
ulrich.bp5 <- as.matrix(ulrich.bp4[,-7])
`

Here I've chosen a colour palette from [I Want Hue](http://tools.medialab.sciences-po.fr/iwanthue/):

``````````````````````
mycols2 <- palette(c(
  rgb(47,75,244, maxColorValue=255),
  rgb(0,129,31, maxColorValue=255),
  rgb(29,0,152, maxColorValue=255),
  rgb(97,105,0, maxColorValue=255),
  rgb(110,0,149, maxColorValue=255),
  rgb(135,215,168, maxColorValue=255),
  rgb(223,0,142, maxColorValue=255),
  rgb(1,211,241, maxColorValue=255),
  rgb(211,90,0, maxColorValue=255),
  rgb(38,149,255, maxColorValue=255),
  rgb(255,126,87, maxColorValue=255),
  rgb(0,59,119, maxColorValue=255),
  rgb(255,78,102, maxColorValue=255),
  rgb(0,106,142, maxColorValue=255),
  rgb(137,0,28, maxColorValue=255),
  rgb(197,160,255, maxColorValue=255),
  rgb(89,53,0, maxColorValue=255),
  rgb(157,0,118, maxColorValue=255),
  rgb(215,188,222, maxColorValue=255),
  rgb(64,23,48, maxColorValue=255),
  rgb(255,127,187, maxColorValue=255)))
``````````````````````

#### 5. Create stacked barplot sorted by mean abundance

`````````
# Export at 5.25 x 6.5
par(mar=c(3.5, 3.5, 1, 10), xpd=NA)
bp <- barplot(ulrich.bp5, col = mycols2, border = NA, width = 2, ylim = c(105,-5), yaxt = "n")
axis(side = 1, at = bp, labels = NA, lwd = 1, tck = -0.02)
axis(side = 2, at = seq(0,100,20), labels = seq(100,0,-20), lwd = 1, tck = -0.02)
mtext("Treatment", side = 1, line = 2.2)
mtext("Relative abundance (%)", side = 2, line = 2.2)
legend("right",inset=c(-0.37,0), fill = mycols2, legend=rownames(ulrich.bp5), border = NA, bty = "n", text.font = 3)
box()
`````````

<img src = "https://user-images.githubusercontent.com/44586553/69078096-bfa6a700-09fc-11ea-9f96-3aad6d425c4b.jpg" width="500" height="600">
