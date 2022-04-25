########For figure 3a, b and extended figure 7c

# Libraries included
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pipeR)
library(devtools)
library(data.table)

library(convert)
library(devtools)

library(RColorBrewer) # Load a package giving more colors
library(pheatmap) # load a package for making heatmaps
library(GOplot)
library(DOSE)
library(ggpubr)
library(rstatix)

organism = "org.Mm.eg.db"
library(clusterProfiler)
library(wordcloud)
library(pathview)
library(organism, character.only = TRUE)
library(sojourner)
library(celltrackR)
library(flowcatchR)
library(ggbeeswarm)
library(lm.beta)

############ Import cell tracks

getwd()
## setwd() is required to set the list of files from the directory to analyze.

tracks <- map(list.files(), function(x){
  fread(x, header = T)})

tracks <- set_names(tracks, file_paths)



############### classify samples and convert to tracks 
gs <- map(tracks, function(x){
  as.data.frame(x) %>%
    slice(4:nrow(x)) %>%
    filter(!is.na(as.numeric(TRACK_ID))) %>%
    dplyr::select(TRACK_ID, POSITION_T, POSITION_X, POSITION_Y)%>%
    arrange(TRACK_ID) %>%
    mutate_if(sapply(., is.character), as.numeric) %>%
    as.tracks(id.column = 1,
              time.column = 2,
              pos.columns = c(3, 4)) })

str(gs)
view(gs)

b_t <- c(gs[[1]], gs[[2]], gs[[11]], gs[[12]], gs[[21]], gs[[22]], gs[[33]], gs[[34]], gs[[35]], gs[[36]])
c_t <- c(gs[[4]], gs[[13]], gs[[14]], gs[[23]], gs[[37]], gs[[38]], gs[[39]], gs[[40]])
d_t <- c(gs[[5]], gs[[6]], gs[[15]], gs[[16]], gs[[24]], gs[[25]], gs[[26]], gs[[41]], gs[[42]], gs[[43]])
e_t <- c(gs[[7]], gs[[8]], gs[[17]], gs[[18]], gs[[27]], gs[[28]], gs[[29]], gs[[44]], gs[[45]], gs[[46]], gs[[47]])

#biological replicates

b_t_1 <- c(gs[[1]], gs[[2]])
b_t_2 <- c(gs[[11]], gs[[12]])
b_t_3 <- c(gs[[21]], gs[[22]])
b_t_4 <- c(gs[[33]], gs[[34]], gs[[35]], gs[[36]])
c_t_1 <- gs[[4]]
c_t_2 <- c(gs[[13]], gs[[14]])
c_t_3 <- gs[[23]]
c_t_4 <- c(gs[[37]], gs[[38]], gs[[39]], gs[[40]])
d_t_1 <- c(gs[[5]], gs[[6]])
d_t_2 <- c(gs[[15]], gs[[16]])
d_t_3 <- c(gs[[24]], gs[[25]], gs[[26]])
d_t_4 <- c(gs[[41]], gs[[42]], gs[[43]])
e_t_1 <- c(gs[[7]], gs[[8]])
e_t_2 <- c(gs[[17]], gs[[18]])
e_t_3 <- c(gs[[27]], gs[[28]], gs[[29]])
e_t_4 <- c(gs[[44]], gs[[45]], gs[[46]], gs[[47]])

#technical replicates

b_t_1_a <- gs[[1]]
b_t_1_b <- gs[[2]]
b_t_2_a <- gs[[11]]
b_t_2_b <- gs[[12]]
b_t_3_a <- gs[[21]]
b_t_3_b <- gs[[22]]
b_t_4_a <- gs[[33]]
b_t_4_b <- gs[[34]]
b_t_4_c <- gs[[35]]
b_t_4_d <- gs[[36]]
c_t_1_a <- gs[[4]]
c_t_2_a <- gs[[13]]
c_t_2_b <- gs[[14]]
c_t_3_a <- gs[[23]]
c_t_4_a <- gs[[37]]
c_t_4_b <- gs[[38]]
c_t_4_c <- gs[[39]]
c_t_4_d <- gs[[40]]
d_t_1_a <- gs[[5]]
d_t_1_b <- gs[[6]]
d_t_2_a <- gs[[15]]
d_t_2_b <- gs[[16]]
d_t_3_a <- gs[[24]]
d_t_3_b <- gs[[25]]
d_t_3_b <- gs[[26]]
d_t_4_a <- gs[[41]]
d_t_4_b <- gs[[42]]
d_t_4_c <- gs[[43]]
e_t_1_a <- gs[[7]]
e_t_1_b <- gs[[8]]
e_t_2_a <- gs[[17]]
e_t_2_b <- gs[[18]]
e_t_3_a <- gs[[27]]
e_t_3_b <- gs[[28]]
e_t_3_c <- gs[[29]]
e_t_4_a <- gs[[44]]
e_t_4_b <- gs[[45]]
e_t_4_c <- gs[[46]]
e_t_4_d <- gs[[47]]

############# filtering based on speed

b_fast.tracks <- selectTracks( b_t, speed, 0.3, Inf )
b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.tracks <- selectTracks( c_t, speed, 0.3, Inf )
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.tracks <- selectTracks( d_t, speed, 0.3, Inf )
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.tracks <- selectTracks( e_t, speed, 0.3, Inf )
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]

b_1_fast.tracks <- selectTracks( b_t_1, speed, 0.3, Inf )
b_2_fast.tracks <- selectTracks( b_t_2, speed, 0.3, Inf )
b_3_fast.tracks <- selectTracks( b_t_3, speed, 0.3, Inf )
b_4_fast.tracks <- selectTracks( b_t_4, speed, 0.3, Inf )
c_1_fast.tracks <- selectTracks( c_t_1, speed, 0.3, Inf )
c_2_fast.tracks <- selectTracks( c_t_2, speed, 0.3, Inf )
c_3_fast.tracks <- selectTracks( c_t_3, speed, 0.3, Inf )
c_4_fast.tracks <- selectTracks( c_t_4, speed, 0.3, Inf )
d_1_fast.tracks <- selectTracks( d_t_1, speed, 0.3, Inf )
d_2_fast.tracks <- selectTracks( d_t_2, speed, 0.3, Inf )
d_3_fast.tracks <- selectTracks( d_t_3, speed, 0.3, Inf )
d_4_fast.tracks <- selectTracks( d_t_4, speed, 0.3, Inf )
e_1_fast.tracks <- selectTracks( e_t_1, speed, 0.3, Inf )
e_2_fast.tracks <- selectTracks( e_t_2, speed, 0.3, Inf )
e_3_fast.tracks <- selectTracks( e_t_3, speed, 0.3, Inf )
e_4_fast.tracks <- selectTracks( e_t_4, speed, 0.3, Inf )

b_1_fast.fast.tracks <- selectTracks( b_t_1, speed, 2, Inf )
b_2_fast.fast.tracks <- selectTracks( b_t_2, speed, 2, Inf )
b_4_fast.fast.tracks <- selectTracks( b_t_4, speed, 2, Inf )
c_1_fast.fast.tracks <- selectTracks( c_t_1, speed, 2, Inf )
c_2_fast.fast.tracks <- selectTracks( c_t_2, speed, 2, Inf )
c_4_fast.fast.tracks <- selectTracks( c_t_4, speed, 2, Inf )
d_1_fast.fast.tracks <- selectTracks( d_t_1, speed, 2, Inf )
d_2_fast.fast.tracks <- selectTracks( d_t_2, speed, 2, Inf )
d_4_fast.fast.tracks <- selectTracks( d_t_4, speed, 2, Inf )
e_1_fast.fast.tracks <- selectTracks( e_t_1, speed, 2, Inf )
e_2_fast.fast.tracks <- selectTracks( e_t_2, speed, 2, Inf )
e_4_fast.fast.tracks <- selectTracks( e_t_4, speed, 2, Inf )

b_1_a_fast.tracks <- selectTracks( b_t_1_a, speed, 0.3, Inf )
b_1_b_fast.tracks <- selectTracks( b_t_1_b, speed, 0.3, Inf )
b_2_a_fast.tracks <- selectTracks( b_t_2_a, speed, 0.3, Inf )
b_2_b_fast.tracks <- selectTracks( b_t_2_b, speed, 0.3, Inf )
#b_3_a_fast.tracks <- selectTracks( b_t_3_a, speed, 0.3, Inf )
#b_3_b_fast.tracks <- selectTracks( b_t_3_b, speed, 0.3, Inf )
b_4_a_fast.tracks <- selectTracks( b_t_4_a, speed, 0.3, Inf )
b_4_b_fast.tracks <- selectTracks( b_t_4_b, speed, 0.3, Inf )
b_4_c_fast.tracks <- selectTracks( b_t_4_c, speed, 0.3, Inf )
b_4_d_fast.tracks <- selectTracks( b_t_4_d, speed, 0.3, Inf )
c_1_a_fast.tracks <- selectTracks( c_t_1_a, speed, 0.3, Inf )
c_2_a_fast.tracks <- selectTracks( c_t_2_a, speed, 0.3, Inf )
c_2_b_fast.tracks <- selectTracks( c_t_2_b, speed, 0.3, Inf )
#c_3_a_fast.tracks <- selectTracks( c_t_3_a, speed, 0.3, Inf )
c_4_a_fast.tracks <- selectTracks( c_t_4_a, speed, 0.3, Inf )
c_4_b_fast.tracks <- selectTracks( c_t_4_b, speed, 0.3, Inf )
c_4_c_fast.tracks <- selectTracks( c_t_4_c, speed, 0.3, Inf )
c_4_d_fast.tracks <- selectTracks( c_t_4_d, speed, 0.3, Inf )
d_1_a_fast.tracks <- selectTracks( d_t_1_a, speed, 0.3, Inf )
d_1_b_fast.tracks <- selectTracks( d_t_1_b, speed, 0.3, Inf )
d_2_a_fast.tracks <- selectTracks( d_t_2_a, speed, 0.3, Inf )
d_2_b_fast.tracks <- selectTracks( d_t_2_b, speed, 0.3, Inf )
#d_3_a_fast.tracks <- selectTracks( d_t_3_a, speed, 0.3, Inf )
#d_3_b_fast.tracks <- selectTracks( d_t_3_b, speed, 0.3, Inf )
d_4_a_fast.tracks <- selectTracks( d_t_4_a, speed, 0.3, Inf )
d_4_b_fast.tracks <- selectTracks( d_t_4_b, speed, 0.3, Inf )
d_4_c_fast.tracks <- selectTracks( d_t_4_c, speed, 0.3, Inf )
e_1_a_fast.tracks <- selectTracks( e_t_1_a, speed, 0.3, Inf )
e_1_b_fast.tracks <- selectTracks( e_t_1_b, speed, 0.3, Inf )
e_2_a_fast.tracks <- selectTracks( e_t_2_a, speed, 0.3, Inf )
e_2_b_fast.tracks <- selectTracks( e_t_2_b, speed, 0.3, Inf )
#e_3_a_fast.tracks <- selectTracks( e_t_3_a, speed, 0.3, Inf )
#e_3_b_fast.tracks <- selectTracks( e_t_3_b, speed, 0.3, Inf )
#e_3_c_fast.tracks <- selectTracks( e_t_3_c, speed, 0.3, Inf )
e_4_a_fast.tracks <- selectTracks( e_t_4_a, speed, 0.3, Inf )
e_4_b_fast.tracks <- selectTracks( e_t_4_b, speed, 0.3, Inf )
e_4_c_fast.tracks <- selectTracks( e_t_4_c, speed, 0.3, Inf )
e_4_d_fast.tracks <- selectTracks( e_t_4_d, speed, 0.3, Inf )

b_1_a_fast.fast.tracks <- selectTracks( b_t_1_a, speed, 2, Inf )
b_1_b_fast.fast.tracks <- selectTracks( b_t_1_b, speed, 2, Inf )
b_2_a_fast.fast.tracks <- selectTracks( b_t_2_a, speed, 2, Inf )
b_2_b_fast.fast.tracks <- selectTracks( b_t_2_b, speed, 2, Inf )
#b_3_a_fast.tracks <- selectTracks( b_t_3_a, speed, 0.3, Inf )
#b_3_b_fast.tracks <- selectTracks( b_t_3_b, speed, 0.3, Inf )
b_4_a_fast.fast.tracks <- selectTracks( b_t_4_a, speed, 2, Inf )
b_4_b_fast.fast.tracks <- selectTracks( b_t_4_b, speed, 2, Inf )
b_4_c_fast.fast.tracks <- selectTracks( b_t_4_c, speed, 2, Inf )
b_4_d_fast.fast.tracks <- selectTracks( b_t_4_d, speed, 2, Inf )
c_1_a_fast.fast.tracks <- selectTracks( c_t_1_a, speed, 2, Inf )
c_2_a_fast.fast.tracks <- selectTracks( c_t_2_a, speed, 2, Inf )
c_2_b_fast.fast.tracks <- selectTracks( c_t_2_b, speed, 2, Inf )
#c_3_a_fast.tracks <- selectTracks( c_t_3_a, speed, 0.3, Inf )
c_4_a_fast.fast.tracks <- selectTracks( c_t_4_a, speed, 2, Inf )
c_4_b_fast.fast.tracks <- selectTracks( c_t_4_b, speed, 2, Inf )
c_4_c_fast.fast.tracks <- selectTracks( c_t_4_c, speed, 2, Inf )
c_4_d_fast.fast.tracks <- selectTracks( c_t_4_d, speed, 2, Inf )
d_1_a_fast.fast.tracks <- selectTracks( d_t_1_a, speed, 2, Inf )
d_1_b_fast.fast.tracks <- selectTracks( d_t_1_b, speed, 2, Inf )
d_2_a_fast.fast.tracks <- selectTracks( d_t_2_a, speed, 2, Inf )
d_2_b_fast.fast.tracks <- selectTracks( d_t_2_b, speed, 2, Inf )
#d_3_a_fast.tracks <- selectTracks( d_t_3_a, speed, 0.3, Inf )
#d_3_b_fast.tracks <- selectTracks( d_t_3_b, speed, 0.3, Inf )
d_4_a_fast.fast.tracks <- selectTracks( d_t_4_a, speed, 2, Inf )
d_4_b_fast.fast.tracks <- selectTracks( d_t_4_b, speed, 2, Inf )
d_4_c_fast.fast.tracks <- selectTracks( d_t_4_c, speed, 2, Inf )
e_1_a_fast.fast.tracks <- selectTracks( e_t_1_a, speed, 2, Inf )
e_1_b_fast.fast.tracks <- selectTracks( e_t_1_b, speed, 2, Inf )
e_2_a_fast.fast.tracks <- selectTracks( e_t_2_a, speed, 2, Inf )
e_2_b_fast.fast.tracks <- selectTracks( e_t_2_b, speed, 2, Inf )
#e_3_a_fast.tracks <- selectTracks( e_t_3_a, speed, 0.3, Inf )
#e_3_b_fast.tracks <- selectTracks( e_t_3_b, speed, 0.3, Inf )
#e_3_c_fast.tracks <- selectTracks( e_t_3_c, speed, 0.3, Inf )
e_4_a_fast.fast.tracks <- selectTracks( e_t_4_a, speed, 2, Inf )
e_4_b_fast.fast.tracks <- selectTracks( e_t_4_b, speed, 2, Inf )
e_4_c_fast.fast.tracks <- selectTracks( e_t_4_c, speed, 2, Inf )
e_4_d_fast.fast.tracks <- selectTracks( e_t_4_d, speed, 2, Inf )


#Sampling from all different replicates

bt1.random <- b_t_1[sample(1:length(b_t_1), 70, replace = F)]
bt2.random <- b_t_2[sample(1:length(b_t_2), 70, replace = F)]
bt4.random <- b_t_4[sample(1:length(b_t_4), 70, replace = F)]
bt.random <- c(bt1.random, bt2.random, bt4.random)

ct1.random <- c_t_1[sample(1:length(c_t_1), 70, replace = F)]
ct2.random <- c_t_2[sample(1:length(c_t_2), 70, replace = F)]
ct4.random <- c_t_4[sample(1:length(c_t_4), 70, replace = F)]
ct.random <- c(ct1.random, ct2.random, ct4.random)

dt1.random <- d_t_1[sample(1:length(d_t_1), 70, replace = F)]
dt2.random <- d_t_2[sample(1:length(d_t_2), 70, replace = F)]
dt4.random <- d_t_4[sample(1:length(d_t_4), 70, replace = F)]
dt.random <- c(dt1.random, dt2.random, dt4.random)

et1.random <- e_t_1[sample(1:length(e_t_1), 70, replace = F)]
et2.random <- e_t_2[sample(1:length(e_t_2), 70, replace = F)]
et4.random <- e_t_4[sample(1:length(e_t_4), 70, replace = F)]
et.random <- c(et1.random, et2.random, et4.random)

b_1_fast.tracks.random <- b_1_fast.tracks[sample(1:length(b_1_fast.tracks), 50, replace = F)]
b_2_fast.tracks.random <- b_2_fast.tracks[sample(1:length(b_2_fast.tracks), 50, replace = F)]
b_4_fast.tracks.random <- b_4_fast.tracks[sample(1:length(b_4_fast.tracks), 50, replace = F)]
bt.fast.random <- c(b_1_fast.tracks.random, b_2_fast.tracks.random, b_4_fast.tracks.random)

c_1_fast.tracks.random <- c_1_fast.tracks[sample(1:length(c_1_fast.tracks), 50, replace = F)]
c_2_fast.tracks.random <- c_2_fast.tracks[sample(1:length(c_2_fast.tracks), 50, replace = F)]
c_4_fast.tracks.random <- c_4_fast.tracks[sample(1:length(c_4_fast.tracks), 50, replace = F)]
ct.fast.random <- c(c_1_fast.tracks.random, c_2_fast.tracks.random, c_4_fast.tracks.random)

d_1_fast.tracks.random <- d_1_fast.tracks[sample(1:length(d_1_fast.tracks), 50, replace = F)]
d_2_fast.tracks.random <- d_2_fast.tracks[sample(1:length(d_2_fast.tracks), 50, replace = F)]
d_4_fast.tracks.random <- d_4_fast.tracks[sample(1:length(d_4_fast.tracks), 50, replace = F)]
dt.fast.random <- c(d_1_fast.tracks.random, d_2_fast.tracks.random, d_4_fast.tracks.random)

e_1_fast.tracks.random <- e_1_fast.tracks[sample(1:length(e_1_fast.tracks), 50, replace = F)]
e_2_fast.tracks.random <- e_2_fast.tracks[sample(1:length(e_2_fast.tracks), 50, replace = F)]
e_4_fast.tracks.random <- e_4_fast.tracks[sample(1:length(e_4_fast.tracks), 50, replace = F)]
et.fast.random <- c(e_1_fast.tracks.random, e_2_fast.tracks.random, e_4_fast.tracks.random)

#saving template tracks

write.csv(normalizeTracks(bt.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/bt.random.csv")
write.csv(normalizeTracks(ct.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/ct.random.csv")
write.csv(normalizeTracks(dt.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/dt.random.csv")
write.csv(normalizeTracks(et.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/et.random.csv")

write.csv(normalizeTracks(bt.fast.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/bt.fast.random.csv")
write.csv(normalizeTracks(ct.fast.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/ct.fast.random.csv")
write.csv(normalizeTracks(dt.fast.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/dt.fast.random.csv")
write.csv(normalizeTracks(et.fast.random),"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/et.fast.random.csv")




#####tracks starplot


par(mfrow = c(2,2))
plot( normalizeTracks(bt.random), main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(ct.random), main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(dt.random), main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(et.random), main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))
dev.off()

par(mfrow = c(2,2))
plot( normalizeTracks(bt.fast.random), main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(ct.fast.random), main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(dt.fast.random), main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(et.fast.random), main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))
dev.off()


par(mfrow = c(2,2))
plot( normalizeTracks(b_fast.tracks.random), main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(c_fast.tracks.random), main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(d_fast.tracks.random), main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(e_fast.tracks.random), main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))
dev.off()





############# template tracks used imported for usage

bt.random.csv <- fread("bt.random.csv", header = T)
ct.random.csv <- fread("ct.random.csv", header = T)
dt.random.csv <- fread("dt.random.csv", header = T)
et.random.csv <- fread("et.random.csv", header = T)



b_t2 <- bt.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))

c_t2 <- ct.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))

d_t2 <- dt.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))

e_t2 <- et.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))


par(mfrow = c(2,2))
plot( b_t2, main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( c_t2, main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( d_t2, main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( e_t2, main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))






bt1.random <- b_t_1[sample(1:length(b_t_1), 70, replace = F)]
bt2.random <- b_t_2[sample(1:length(b_t_2), 70, replace = F)]
bt4.random <- b_t_4[sample(1:length(b_t_4), 70, replace = F)]
bt.random <- c(bt1.random, bt2.random, bt4.random)

ct1.random <- c_t_1[sample(1:length(c_t_1), 70, replace = F)]
ct2.random <- c_t_2[sample(1:length(c_t_2), 70, replace = F)]
ct4.random <- c_t_4[sample(1:length(c_t_4), 70, replace = F)]
ct.random <- c(ct1.random, ct2.random, ct4.random)

dt1.random <- d_t_1[sample(1:length(d_t_1), 70, replace = F)]
dt2.random <- d_t_2[sample(1:length(d_t_2), 70, replace = F)]
dt4.random <- d_t_4[sample(1:length(d_t_4), 70, replace = F)]
dt.random <- c(dt1.random, dt2.random, dt4.random)

et1.random <- e_t_1[sample(1:length(e_t_1), 70, replace = F)]
et2.random <- e_t_2[sample(1:length(e_t_2), 70, replace = F)]
et4.random <- e_t_4[sample(1:length(e_t_4), 70, replace = F)]
et.random <- c(et1.random, et2.random, et4.random)

b_1_fast.tracks.random <- b_1_fast.tracks[sample(1:length(b_1_fast.tracks), 50, replace = F)]
b_2_fast.tracks.random <- b_2_fast.tracks[sample(1:length(b_2_fast.tracks), 50, replace = F)]
b_4_fast.tracks.random <- b_4_fast.tracks[sample(1:length(b_4_fast.tracks), 50, replace = F)]
bt.fast.random <- c(b_1_fast.tracks.random, b_2_fast.tracks.random, b_4_fast.tracks.random)

c_1_fast.tracks.random <- c_1_fast.tracks[sample(1:length(c_1_fast.tracks), 50, replace = F)]
c_2_fast.tracks.random <- c_2_fast.tracks[sample(1:length(c_2_fast.tracks), 50, replace = F)]
c_4_fast.tracks.random <- c_4_fast.tracks[sample(1:length(c_4_fast.tracks), 50, replace = F)]
ct.fast.random <- c(c_1_fast.tracks.random, c_2_fast.tracks.random, c_4_fast.tracks.random)

d_1_fast.tracks.random <- d_1_fast.tracks[sample(1:length(d_1_fast.tracks), 50, replace = F)]
d_2_fast.tracks.random <- d_2_fast.tracks[sample(1:length(d_2_fast.tracks), 50, replace = F)]
d_4_fast.tracks.random <- d_4_fast.tracks[sample(1:length(d_4_fast.tracks), 50, replace = F)]
dt.fast.random <- c(d_1_fast.tracks.random, d_2_fast.tracks.random, d_4_fast.tracks.random)

e_1_fast.tracks.random <- e_1_fast.tracks[sample(1:length(e_1_fast.tracks), 50, replace = F)]
e_2_fast.tracks.random <- e_2_fast.tracks[sample(1:length(e_2_fast.tracks), 50, replace = F)]
e_4_fast.tracks.random <- e_4_fast.tracks[sample(1:length(e_4_fast.tracks), 50, replace = F)]
et.fast.random <- c(e_1_fast.tracks.random, e_2_fast.tracks.random, e_4_fast.tracks.random)


view(c_fast.tracks.random)
write.csv(b_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/b_fast.tracks.random.csv")
write.csv(c_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/c_fast.tracks.random.csv")
write.csv(d_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/d_fast.tracks.random.csv")
write.csv(e_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/e_fast.tracks.random.csv")

write.csv(bt.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/bt.random.csv")
write.csv(ct.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/ct.random.csv")
write.csv(dt.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/dt.random.csv")
write.csv(et.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/et.random.csv")

write.csv(bt.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/bt.fast.random.csv")
write.csv(ct.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/ct.fast.random.csv")
write.csv(dt.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/dt.fast.random.csv")
write.csv(et.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_3/et.fast.random.csv")




########################### msd

msddata <- rbind(b_1.msd, b_2.msd, b_4.msd, c_1.msd, c_2.msd, c_4.msd, d_1.msd, d_2.msd, d_4.msd,
                 e_1.msd, e_2.msd, e_4.msd)


msddata <- rbind( b_t_1.msd, b_t_2.msd, b_t_3.msd, b_t_4.msd, c_t_1.msd, c_t_2.msd, c_t_3.msd, c_t_4.msd, d_t_1.msd, d_t_2.msd, d_t_3.msd, d_t_4.msd, e_t_1.msd, e_t_2.msd, e_t_3.msd, e_t_4.msd)
msddata <- rbind( c_t_1.msd, c_t_2.msd, c_t_3.msd, c_t_4.msd)



msdata1 <- rbind(b_1.msd, c_1.msd, d_1.msd, e_1.msd)

msdata2 <- rbind(b_2.msd, c_2.msd, d_2.msd, e_2.msd)

msdata4 <- rbind(b_4.msd, c_4.msd, d_4.msd, e_4.msd)

msdata <- rbind(b_1.msd, b_2.msd, b_4.msd, c_1.msd, c_2.msd, c_4.msd, d_1.msd, d_2.msd, d_4.msd,
                e_1.msd, e_2.msd, e_4.msd)

msdata <- rbind(b_1.msd, b_2.msd, b_4.msd, c_1.msd, c_2.msd, c_4.msd, d_1.msd, d_2.msd, d_4.msd,
                e_1.msd, e_2.msd, e_4.msd)

msdatat <- rbind(bt.random.msd, ct.random.msd, dt.random.msd, et.random.msd)

msdatatf <- rbind(bt.fast.random.msd, ct.fast.random.msd, dt.fast.random.msd, et.fast.random.msd)



p1 <- ggplot(msdatatf, aes( x = dt*15 , y = mean, color = condition, fill = condition)) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line() +
  labs( x = expression( paste(Delta,"t (min)") ),
        y = "mean square displacement") 
#  stat_compare_means(aes(group = condition), label = "p.signif", 
#                     label.y = c(16, 25, 29))

p1 + scale_x_continuous(limits = c(30, 290)) + theme_classic()




# Speed

bt.random.speeds <- sapply( bt.random, speed )
str( bt.random.speeds )
summary( bt.random.speeds )

ct.random.speeds <- sapply( ct.random, speed )
str( ct.random.speeds )
summary( ct.random.speeds )

dt.random.speeds <- sapply( dt.random, speed )
str( dt.random.speeds )
summary( dt.random.speeds )

et.random.speeds <- sapply( et.random, speed )
str( et.random.speeds )
summary( et.random.speeds )

bt.fast.random.speeds <- sapply( bt.fast.random, speed )
str( bt.fast.random.speeds )
summary( bt.fast.random.speeds )

ct.fast.random.speeds <- sapply( ct.fast.random, speed )
str( ct.fast.random.speeds )
summary( ct.fast.random.speeds )

dt.fast.random.speeds <- sapply( dt.fast.random, speed )
str( dt.fast.random.speeds )
summary( dt.fast.random.speeds )

et.fast.random.speeds <- sapply( et.fast.random, speed )
str( et.fast.random.speeds )
summary( et.fast.random.speeds )


db_t <- data.frame(speed = bt.random.speeds, condition = "SC")
dc_t <- data.frame(speed = ct.random.speeds, condition = "SC+Wnt")
dd_t <- data.frame(speed = dt.random.speeds, condition = "SC+PC")
de_t <- data.frame(speed = et.random.speeds, condition = "SC+PC+i")


#####################Speed fast

db_t_f <- data.frame(speed = bt.fast.random.speeds, condition = "SC")
dc_t_f <- data.frame(speed = ct.fast.random.speeds, condition = "SC+Wnt")
dd_t_f <- data.frame(speed = dt.fast.random.speeds, condition = "SC+PC")
de_t_f <- data.frame(speed = et.fast.random.speeds, condition = "SC+PC+i")


d <- rbind(db_t, dc_t, dd_t,de_t)
d_f <- rbind(db_t_f, dc_t_f, dd_t_f,de_t_f)


b_mean.speed <- mean(bt.random.speeds)
c_mean.speed <- mean(ct.random.speeds)
d_mean.speed <- mean(dt.random.speeds)
e_mean.speed <- mean(et.random.speeds)

b_mean.speed.fast <- mean(bt.fast.random.speeds)
c_mean.speed.fast<- mean(ct.fast.random.speeds)
d_mean.speed.fast <- mean(dt.fast.random.speeds)
e_mean.speed.fast <- mean(et.fast.random.speeds)

b_median.speed <- median(bt.random.speeds)
c_median.speed <- median(ct.random.speeds)
d_median.speed <- median(dt.random.speeds)
e_median.speed <- median(et.random.speeds)

b_median.speed.fast <- median(bt.fast.random.speeds)
c_median.speed.fast<- median(ct.fast.random.speeds)
d_median.speed.fast <- median(dt.fast.random.speeds)
e_median.speed.fast <- median(et.fast.random.speeds)

b_std <- sd(bt.random.speeds)/sqrt(length(bt.random.speeds))
c_std <- sd(ct.random.speeds)/sqrt(length(ct.random.speeds))
d_std <- sd(dt.random.speeds)/sqrt(length(dt.random.speeds))
e_std <- sd(et.random.speeds)/sqrt(length(et.random.speeds))

b_std.fast <- sd(bt.fast.random.speeds)/sqrt(length(bt.fast.random.speeds))
c_std.fast <- sd(ct.fast.random.speeds)/sqrt(length(ct.fast.random.speeds))
d_std.fast <- sd(dt.fast.random.speeds)/sqrt(length(dt.fast.random.speeds))
e_std.fast <- sd(et.fast.random.speeds)/sqrt(length(et.fast.random.speeds))



bt1.random.speeds <- sapply( bt1.random, speed )
str( bt1.random.speeds )
summary( bt1.random.speeds )

bt2.random.speeds <- sapply( bt2.random, speed )
str( bt2.random.speeds )
summary( bt2.random.speeds )

bt4.random.speeds <- sapply( bt4.random, speed )
str( bt4.random.speeds )
summary( bt4.random.speeds )

ct1.random.speeds <- sapply( ct1.random, speed )
str( ct1.random.speeds )
summary( ct1.random.speeds )

ct2.random.speeds <- sapply( ct2.random, speed )
str( ct2.random.speeds )
summary( ct2.random.speeds )

ct4.random.speeds <- sapply( ct4.random, speed )
str( ct4.random.speeds )
summary( ct4.random.speeds )

dt1.random.speeds <- sapply( dt1.random, speed )
str( dt1.random.speeds )
summary( dt1.random.speeds )

dt2.random.speeds <- sapply( dt2.random, speed )
str( dt2.random.speeds )
summary( dt2.random.speeds )

dt4.random.speeds <- sapply( dt4.random, speed )
str( dt4.random.speeds )
summary( dt4.random.speeds )

et1.random.speeds <- sapply( et1.random, speed )
str( et1.random.speeds )
summary( et1.random.speeds )

et2.random.speeds <- sapply( et2.random, speed )
str( et2.random.speeds )
summary( et2.random.speeds )

et4.random.speeds <- sapply( et4.random, speed )
str( et4.random.speeds )
summary( et4.random.speeds )

b_1_fast.tracks.random.speeds <- sapply( b_1_fast.tracks.random, speed )
str( b_1_fast.tracks.random.speeds )
summary( b_1_fast.tracks.random.speeds )

b_2_fast.tracks.random.speeds <- sapply( b_2_fast.tracks.random, speed )
str( b_2_fast.tracks.random.speeds )
summary( b_2_fast.tracks.random.speeds )

b_4_fast.tracks.random.speeds <- sapply( b_4_fast.tracks.random, speed )
str( b_4_fast.tracks.random.speeds )
summary( b_4_fast.tracks.random.speeds )

c_1_fast.tracks.random.speeds <- sapply( c_1_fast.tracks.random, speed )
str( c_1_fast.tracks.random.speeds )
summary( c_1_fast.tracks.random.speeds )

c_2_fast.tracks.random.speeds <- sapply( c_2_fast.tracks.random, speed )
str( c_2_fast.tracks.random.speeds )
summary( c_2_fast.tracks.random.speeds )

c_4_fast.tracks.random.speeds <- sapply( c_4_fast.tracks.random, speed )
str( c_4_fast.tracks.random.speeds )
summary( c_4_fast.tracks.random.speeds )

d_1_fast.tracks.random.speeds <- sapply( d_1_fast.tracks.random, speed )
str( d_1_fast.tracks.random.speeds )
summary( d_1_fast.tracks.random.speeds )

d_2_fast.tracks.random.speeds <- sapply( d_2_fast.tracks.random, speed )
str( d_2_fast.tracks.random.speeds )
summary( d_2_fast.tracks.random.speeds )

d_4_fast.tracks.random.speeds <- sapply( d_4_fast.tracks.random, speed )
str( d_4_fast.tracks.random.speeds )
summary( d_4_fast.tracks.random.speeds )

e_1_fast.tracks.random.speeds <- sapply( e_1_fast.tracks.random, speed )
str( e_1_fast.tracks.random.speeds )
summary( e_1_fast.tracks.random.speeds )

e_2_fast.tracks.random.speeds <- sapply( e_2_fast.tracks.random, speed )
str( e_2_fast.tracks.random.speeds )
summary( e_2_fast.tracks.random.speeds )

e_4_fast.tracks.random.speeds <- sapply( e_4_fast.tracks.random, speed )
str( e_4_fast.tracks.random.speeds )
summary( e_4_fast.tracks.random.speeds )

b1_mean.speed <- mean(bt1.random.speeds)
b2_mean.speed <- mean(bt2.random.speeds)
b4_mean.speed <- mean(bt4.random.speeds)

c1_mean.speed <- mean(ct1.random.speeds)
c2_mean.speed <- mean(ct2.random.speeds)
c4_mean.speed <- mean(ct4.random.speeds)

d1_mean.speed <- mean(dt1.random.speeds)
d2_mean.speed <- mean(dt2.random.speeds)
d4_mean.speed <- mean(dt4.random.speeds)

e1_mean.speed <- mean(et1.random.speeds)
e2_mean.speed <- mean(et2.random.speeds)
e4_mean.speed <- mean(et4.random.speeds)

b1_mean.speed.fast <- mean(b_1_fast.tracks.random.speeds)
b2_mean.speed.fast <- mean(b_2_fast.tracks.random.speeds)
b4_mean.speed.fast <- mean(b_4_fast.tracks.random.speeds)

c1_mean.speed.fast <- mean(c_1_fast.tracks.random.speeds)
c2_mean.speed.fast <- mean(c_2_fast.tracks.random.speeds)
c4_mean.speed.fast <- mean(c_4_fast.tracks.random.speeds)

d1_mean.speed.fast <- mean(d_1_fast.tracks.random.speeds)
d2_mean.speed.fast <- mean(d_2_fast.tracks.random.speeds)
d4_mean.speed.fast <- mean(d_4_fast.tracks.random.speeds)

e1_mean.speed.fast <- mean(e_1_fast.tracks.random.speeds)
e2_mean.speed.fast <- mean(e_2_fast.tracks.random.speeds)
e4_mean.speed.fast <- mean(e_4_fast.tracks.random.speeds)

b1_median.speed <- median(bt1.random.speeds)
b2_median.speed <- median(bt2.random.speeds)
b4_median.speed <- median(bt4.random.speeds)

c1_median.speed <- median(ct1.random.speeds)
c2_median.speed <- median(ct2.random.speeds)
c4_median.speed <- median(ct4.random.speeds)

d1_median.speed <- median(dt1.random.speeds)
d2_median.speed <- median(dt2.random.speeds)
d4_median.speed <- median(dt4.random.speeds)

e1_median.speed <- median(et1.random.speeds)
e2_median.speed <- median(et2.random.speeds)
e4_median.speed <- median(et4.random.speeds)

b1_median.speed.fast <- median(b_1_fast.tracks.random.speeds)
b2_median.speed.fast <- median(b_2_fast.tracks.random.speeds)
b4_median.speed.fast <- median(b_4_fast.tracks.random.speeds)

c1_median.speed.fast <- median(c_1_fast.tracks.random.speeds)
c2_median.speed.fast <- median(c_2_fast.tracks.random.speeds)
c4_median.speed.fast <- median(c_4_fast.tracks.random.speeds)

d1_median.speed.fast <- median(d_1_fast.tracks.random.speeds)
d2_median.speed.fast <- median(d_2_fast.tracks.random.speeds)
d4_median.speed.fast <- median(d_4_fast.tracks.random.speeds)

e1_median.speed.fast <- median(e_1_fast.tracks.random.speeds)
e2_median.speed.fast <- median(e_2_fast.tracks.random.speeds)
e4_median.speed.fast <- median(e_4_fast.tracks.random.speeds)

b1_std <- sd(bt1.random.speeds)/sqrt(length(bt1.random.speeds))
b2_std <- sd(bt2.random.speeds)/sqrt(length(bt2.random.speeds))
b4_std <- sd(bt4.random.speeds)/sqrt(length(bt4.random.speeds))

c1_std <- sd(ct1.random.speeds)/sqrt(length(ct1.random.speeds))
c2_std <- sd(ct2.random.speeds)/sqrt(length(ct2.random.speeds))
c4_std <- sd(ct4.random.speeds)/sqrt(length(ct4.random.speeds))

d1_std <- sd(dt1.random.speeds)/sqrt(length(dt1.random.speeds))
d2_std <- sd(dt2.random.speeds)/sqrt(length(dt2.random.speeds))
d4_std <- sd(dt4.random.speeds)/sqrt(length(dt4.random.speeds))

e1_std <- sd(et1.random.speeds)/sqrt(length(et1.random.speeds))
e2_std <- sd(et2.random.speeds)/sqrt(length(et2.random.speeds))
e4_std <- sd(et4.random.speeds)/sqrt(length(et4.random.speeds))

b1_std.fast <- sd(b_1_fast.tracks.random.speeds)/sqrt(length(b_1_fast.tracks.random.speeds))
b2_std.fast <- sd(b_2_fast.tracks.random.speeds)/sqrt(length(b_2_fast.tracks.random.speeds))
b4_std.fast <- sd(b_4_fast.tracks.random.speeds)/sqrt(length(b_4_fast.tracks.random.speeds))

c1_std.fast <- sd(c_1_fast.tracks.random.speeds)/sqrt(length(c_1_fast.tracks.random.speeds))
c2_std.fast <- sd(c_2_fast.tracks.random.speeds)/sqrt(length(c_2_fast.tracks.random.speeds))
c4_std.fast <- sd(c_4_fast.tracks.random.speeds)/sqrt(length(c_4_fast.tracks.random.speeds))

d1_std.fast <- sd(d_1_fast.tracks.random.speeds)/sqrt(length(d_1_fast.tracks.random.speeds))
d2_std.fast <- sd(d_2_fast.tracks.random.speeds)/sqrt(length(d_2_fast.tracks.random.speeds))
d4_std.fast <- sd(d_4_fast.tracks.random.speeds)/sqrt(length(d_4_fast.tracks.random.speeds))

e1_std.fast <- sd(e_1_fast.tracks.random.speeds)/sqrt(length(e_1_fast.tracks.random.speeds))
e2_std.fast <- sd(e_2_fast.tracks.random.speeds)/sqrt(length(e_2_fast.tracks.random.speeds))
e4_std.fast <- sd(e_4_fast.tracks.random.speeds)/sqrt(length(e_4_fast.tracks.random.speeds))

dmean <- data.frame(mean = c(b1_mean.speed, b2_mean.speed, b4_mean.speed, c1_mean.speed, c2_mean.speed, c4_mean.speed,
                             d1_mean.speed, d2_mean.speed, d4_mean.speed, e1_mean.speed, e2_mean.speed, e4_mean.speed), 
                    condition = c("SC", "SC", "SC", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+PC", "SC+PC", "SC+PC", 
                                  "SC+PC+i", "SC+PC+i", "SC+PC+i"), replicate = c("SC_1", "SC_2", "SC_3", "SC+Wnt_1", "SC+Wnt_2", "SC+Wnt_4", "SC+PC_1", "SC+PC_2", "SC+PC_4", 
                                                                                  "SC+PC+i_1", "SC+PC+i_2", "SC+PC+i_4"))

dmedian <- data.frame(median = c(b1_median.speed, b2_median.speed, b4_median.speed, c1_median.speed, c2_median.speed, c4_median.speed,
                                 d1_median.speed, d2_median.speed, d4_median.speed, e1_median.speed, e2_median.speed, e4_median.speed), 
                      condition = c("SC", "SC", "SC", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+PC", "SC+PC", "SC+PC", 
                                    "SC+PC+i", "SC+PC+i", "SC+PC+i"), replicate = c("SC_1", "SC_2", "SC_3", "SC+Wnt_1", "SC+Wnt_2", "SC+Wnt_4", "SC+PC_1", "SC+PC_2", "SC+PC_4", 
                                                                                    "SC+PC+i_1", "SC+PC+i_2", "SC+PC+i_4"))

dmeanf <- data.frame(mean = c(b1_mean.speed.fast, b2_mean.speed.fast, b4_mean.speed.fast, c1_mean.speed.fast, c2_mean.speed.fast, c4_mean.speed.fast,
                              d1_mean.speed.fast, d2_mean.speed.fast, d4_mean.speed.fast, e1_mean.speed.fast, e2_mean.speed.fast, e4_mean.speed.fast), 
                     condition = c("SC", "SC", "SC", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+PC", "SC+PC", "SC+PC", 
                                   "SC+PC+i", "SC+PC+i", "SC+PC+i"), replicate = c("SC_1", "SC_2", "SC_3", "SC+Wnt_1", "SC+Wnt_2", "SC+Wnt_4", "SC+PC_1", "SC+PC_2", "SC+PC_4", 
                                                                                   "SC+PC+i_1", "SC+PC+i_2", "SC+PC+i_4"))

dmedianf <- data.frame(median = c(b1_median.speed.fast, b2_median.speed.fast, b4_median.speed.fast, c1_median.speed.fast, c2_median.speed.fast, c4_median.speed.fast,
                                  d1_median.speed.fast, d2_median.speed.fast, d4_median.speed.fast, e1_median.speed.fast, e2_median.speed.fast, e4_median.speed.fast), 
                       condition = c("SC", "SC", "SC", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+PC", "SC+PC", "SC+PC", 
                                     "SC+PC+i", "SC+PC+i", "SC+PC+i"), replicate = c("SC_1", "SC_2", "SC_3", "SC+Wnt_1", "SC+Wnt_2", "SC+Wnt_4", "SC+PC_1", "SC+PC_2", "SC+PC_4", 
                                                                                     "SC+PC+i_1", "SC+PC+i_2", "SC+PC+i_4"))
dc_t <- data.frame(speed = ct.random.speeds, condition = "SC+Wnt")
dd_t <- data.frame(speed = dt.random.speeds, condition = "SC+PC")
de_t <- data.frame(speed = et.random.speeds, condition = "SC+PC+i")

mmean <- data.frame(mean=c(b_mean.speed, c_mean.speed, d_mean.speed, e_mean.speed), 
                    sem=c(b_std, c_std, d_std, e_std), condition= c("SC","SC+Wnt", "SC+PC", "SC+PC+i","SC", "SC+Wnt", "SC+PC", "SC+PC+i"))

mmedian <- data.frame(median=c(b_median.speed, c_median.speed, d_median.speed, e_median.speed), 
                      sem=c(b_std, c_std, d_std, e_std), condition= c("SC","SC+Wnt", "SC+PC", "SC+PC+i","SC", "SC+Wnt", "SC+PC", "SC+PC+i"))

mmeanf <- data.frame(mean=c(b_mean.speed.fast, c_mean.speed.fast, d_mean.speed.fast, e_mean.speed.fast), 
                     sem=c(b_std.fast, c_std.fast, d_std.fast, e_std.fast), condition= c("SC","SC+Wnt", "SC+PC", "SC+PC+i","SC", "SC+Wnt", "SC+PC", "SC+PC+i"))

mmedianf <- data.frame(median=c(b_median.speed.fast, c_median.speed.fast, d_median.speed.fast, e_median.speed.fast), 
                       sem=c(b_std.fast, c_std.fast, d_std.fast, e_std.fast), condition= c("SC","SC+Wnt", "SC+PC", "SC+PC+i","SC", "SC+Wnt", "SC+PC", "SC+PC+i"))


mp <- ggplot(mmeanf[1:4,], aes(x=condition, y= mean, fill=condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sem), width=.2, position = position_dodge(.9)) +
  geom_quasirandom(dmeanf, mapping = aes(x=condition, y=mean, fill=condition), method= "quasirandom", dodge.width=2,size = 2,alpha=.08, color="black") +
  theme_classic() +
  ylim(0,3)

mp

mp <- ggplot(mmedianf[1:4,], aes(x=condition, y= median, fill=condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=median, ymax=median+sem), width=.2, position = position_dodge(.9)) +
  geom_quasirandom(dmedianf, mapping = aes(x=condition, y=median, fill=condition), method= "quasirandom", dodge.width=2,size = 2,alpha=.08, color="black") +
  theme_classic() +
  ylim(0,2.5)

mp




#directionality
#based on DiPer

b_dir <- data.frame(directionality_ratio = c(1, 0.662080078, 0.549937134, 0.477047374, 0.421496178, 0.391436435, 0.377909197, 0.354733952, 0.340256607,
                                             0.322398612, 0.311984215, 0.30319092, 0.299231629, 0.28981961, 0.274961877, 0.273949291, 0.274121357,
                                             0.256631378, 0.252964107), 
                    SEM = c(0,0.024008245, 0.02202398, 0.020348071, 0.019254861, 0.017949612, 0.017461004, 0.017635317,
                            0.018067597, 0.017432336, 0.017316941, 0.017097305, 0.018026723, 0.017113889, 0.017269495,
                            0.018356142, 0.018208373, 0.018010356, 0.017760944),
                    time = c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                             195, 210, 225, 240, 255, 270), 
                    condition = c("SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC",
                                  "SC", "SC", "SC", "SC", "SC", "SC", "SC"))

c_dir <- data.frame(directionality_ratio = c(1, 0.713435707, 0.566279665, 0.538154148, 0.496925817, 0.471110265, 0.442941256,
                                             0.434993633, 0.418104729, 0.396264801, 0.381229457, 0.376790935, 0.362815537,
                                             0.35512762, 0.346037719, 0.34024185, 0.331745227, 0.320669978, 0.315035353),
                    SEM = c(0, 0.021347478, 0.022899443, 0.021893858, 0.020800602, 0.020637294, 0.02039385, 0.02038432,
                            0.01983696, 0.019485707, 0.019665489, 0.019607324, 0.019569159, 0.019337261, 0.019264371,
                            0.019436864, 0.019712482, 0.01965481, 0.019165852),
                    time= c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                            195, 210, 225, 240, 255, 270), 
                    condition = c("SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt",
                                  "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt"))

d_dir <- data.frame(directionality_ratio = c(1, 0.689662153, 0.531030501, 0.485461547, 0.444981635, 0.422936513, 0.407079627,
                                             0.393273055, 0.366400773, 0.342393511, 0.322826438, 0.328674614, 0.314934903,
                                             0.304794333, 0.301269991, 0.280817085, 0.273382941, 0.277379593, 0.270572162),
                    SEM = c(0, 0.023200034, 0.020664055, 0.019870195, 0.019617406, 0.01908628, 0.018402505, 0.018546908,
                            0.018822209, 0.018450514, 0.019240077, 0.01870074, 0.019048253, 0.019257085, 0.019383445,
                            0.01846525, 0.018858018, 0.018719036, 0.018591558),
                    time= c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                            195, 210, 225, 240, 255, 270), 
                    condition = c("SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", 
                                  "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC"))

e_dir <- data.frame(directionality_ratio = c(1, 0.701667452, 0.537200346, 0.478979771, 0.426837702, 0.380889773, 0.37340414,
                                             0.355480521, 0.339227028, 0.327367354, 0.314174523, 0.293544139, 0.291157575,
                                             0.283808345, 0.269883369, 0.263644404, 0.262993466, 0.255346684, 0.25420797),
                    SEM = c(0, 0.021548833, 0.021612329, 0.019647063, 0.018738367, 0.017977935, 0.018052047, 0.01790108,
                            0.016714176, 0.015967777, 0.016481424, 0.016338921, 0.016088542, 0.016657846, 0.016784674,
                            0.017088774, 0.01685693, 0.017139392, 0.017512149), 
                    time= c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                            195, 210, 225, 240, 255, 270), 
                    condition = c("SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i",
                                  "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i"))


dir <- rbind(b_dir, c_dir, d_dir, e_dir)


p1 <- ggplot( dir, aes( x = time , y = directionality_ratio, color = condition, fill = condition ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = directionality_ratio - SEM, ymax = directionality_ratio + SEM) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = "time(min)",
        y = "directionality ratio" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )
