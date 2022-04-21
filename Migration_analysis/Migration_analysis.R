########For figure 3a, b and extended figure 7c

############ Import cell tracks

getwd()
setwd("I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new")

setwd("I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3")

tracks_koen <- map(list.files(), function(x){
  fread(x, header = T)})

file_paths <- fs::dir_ls("I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_4")

tracks_koen <- set_names(tracks_koen, file_paths)



############### classify samples and convert to tracks 
str(gs)
gs <- map(tracks_koen, function(x){
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


############# filtering based on speed


b_fast.tracks <- selectTracks( b_t, speed, 0.3, Inf )
b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.tracks <- selectTracks( c_t, speed, 0.3, Inf )
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.tracks <- selectTracks( d_t, speed, 0.3, Inf )
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.tracks <- selectTracks( e_t, speed, 0.3, Inf )
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]




b_fast.tracks <- selectTracks( b_t, speed, 0.3, Inf )
b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.tracks <- selectTracks( c_t, speed, 0.3, Inf )
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.tracks <- selectTracks( d_t, speed, 0.3, Inf )
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.tracks <- selectTracks( e_t, speed, 0.3, Inf )
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]


b_fast.fast.tracks <- selectTracks( b_t, speed, 2, Inf )
b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.fast.tracks <- selectTracks( c_t, speed, 2, Inf )
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.fast.tracks <- selectTracks( d_t, speed, 2, Inf )
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.fast.tracks <- selectTracks( e_t, speed, 2, Inf )
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]



par(mfrow = c(2,2))
plot( normalizeTracks(b_fast.tracks), main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(c_fast.tracks), main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(d_fast.tracks), main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(e_fast.tracks), main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))
dev.off()


par(mfrow = c(2,2))
plot( normalizeTracks(b_fast.fast.tracks), main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(c_fast.fast.tracks), main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(d_fast.fast.tracks), main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(e_fast.fast.tracks), main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))
dev.off()


############# template tracks used imported for usage

bt.random.csv <- fread("bt.random.csv", header = T)
ct.random.csv <- fread("ct.random.csv", header = T)
dt.random.csv <- fread("dt.random.csv", header = T)
et.random.csv <- fread("et.random.csv", header = T)



b_t <- bt.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))

c_t <- ct.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))

d_t <- dt.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))

e_t <- et.random.csv %>%
  as.tracks(id.column = 2,
            time.column = 3,
            pos.columns = c(4, 5))






#####tracks

par(mfrow = c(2,2))
plot( normalizeTracks(b_t)), main = "SC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(c_t), main = "SC+Wnt", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(d_t), main = "SC+PC", xlim=c(-80,80), ylim=c(-80,80))
plot( normalizeTracks(e_t), main = "SC+PC+i", xlim=c(-80,80), ylim=c(-80,80))




b_t <- c(gs[[1]], gs[[2]], gs[[11]], gs[[12]], gs[[21]], gs[[22]], gs[[33]], gs[[34]], gs[[35]], gs[[36]])
c_t <- c(gs[[4]], gs[[13]], gs[[14]], gs[[23]], gs[[37]], gs[[38]], gs[[39]], gs[[40]])
d_t <- c(gs[[5]], gs[[6]], gs[[15]], gs[[16]], gs[[24]], gs[[25]], gs[[26]], gs[[41]], gs[[42]], gs[[43]])
e_t <- c(gs[[7]], gs[[8]], gs[[17]], gs[[18]], gs[[27]], gs[[28]], gs[[29]], gs[[44]], gs[[45]], gs[[46]], gs[[47]])

############### subtract samples to velocity and same number

wer <- sample(1:length(b_fast.tracks), 300, replace = F)

length(b_t)
length(b_fast.tracks)
length(b_fast.fast.tracks)
length(b_fast.tracks) - length(b_fast.fast.tracks)
length(b_t) - length(b_fast.tracks)

length(c_t)
length(c_fast.tracks)
length(c_fast.fast.tracks)
length(c_fast.tracks) - length(c_fast.fast.tracks)
length(c_t) - length(c_fast.tracks)

length(d_t)
length(d_fast.tracks)
length(d_fast.fast.tracks)
length(d_fast.tracks) - length(d_fast.fast.tracks)
length(d_t) - length(d_fast.tracks)


length(e_t)
length(e_fast.tracks)
length(e_fast.fast.tracks)
length(e_fast.tracks) - length(e_fast.fast.tracks)
length(e_t) - length(e_fast.tracks)


b_fast.tracks <- selectTracks( b_t, speed, 0.3, Inf )
b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.tracks <- selectTracks( c_t, speed, 0.3, Inf )
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.tracks <- selectTracks( d_t, speed, 0.3, Inf )
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.tracks <- selectTracks( e_t, speed, 0.3, Inf )
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]


b_fast.fast.tracks <- selectTracks( b_t, speed, 2, Inf )
b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.fast.tracks <- selectTracks( c_t, speed, 2, Inf )
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.fast.tracks <- selectTracks( d_t, speed, 2, Inf )
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.fast.tracks <- selectTracks( e_t, speed, 2, Inf )
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]

b_fast.tracks.random <- b_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
c_fast.tracks.random <- c_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
d_fast.tracks.random <- d_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]
e_fast.tracks.random <- e_fast.tracks[sample(1:length(b_fast.tracks), 300, replace = F)]

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
write.csv(b_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/b_fast.tracks.random.csv")
write.csv(c_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/c_fast.tracks.random.csv")
write.csv(d_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/d_fast.tracks.random.csv")
write.csv(e_fast.tracks.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/e_fast.tracks.random.csv")

write.csv(bt.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/bt.random.csv")
write.csv(ct.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/ct.random.csv")
write.csv(dt.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/dt.random.csv")
write.csv(et.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/et.random.csv")

write.csv(bt.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/bt.fast.random.csv")
write.csv(ct.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/ct.fast.random.csv")
write.csv(dt.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/dt.fast.random.csv")
write.csv(et.fast.random,"I:/Group Rheenen/ExpDATA/DL/Analysis/SIvLI paper/migration/new/tracks_koen_3/et.fast.random.csv")




########################### msd





msddata <- fread("msddata.csv")
msdata <- rbind(b_1.msd, b_2.msd, b_4.msd, c_1.msd, c_2.msd, c_4.msd, d_1.msd, d_2.msd, d_4.msd,
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



view(msddata)
####Mean square displacement
p1 <- ggplot(msdatatf, aes( x = dt*15 , y = mean, color = condition, fill = condition)) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line() +
  labs( x = expression( paste(Delta,"t (min)") ),
        y = "mean square displacement") 
#  stat_compare_means(aes(group = condition), label = "p.signif", 
#                     label.y = c(16, 25, 29))

p1 + scale_x_continuous(limits = c(30, 290)) + theme_classic()



p1 <- ggplot( msddata, aes( x = dt*15 , y = mean, color = replicate, fill = replicate)) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line() +
  labs( x = expression( paste(Delta,"t (min)") ),
        y = "mean square displacement") 
#  stat_compare_means(aes(group = condition), label = "p.signif", 
#                     label.y = c(16, 25, 29))

p1 + scale_x_continuous(limits = c(30, 290)) + theme_classic()

#Persistence

b_t.acov <- aggregate( b_fast.tracks.random, overallDot, FUN = "mean.se" )
b_t.acov$condition <- "SC"
b_t.acov$dt <- b_t.acov$i * timeStep( b_t ) * 10
b_t.acov[,2:4] <- b_t.acov[,2:4] / b_t.acov$mean[1]

c_t.acov <- aggregate( c_fast.tracks.random, overallDot, FUN = "mean.se" )
c_t.acov$condition <- "SC+Wnt"
c_t.acov$dt <- c_t.acov$i * timeStep( c_t ) * 10
c_t.acov[,2:4] <- c_t.acov[,2:4] / c_t.acov$mean[1]

d_t.acov <- aggregate( d_fast.tracks.random, overallDot, FUN = "mean.se" )
d_t.acov$condition <- "SC+PC"
d_t.acov$dt <- d_t.acov$i * timeStep( d_t ) * 10
d_t.acov[,2:4] <- d_t.acov[,2:4] / d_t.acov$mean[1]

e_t.acov <- aggregate( e_fast.tracks.random, overallDot, FUN = "mean.se" )
e_t.acov$condition <- "SC+PC+i"
e_t.acov$dt <- e_t.acov$i * timeStep( e_t ) * 10
e_t.acov[,2:4] <- e_t.acov[,2:4] / e_t.acov$mean[1]


acovdata <- rbind( b_t.acov,c_t.acov,d_t.acov,e_t.acov)


## plot dir without folders

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


### plot MSD without folders


#MSD









b_msd <- data.frame(msd = c(8.475169357,19.2841133, 29.31538364, 40.88845524, 51.81546821, 61.65320753, 71.43896987, 80.35781459, 88.28840664,
                            95.66997578, 99.98160743, 100.3347125, 108.2288975, 110.3929221, 107.1870041, 97.73607458, 97.94967399, 70.95343922,79.93133386), 
                    SEM = c(3.903062696, 9.994340678, 15.18781093, 21.85781445, 28.32026294, 33.89550416, 39.28289754, 43.26508639, 46.27301281, 49.15042418, 48.9829859,
                            45.85448079, 49.85392837, 47.96250779, 42.30017044, 34.88323407, 36.81480378, 44.25571059, 51.73180995),
                    time = c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                             195, 210, 225, 240, 255, 270), 
                    condition = c("SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC", "SC",
                                  "SC", "SC", "SC", "SC", "SC", "SC", "SC"))

c_msd <- data.frame(msd = c(12.3193082, 26.40357665, 42.05931289, 59.80885697, 78.34049791,  96.74330326, 115.2226528, 134.6304687, 155.4552044, 177.9378179,
                            203.2753127, 228.1502333, 254.6238313, 281.9283697, 304.4132198, 318.5439594, 328.4465862, 354.9188846, 384.7323406),
                    SEM = c(0, 0.021347478, 0.022899443, 0.021893858, 0.020800602, 0.020637294, 0.02039385, 0.02038432,
                            0.01983696, 0.019485707, 0.019665489, 0.019607324, 0.019569159, 0.019337261, 0.019264371,
                            0.019436864, 0.019712482, 0.01965481, 0.019165852),
                    time= c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                            195, 210, 225, 240, 255, 270), 
                    condition = c("SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt",
                                  "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt", "SC+Wnt"))

d_msd <- data.frame(msd = c(1, 0.689662153, 0.531030501, 0.485461547, 0.444981635, 0.422936513, 0.407079627,
                            0.393273055, 0.366400773, 0.342393511, 0.322826438, 0.328674614, 0.314934903,
                            0.304794333, 0.301269991, 0.280817085, 0.273382941, 0.277379593, 0.270572162),
                    SEM = c(0, 0.023200034, 0.020664055, 0.019870195, 0.019617406, 0.01908628, 0.018402505, 0.018546908,
                            0.018822209, 0.018450514, 0.019240077, 0.01870074, 0.019048253, 0.019257085, 0.019383445,
                            0.01846525, 0.018858018, 0.018719036, 0.018591558),
                    time= c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                            195, 210, 225, 240, 255, 270), 
                    condition = c("SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", 
                                  "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC", "SC+PC"))

e_msd <- data.frame(msd = c(1, 0.701667452, 0.537200346, 0.478979771, 0.426837702, 0.380889773, 0.37340414,
                            0.355480521, 0.339227028, 0.327367354, 0.314174523, 0.293544139, 0.291157575,
                            0.283808345, 0.269883369, 0.263644404, 0.262993466, 0.255346684, 0.25420797),
                    SEM = c(0, 0.021548833, 0.021612329, 0.019647063, 0.018738367, 0.017977935, 0.018052047, 0.01790108,
                            0.016714176, 0.015967777, 0.016481424, 0.016338921, 0.016088542, 0.016657846, 0.016784674,
                            0.017088774, 0.01685693, 0.017139392, 0.017512149), 
                    time= c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180,
                            195, 210, 225, 240, 255, 270), 
                    condition = c("SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i",
                                  "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i", "SC+PC+i"))

rmsd <- rbind(b_msd, c_msd, d_msd, e_msd)


p1 <- ggplot( rmsd, aes( x = time , y = msd, color = condition, fill = condition ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = msd - SEM, ymax = msd + SEM) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = "time(min)",
        y = "msd" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )


