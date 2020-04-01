############################
### Install the packages ### (need to be done everytime R version is updated)
############################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#     BiocManager::install()

#  install.packages("tidyverse")
#  install.packages("splancs")
#  install.packages("caroline")
#  install.packages("viridis")
#  BiocManager::install.packages("flowCore")


### Load packages
library(tidyverse)
library(splancs)
library(caroline)
library(viridis)
library(flowCore)
library(grid)
library(gridExtra)

# Population-specific colors
group_colors <- c(unknown="grey", beads="red3", prochloro=viridis::viridis(4)[1],synecho=viridis::viridis(4)[2],picoeuk=viridis::viridis(4)[3], croco=viridis::viridis(4)[4])

# functions for plotting cytograms

plot_cytogram <- function (evtopp, para.x = "X.FSC.2", para.y = "X.692.40", ...){
    cols <- colorRampPalette(viridis(7))
    par(pty = "s")
    plot(evtopp[, c(para.x, para.y)], pch = 16, cex = 0.3, col = densCols(log10(evtopp[, c(para.x, para.y)]),
          colramp = cols))
} # end plot.cytogram

plot_vct_cytogram <- function(opp, para.x = 'X.FSC.2', para.y = 'X.692.40', transform = T, ...){
        p <- ggplot2::ggplot(opp) + ggplot2::stat_bin_2d(ggplot2::aes_string(para.x,
                para.y, fill = "pop", alpha = quote(..count..)), colour = NA,
                bins = 100, show.legend = T) + ggplot2::theme_bw() +
                ggplot2::coord_fixed() + ggplot2::stat_density_2d(ggplot2::aes_string(para.x,
                para.y, color = "pop"), bins = 5, show.legend = F) +
                ggplot2::scale_fill_manual(values = group_colors) + ggplot2::scale_alpha_continuous(range = c(0.3,1)) +
                ggplot2::scale_color_manual(values = group_colors) +
                ggplot2::guides(color = "none", alpha = "none", fill = ggplot2::guide_legend(override.aes = list(size = 2,
                    alpha = 0.5), title = "population")) +
                theme(aspect.ratio = 1)
            if (transform)
                p <- p + ggplot2::scale_y_continuous(trans = "log10") +
                    ggplot2::scale_x_continuous(trans = "log10")
            return(p)
}

concatenateFCS <- function(fcs_list, n=100000, min.fsc=0, min.pe=0, min.chl=0, transform=TRUE,...){
    n <- as.numeric(n)
    DF <- NULL
    i <- 0
    for (file in fcs_list){
        message(round(100*i/length(fcs_list)), "% completed \r", appendLF=FALSE)
        tryCatch({
            df <- read.FCS(file, transformation = F, emptyValue=F) # create a flowFrame WITHOUT transformation
            df <- caroline::tab2df(exprs(df))
            df <- subset(df, X.FSC.2 > min.fsc & X.580.30 > min.pe & X.692.40 > min.chl)
            df <- df[round(seq(1,nrow(df), length.out=round(n/length(fcs_list)))),]

            if(any(is.na(df))) next
                DF <- rbind(DF, df)
        }, error = function(e) {
            cat(paste("Error with file", file, ":", e))
        })
        i <- i + 1
        flush.console()
    }
    return(DF)
} # end of concatenateFCS

#################
## SOM BUCKETS ##
#################

experiment <- 'NH_incubations'
path <- '~/Documents/Influx/Gradients_3/NH_incubations/'    # Set path according to your directory names
fcs_path <- paste0(path, 'FCS/')

bucket_file <- paste0(path, experiment, '_buckets_12.csv')  # Result of running Matthias' self-organizing map script
bucket <- read.table(bucket_file, sep = '\t', header = T)
bmu <- sort(unique(bucket$bmu))

# Group files by SOM bucket
bucket_list <- list()
for (type in bmu){
    print(paste0('bucket = ', as.character(type)))
    these <- subset(bucket, bmu == type)    # Grab part of data frame for this bucket
    opp_list <- paste0(fcs_path, these$file)
    print(paste0(as.character(length(opp_list)), ' files'))
    bucket_list[[type]] <- opp_list
}

summary_table <- NULL   # initialize the final summary table

#################################
### LOOP AND GATE EACH BUCKET ###
#################################


#for (i in seq(1, length(bucket_list))) {    # step through each SOM bucket
    i <- 1
    file_list <- bucket_list[[i]]
    print(paste("processing bucket:", i))

    # Concatenate files in this bucket to gate
    OPP <- concatenateFCS(file_list)
    OPP$pop <- "unknown"
    OPP <- OPP[, c("X.FSC.2", "X.SSC", "X.530.40", "X.580.30", "X.610.20", "X.692.40", "pop")]

    ##################
    ### GATE BEADS ###
    ##################

    param <- c("X.692.40", "X.580.30")

    # draw gates
    print("Gating Beads")
    plot_cytogram(OPP, param[1], param[2], main="BEADS")
    poly_beads <- getpoly(quiet=TRUE)
    polygon(poly_beads, lwd = 2, border = group_colors[2])

    # filter particles
    beads <- subset(OPP, inout(OPP[, param], poly = poly_beads, bound = TRUE, quiet = TRUE))
    OPP[row.names(beads), 'pop'] <- "beads"
    #pop <- OPP$pop  # Pull this out so we can add it back in later.  Concatenated dataframes are tricky.
    #beads <- beads[, c("X.FSC.2", "X.SSC", "X.530.40", "X.580.30", "X.610.20", "X.692.40")]
    #norm <- beads %>%
    #  summarise_all(median)   # Find median for all bead parameters

    #norm_opp <- data.frame(mapply("/", OPP, norm))  # normalize to bead median
    #norm_opp$pop <- pop     # put the labels back

############################
### GATING SYNECHOCOCCUS ###
############################

# Select parameters for gating
param <- c("X.FSC.2", "X.580.30")

# draw gates
x <- subset(OPP, pop=="unknown")

print("Gating Synecho")
par(mfrow = c(1, 1))
    plot_cytogram(x, param[1], param[2], main="SYN")
    points(subset(OPP, pop == "beads")[, param], col = group_colors[2], pch = 16, cex = 0.4)  # beads for reference

poly_syn <- getpoly(quiet = TRUE)
polygon(poly_syn, lwd = 2,  border = group_colors[4])

# label particles
syn <- subset(x, inout(x[, param], poly = poly_syn, bound = TRUE, quiet = TRUE))
OPP[row.names(syn), 'pop'] <- "synecho"

############################
### GATE PROCHLOROCOCCUS ###
############################

# Note that for some stations, there may be no Pro present

# Select parameters for gating
param <- c("X.FSC.2", "X.692.40")

# draw gates
x <- subset(OPP, pop=="unknown")

print("gating PRO")
par(mfrow = c(1, 1))
    plot_cytogram(x, param[1], param[2], main="PRO then PicoEuk"); abline(v=1, h=1, col='grey',lty=2)
    points(subset(OPP, pop == "beads")[, param], col = group_colors[2], pch = 16, cex = 0.4)  # beads for reference
    points(syn[,param], col = group_colors[4], pch = 16, cex = 0.4)   # syn for reference

poly_pro <- getpoly(quiet = TRUE)
polygon(poly_pro, lwd = 2,  border = group_colors[3])

# label particles
pro <- subset(x, inout(x[, param], poly = poly_pro, bound = TRUE, quiet = TRUE))
OPP[row.names(pro), 'pop'] <- "prochloro"

###########################
### GATE PICOEUKARYOTES ###
###########################
#  pico: everything left that is bigger and has more chlorophyll fluo than Pro and Syn

# Select parameters for gating
param <- c("X.FSC.2", "X.692.40")

# draw gates
x <- subset(OPP, pop=="unknown")

print("gating PICOEUKS")
par(mfrow = c(1, 1))
    plot_cytogram(x, param[1], param[2], main="PRO then PicoEuk"); abline(v=1, h=1, col='grey',lty=2)
    points(subset(OPP, pop == "beads")[, param], col = group_colors[2], pch = 16, cex = 0.4)  # beads for reference
    points(syn[,param], col = group_colors[4], pch = 16, cex = 0.4)   # syn for reference
    points(pro[,param], col = group_colors[3], pch = 16, cex = 0.4)   # pro for reference

poly_pico <- getpoly(quiet=TRUE)
polygon(poly_pico, lwd = 2,  border=group_colors[5])

# label particles

pico <- subset(x, inout(x[, param], poly = poly_pico, bound = TRUE, quiet = TRUE))
    ### PLOT TO CHECK GATES ###

    p <- list()
        p[[1]] <- plot_vct_cytogram(OPP, "X.FSC.2","X.692.40", transform = FALSE)
        p[[2]] <- plot_vct_cytogram(OPP, "X.FSC.2","X.580.30", transform = FALSE)
        p[[3]] <- plot_vct_cytogram(OPP, "X.692.40","X.580.30", transform = FALSE)
        p[[4]] <- plot_vct_cytogram(OPP, "X.SSC","X.692.40", transform = FALSE)
        tx <- textGrob(paste0('Bucket ', as.character(i)))
    grid.arrange(grobs = p, top = tx, nrow = 2)


    # save gating parameters for this bucket

    gfile <- paste0(path, experiment, '_gates_Bucket_', i, '.RData')
    save(list = c('poly_beads', 'poly_syn', 'poly_pro', 'poly_pico'), file = gfile)

    remove("OPP")

    ########################################
    ### LOOP THROUGH FILES & APPLY GATES ###
    ########################################

    for (file in file_list){    # step through files in current bucket
    #file <- file_list[3]

    fcs <- read.FCS(file, transformation = T, emptyValue = F) # create a flowFrame
    opp <- tab2df(exprs(fcs)) # turn flowFrame into a table

    plot_cytogram(opp)

    ########################
    ### Apply bead gates ###
    ########################

    beads <- subset(opp, inout(opp[, param], poly = poly_beads, bound = TRUE, quiet = TRUE))

    # label particles
    opp[row.names(beads), 'pop'] <- "beads"

    #####################
    ### PHYTOPLANKTON ###
    #####################

    opp$pop <- "unknown"

    ### SYNECHOCOCCUS ###

    # Select parameters
    param <- c("X.FSC.2", "X.580.30")

    x <- subset(opp, pop=="unknown")
    syn <- subset(x, inout(x[, param], poly = poly_syn, bound = TRUE, quiet = TRUE))
    opp[row.names(syn), 'pop'] <- "synecho"

    ### PRO ###
    # Select parameters

    param <- c("X.FSC.2", "X.692.40")
    x <- subset(opp, pop == "unknown")
    pro <- subset(x, inout(x[, param], poly = poly_pro, bound = TRUE, quiet = TRUE))
    opp[row.names(pro), 'pop'] <- "prochloro"

    ### PicoEUK ###
    x <- subset(opp, pop == "unknown")
    pico <- subset(x, inout(x[, param], poly = poly_pico, bound = TRUE, quiet = TRUE))
    opp[row.names(pico), 'pop'] <- "picoeuk"

    #################
    ### SAVE PLOT ###
    #################

    fig_file <- paste0(path, "Cytograms/", file, ".png")
    graphics.off()  # close any quartz windows so png will save a figure

    png(fig_file, width=9, height=9, unit='in', res=300)
    p <- list()
        p[[1]] <- plot_vct_cytogram(opp, "X.FSC.2","X.692.40", transform = FALSE)
        p[[2]] <- plot_vct_cytogram(opp, "X.FSC.2","X.580.30", transform = FALSE)
        p[[3]] <- plot_vct_cytogram(opp, "X.692.40","X.580.30", transform = FALSE)
        p[[4]] <- plot_vct_cytogram(opp, "X.SSC","X.692.40", transform = FALSE)
    tx <- textGrob(basename(file))
    grid.arrange(grobs = p, top = tx, nrow = 2)

    dev.off()

    #####################
    ### SUMMARY STATS ###
    #####################

    stat.table <- NULL
    for(i in unique(opp$pop)){
        print(i)
        if(i == 'unknown' | i == 'beads') next
        p <- subset(opp, pop == i)
        n <- nrow(p)
        if(n ==0) {
            fsc <- 0
            chl <- 0
            }else{
                fsc <- round(median(p$X.FSC.2),4)
                chl <- round(median(p$X.692.40),4)
                ssc <- round(median(p$X.SSC),4)
                pe <- round(median(p$X.580.30),4)
            }

        var <- cbind(i,n,fsc,chl,ssc,pe)
        stat.table <- rbind(stat.table, var)

    }   #end stat loop


    remove("opp")  # clear variable for the next run


    table <- data.frame(cbind(stat.table, file=basename(file)))
    summary.table <- rbind(summary.table, table)
    } # end file_list loop

#}   # end of bucket loop

file_summary <- paste0(path, experiment, '_summary.csv')
write.csv(summary.table, file = file_summary, row.names=FALSE)


# ####################
# ### ADD METADATA ###
# ####################
# meta <- read.csv("meta.txt")
# meta$file <- paste0(meta$filename, ".fcs")

# summary.table <- read.csv(file=paste(savepath,"/summary.csv", sep=""))
# merged <- merge(summary.table, meta, by="file", all.x=TRUE) # add metadata to summary table
# merged$abundance <- round(1000*merged$n / merged$volume_ul) #calculate cell abundance




# ################
# ### PLOTTING ###
# ################
# library(tidyverse)
# library(viridis)
# merged %>%
#     group_by(station, depth, i) %>%
#     summarise(abundance=mean(abundance)) %>%
#     ggplot() +
#     geom_point(aes(station, depth, colour=abundance), pch=16, size=4) +
#     scale_color_viridis() +
#     scale_y_reverse() +
#     theme_bw() +
#     facet_grid(~ i , scales='free')


# library(MBA)
# df <- subset(merged, i == 'prochloro')
# df <- subset(merged, i == 'synecho')
# df <- subset(merged, i == 'picoeuk')

# mba <- mba.surf(df[,c('station', 'depth', 'abundance')], 50, 50)
# dimnames(mba$xyz.est$z) <- list(mba$xyz.est$x, mba$xyz.est$y)
# df2 <- reshape::melt(mba$xyz.est$z, varnames = c('station', 'depth'), value.name = 'abundance')

# ggplot(data=df2, aes(station, depth))+
#        geom_raster(aes(fill = value), interpolate = F, hjust = 0.5, vjust = 0.5) +
#        geom_point(data = df, aes(station, depth), colour = 'white') +
#        scale_y_reverse() +
#        scale_fill_gradientn(colours = viridis(7)) +
#        theme_bw()
