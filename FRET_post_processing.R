#INFO------------------------------------------------------------------------------
# Name: F-Lab FRET Post processing Tool
# Desc: Tool to calculate FRET Efficiencies and create standardized graphs to increase reproducability and comparability
#
# In:   Results file from ImageJ (###format?), nd-file from Microscope
# Out:  Tables and plots after the following operations:
#       - backgrund correction
#       - baseline correction
#       - ratio calculation
#       - leakage correction (using spectral decomposition, quantum yields and stochimetrics)
#       User can: 
#       - group ROIs (creates averages and standard deviations)
#       - select data to be shown in report
#       - export table (csv,tsv) and plot (svg,png)
#
# '###' marks code for attention (possible sources of error)
# written in R version 3.2.5 (2016-04-14)

#LIBRARIES-------------------------------------------------------------------------
#library(ggplot2) ###Version?
#library(shiny) ###Version? 
#library(baseline) http://crantastic.org/packages/baseline

#FUNCTIONS-------------------------------------------------------------------------

read.tcsv = function(file, header, sep, ncol, nrows) {
  #desc: read function for horizontal csv, transposes the text before calling read.csv
  #in:   like read.csv, csv-file in horizontal layout, ncol & nrows refer to the input
  #out:  dataframe
  # source: http://stackoverflow.com/questions/17288197/reading-a-csv-file-organized-horizontally
  
  x = readLines(file, n=nrows)
  
  .splitvar = function(x, sep, ncol) {
    var = unlist(strsplit(x, split=sep))
    length(var) = ncol
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=ncol))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header)
  return(out)
} #read transposed csv

#UI----------
#shiny app
  #setup tab - 
    #raw data source (from ImageJ) [else: call ImageJ from R?? Bio7]
    #settings; import from nd_file, editable settings, option to save as nd_file
  #show graph & tables (allow export)  constantly updated depending on active corrections :)
    #bg-correction: checkbox: activates select background ROI
    #baseline-correction: checkbox: activates select baseline area tool (brush) & dropdown basleine fitting type
    #donor-acceptor bleed-through-correction: checkbox: activates textbox: input factor 
    #channel select & grouping
    #events toggle show/hide / add

#DATA INPUT------------------------------------------------------------------------
#nd file (experiment parameters and event log)
### add prompt message & set directory to last used ###if no nd file exists: ignore / create via UI
file_nd <- file.choose()
nd <- readLines(file_nd)

#intensity files  
### add prompt message, ### file selector often apears in background
intensities_d  <- read.table(file.choose(), header=TRUE, sep = "\t") ### determine nrows if file size causes performance inssues 

### add other channels, ### how to deal with multiple donors/acceptors? (in 4 channel FRET): use directory that contains all data?

#additional info from user: ###maybe add to nd-file?
bl_startFr <- 2#readline("baseline start frame: ") #default: 2 or select values from diagram brush (multiple sections possible??)
bl_endFr <- 30#readline("baseline end frame: ") #default?
bg_col <- 2#readline("Which column countains the background values? ") #only numbers or start/end (use function?)
loop1_delay <- 1#readline("loop1 delay in frames: ")
loop2_delay <- 2#readline("loop2 delay in frames: ")
clockfactor <- 1.5#readline("clockfactor: ")

#DATA WRANGLING-------------------------------------------------------------------

#.ND FILE---
#This file is composed of two sections: parameters (vertical key-value pairs, csv) and events (table without headers, csv)

#parsing 
NParameters <- grep("Event1", nd, value = FALSE) - 1 #number of parameters? (this line separates both sections) 
parameters <- read.tcsv(file_nd, header = TRUE, sep = ",", ncol= 2 , nrows = NParameters) #custom function to transpose
events <- read.csv(file_nd, header = FALSE, sep=",", skip = NParameters, nrows = parameters$NEvents)
  
#add headers for events-dataframe
names(events) <- c("ID", "description", "NA_NA" , "time", "NA_num1", "frame", "NA_num2", "NA_log", "color_num") ### add missing headers

#correct data classes
# parameters
parameters$StartTime1 <- strptime(parameters$StartTime1,"%Y%m%d %H:%M:%OS") #use op <- options(digits.secs = 3) to show ms, options(op) to reset
p_logicals <- grep('^Do', names(parameters), value = TRUE)  #all headers starting with "Do" (this seems to indicate all logical parameters)
#p_logicals <- add WaveInFileName ###
parameters[p_logicals] <- lapply(parameters[p_logicals], function(y) as.logical(gsub(" ", "", y))) #remove whitespaces and convert to logical 

# events
events$time <- strptime(events$time,"%Y%m%d %H:%M:%OS") #use op <- options(digits.secs = 3) to show ms, options(op) to reset
events$ID <- as.character(events$ID)
events$description <- as.character(events$description)
events$NA_log <- lapply(events$NA_log, function(y) as.logical(gsub(" ", "", y)))

#INTENSITY FILES---

##RENAME COLUMN HEADERS---
names(intensities_d)[1] <- "frame"
nROIs <- ncol(intensities_d) - 1
names(intensities_d)[2:(nROIs+1)] <- paste("ROI", 1:(nROIs), sep="")
names(intensities_d)[bg_col] <- "ROI_bg"

#CALCULATIONS---------------------------------------------------------------------

#ADD TIME COLUMN---
intensities_d$time_s <- intensities_d$frame * clockfactor ##convert frames to time in sec ### use time class for min:s.ms?

#BACKGROUND CORRECTION---
ROIs <- grep("ROI", names(intensities_d), value = TRUE) #list of all ROIs
intensities_d[ROIs]  <- intensities_d[ROIs] - intensities_d$ROI_bg #subtract background (ROI_bg) from all ROIs

# #BASELINE CORRECTION---
#if true
#plot graph of all rois (wo bg) of intensities and use brush to mark baseline areas (click checkbox to add, else replace), allow zoom? brush in plotly?
#intensities_d$ROIavg <- avg(intensities_d[ROIs, -ROI_bg])

#View(donorInt)
# # show graphs (overlapping or average) in new window and allow UI entry of the following parameters: 


# # input validation: correct datatype & within range
# 

# #---
# 
# #DRUG ADDITION PROTOCOL
# # array of name, conc, start and stop used to make marks on graphs
# # include loop delay
# # use numbers for shadow in graph
# #---
#  
# #MERGING SELECTED ROIS INTO GROUPS
# # make selection, add group as new column in table (average, stdv), optional: show stdv shadow in graph
# #--
# 
# #EXPORT DATA---
# # select options: donor, acceptor, ratios
# #                 single channels, groups
# #                 order by ROI, order by value type
# file_exp <- file.choose(new = TRUE) #default: same folder as input files
# #---

#PLOT---
# ggplot2 or plotly (zoomable)
# geom line donor
# geom line acceptor
# geom line ratio
# export svg (Cairo for smoothing?)
# svg(filename="Std_SVG.svg",
#     width=5,
#     height=4,
#     pointsize=12)
# my_sc_plot(data)
# dev.off()
# 
# png(filename="Std_PNG_cairo.png",
#     type="cairo",
#     units="in",
#     width=5,
#     height=4,
#     pointsize=12,
#     res=96)
# my_sc_plot(data)
# dev.off()
#---

#OUTPUT
print(parameters)
print(events)
print(intensities_d)
#write.csv(df, file="somedf.csv")

## UPGRADES
## interactive ui: shiny
## which channels to show, which to merge
## automatic roi selection.
## ratio heatmap per frame (picture, if ImageJ controlable anddoesnt take too long to compute) -> video contol: regular and by graph brush
#ratiometric FRET -> FRET Efficiency!
# at least 2 corrections:
# contribution from direct A excitation to IA
# ratio between the D and A fluorescence quantum yields 
# If IA and ID were not determined using spectral decomposition two additional factors must be accounted for:
# leakage" of D emission into the value of  IA
# leakage of A emission into the value of ID.


# set up FRET system:
#choose donor(s) and acceptor(s), excitation and emission wavelengths, concentration ratios & quantum yields
#spectra: http://www.tsienlab.ucsd.edu/Documents.htm -> xls file -> text file (how to parse synonymous fluorophor names?)
#http://www.fluortools.com/software/ae/documentation/tools/FRET
#this is done via direct measurement not calculation
