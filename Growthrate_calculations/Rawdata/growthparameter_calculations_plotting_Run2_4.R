# Script to calculate growth rates for plate reader experiments

# This script processes raw optical density (OD) data from plate reader experiments, 
# calculates growth rates, yields, and lag times for each well, and outputs results for further analysis.

# PARAMETERS TO CUSTOMIZE:

# Minimum number of data points required for calculating growth rates.
# Ideally, this should cover at least two doubling times of the organism.
# Every datapoint is 30 sec, so 120 points is 1 hour (still not too much)
DatapointsinH <- 240

# Set the working directory for accessing input and output files.
# Uncomment and update the following line to your directory path if needed:
# setwd("/path/to/your/directory")

# Path to the raw data file (comma-separated). 
# The file should have:
# 1) Time in the first column.
# 2) Well names in the first row (e.g., A1, H11).
# 3) OD measurements in subsequent columns.
rawdata_file <- "../Rawdata/SM Run 2_4.csv"

# Path to the metadata file. This file must include:
# 1) A "well" column matching the well names in the raw data file.
# 2) A "type" column with well classifications:
#    - "s": sample wells
#    - "n": negative controls
#    - "p": positive controls
#    - "B": background wells
#    - 0: wells to discard
metadata_file <- "../Rawdata/metadata run 2CSV.csv"

# Create filenames for output
csvname <- sub("^.*/", "", rawdata_file)
basename <- substr(csvname, 1, nchar(csvname) - 4)
# Replace .csv at end with .png for image output
figname <- sub("\\.csv$", ".png", csvname)
# Replace "SM Run " with results_r for results output
resultsname <- sub("SM Run ", "results_r", csvname)

# Folder to save output files.
output_folder <- "../Rawdata/output/"

# Set to FALSE if the time in the raw data file is in hours.
time_seconds <- TRUE

# Set to FALSE if the input files are semicolon-separated instead of comma-separated.
comma_seperated <- TRUE

# LIBRARY INSTALLATION (if not installed already):
# Uncomment and run the following lines to install required packages:
# install.packages("tidyverse")    # Includes ggplot2, dplyr, tidyr, and more
# install.packages("growthrates") # For growth rate fitting
# install.packages("dplyr")

# Load necessary libraries
library(tidyverse)
library(growthrates)
library(dplyr)

# STEP 1: Load raw data and metadata
if (comma_seperated) {
  rawdata <- read.csv(rawdata_file, header = TRUE, check.names = FALSE)
  metadata <- read.csv(metadata_file, header = TRUE)
} else {
  rawdata <- read.csv2(rawdata_file, header = TRUE, check.names = FALSE)
  metadata <- read.csv2(metadata_file, header = TRUE)
}

# STEP 2: Preprocess raw data
# Remove rows and columns with only NAs or empty values
rawdata <- rawdata[!apply(is.na(rawdata) | rawdata == "", 1, all),]

# Rename the first column to "Time" for clarity
names(rawdata)[1] <- "Time"

# Convert all columns to numeric for calculations
rawdata <- mutate_all(rawdata, function(x) as.numeric(as.character(x)))

# Convert time to hours if needed
if (time_seconds) {
  rawdata$Time <- (rawdata$Time - rawdata$Time[1]) / 3600
}

#create long format for plotting
plot_data <- pivot_longer(rawdata, !Time, names_to = "well", values_to = "OD")
plot_data$plate_col<- as.integer(substring(plot_data$well, 2))
plot_data$plate_row <- substr(plot_data$well,0,1)

#adjust the row and col with the bacteria and medium from the metadata file
#first split the well in metadata in col and row
metadata$plate_col <- as.integer(substring(metadata$well, 2))
metadata$plate_row <- substr(metadata$well,0,1)

#Link the row to the medium (correct except for row 12, so add / <bacteria in row 12>)
#Link the column to bacteria (except for row 12, there add MM)
metadata <- metadata %>%
  group_by(plate_row) %>%
  mutate(
    mm_medium = medium[medium != "MM"][1],              # first non-MM medium
    mm_bacteria = bacteria[medium == "MM"][1],          # bacteria from MM row
    corrected_medium = ifelse(
      medium == "MM",
      paste(mm_medium, ifelse(bacteria == "None", "", bacteria), sep = "\n"),
      paste(medium, ifelse(mm_bacteria == "None", "", mm_bacteria), sep = "\n")
    )
  ) %>%
  ungroup() %>%
  select(-mm_medium, -mm_bacteria)  # optional: clean up helper columns

plot_data_adjusted <- merge(plot_data, metadata, by=c('plate_col', 'plate_row'))
plot_data_adjusted$row_medium <- paste(plot_data_adjusted$plate_row, plot_data_adjusted$corrected_medium, sep = " ")
plot_data_adjusted$col_bacteria <- ifelse(
  plot_data_adjusted$plate_col == 12,
  paste(plot_data_adjusted$plate_col+10, "MM"),
  paste(plot_data_adjusted$plate_col+10, plot_data_adjusted$bacteria)
)

#code for plotting in well
#split data column and rows of 96well plate
ggplot(plot_data_adjusted, aes(x=Time, y=OD, )) +
  facet_grid(row_medium ~ col_bacteria) +
  geom_line() +
  theme_bw() +
  labs(title=paste('Timecourse of plate', basename)) +
  theme_bw() +
  theme(axis.text.x = element_blank())
#Export the graph
ggsave(paste(output_folder, figname, sep = ""))

# STEP 3: Subtract background OD
background_wells <- metadata[metadata$type == "B",]$well
if (basename == "SM Run 2_1" | basename == "SM Run 2_2" | basename == "SM Run 2_3"| basename == "SM Run 2_4") {
  background_wells <- background_wells[!grepl("12$", background_wells)]
}
average_background <- rowMeans(rawdata[, background_wells], na.rm = TRUE)

# Plot the background OD over time for visualization
ggplot(data.frame(bg = average_background, time = rawdata$Time), aes(x = time, y = bg)) +
  geom_line() + 
  xlab("Time (hours)") +
  ggtitle(paste("Average Background OD Over Time of plate ", basename))
ggsave(paste(output_folder, "background ", figname, sep = ""))

# Subtract average background OD from all wells
rawdata_bg <- cbind(rawdata[1], rawdata[,-1] - average_background)
cat("Lowest value after subtracting background: ", min(rawdata_bg[,-1]))

# Replace values <= 0 with a small positive constant
rawdata_bg[,-1] <- apply(rawdata_bg[,-1], 2, function(x) ifelse(x <= 0, 1e-9, x))

# STEP 4: Create a long-format dataset for plotting
plot_data <- pivot_longer(rawdata_bg, !Time, names_to = "well", values_to = "OD")
plot_data$col <- as.integer(substring(plot_data$well, 2))
plot_data$row <- substr(plot_data$well, 1, 1)

# STEP 5: Filter metadata to include only samples, positive, and negative controls
metadata <- subset(metadata, type == "s" | type == "p" | type == "n")
growth_data <- rawdata[c("Time", metadata$well)]
rawdata_bg <- rawdata_bg[c("Time", metadata$well)]

# STEP 5.5: Filter out the columns in which the last value is not at least 0.05 OD points higher 
        # than the first value (ergo, the bacteria has not grown) - value can be changed if large background
        # If values are needed for downstream calculations this can be commented out.
# diffs <- growth_data[nrow(growth_data), ] - growth_data[1, ]
# growth_data_filtered <- growth_data[, diffs >= 0.05]
# 
# #Test what has been removed by filtering
# filtered_data <- growth_data[, c(1, which(diffs < 0.05))]
# filtered_data_plot <- filtered_data %>%
#   pivot_longer(
#     cols = -1,                         # All columns except the first (x-axis)
#     names_to = "sample",              # New column with sample names
#     values_to = "value"               # New column with measurement values
#   )
# # Plot with ggplot2
# ggplot(filtered_data_plot, aes(x =Time, y = value, color = sample)) +
#   geom_line() +
#   labs(x = names(filtered_data_plot)[1], y = "Value") +
#   ggtitle(paste("Removed rimecourses of plate ", basename)) +
#   theme_minimal()
# ggsave(paste(output_folder, "removed ", figname, sep = ""))
growth_data_filtered <- growth_data

# STEP 6: Calculate growth properties
# Calculate yield as the maximum OD for each well
yield <- apply(select(rawdata_bg, -Time), 2, max, na.rm = TRUE)

# Prepare growth data by capping OD values to avoid artifacts (like double growthcurves)
# Edit this value if you still see the "double growthcurve" effect- consider that in our 
# experiment we are only interested in the first part of the growthcurve since this is the 
# initial growth that the bacteria will also have when grown together with another strain
# 
# Note: this is now left out
#
# growth_data[-1][growth_data[-1] > 0.6] <- 0.6

# Fit growth rates for each well using the `growthrates` packag
##ALSO PLOT IMAGES TO CHECK IF THE GROWTHRATE ESTIMATE IS CORRECT
all_fits <- list()
for (i in colnames(growth_data_filtered)[-1]) {
  fit <- fit_easylinear(growth_data_filtered$Time, growth_data_filtered[[i]], h = DatapointsinH)
  # Display the plot in the R window
  plot(fit, log = "y", main = paste("Fit for plate ", basename, " Well ", i, " gr= ",round(coef(fit)[3], 3)," y= ", round(yield[i], 3) , sep = ""))
  # Add a horizontal dashed line at yield[i]
  abline(h = yield[i]+mean(average_background), lty = "dashed", col = "red")
  # Save the plot as a JPEG file
  jpeg(file = paste(output_folder, basename, " fit_well_", i, ".jpeg", sep = ""))
  plot(fit, log = "y", main = paste("Fit for plate ", basename, " Well ", i, " gr= ",round(coef(fit)[3], 3)," y= ", round(yield[i], 3) , sep = ""))
  # Add a horizontal dashed line at yield[i]
  abline(h = yield[i]+mean(average_background), lty = "dashed", col = "red")
  dev.off()
  all_fits <- c(all_fits, fit)
}

# Collect growth rate and lag time estimates
growthrates <- unlist(lapply(all_fits, function(x) coef(x)[3][[1]]))
lags <- unlist(lapply(all_fits, function(x) coef(x)[4][[1]]))

# STEP 7: Combine results and export
results <- metadata
results$growthrate <- growthrates
results$yield <- yield
results$lag <- lags
results$medium <- substr(results$medium, 1, nchar(results$medium) - 3)

# Save results to a CSV file
write_csv2(results, paste(output_folder, resultsname, sep = ""))

