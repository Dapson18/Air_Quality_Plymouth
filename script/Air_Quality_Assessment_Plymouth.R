
library(xts)
library(lubridate)
library(tidyverse)
library(dplyr)
library(zoo)
library(scales)
library(gridExtra)
library(gridtext)
library(grid)
library(showtext)
library(tidyr)


# Load the Aeroqual sensor data for three transects
data <- read.csv(file.choose())
data1 <- read.csv(file.choose())
data2 <- read.csv(file.choose())

#Select only necessary columns
data <- data[, c("Date.Time", "NO2.ppm.")]
data1 <- data1[, c("Date.Time", "NO2.ppm.")]
data2 <- data2[, c("Date.Time", "NO2.ppm.")]


# Function to convert and clean date columns to required format as some are in character format
convert_dates <- function(data, date_col, date_format) {
  data[[date_col]] <- strptime(data[[date_col]], format = date_format)
  data[[date_col]] <- ymd_hms(data[[date_col]])
  return(data)
}

# Convert Date.Time to proper format
data <- convert_dates(data, "Date.Time", "%m/%d/%Y %H:%M")
data1 <- convert_dates(data1, "Date.Time", "%d %b %Y %H:%M")
data2 <- convert_dates(data2, "Date.Time", "%d %b %Y %H:%M")

# Create a factor to segregate the data by time of day
data$time_of_day <- 'Morning'
data1$time_of_day <- 'Afternoon'
data2$time_of_day <- 'Evening'

#########################################################################################

#convert NO2 for all three transects from ppm to ugm-3
data$NO2_ugm3 <- data$NO2.ppm. * 1 / 0.001
data1$NO2_ugm3 <- data1$NO2.ppm. * 1 / 0.001
data2$NO2_ugm3 <- data2$NO2.ppm. * 1 / 0.001

# Filter out negative or zero NO2 values
data <- filter(data, NO2_ugm3 > 0)
data1 <- filter(data1, NO2_ugm3 > 0)
data2 <- filter(data2, NO2_ugm3 > 0)

######################################################################################

# Function to perform background correction using 5th percentile
correct_background <- function(data, column) {
  background <- quantile(data[[column]], 0.05)
  # Remove '_ugm3' if it's in the column name and then add '_corrected'
  new_column_name <- gsub("_ugm3", "", column)
  data[[paste0(new_column_name, "_corrected")]] <- data[[column]] - background
  return(data)
}

# Apply background correction for each transect
data <- correct_background(data, "NO2_ugm3")
data1 <- correct_background(data1, "NO2_ugm3")
data2 <- correct_background(data2, "NO2_ugm3")

# Time frames for each transect
timeframes <- list(
  morning = list(start = '11:21:00', end = '12:20:00'),
  afternoon = list(start = '13:50:00', end = '14:57:00'),
  evening = list(start = '16:35:00', end = '17:38:00')
)

# Filter data based on the selected timeframes
filter_data <- function(data, start_time, end_time) {
  data <- filter(data, Date.Time >= as.POSIXct(paste('2024-02-20', start_time, sep = " ")),
                 Date.Time <= as.POSIXct(paste('2024-02-20', end_time, sep = " ")))
  return(data)
}

data <- filter_data(data, timeframes$morning$start, timeframes$morning$end)
data1 <- filter_data(data1, timeframes$afternoon$start, timeframes$afternoon$end)
data2 <- filter_data(data2, timeframes$evening$start, timeframes$evening$end)

############################################################################

# Combine all data into a single data frame
df <- bind_rows(data, data1, data2)

# Plot NO2 levels for different times of day
ggplot(df, aes(x = time_of_day, y = NO2_corrected, colour = time_of_day)) +
  geom_boxplot(varwidth = TRUE) +
  labs(title = "Boxplot of Corrected NO2 Levels at Different Times of Day", 
       x = 'Time of Day', y = 'NO2 Corrected (µg/m³)') +
  stat_summary(fun = mean, colour = "darkblue", geom = "point", shape = 18, size = 3, show.legend = FALSE) +
  stat_summary(fun = mean, colour = "darkblue", geom = "text", show.legend = FALSE, vjust = -0.7, aes(label = round(..y.., digits = 3))) +
  theme(legend.position = "none")

#############################################################################

# Handling environmental data (Temperature, Relative Humidity, Wind Speed, etc.)
kes <- read.csv(file.choose(), skip = 9)
kes <- kes[-1, ]  # Remove empty first row
kes <- kes[, -17]  # Remove column with NA values
kes$Time <- as.POSIXct(kes$Time, format = "%Y-%m-%d %H:%M")

# Convert all columns except Time to numeric
kes <- kes %>% mutate_at(vars(-Time), as.numeric)

# Resample data to 1-minute intervals
kes22 <- kes %>%
  group_by(Date.Time = floor_date(Time, "1 min")) %>%
  summarise(across(c(Temp, Rel..Hum., Wind.Speed, Headwind), mean))

# Filter environmental data for the first transect
kes22 <- filter(kes22, Date.Time >= as.POSIXct(paste('2024-02-20', '11:10:00', sep = " ")),
                Date.Time <= as.POSIXct(paste('2024-02-20', '12:20:00', sep = " ")))

# Background temperature data

T_initial_O <- 11.0533  # initial observed temperature
T_final_O <- 13.98    # final observed temperature
t_total <- 70      # total time of your measurements

# Function to apply background correction
correct_temp <- function(measured_temp, elapsed_time, T_initial_O, T_final_O, t_total) {
  T_corrected <- measured_temp - ((T_final_O - T_initial_O) * (elapsed_time / t_total))
  return(T_corrected)
}

# Calculate elapsed time for each measurement from the start
kes22$elapsed_time <- as.numeric(difftime(kes22$Date.Time, min(kes22$Date.Time), units = "mins"))

# Apply the correction to each temperature reading
kes22$Temp_Corrected <- mapply(correct_temp, kes22$Temp, kes22$elapsed_time, MoreArgs = list(T_initial_O = T_initial_O, T_final_O = T_final_O, t_total = t_total))

########################
##Resolve background concentration for relative humidity

# Let us resolve and correct Humidity also
RH_initial_O <- 87.0533  # initial observed R.H.
RH_final_O <- 83.03    # final observed R.H.
RH_total <- 70      # total time of your measurements

# Function to apply background correction
correct_HUM <- function(measured_rh, elapsed_time, RH_initial_O, RH_final_O, RH_total) {
  RH_corrected <- measured_rh - ((RH_final_O - RH_initial_O) * (elapsed_time / RH_total))
  return(RH_corrected)
}

# Calculate elapsed time for each measurement from the start
kes22$elapsed_time_rh <- as.numeric(difftime(kes22$Date.Time, min(kes22$Date.Time), units = "mins"))

# Apply the correction to each Relative Humidity reading
kes22$Rel.Hum_Corrected <- mapply(correct_HUM, kes22$Rel..Hum., kes22$elapsed_time_rh, MoreArgs = list(RH_initial_O = RH_initial_O, RH_final_O = RH_final_O, RH_total = RH_total))


# Take the time frame we want i.e. the one for the first transect when GPS Started taking reading
start_t <- '11:21:00'  # Adjust format if necessary
end_t <- '12:20:00'


filtered_kes22 <- filter(kes22, Date.Time >= as.POSIXct(paste('2024-02-20', start_t, sep = " ")),
                         Date.Time <= as.POSIXct(paste('2024-02-20', end_t, sep = " ")))


######################################################

# Merge environmental data with NO2 data
merged_data <- merge(filtered_kes22, df, by = 'Date.Time')

# Linear model for predicting NO2 based on environmental factors
lm_model <- lm(NO2_corrected ~ Temp_Corrected + Rel.Hum_Corrected + Wind.Speed + Headwind, data = merged_data)
summary(lm_model)

# Plot residuals vs predicted values for linearity test
plot(predict(lm_model), residuals(lm_model))
title("Residuals vs Predicted Values")

# Normality of residuals using Q-Q plot
qqnorm(residuals(lm_model))
qqline(residuals(lm_model))



# Get regression summary
full_morning_summary <- summary(lm_model)

# Save summary to CSV
write.csv(full_morning_summary$coefficients, "full_morning_regression_summary.csv")




####################################################################################

step_model <- step(lm(NO2_corrected ~ Temp_Corrected + Rel.Hum_Corrected + Wind.Speed + Headwind, data = merged_data))
summary(step_model)
##########################################

# Plot residuals vs predicted values for linearity test
plot(predict(step_model), residuals(step_model))
title("Residuals vs Predicted Values")

# Normality of residuals using Q-Q plot
qqnorm(residuals(step_model))
qqline(residuals(step_model))


# Get regression summary
stepwise_morning_summary <- summary(step_model)

# Save summary to CSV
write.csv(stepwise_morning_summary$coefficients, "stepwise_morning_regression_summary.csv")

#####################################

ggplot(merged_data, aes(x = Temp_Corrected, y = NO2_corrected)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Linear fit line
  labs(title = "NO2 vs Temperature (11 AM to 1 PM)", x = "Temperature (°C)", y = "NO2 (µg/m³)")






######################################################
# Let us do the same for the third transect
kes3 <- read.csv(file.choose(), skip = 9)

# Remove the empty first row
kes3 <- kes3[-1, ]

# Delete column X as it only has NA as values
kes3 <- kes3[, -17]

# Convert to datetime or POSIXct, and the other variables to numeric
kes3$Time <- as.POSIXct(kes3$Time, format = "%Y-%m-%d %H:%M")



#############
# Convert every column in kes to numeric apart from Time column
kes3 <- kes3 %>% mutate_at(vars(-Time), as.numeric)

## Resample data to 1-minute intervals
kes33 <- kes3 %>% group_by(Date.Time = floor_date(Time, "1 min")) %>% 
  summarise(across(c(Temp, Wet.Bulb.Temp., Rel..Hum., Baro., Altitude, Station.P., Wind.Speed, 
                     Heat.Index, Dew.Point, Dens..Alt., Crosswind, Headwind, Mag..Dir.,
                     True.Dir., Wind.Chill), mean))

kes33 <- data.frame(kes33)

############################################################

##############################################################

# Let us select only 4 predictors

kes33 <-  kes33[, c("Date.Time", "Temp", "Rel..Hum.", 'Headwind', 'Wind.Speed')]

#### Pick the start and end time for temperature i.e. start to finish of transect
s_t <- '16:35:00'  # Adjust format if necessary
e_t <- '17:38:00'


kes33 <- filter(kes33, Date.Time >= as.POSIXct(paste('2024-02-20', s_t, sep = " ")),
                Date.Time <= as.POSIXct(paste('2024-02-20', e_t, sep = " ")))

####### resolve and correct the temperature
######
# Background temperature data

T_initial_O <- 14.0133
T_final_O <- 13.01  
t_total <- 63     

# Function to apply background correction
correct_temp <- function(measured_temp, elapsed_time, T_initial_O, T_final_O, t_total) {
  T_corrected <- measured_temp - ((T_final_O - T_initial_O) * (elapsed_time / t_total))
  return(T_corrected)
}

# Calculate elapsed time for each measurement from the start
kes33$elapsed_time <- as.numeric(difftime(kes33$Date.Time, min(kes33$Date.Time), units = "mins"))

# Apply the correction to each temperature reading
kes33$Temp_Corrected <- mapply(correct_temp, kes33$Temp, kes33$elapsed_time, MoreArgs = list(T_initial_O = T_initial_O, T_final_O = T_final_O, t_total = t_total))


################## Let us resolve and correct Humidity also
###### ##Resolve background concentration for relative humidity

# Let us resolve and correct Humidity also
RH_initial_O <- 65.64667
RH_final_O <- 80.20333
RH_total <- 63

# Function to apply background correction
correct_HUM <- function(measured_rh, elapsed_time, RH_initial_O, RH_final_O, RH_total) {
  RH_corrected <- measured_rh - ((RH_final_O - RH_initial_O) * (elapsed_time / RH_total))
  return(RH_corrected)
}

# Calculate elapsed time for each measurement from the start
kes33$elapsed_time_rh <- as.numeric(difftime(kes33$Date.Time, min(kes33$Date.Time), units = "mins"))

# Apply the correction to each temperature reading
kes33$Rel.Hum_Corrected <- mapply(correct_HUM, kes33$Rel..Hum., kes33$elapsed_time_rh, MoreArgs = list(RH_initial_O = RH_initial_O, RH_final_O = RH_final_O, RH_total = RH_total))

# Take the time frame we want i.e. the one for the third transect
start_ <- '16:35:00'  # Adjust format if necessary
end_ <- '17:38:00'


filtered_kes33 <- filter(kes33, Date.Time >= as.POSIXct(paste('2024-02-20', start_, sep = " ")),
                         Date.Time <= as.POSIXct(paste('2024-02-20', end_, sep = " ")))


###################################################################
merged_data2 <- merge(filtered_kes33, data2, by = 'Date.Time')

lm_2 <- lm(NO2_corrected ~ Temp_Corrected + Rel.Hum_Corrected + Headwind, data = merged_data2)

summary(lm_2)


############################ Using Step Wise


step_model2 <- step(lm(NO2_corrected ~ Temp_Corrected + Rel.Hum_Corrected + Wind.Speed + Headwind, data = merged_data2))
summary(step_model2)


#############################

plot(predict(step_model2), residuals(step_model2))
qqnorm(residuals(step_model2))
qqline(residuals(step_model2))





######################################################
########## TRANSECT 1 (Morning)

######### Plots for Temperature

# Summary statistics for each variable
summary(merged_data)

# Calculate additional descriptive statistics
summary_stat_morning <- merged_data %>%
  summarise(
    mean_NO2 = mean(NO2_corrected, na.rm = TRUE),
    sd_NO2 = sd(NO2_corrected, na.rm = TRUE),
    var_NO2 = var(NO2_corrected, na.rm = TRUE),
    median_NO2 = median(NO2_corrected, na.rm = TRUE),
    IQR_NO2 = IQR(NO2_corrected, na.rm = TRUE),
    mean_temp = mean(Temp_Corrected, na.rm = TRUE),
    sd_temp = sd(Temp_Corrected, na.rm = TRUE),
    var_temp = var(Temp_Corrected, na.rm = TRUE),
    median_temp = median(Temp_Corrected, na.rm = TRUE),
    IQR_temp = IQR(Temp_Corrected, na.rm = TRUE),
    mean_RelHum = mean(Rel.Hum_Corrected, na.rm = TRUE),
    sd_RelHum = sd(Rel.Hum_Corrected, na.rm = TRUE),
    var_RelHum = var(Rel.Hum_Corrected, na.rm = TRUE),
    median_RelHum = median(Rel.Hum_Corrected, na.rm = TRUE),
    IQR_RelHum = IQR(Rel.Hum_Corrected, na.rm = TRUE),
    mean_hw = mean(Headwind, na.rm = TRUE),
    sd_hw = sd(Headwind, na.rm = TRUE),
    var_hw = var(Headwind, na.rm = TRUE),
    median_hw = median(Headwind, na.rm = TRUE),
    IQR_hw = IQR(Headwind, na.rm = TRUE),
    mean_hw = mean(Wind.Speed, na.rm = TRUE),
    sd_hw = sd(Wind.Speed, na.rm = TRUE),
    var_hw = var(Wind.Speed, na.rm = TRUE),
    median_hw = median(Wind.Speed, na.rm = TRUE),
    IQR_hw = IQR(Wind.Speed, na.rm = TRUE),
  )



# Reshape the data to a long format for easier plotting with facets
t1_long <- pivot_longer(merged_data, cols = c(Temp_Corrected, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t1_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Temp_Corrected" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Corrected Temperature (°C) and NO2 Levels (µg/m³) Over Time (Morning)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')





######### Plots for Relative Humidity

t2_long <- pivot_longer(merged_data, cols = c(Rel.Hum_Corrected, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t2_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Rel.Hum_Corrected" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Corrected Relative Humidity (%) and NO2 Levels (µg/m³) Over Time (Morning)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')


#####

##########
######### Plots for Headwind

t3_long <- pivot_longer(merged_data, cols = c(Headwind, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t3_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Headwind" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Graph of Headwind (m/s) and NO2 Levels (µg/m³) Over Time (Morning)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')



##################################################

###### Wind Speed


t4_long <- pivot_longer(merged_data, cols = c(Wind.Speed, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t4_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Wind.Speed" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Corrected Windspeed (m/s) and NO2 Levels (µg/m³) Over Time (Morning)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')



#######################################################################
#######################################

###################################################
#### THIRD TRANSECT

# Summary statistics for each variable
summary(merged_data2)

# Calculate additional descriptive statistics
summary_stat_evening <- merged_data2 %>%
  summarise(
    mean_NO2 = mean(NO2_corrected, na.rm = TRUE),
    sd_NO2 = sd(NO2_corrected, na.rm = TRUE),
    var_NO2 = var(NO2_corrected, na.rm = TRUE),
    median_NO2 = median(NO2_corrected, na.rm = TRUE),
    IQR_NO2 = IQR(NO2_corrected, na.rm = TRUE),
    mean_temp = mean(Temp_Corrected, na.rm = TRUE),
    sd_temp = sd(Temp_Corrected, na.rm = TRUE),
    var_temp = var(Temp_Corrected, na.rm = TRUE),
    median_temp = median(Temp_Corrected, na.rm = TRUE),
    IQR_temp = IQR(Temp_Corrected, na.rm = TRUE),
    mean_RelHum = mean(Rel.Hum_Corrected, na.rm = TRUE),
    sd_RelHum = sd(Rel.Hum_Corrected, na.rm = TRUE),
    var_RelHum = var(Rel.Hum_Corrected, na.rm = TRUE),
    median_RelHum = median(Rel.Hum_Corrected, na.rm = TRUE),
    IQR_RelHum = IQR(Rel.Hum_Corrected, na.rm = TRUE),
    mean_hw = mean(Headwind, na.rm = TRUE),
    sd_hw = sd(Headwind, na.rm = TRUE),
    var_hw = var(Headwind, na.rm = TRUE),
    median_hw = median(Headwind, na.rm = TRUE),
    IQR_hw = IQR(Headwind, na.rm = TRUE)
  )



#### Plotting with Facet for third transect

# Reshape the data to a long format for easier plotting with facets
t5_long <- pivot_longer(merged_data2, cols = c(Temp_Corrected, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t5_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Temp_Corrected" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Corrected Temperature (°C) and NO2 Levels (µg/m³) Over Time (Evening)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')



######### Plots for Relative Humidity

t6_long <- pivot_longer(merged_data2, cols = c(Rel.Hum_Corrected, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t6_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Rel.Hum_Corrected" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Corrected Relative Humidity (%) and NO2 Levels (µg/m³) Over Time (Evening)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')


#####

##########
######### Plots for Headwind

t7_long <- pivot_longer(merged_data2, cols = c(Headwind, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t7_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Headwind" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Graph of Headwind (m/s) and NO2 Levels (µg/m³) Over Time (Evening))") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')



##################################################

###### Wind Speed


t8_long <- pivot_longer(merged_data2, cols = c(Wind.Speed, NO2_corrected), 
                        names_to = "variable", values_to = "value",
                        names_repair = "unique")

# Enhanced Combined Plot with Improved Aesthetics
ggplot(t8_long, aes(x = Date.Time, y = value, color = variable)) +
  geom_line(size = 1.5) + # Thicker lines for better visibility
  geom_point(size = 2, shape = 21, fill = "white") + # Points with white fill
  scale_color_manual(values = c("Wind.Speed" = "#377EB8", "NO2_corrected" = "#E41A1C")) +
  labs(y = "Measurement", x = 'Time', 
       title = "Graph of Windspeed (m/s) and NO2 Levels (µg/m³) Over Time (Evening)") +
  theme_minimal(base_size = 13) + # Larger base text size for readability
  theme(legend.position = "bottom", # Move legend to bottom
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14), # Centered and larger title
        axis.text.x = element_text(), # Rotate x-axis text
        panel.background = element_rect(fill = "white"), # Light gray background
        panel.grid.major = element_line(color = "gray70"), # Custom major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  facet_grid(variable ~ ., scales = 'free_y')



##########################################################################################
################ FULL TRANSECT #####################################

############# PLOT ALL

Transect1 <- merged_data



Transect3 <- merged_data2

#############
#Join the two dataframes

full_transect <- rbind(Transect1, Transect3)

summary(full_transect)

####################### Model
lm_full <- lm(NO2_corrected ~ Temp_Corrected + Rel.Hum_Corrected + Wind.Speed + Headwind, data = full_transect)

summary(lm_full)

##############################################################
##### Step Wise model for Full Transect ####################################


step_full_transect <- step(lm(NO2_corrected ~ Temp_Corrected + Rel.Hum_Corrected + Wind.Speed + Headwind, data = full_transect))
summary(step_full_transect)


#############################

plot(predict(step_full_transect), residuals(step_full_transect))
qqnorm(residuals(step_full_transect))
qqline(residuals(step_full_transect))


###  the plots indicate Linearity and are normally distribured
###############################


# Create the scatter plot with ggplot2 for full transect

#####Temperature

ggplot(full_transect, aes(x = Temp_Corrected, y = NO2_corrected)) +
  geom_point(aes(color = time_of_day, shape = time_of_day), size = 2.2) +  # Add points
  geom_smooth(method = "lm", colour = "red", se = FALSE) + # Add regression line
  labs(title = "Scatter Plot of NO2 vs Temperature with a Fitted Regression Line",
       x = "Temperature (°C)",
       y = "NO2_Corrected (µg/m³)", color ='Time of Day', shape='Time of Day') +
  scale_color_manual(values = c("Morning" = "blue", "Evening" = "orange")) + # Customize colors
  scale_shape_manual(values = c("Morning" = 16, "Evening" = 17)) +
  theme(axis.title.x = element_text(vjust = 0, size = 13, face='italic'),
        axis.title.y = element_text(vjust = 2, size = 13, face='italic'),
        axis.text = element_text(size = 12),
        plot.title = element_text(face='bold.italic', size = 13),
        panel.background = element_rect(fill = "white"),  # Set panel background to white
        #panel.grid.major = element_line(color = "grey", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5)) 


## RELATIVE HUMIDITY

ggplot(full_transect, aes(x = Rel.Hum_Corrected, y = NO2_corrected)) +
  geom_point(aes(color = time_of_day, shape = time_of_day), size = 2.2) +  # Add points
  geom_smooth(method = "lm", colour = "red", se = FALSE) + # Add regression line
  labs(title = "Scatter Plot of NO2 vs Relative Humidity with a Fitted Regression Line",
       x = "Relative Humidity (%)",
       y = "NO2_Corrected (µg/m³)", color ='Time of Day', shape='Time of Day') +
  scale_color_manual(values = c("Morning" = "blue", "Evening" = "orange")) + # Customize colors
  scale_shape_manual(values = c("Morning" = 16, "Evening" = 17)) +
  theme(axis.title.x = element_text(vjust = 0, size = 13, face='italic'),
        axis.title.y = element_text(vjust = 2, size = 13, face='italic'),
        axis.text = element_text(size = 12),
        plot.title = element_text(face='bold.italic', size = 13),
        panel.background = element_rect(fill = "white"),  # Set panel background to white
        #panel.grid.major = element_line(color = "grey", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5)) 



###### Headwind

ggplot(full_transect, aes(x = Headwind, y = NO2_corrected)) +
  geom_point(aes(color = time_of_day, shape = time_of_day), size = 2.2) +  # Add points
  geom_smooth(method = "lm", colour = "red", se = FALSE) + # Add regression line
  labs(title = "Scatter Plot of NO2 vs Headwind with a Fitted Regression Line",
       x = "Headwind (m/s)",
       y = "NO2_Corrected (µg/m³)", color ='Time of Day', shape='Time of Day') +
  scale_color_manual(values = c("Morning" = "blue", "Evening" = "orange")) + # Customize colors
  scale_shape_manual(values = c("Morning" = 16, "Evening" = 17)) +
  theme(axis.title.x = element_text(vjust = 0, size = 13, face='italic'),
        axis.title.y = element_text(vjust = 2, size = 13, face='italic'),
        axis.text = element_text(size = 12),
        plot.title = element_text(face='bold.italic', size = 13),
        panel.background = element_rect(fill = "white"),  # Set panel background to white
        #panel.grid.major = element_line(color = "grey", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5)) 


####################  WINDSPEED

ggplot(full_transect, aes(x = Wind.Speed, y = NO2_corrected)) +
  geom_point(aes(color = time_of_day, shape = time_of_day), size = 2.2) +  # Add points
  geom_smooth(method = "lm", colour = "red", se = FALSE) + # Add regression line
  labs(title = "Scatter Plot of NO2 vs Windspeed with a Fitted Regression Line",
       x = "Wind Speed (m/s)",
       y = "NO2_Corrected (µg/m³)", color ='Time of Day', shape='Time of Day') +
  scale_color_manual(values = c("Morning" = "blue", "Evening" = "orange")) + # Customize colors
  scale_shape_manual(values = c("Morning" = 16, "Evening" = 17)) +
  theme(axis.title.x = element_text(vjust = 0, size = 13, face='italic'),
        axis.title.y = element_text(vjust = 2, size = 13, face='italic'),
        axis.text = element_text(size = 12),
        plot.title = element_text(face='bold.italic', size = 13),
        panel.background = element_rect(fill = "white"),  # Set panel background to white
        #panel.grid.major = element_line(color = "grey", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5)) 



###############################################################################
#######################################################################


