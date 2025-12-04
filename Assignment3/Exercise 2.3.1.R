# install.packages(c("ctsmTMB","data.table","ggplot2","gridExtra"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(patchwork) 

df <- read_csv("ex3_largecase.csv")
df$Timestamp <- ymd_hms(df$Timestamp, tz = "UTC")
unique_events <- unique(df$Event_ID)

############################
#2.3.1. Plot all the events 
############################

library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(ctsmTMB)
library(gridExtra)
library(cowplot)

data <- read.csv("ex3_largecase.csv")
data$Timestamp <- as.POSIXct(data$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
events <- unique(data$Event_ID)

for (event in events) {
  event_data <- data %>% filter(Event_ID == event)
  
  # Plot Rainfall
  rainfall_plot <- ggplot(event_data, aes(x = Timestamp, y = Rainfall)) +
    geom_line(color = "blue") +
    scale_x_datetime(date_breaks = "6 hours", date_labels = "%H:%M") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(title = paste("Rainfall - Event", event), x = "Time", y = "Rainfall (µm/min)") +
    theme_minimal()
  
  # Plot Pumpflow
  pumpflow_plot <- ggplot(event_data, aes(x = Timestamp, y = Pumpflow)) +
    geom_line(color = "red") +
    scale_x_datetime(date_breaks = "6 hours", date_labels = "%H:%M") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(title = paste("Pumpflow - Event", event), x = "Time", y = "Pumpflow (m^3/min)") +
    theme_minimal()
  
  # Plot Volume
  volume_plot <- ggplot(event_data, aes(x = Timestamp, y = Volume)) +
    geom_line(color = "darkgreen") +
    scale_x_datetime(date_breaks = "6 hours", date_labels = "%H:%M") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(title = paste("Volume - Event", event), x = "Time", y = "Volume (m^3)") +
    theme_minimal()
  
  all_plots <- append(all_plots, list(rainfall_plot, pumpflow_plot, volume_plot))
}

combined_plot <- plot_grid(
  plotlist = all_plots,
  ncol = 3,
  align = "v", 
  axis = "l"   
)

print(combined_plot)
ggsave("C:/Users/alba/Desarrollo_DTU/2025_Fall/02427_Advanced_Time_Series/Computer_Exercise3/combined_plot_2_3_1.png",
       plot = combined_plot, width = 20.69,
       height = 8.27, units = "in", dpi = 300)

