# Analyze and plot externalization pilot data for CherISH WP1

library(tidyverse)
library(scales)

# Read data
e <- read.csv("\\\\KFS.oeaw.ac.at\\Fileserver\\ProjektDaten\\CherISH\\data\\wp-1\\Behav\\extern_task\\externData.csv")

# Convert to factor where necessary
e <- e %>% mutate_at(c('participant','sourceInt','trajectory','side','respType'),
                     as.factor)
e <- e %>% mutate_at(c('onsetDistance','offsetDistance'),
                     as.numeric)
e <- e %>% mutate_at(c('extRating'),as.integer)

# Rescale externalization rating
e <- e %>% mutate(extRatingRescaled = 
                    rescale(e$extRating, to=c(-10,10), from=range(e$extRating)))

# Use rescale to find out the head and PPS boundaries, respectively
# rescale(c(0,680,815,1500), to = c(-10,10), from = range(c(0,1500)))
intern <- rescale(c(0,709,791,1500), to = c(-10,10), from = range(c(0,1500)))
PPS <- rescale(c(0,449,1051,1500), to = c(-10,10), from = range(c(0,1500)))
intern <- abs(intern[2])
PPS <- abs(PPS[2])

# Calculate correlation
corr_distVrating <- cor(e$distance,abs(e$extRatingRescaled), method = "pearson")

# Make verbose variables for plotting
e <- e %>% mutate(source_intensity = case_when(
  e$sourceInt == 0 ~ "low_intensity",
  e$sourceInt == 1 ~ "high_intensity"),
  sound_trajectory = case_when(
    e$trajectory == 0 ~ "stationary",
    e$trajectory == 1 ~ "looming",
    e$trajectory == 2 ~ "receding",
    e$trajectory == 3 ~ "rotate_PPS",
    e$trajectory == 4 ~ "rotate_EPS"),
  response_type = case_when(
    e$respType == 1 ~ "stationary",
    e$respType == 2 ~ "onset",
    e$respType == 3 ~ "offset",
  ))
e <- e %>% mutate_at(c('source_intensity','sound_trajectory','response_type'),as.factor)

# Filter for the three types of response
stat <- e %>% filter(respType == 1)
onset <- e %>% filter(respType == 2)
offset <- e %>% filter(respType == 3)

### Plots grouped by intensity
# Plot stationary responses
ggplot(stat, aes(x=onsetDistance,y=abs(extRatingRescaled),fill=source_intensity)) +
  geom_boxplot(aes(group=paste(onsetDistance,source_intensity))) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 2, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 2, yend = intern), color = "green",
               linetype = "dashed") +
  labs(title="Responses to stationary sounds")

# Plot onset responses
ggplot(onset, aes(x=onsetDistance,y=abs(extRatingRescaled),fill=source_intensity)) +
  geom_boxplot(aes(group=paste(onsetDistance,source_intensity))) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 2, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 2, yend = intern), color = "green",
               linetype = "dashed") +
  labs(title="Responses to the onset of moving sounds")

# Plot offset responses
ggplot(offset, aes(x=offsetDistance,y=abs(extRatingRescaled),fill=source_intensity)) +
  geom_boxplot(aes(group=paste(offsetDistance,source_intensity))) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 2, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 2, yend = intern), color = "green",
               linetype = "dashed") +
  labs(title="Responses to the offset of moving sounds")

### Plots grouped by trajectory

# Filter for trajectories
loom <- e %>% filter(trajectory==1)
rec <- e %>% filter(trajectory==2)
pps <- e %>% filter(trajectory==3)
eps <- e %>% filter(trajectory==4)

# Plot onset and offset responses by trajectory
ggplot(e, aes(x=sound_trajectory,y=abs(extRatingRescaled),fill=factor(response_type, level = c("onset","offset")))) +
  geom_boxplot(aes(group=paste(sound_trajectory,response_type))) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 5, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 5, yend = intern), color = "green",
               linetype = "dashed") +
  scale_x_discrete(limits = c("looming","receding","rotate_PPS","rotate_EPS")) +
  scale_fill_discrete(limits = c("onset","offset")) +
  labs(title="Offset and onset responses",
       x = "Sound trajectory",
       y = "Distance rating [a.u.]",
       fill = "Response type")


