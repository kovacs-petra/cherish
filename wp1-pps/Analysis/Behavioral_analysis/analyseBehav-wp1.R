# Analyze and plot externalization pilot data for CherISH WP1

library(tidyverse)
library(lme4)
library(broom.mixed)

# Read data
d <- read.csv("\\\\kfs.oeaw.ac.at\\Fileserver\\ProjektDaten\\CherISH\\data\\wp-1\\Behav\\main_task\\mainTaskData.csv")

# Convert to factor where necessary
d <- d %>% mutate_at(c('subNum','onsetDistance','offsetDistance','direction',
                       'trajectory','offsetAzimuth','targetTrial','congruence',
                       'targetAzimuth', 'sourceInt','accuracy', 'trigger'), as.factor)
d <- d %>% mutate_at(c('blockNo','trialNo','frequency','totalDur',
                       'durStatOnset','respTime','promptness'), as.numeric)

# Add verbose variables for plotting
d <- d %>% mutate(
  source_intensity = case_when(
    d$sourceInt == 0 ~ "low_intensity",
    d$sourceInt == 1 ~ "high_intensity"),
  distance_at_offset = case_when(
    d$offsetDistance == 0.2 ~ "PPS",
    d$offsetDistance == 2 ~ "EPS"),
  sound_direction = case_when(
    d$direction == 1 ~ "radial",
    d$direction == 2 ~ "angular"),
  sound_trajectory = case_when(
    d$trajectory == 1 ~ "looming",
    d$trajectory == 2 ~ "receding",
    d$trajectory == 3 ~ "rotate_PPS",
    d$trajectory == 4 ~ "rotate_EPS"))

d <- d %>% mutate_at(c('source_intensity','distance_at_offset',
                       'sound_direction','sound_trajectory'),as.factor)

# Filter out too fast and slow values (this also gets rid of NaN values)
d_rt <- d %>% filter(respTime >= 100)
d_rt <- d_rt %>% filter(respTime <= 1000)

# Filter for accurate responses
d_acc <- d_rt %>% filter(accuracy == 1)

# Plot promptness vs. trajectory
ggplot(d_acc, aes(sound_trajectory, promptness)) +
  geom_boxplot(width=.5) +
  # Reorder trajectory factors
  scale_x_discrete(limits = c("looming","receding","rotate_PPS","rotate_EPS")) +
  theme_light()

# Plot also RT vs. trajectory for good measure
ggplot(d, aes(sound_trajectory, respTime)) +
  geom_boxplot()

# Group by source intensity
ggplot(d_acc, aes(sound_trajectory, promptness,fill=source_intensity)) +
  geom_boxplot(aes(group=paste(sound_trajectory,source_intensity)),width=.3) +
  # Reorder trajectory factors
  scale_x_discrete(limits = c("looming","receding","rotate_PPS","rotate_EPS"))

# Fit a model
lm1 <- lm(promptness ~ sound_direction + distance_at_offset + 
              source_intensity + sound_direction:distance_at_offset,
            data = d)

tidy(lm1, conf.int = TRUE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Resources
ggplot(d, aes(x=distance,y=abs(extRatingRescaled),fill=source_intensity)) +
  geom_boxplot(aes(group=paste(distance,source_intensity))) +
  # geom_point(aes(group = source_intensity), position=position_dodge(width=.15),
  # shape="x") +
  facet_wrap(vars(stimulus_type)) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 2, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 2, yend = intern), color = "green",
               linetype = "dashed")

# # Now facet by intensity, and group by stim type (so, the other way around)
ggplot(d, aes(x=distance,y=abs(extRatingRescaled),fill=stimulus_type)) +
  geom_boxplot(aes(group=paste(distance,stimulus_type))) +
  facet_wrap(vars(source_intensity)) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 2, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 2, yend = intern), color = "green",
               linetype = "dashed")

# Facet by subject
ggplot(d, aes(x=distance,y=abs(extRatingRescaled),fill=source_intensity)) +
  geom_boxplot(aes(group=paste(distance,source_intensity,stimulus_type))) +
  geom_point(aes(group = source_intensity), position=position_dodge(width=.15),
             shape = "x") +
  facet_grid(rows=vars(stimulus_type),cols=vars(subNum)) +
  # Mark end of arms with dashed horizontal line:
  geom_segment(aes(x = 0, y = PPS, xend = 2, yend = PPS), color = "red",
               linetype = "dashed") +
  # Mark out-of-head border with dashed horizontal line:
  geom_segment(aes(x = 0, y = intern, xend = 2, yend = intern), color = "green",
               linetype = "dashed")

