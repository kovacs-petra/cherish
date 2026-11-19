# Analyze and plot externalization task data for CherISH WP1

library(tidyverse)
library(lme4)
library(broom.mixed)
library(modelsummary)
library(scales)

# Set colors
loom_color = "#dc3220"
rec_color = "#1A85FF"
pps_color = "#4B0092"
eps_color = "#1AFF1A"

# Read data
path_to_data <- "\\promptness_ERP_artrej.csv"
d <- path_to_data

# Rename some variables
d <- rename(d, motion_condition = trajectory, motion_trajectory = direction)

# Convert to factor where necessary
d <- d %>% mutate_at(c('subNum','onsetDistance','offsetDistance','motion_trajectory',
                       'motion_condition','offsetAzimuth','targetTrial','congruence',
                       'targetAzimuth', 'sourceInt','accuracy', 'trigger'), as.factor)
d <- d %>% mutate_at(c('blockNo','trialNo','frequency','totalDur',
                       'durStatOnset','respTime','promptness','dbFS','meanERP'), as.numeric)

# Add verbose variables for plotting
d <- d %>% mutate(
  source_intensity = case_when(
    d$sourceInt == 0 ~ "low_intensity",
    d$sourceInt == 1 ~ "high_intensity"),
  distance_at_offset = case_when(
    d$offsetDistance == 0.2 ~ "PPS",
    d$offsetDistance == 2 ~ "EPS"),
  motion_trajectory = case_when(
    d$motion_trajectory == 1 ~ "radial",
    d$motion_trajectory == 2 ~ "angular"),
  motion_condition = case_when(
    d$motion_condition == 1 ~ "Looming",
    d$motion_condition == 2 ~ "Receding",
    d$motion_condition == 3 ~ "Rotating in PPS",
    d$motion_condition == 4 ~ "Rotating in EPS"))

d <- d %>% mutate_at(c('source_intensity','distance_at_offset',
                       'motion_trajectory','motion_condition'),as.factor)

# Filter out too fast and slow values (this also gets rid of NaN values)
d_rt <- d %>% filter(respTime >= 100)
d_rt <- d_rt %>% filter(respTime <= 1000)

# Filter for accurate responses
d_acc <- d_rt %>% filter(accuracy == 1) %>% 
  drop_na(meanERP)

# # Plot promptness vs. trajectory
ggplot(d_acc, aes(motion_condition, promptness)) +
  geom_boxplot(fill = c(loom_color,rec_color,pps_color,eps_color),width=.2) +
  # Reorder trajectory factors
  scale_x_discrete(limits = c("Looming","Receding","Rotating in PPS","Rotating in EPS")) +
  labs(x = "Cue motion condition",
       y = "Promptness (1/s)") +
  # scale_color_manual(values=c(loom_color,rec_color,pps_color,eps_color)) +
  scale_y_continuous(labels = label_number(scale = 1e3)) +
  annotate("text", x = -Inf,y = Inf,label = "1e-3", hjust = -0.1,
           vjust = 1.1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 10, 10, 30)) +
  theme_classic(base_size = 20)

# Plot also RT vs. trajectory for good measure
# ggplot(d, aes(sound_trajectory, respTime)) +
#   geom_boxplot()

# Group by source intensity
# ggplot(d_acc, aes(sound_trajectory, promptness,fill=source_intensity)) +
#   geom_boxplot(aes(group=paste(sound_trajectory,source_intensity)),width=.3) +
#   # Reorder trajectory factors
#   scale_x_discrete(limits = c("looming","receding","rotate_PPS","rotate_EPS"))

# Plot by dbFS
ggplot(d_acc,aes(dbFS,promptness,color=motion_condition)) +
  geom_smooth(aes(linetype = motion_trajectory,fill=distance_at_offset), se = TRUE,
              method = "lm", alpha = .1, linewidth = 2) +
  labs(x = "Source intensity (dB FS)",
       y = "Promptness (1/s)",
       linetype = "Motion trajectory",
       fill = "Distance at offset",
       color = "Motion condition") +
  scale_color_manual(values=c(loom_color,rec_color,pps_color,eps_color)) +
  scale_y_continuous(labels = label_number(scale = 1e3)) +
  annotate("text", x = -Inf,y = Inf,label = "1e-3", hjust = -0.1,
           vjust = 1.1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 10, 10, 30)) +
  theme_classic(base_size = 15)

# Fit a model using the numeric dbFS values as intensity predictors
lm_behav <- lmer(promptness ~ motion_trajectory + distance_at_offset + dbFS +
                   motion_trajectory:distance_at_offset + 
                    (1|subNum), data = d_acc)

lm_behav_tidy <- tidy(lm_behav, conf.int = TRUE)

modelplot(lm_behav,coef_map = c("motion_trajectorradial" = "Motion trajectory (radial)",
                                "distance_at_offsetPPS" = "Distance at offset (PPS)",
                                "dbFS" = "dbFS",
                                "motion_trajectoryradial:distance_at_offsetPPS"="Motion trajectory (radial) : Distance at offset (PPS)"))


# Plot by the three factors examined in the model
ggplot(d_acc, aes(motion_trajectory, promptness,fill=distance_at_offset)) +
  geom_boxplot(aes(group=paste(motion_trajectory,distance_at_offset)),width=.3) +
  facet_wrap(vars(source_intensity),nrow=1,ncol=2) 

# Fit a neurobehavioral model including mean ERP as a predictor
lm_neurobehav <- lmer(promptness ~ 
                        motion_trajectory + distance_at_offset + dbFS + meanERP +
                        motion_trajectory:distance_at_offset + 
                        (1|subNum), data = d_acc)

lm_neurobehav_tidy <- tidy(lm_neurobehav, conf.int = TRUE)

modelplot(lm_neurobehav,coef_map = c("motion_trajectorradial" = "Motion trajectory (radial)",
                                "distance_at_offsetPPS" = "Distance at offset (PPS)",
                                "dbFS" = "dbFS",
                                "motion_trajectoryradial:distance_at_offsetPPS"="Motion trajectory (radial) : Distance at offset (PPS)",
                                "meanERP" = "Mean ERP amplitude"))



