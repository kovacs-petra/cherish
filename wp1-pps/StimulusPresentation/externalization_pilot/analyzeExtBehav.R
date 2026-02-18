# Analyze and plot externalization pilot data for CherISH WP1

library(tidyverse)
library(scales)

# Read data
d <- read.csv("C:\\Users\\pkovacs\\Documents\\GitHub\\cherish\\wp1-pps\\StimulusPresentation\\externalization_pilot\\externData.csv")

# Convert to factor where necessary
d <- d %>% mutate_at(c('subNum','sourceInt','f0cond','f0',
                  'azimuth'),as.factor)

# Rescale externalization rating
d <- d %>% mutate(extRatingRescaled = 
                    rescale(d$extRating, to=c(-10,10), from=range(d$extRating)))

# Plot distance vs. absolute rating
ggplot(d, aes(distance,abs(extRatingRescaled))) +
  geom_boxplot(aes(group = distance, color = sourceInt)) +
  facet_grid(rows = vars(sourceInt),cols= vars(f0cond)) +
  geom_smooth(method = "lm", se = FALSE)

# Calculate correlation
corr_distVrating <- cor(d$distance,abs(d$extRatingRescaled), method = "pearson")
