
#create cadmium plot for report

# Copyright 2017 Province of British Columbia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.


library (ggplot2)
library(readr)
library(dplyr) # JS: For arrange()
library(magrittr) # JS: for %<>%)

Cd<- read_csv("C:/R-repositories/SSD-fitting/Cd-example/Cd_BC_Multi.csv")
View(Cd)

#specificy factor for later sorting
Cd$Group <- factor(Cd$Group, levels = c("Amphibian","Fish (non-salmonid)",
                                        "Fish (salmonid)","Invertebrate",
                                        "Plant"))

# Arrange by Group (factor) and species (in alphabetic order)
Cd %<>% arrange(Group, Species)

# Create a data frame for the plotting position of each species
PlotPos <- select(Cd, Group, Species) %>%
  arrange(Group, Species) %>%
  unique()

# Create a plotting position for each species
PlotPos %<>% mutate(Position = 1:nrow(PlotPos))

# Update Cd with plotting positions
Cd %<>% left_join(PlotPos, by = c("Species", "Group"))

# # Force ggplot to arrange by factor, created above
Cdplot<-ggplot(data = Cd, aes(x = Conc, y = Position, shape = Group)) +
  # JS: moved size argument up to geom_point
  geom_point(position = "jitter", size = 2) +
  scale_x_log10(breaks = 10 ^ (-2:3)) +
  # JS: Reverse the plotting position (i.e., 1 at the top), create a series of 
  # breaks (number of species), and use the linked species names as the labels
  scale_y_continuous(trans = "reverse", breaks = seq(1, max(Cd$Position)),
                                                   labels = unique(Cd$Species)) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4)) +
  xlab("Log Concentration") +
  # JS: remove 'Position' as the label
  ylab("") +
  #JS:  Remove minor grid lines; too busy
  theme(panel.grid.minor = element_blank())

Cdplot

# JS: Changed plot dimensions to account for legend; looks a bit better
ggsave("Plots/Cdplot.png", Cdplot, width = 8.80, height = 5.57)
