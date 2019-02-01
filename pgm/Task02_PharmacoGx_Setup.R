library(PharmacoGx)
library(dplyr)
a = availablePSets()
a.keep = a %>% filter(Dataset.Type == "sensitivity")

for (name in a.keep$PSet.Name) {
  print(name)
  downloadPSet(as.character(name))
}