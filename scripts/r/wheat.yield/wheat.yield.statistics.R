yield1 = read.csv("rotation1scripts_v4/original_data/yield.data/FAOSTAT_data_9-25-2019.csv", stringsAsFactors = F)
yield2 = read.csv("rotation1scripts_v4/original_data/yield.data/FAOSTAT_data_9-25-2019 (1).csv", stringsAsFactors = F)
yield3 = read.csv("rotation1scripts_v4/original_data/yield.data/FAOSTAT_data_9-25-2019 (2).csv", stringsAsFactors = F)
plot(diff(yield1$Value))
library(ggplot2)

yield3 = yield3[which(yield3$Unit == "tonnes"), ]
ggplot(data = yield3, aes(x = Year, y = Value, group = Area, color = Area)) + geom_line() + theme(legend.position = "none")


