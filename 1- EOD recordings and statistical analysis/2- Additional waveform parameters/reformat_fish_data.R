fish_data_path<-file.path("~/Desktop/mau_eod_analysis/input_data/fish_data.txt")
fish_data<-read.csv(fish_data_path,header=T,sep = "\t")


fish_data <- fish_data %>%
  mutate(Date = as.Date(Date, format = "%d-%b"))

# Calculate n_days
fish_data <- fish_data %>%
  group_by(Individual, Period) %>%
  mutate(n_days = as.numeric(Date - min(Date)))

# View the result
print(fish_data)

write.table(fish_data,"~/Desktop/mau_eod_analysis/input_data/fish_data.updated.txt",quote = F,sep = "\t")
