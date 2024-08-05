setwd("E:/youyi_fucha/first_test/gradient_voxel/Step_01_subject_data")
library(bruceR)
data_orgin <- import("all_first_subjects.xlsx")
data_bids <- import("bids_data_subjects.xlsx")

data_orgin %>% inner_join(data_bids) -> all_data

all_data %>% filter(mean_fd<=0.2,fd_higer_percent<0.25) -> motion_valid_data
export(motion_valid_data,"motion_valid.xlsx")



all_data %>% filter(mean_fd>0.2,fd_higer_percent>0.25)