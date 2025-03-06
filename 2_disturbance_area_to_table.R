rm(list=ls())
gc()

rnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etel%C3%A4-Karjala","Keski-Suomi",
            "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-H%C3%A4me","Pohjanmaa","Varsinais-Suomi",
            "Etel%C3%A4-Pohjanmaa","P%C3%A4ij%C3%A4t-H%C3%A4me","Satakunta","Kymenlaakso",
            "Kainuu","Etel%C3%A4-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
regnames <- c("Uusimaa","Ahvenanmaa","Keski-Pohjanmaa","Pirkanmaa","Etela-Karjala","Keski-Suomi",
              "Pohjois-Savo","Lappi_E","Lappi_P","Kanta-Hame","Pohjanmaa","Varsinais-Suomi",
              "Etela-Pohjanmaa","Paijat-Hame","Satakunta","Kymenlaakso",
              "Kainuu","Etela-Savo","Pohjois-Karjala","Pohjois-Pohjanmaa")
rnos <- c(1:8,8:19)
rids <- c(1,3:length(rnames))
r_noi <- 1
load(file = paste0("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Results/validation_stats_rno1.rdata"))
damage_areas <- array(0,c(length(rids),nrow(out[[1]]),8),
                      dimnames = list(regnames[rids],rownames(out[[1]]),
                                      c("bb_stats","bb_sim_harv","bb_sim_all",
                                        "w_stats","w_sim_harv","w_sim_all",
                                        "bb_dam_sim_med/dam_stat_med",
                                        "w_dam_sim_med/dam_stat_med")))

for(ri in 1:length(rids)){
  r_noi <- rids[ri]
  r_no <- rnos[r_noi]
  print(paste0("Load validation_stats_rno",r_no,".rdata"))
  load(file = paste0("/scratch/project_2000994/PREBASruns/adaptFirst/Rsrc/Results/validation_stats_rno",r_no,".rdata"))
  damage_areas[ri,,1] <- out[[1]][,1]
  damage_areas[ri,,2] <- out[[1]][,2]
  damage_areas[ri,,3] <- out[[1]][,3]
  damage_areas[ri,,4] <- out[[2]][,1]
  damage_areas[ri,,5] <- out[[2]][,2]
  damage_areas[ri,,6] <- out[[2]][,3]
  damage_areas[ri,,7] <- out[[3]][,"pbb_median_decl/median_all"]
  damage_areas[ri,,8] <- out[[3]][,"pw_median_decl/median_all"]
}
