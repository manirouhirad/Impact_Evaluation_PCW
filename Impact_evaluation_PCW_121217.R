#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Description:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Created by Mani Rouhi Rad on 11/13/2017
# 
# This version includes early diff in diff results of PIEA land use changes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         To Do's:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Parallel trends assumption: Different combinations of treatment and control to check the parallel trends assumption.
# 2) Add land use for each year as control. Provide summary statistics for land use change. Blair is taking care of that.\
# 3) Check how these land use changes compare to a broader range of incentive based land use change programs
# 4) how much increase do we expect? we have hectares that have changed for each region
# 5) try the regressions with net forest loss
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Questions:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Important mentions:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data_reg = data_reg[year > 2004] line 276-ish. If we get better rainfall data, we can use data prior to 2005 as well.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########                                  Load Packages    ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rgeos)
library(rgdal)
library(plyr)
library(reshape2)
library(dplyr)
library(data.table)
library(ggplot2)
library(scales)
library(R.matlab)
library(sp)
library(sf)
library(maptools)
library(broom)
library(quantreg)
library(foreign)
library(raster)
library(plm)
library(lmtest)
library(multiwayvcov)
library(RColorBrewer)
library(readstata13)
library(parallel)
library(here)


# setwd("/Users/Mani/Dropbox/Panama S2/Impact Evaluation") # mac
# setwd("C:/Users/mr2284/Dropbox/Panama S2/Impact Evaluation") # windows

simpleCap <- function(x) {                                          # 1) Function defined to make all caps to propper
  paste(substring(x, 1,1), tolower(substring(x, 2)),                # 2)
        sep="")                       }                             # 3)

FN_area_sp = function(input_sp, scenario){                          # 1) function defined to find the area
  foo = st_as_sf(input_sp)
  a = as.data.table(input_sp@data)
  b = as.data.table(st_area(foo))
  b[, paste(scenario, "ha", sep="_"):=as.numeric(x*0.0001)]
  b=b[, 2]
  output_data=data.table(a, b)
  return(output_data)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########                                  Input data    ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   I. River discharge
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# river_discharge=data.table(read.csv("./Data/Map Area_update.csv"))
river_discharge=data.table(read.csv("./Data/Map Area2.csv"))
river_discharge=melt(river_discharge, id=c("River", "Year", "Day"), variable.name = "month", value.name = "discharge_m3_s")
river_discharge=river_discharge[complete.cases(Day)]
river_discharge[ , month:=as.character(month)]
river_discharge[ , month:=match(month,month.abb)]

# river_discharge[, foo_date:=as.character(paste(Year, month, Day, sep = "/"))]
# river_discharge[, date_read := as.Date(foo_date, "%Y/%b/%d")]
river_discharge=river_discharge[,.(River, year=Year, month, day=Day, discharge_m3_s)]
river_discharge[, discharge_m3_s:=as.numeric(discharge_m3_s)]

river_discharge[, River:=gsub("Cano", "Cano Quebrado", River)]
river_discharge[, River:=gsub("Ciri", "Ciri grande", River)]
river_discharge[, unique(River)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   II. River - watershed connection
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
river_names              = readOGR(dsn="./Data/Rivers",layer="Rivers_PCWCopy")
river_names=data.table(river_names@data)
river_names=river_names[,.(NOMBRE, SUBCUENC_2)]

river_names[, River:=simpleCap(NOMBRE)]
spanish_cano=river_names[, unique(River)][3]
river_names[, River:=gsub(spanish_cano, "Cano Quebrado", River)]
river_names[, unique(River)]

setkey(river_names,     River)
setkey(river_discharge, River)

river_discharge=river_names[river_discharge]
river_discharge=river_discharge[complete.cases(NOMBRE)]
river_discharge=river_discharge[,.(River, SUBCUENC_2, year, month, day, discharge_m3_s)]

weird_char = unique(river_discharge$SUBCUENC_2)[1]
weird_char = substr(weird_char, 1, 4) # win
# weird_char = substr(weird_char, 1, 3) # mac

river_discharge[, SUBCUENC_2 := gsub(weird_char, "Rio", SUBCUENC_2)]
river_discharge[, SUBCUENC_2 := gsub("BoquerÃ³n", "Boqueron", SUBCUENC_2)]
river_discharge[, SUBCUENC_2 := gsub("CaÃ±o", "Cano", SUBCUENC_2)]

weird_char_ciri=unique(river_discharge$SUBCUENC_2)[4]
weird_char_ciri=substr(weird_char_ciri, 5, 8)
weird_char_pequeni=unique(river_discharge$SUBCUENC_2)[6]
weird_char_pequeni=substr(weird_char_pequeni, 5, 13)
river_discharge[, SUBCUENC_2 := gsub(weird_char_ciri,    "Ciri",    SUBCUENC_2)] # win
river_discharge[, SUBCUENC_2 := gsub(weird_char_pequeni, "Pequeni", SUBCUENC_2)] # win
river_discharge[, SUBCUENC_2 := gsub("GatÃºn", "Gatun", SUBCUENC_2)]
river_discharge[, unique(SUBCUENC_2)]

river_discharge[year==2016 & River=="Gatun",.N, by="month"]
river_discharge[year==2016 & River=="Chagres",.N, by="month"]
river_discharge[,.N, by=River]

# Except Indio, which will be thrown out, the panel is almost balanced. Only Chagres has more data which may mean there is "extra" data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   II. Hectares in each region-year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PIEA                 = readOGR(dsn="./Data/PIEA",       layer="PIEA_polygons")
watershed_original   = readOGR(dsn="./Data/Watersheds", layer="Watersheds_PCW")
watershed_original   = gBuffer(watershed_original,    byid = TRUE, width = 0)
watershed_PIEA_final = raster::intersect(PIEA, watershed_original)
# writeOGR(watershed_PIEA_final, 
#          dsn="./Data/PIEA", layer="PIEA_by_watershed", driver="ESRI Shapefile")


watershed_PIEA_final=FN_area_sp(watershed_PIEA_final, "area_PIEA_watershed")
watershed_PIEA_final=watershed_PIEA_final[,.(Uso, SUBCUENC_2, area_PIEA_watershed_ha, area_watershed_ha=area_Water)]
watershed_PIEA_final[, Uso:=gsub(": Programa PIEA", "", Uso)]
watershed_PIEA_final_by_land_use = watershed_PIEA_final[,   sum(area_PIEA_watershed_ha),                          by="Uso"]
watershed_PIEA_final_by_watershed= watershed_PIEA_final[, .(sum(area_PIEA_watershed_ha), sum(area_watershed_ha)), by="SUBCUENC_2"]
watershed_PIEA_final_by_watershed[, percent_area:=V1/V2*100]
setkey(watershed_PIEA_final_by_watershed, percent_area)

watershed_PIEA_final[, SUBCUENC_2 := gsub(weird_char, "Rio", SUBCUENC_2)]
watershed_PIEA_final[, SUBCUENC_2 := gsub("BoquerÃ³n", "Boqueron", SUBCUENC_2)]
watershed_PIEA_final[, SUBCUENC_2 := gsub("CaÃ±o", "Cano", SUBCUENC_2)]
watershed_PIEA_final[, SUBCUENC_2 := gsub(weird_char_ciri, "Ciri", SUBCUENC_2)]
watershed_PIEA_final[, SUBCUENC_2 := gsub(weird_char_pequeni, "Pequeni", SUBCUENC_2)] # mac
watershed_PIEA_final[, SUBCUENC_2 := gsub("GatÃºn", "Gatun", SUBCUENC_2)]
watershed_PIEA_final[, unique(SUBCUENC_2)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   III. Rainfall (Other climatic variables?)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################
# rainfall_PIEA=read.dta13("./Data/Rainfall/acp_15Minstorms_Steve_v18112016.dta")
# rainfall_PIEA=data.table(rainfall_PIEA)
# rainfall_PIEA[,month:=substr(date, 4, 5)]
# rainfall_PIEA[,day  :=substr(date, 1, 2)]
# rainfall_PIEA[,total_rain_day:=sum(Totalrain), by=c("site", "Year", "month", "day")]
# rainfall_PIEA=unique(rainfall_PIEA, by=c("site", "Year", "month", "day"))
# rainfall_PIEA[,month:=as.numeric(month)]
# rainfall_PIEA[,Year :=as.numeric(Year)]
# rainfall_PIEA[,day  :=as.numeric(day)]

################
rainfall_PIEA2=read.csv("./Data/Rainfall/ACP_Rainfall_hourly.csv")
rainfall_PIEA2=data.table(rainfall_PIEA2)
rainfall_PIEA2[,DATETIME.dd.mm.yyyy.hh.mm.:=as.character(DATETIME.dd.mm.yyyy.hh.mm.)]
rainfall_PIEA2[, day:=substr(DATETIME.dd.mm.yyyy.hh.mm., 1, regexpr('/', DATETIME.dd.mm.yyyy.hh.mm.)-1)]
rainfall_PIEA2[,first_foo:=substr(DATETIME.dd.mm.yyyy.hh.mm., regexpr('/', DATETIME.dd.mm.yyyy.hh.mm.)+1, 20)]
rainfall_PIEA2[, month:=substr(first_foo, 1, regexpr('/', first_foo)-1)]
rainfall_PIEA2[,first_foo:=substr(first_foo, regexpr('/', first_foo)+1, 20)]
rainfall_PIEA2[, year:=substr(first_foo, 1, regexpr(' ', first_foo)-1)]

rainfall_PIEA2[, day:=  as.numeric(as.character(day))]
rainfall_PIEA2[, month:=as.numeric(as.character(month))]
rainfall_PIEA2[, year:= as.numeric(as.character(year))]
rainfall_PIEA2=rainfall_PIEA2[, .(Station, year, month, day, total_rain_day=Rain..mm.)]
rainfall_PIEA2[, total_rain_day:=sum(total_rain_day), by=c("Station", "year", "month", "day")]
rainfall_PIEA2=unique(rainfall_PIEA2, by=c("Station", "year", "month", "day"))

####### check these two for rainfall before 2004
rainfall_PIEA2[,unique(Station)]
rainfall_PIEA[, unique(site)]

rainfall_PIEA2[Station== "CHORRO", .N, by="year"]
rainfall_PIEA[ site == "CHORRO"  , .N, by="Year"]


# for (i in 1:nrow(rivers_used)) {
#   temp=rainfall_PIEA2[Station== rivers_used[i], .N, by="year"]
#   print(rivers_used[i])
#   print(temp)
#   readline(prompt="Press [enter] to continue")
# }
######## only SANMIGUEL and CHAGRECITO don't have enough data. Unfortunately, ZANGUENA only has data starting 2004. This is in Cano Quebrado catchment so would have been helpful.
setnames(rainfall_PIEA2, "Station", "site")
###############


rainfall_stations = readOGR(dsn="./Data/Rainfall/stations", layer="stations")
# foo1 = watershed_original[!is.na(over(rainfall_stations, as(watershed_original,"SpatialPolygons"))),]

bar_1=over(rainfall_stations, as(watershed_original,"SpatialPolygons"))
bar_1=data.table(bar_1)
setnames(bar_1, "bar_1", "over_id")
bar_2=data.table(bar_1, rainfall_stations@data)
bar_2=bar_2[complete.cases(over_id)]
bar_3=data.table(watershed_original@data)

setkey(bar_2, over_id)
setkey(bar_3, OBJECTID_1)
bar= bar_2[bar_3]
bar=bar[complete.cases(CODENAME)]
setkey(bar, CODENAME)
rainfall_stations=bar[,.(OBJECTID, over_id, CODENAME, SUBCUENC_2)]
rainfall_stations[, CODENAME:=gsub(" ", "", CODENAME, fixed = TRUE)]
setkey(rainfall_stations, CODENAME)
setkey(rainfall_PIEA2, site)
rainfall_PIEA2=rainfall_stations[rainfall_PIEA2]
rainfall_PIEA2=rainfall_PIEA2[complete.cases(CODENAME) & complete.cases(SUBCUENC_2)]
rainfall_PIEA2=rainfall_PIEA2[,.(Station= CODENAME, SUBCUENC_2, date=as.Date(paste(year, month, day, sep = "-")), year, month, day, total_rain_day)]
# rainfall_PIEA2[, date_read:=as.Date(date, "dd/mm/yyyy")]
# rainfall_PIEA2[, date_read:=as.Date(date, "%d/%m/%Y")]

# 
# ggplot(data = rainfall_PIEA2[Station=="CANONES"], aes(x=date_read, y=total_rain_day))+
#   geom_line(size=1.1)+
#   # coord_cartesian(ylim=c(0, 400)) +
#   scale_x_date(labels = date_format("%m-%Y")) #+ 
# # scale_y_discrete(breaks=seq(0, 400, 100))  # Ticks from 0-10, every .25
# rainfall_PIEA2[, unique(Station)]

rainfall_PIEA2[,        unique(SUBCUENC_2)]
# watershed_PIEA_final[, unique(SUBCUENC_2)]
# river_discharge[, unique(SUBCUENC_2)]
# str(river_discharge)
# str(rainfall_PIEA2)

rainfall_PIEA2[SUBCUENC_2== rainfall_PIEA2[,        unique(SUBCUENC_2)][6], unique(Station)]
rainfall_PIEA2[Station   == rainfall_PIEA2[SUBCUENC_2==rainfall_PIEA2[,        unique(SUBCUENC_2)][6], unique(Station)]]
# rainfall_PIEA[Station  == rainfall_PIEA2[SUBCUENC_2==rainfall_PIEA2[,        unique(SUBCUENC_2)][6], unique(Station)]]
### looks good \/

rainfall_PIEA2[, SUBCUENC_2 := gsub(weird_char, "Rio", SUBCUENC_2)]
rainfall_PIEA2[, SUBCUENC_2 := gsub("BoquerÃ³n", "Boqueron", SUBCUENC_2)]
rainfall_PIEA2[, SUBCUENC_2 := gsub("CaÃ±o", "Cano", SUBCUENC_2)]
rainfall_PIEA2[, SUBCUENC_2 := gsub(weird_char_ciri, "Ciri", SUBCUENC_2)] # win
rainfall_PIEA2[, SUBCUENC_2 := gsub(weird_char_pequeni, "Pequeni", SUBCUENC_2)] # win
rainfall_PIEA2[, SUBCUENC_2 := gsub("GatÃºn", "Gatun", SUBCUENC_2)]
rainfall_PIEA2[, unique(SUBCUENC_2)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   IV. Seasons
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seasons=read.csv("./Data/seasons.csv")
seasons=data.table(seasons)
seasons[, Dry_season_start_date := as.Date(Dry_season_start, "%d-%b-%y")]
seasons[, Dry_season_end_date   := as.Date(Dry_season_end,   "%d-%b-%y")]
seasons=seasons[,.(YEAR, Dry_season_start_date, Dry_season_end_date)]
seasons[, Dry_season_length := Dry_season_end_date - Dry_season_start_date]

setkey(seasons, YEAR)
setkey(river_discharge, year)


river_discharge=river_discharge[seasons]
river_discharge=river_discharge[complete.cases(River)]
river_discharge[, date_disch:=as.Date(paste(year, month, day, sep = "/"))]
river_discharge=river_discharge[complete.cases(date_disch)]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   V. Merge
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setkey(rainfall_PIEA2,  SUBCUENC_2, year, month, day)
setkey(river_discharge, SUBCUENC_2, year, month, day)

unique(rainfall_PIEA2,  by="SUBCUENC_2")
unique(river_discharge, by="SUBCUENC_2")

# setkey(watershed_PIEA_final, SUBCUENC_2, year, month, day)
# setkey(watershed_PIEA_final, SUBCUENC_2, Uso)

# data_reg = river_discharge[rainfall_PIEA2, allow.cartesian=T, all.x=T]
# data_reg = rainfall_PIEA2[river_discharge, allow.cartesian=T]
# data_reg = data_reg[complete.cases(date_read)]
# data_reg = data_reg[complete.cases(total_rain_day) & complete.cases(discharge_m3_s)]

data_reg = merge(river_discharge, rainfall_PIEA2, by.x=c("SUBCUENC_2", "year", "month", "day"), by.y=c("SUBCUENC_2", "year", "month", "day") , all.x=T)
# data_reg = data_reg[year > 2004]
data_reg = data_reg[year > 2003]
data_reg = data_reg[,.(SUBCUENC_2, River, year, month, day, Dry_season_start_date, Dry_season_end_date, discharge_m3_s, Station, total_rain_day)]
data_reg[, date_read := as.Date(paste(year, month, day, sep = "-"))]
data_reg=data_reg[complete.cases(date_read)]
data_reg[is.na(total_rain_day), total_rain_day:=0]

# rivers_used=data_reg[, unique(Station), by="SUBCUENC_2"]
# rivers_used=rivers_used[,2]
data_reg[, day_of_dry_season := date_read - Dry_season_start_date]
data_reg[, day_of_wet_season := date_read - (Dry_season_end_date+1)]
data_reg[, day_of_season     := ifelse(date_read < Dry_season_end_date+1, day_of_dry_season, day_of_wet_season)]
data_reg[, dry_dummy         := ifelse(date_read < Dry_season_end_date+1, 1, 0)]
data_reg=data_reg[year<2016]

######## I think merging discharge and rainfall should be done more carefully




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   VI. land use data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are the areas that Blair estimated. I double chekced Trinidad for 2012 and the estimates look good. So I am using the estimates.
# 
land_use=read.csv("./Data/Land Covers/land_use_over_time_R_input.csv")
land_use=data.table(land_use)
land_use[, unique(River)]
land_use=land_use[River!="Las Cascadas"]

# col_land_use=colnames(land_use)[3:5]
# land_use[, (col_land_use):=lapply(.SD, sum), by="River", .SDcols=col_land_use]
# land_use=unique(land_use, by="River")

# land_use=land_use[River=="Trinidad" | River== "Ciri grande" | River=="Cano Quebrado"]
land_use=land_use[,foo:=ifelse(X2003==0 & X2008==0 & X2012==0, 1, 0)]
land_use=land_use[foo==0]
land_use=land_use[Land.use!="Water"]

land_use[River=="Trinidad"]


forest=land_use[Land.use=="Mature Forests" | Land.use=="Secondary Forests" | Land.use=="Forest Plantations"]
cols=colnames(forest)[3:5]
forest[, (cols):=lapply(.SD, sum), by="River", .SDcols=cols] 
forest=unique(forest, by="River")
forest[, forest_2012_2008:=X2012 - X2008]
forest[, forest_2008_2003:=X2008 - X2003]
forest=forest[,.(River, forest_2012_2008, forest_2008_2003)]

pasture=land_use[Land.use=="Pastures"]
pasture_fig=copy(pasture)
pasture[, pasture_2012_2008:=X2012 - X2008]
pasture[, pasture_2008_2003:=X2008 - X2003]
pasture=pasture[,.(River, pasture_2012_2008, pasture_2008_2003)]

convertible=land_use[Land.use=="Pastures" | Land.use=="Crops" | Land.use=="Scrub" | Land.use=="Kans Grass"]
cols=colnames(convertible)[3:5]
convertible[, (cols):=lapply(.SD, sum), by="River", .SDcols=cols] 
convertible=unique(convertible, by="River")
convertible[, convertible_2012_2008:=X2012 - X2008]
convertible[, convertible_2008_2003:=X2008 - X2003]
convertible=convertible[,.(River, convertible_2012_2008, convertible_2008_2003)]

urban=land_use[Land.use=="Towns"]
urban[, urban_2012_2008:=X2012 - X2008]
urban[, urban_2008_2003:=X2008 - X2003]
urban=urban[,.(River, urban_2012_2008, urban_2008_2003)]


non_forest=land_use[Land.use=="Pastures" | Land.use=="Crops" | Land.use=="Scrub" | Land.use=="Kans Grass" | Land.use=="Towns" | Land.use=="Barren Soil"]
cols=colnames(non_forest)[3:5]
non_forest[, (cols):=lapply(.SD, sum), by="River", .SDcols=cols] 
non_forest=unique(non_forest, by="River")
non_forest[, non_forest_2012_2008:=X2012 - X2008]
non_forest[, non_forest_2008_2003:=X2008 - X2003]
non_forest=non_forest[,.(River, non_forest_2012_2008, non_forest_2008_2003)]


setkey(forest, River)
setkey(pasture, River)
setkey(convertible, River)
setkey(non_forest, River)

land_use_2008_2012=forest[pasture]
land_use_2008_2012=land_use_2008_2012[convertible]
land_use_2008_2012=land_use_2008_2012[non_forest]
land_use_2008_2012[, net_2008_2012:=forest_2012_2008 - convertible_2012_2008]

# land_use_fig[, type:=ifelse(Land.use=="Mature Forests" | Land.use=="Secondary Forests" | Land.use=="Forest Plantations", "forest", 
#                             ifelse(Land.use=="Pastures" | Land.use=="Crops" | Land.use=="Scrub" | Land.use=="Kans Grass", "convertible", 
#                                    ifelse(Land.use=="Towns", "towns", "other")))]




forest=land_use[Land.use=="Mature Forests" | Land.use=="Secondary Forests" | Land.use=="Forest Plantations"]
cols=colnames(forest)[3:5]
forest[, (cols):=lapply(.SD, sum), by="River", .SDcols=cols] 
forest=unique(forest, by="River")
forest=forest[,.(River, f_2003=X2003, f_2008=X2008, f_2012=X2012)]


non_forest=land_use[Land.use=="Pastures" | Land.use=="Crops" | Land.use=="Scrub" | Land.use=="Kans Grass" | Land.use=="Towns" | Land.use=="Barren Soil"]
cols=colnames(non_forest)[3:5]
non_forest[, (cols):=lapply(.SD, sum), by="River", .SDcols=cols] 
non_forest=unique(non_forest, by="River")
non_forest=non_forest[,.(River, n_2003=X2003, n_2008=X2008, n_2012=X2012)]

setkey(forest, River)
setkey(non_forest, River)
net=forest[non_forest]
net[, X2003:=f_2003 - n_2003]
net[, X2008:=f_2008 - n_2008]
net[, X2012:=f_2012 - n_2012]
net=net[,.(River, type="net", X2003, X2008, X2012)]

land_use_fig = copy(land_use)
land_use_fig[, type:=ifelse(Land.use=="Mature Forests" | Land.use=="Secondary Forests" | Land.use=="Forest Plantations", "forest", 
                            ifelse(Land.use=="Pastures" | Land.use=="Crops" | Land.use=="Scrub" | Land.use=="Kans Grass" | Land.use=="Towns" | Land.use=="Barren Soil", "non_forest" 
                              , "other"))]

land_use_fig=land_use_fig[type=="forest" | type== "non_forest"]

col_sum=colnames(land_use_fig)[3:5]
land_use_fig[,(col_sum):=lapply(.SD, sum), by=c("River", "type"), .SDcols=col_sum]
land_use_fig=unique(land_use_fig, by=c("type", "River"))

land_use_fig=land_use_fig[,.(River, type, X2003, X2008, X2012)]
land_use_fig=rbind(land_use_fig, net)

setkey(land_use_fig, River)
land_use_fig=land_use_fig[,.(River, type, X2003, X2008, X2012)]
land_use_fig=melt(land_use_fig, id=c("River", "type"))
land_use_fig[, variable:=ifelse(variable=="X2003", 2003, ifelse(variable=="X2008", 2008, 2012))]
foo_d=land_use_fig[variable==2003, value]
land_use_fig[, foo:=foo_d]
land_use_fig[, foo:=value - foo]

ggplot(land_use_fig[type=="forest"], aes(x=factor(variable), y=foo, group=River, col=River))+
  geom_point(size=2)+
  geom_line(size=2)+
  xlab("Year")+
  ylab("Change in forest cover from 2003 (ha)")

###### >>>>> I need to get back to this and add them to the figure with ESE estimates as bars <<<<<<<<<

library(gridExtra)
library(grid)

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
  # heights=unit(c(2, 2, 2), c("in", "in", "in")))
  
}


library(directlabels)


lu_fig1<- ggplot(land_use_fig[type=="forest"], aes(x=factor(variable), y=foo, group=River, col=River))+
  geom_point(size=2)+
  geom_line(size=1.2)+
  theme_bw()+
  xlab("Year")+
  ylab("Ha")+
  ggtitle(paste("forest land"))+
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = River), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8))


lu_fig2<- ggplot(land_use_fig[type=="net"], aes(x=factor(variable), y=foo, group=River, col=River))+
  geom_point(size=2)+
  geom_line(size=1.2)+
  theme_bw()+
  xlab("Year")+
  ylab("Ha")+
  ggtitle(paste("net forest cover"))+
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = River), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8))



lu_fig3 <- ggplot(land_use_fig[type=="non_forest"], aes(x=factor(variable), y=foo, group=River, col=River))+
  geom_point(size=2)+
  geom_line(size=1.2)+
  theme_bw()+
  xlab("Year")+
  ylab("Ha")+
  ggtitle(paste("non forest land"))+
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = River), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8))





multiplot(lu_fig1, lu_fig2, lu_fig3, cols = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
########                                  Data Analysis           ###########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trin     = watershed_PIEA_final_by_watershed[,SUBCUENC_2][25]
ciri     = watershed_PIEA_final_by_watershed[,SUBCUENC_2][31]
cano     = watershed_PIEA_final_by_watershed[,SUBCUENC_2][3]
Boqueron = watershed_PIEA_final_by_watershed[,SUBCUENC_2][4]
Gatun    = watershed_PIEA_final_by_watershed[,SUBCUENC_2][8]

watershed_PIEA_west = watershed_PIEA_final_by_watershed[SUBCUENC_2==trin | SUBCUENC_2==ciri | SUBCUENC_2==cano | SUBCUENC_2==Boqueron | SUBCUENC_2==Gatun]

watershed_PIEA_west = watershed_PIEA_west[ SUBCUENC_2==trin,     SUBCUENC_2:="Trinidad"]
watershed_PIEA_west = watershed_PIEA_west[ SUBCUENC_2==ciri,     SUBCUENC_2:="Ciri grande"]
watershed_PIEA_west = watershed_PIEA_west[ SUBCUENC_2==cano,     SUBCUENC_2:="Cano Quebrado"]
watershed_PIEA_west = watershed_PIEA_west[ SUBCUENC_2==Gatun,    SUBCUENC_2:="Gatun"]
watershed_PIEA_west = watershed_PIEA_west[ SUBCUENC_2==Boqueron, SUBCUENC_2:="Boqueron"]

river_discharge_west= river_discharge[River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado" | River=="Boqueron" | River=="Gatun"]
river_discharge_west=river_discharge_west[date_disch < Dry_season_end_date]
river_discharge_west[, tot_discharge := discharge_m3_s*24*3600]  # M3

setkey(watershed_PIEA_west, SUBCUENC_2)
setkey(river_discharge_west, River)

watershed_PIEA_west=watershed_PIEA_west[river_discharge_west]

summary_PIEA=watershed_PIEA_west[,.(Watershed=SUBCUENC_2, PIEA_hecatres=V1, Total_hectares=V2, percent_area_in_PIEA=percent_area)]
summary_PIEA=unique(summary_PIEA, by="Watershed")
summary_PIEA=rbind(summary_PIEA, list("Chagres", 0, 0, 0), list("Pequeni", 0, 0, 0))

data_pre=watershed_PIEA_west[year < 2008]

data_pre[, discharge_per_ha := tot_discharge/V2] # m3 per day per hectare
data_pre[, extra_discharge  := discharge_per_ha * V1 * .2] # m3 per day for the entire catchment under 20% ESE
data_pre[, mean_extra_discharge :=mean(extra_discharge), by=c("SUBCUENC_2", "month")]
data_pre[, sd_extra_discharge   :=sd(extra_discharge), by=c("SUBCUENC_2", "month")]
data_pre=unique(data_pre, by=c("SUBCUENC_2", "month"))

data_pre=data_pre[,.(SUBCUENC_2, month, mean_extra_discharge, sd_extra_discharge)]

# data_pre[, percent_extra_discharge := extra_discharge/(tot_discharge * 10^6)*100]
# watershed_PIEA_final_by_watershed[, water_increase_treatment_ESE_15:=V1*.15]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
########                                  Regressions           ###########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   0. Teeatment before and after ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#################

data_reg_copy_daily = copy(data_reg)
data_reg_copy_daily[,total_rain_day:=mean(total_rain_day), by=c("River", "year", "month", "day")]
data_reg_copy_daily=unique(data_reg_copy_daily, by=c("River", "year", "month", "day"))
data_reg_copy_daily[, tot_discharge  := discharge_m3_s*24*3600] # total daily Million m3
data_reg_copy_daily[, lag.rain_day   := c(NA, total_rain_day[-.N]), by=c("SUBCUENC_2", "year")]
data_reg_copy_daily[, lag.rain_day_2 := c(NA, lag.rain_day[-.N]),   by=c("SUBCUENC_2", "year")]
data_reg_copy_daily[, rain_past_2_days := lag.rain_day + lag.rain_day_2]
data_reg_copy_daily[, lag.discharge  := c(NA, tot_discharge[-.N]),  by=c("SUBCUENC_2", "year")]
data_reg_daily_dry = data_reg_copy_daily[dry_dummy==1]
data_reg_daily_wet = data_reg_copy_daily[dry_dummy==0]

# sum_stats=data_reg_copy_daily[,.(River, )]


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# 
# FN_before_after_by_river     <- function(year_pre, year_post, River_tr) {
#   # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
#   dt=data_reg_copy_daily[(River==River_tr)]
#   dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
#   coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
#   for(i in 1:12){
#     DiD_prepost = lm(tot_discharge ~ year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
#                      data = dt[ month == i])
#     beta_DiD = summary(DiD_prepost)$coefficients[2,1]
#     std_error= summary(DiD_prepost)$coefficients[2,2]
#     coefs_DiD[i, 2] = beta_DiD
#     coefs_DiD[i, 3] = std_error
#   }
#   dt_ESE=data_pre[SUBCUENC_2==River_tr]
#   
#   l=ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
#     geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
#     geom_line() +
#     geom_point() +
#     geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge), 
#                   width=.2, col="red") +
#     xlab("month")+
#     ylab("Diff in Diff Coefficient")+
#     ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
#   return(l)
# }
# 




FN_before_after_by_river     <- function(year_pre, year_post, River_tr) {
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=data_reg_copy_daily[River==River_tr & month < 5]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  dt_pre=dt[year_treat==0 ,.(year, month, day, tot_discharge)]
  dt_pos=dt[year_treat==1 ,.(year, month, day, tot_discharge)]
  
  dt_pre[,mean_disch_pre:=mean(tot_discharge), by="month"]
  dt_pre[,foo:=.N, by="month"]
  dt_pre[,sd_disch_pre  :=  sd(tot_discharge), by= "month"]
  dt_pre[,sd_disch_pre  :=  (sd_disch_pre)^2/foo]
  dt_pre=unique(dt_pre, by="month")
  dt_pre=dt_pre[,.(month, mean_disch_pre, sd_disch_pre)]
  
  dt_pos[,mean_disch_pos:=mean(tot_discharge), by="month"]
  dt_pos[,foo:=.N, by="month"]
  dt_pos[,sd_disch_pos  :=  sd(tot_discharge), by= "month"]
  dt_pos[,sd_disch_pos  :=  (sd_disch_pos)^2/foo]
  dt_pos=unique(dt_pos, by="month")
  dt_pos=dt_pos[,.(month, mean_disch_pos, sd_disch_pos)]
  
  setkey(dt_pre, month)
  setkey(dt_pos, month)
  dt=dt_pre[dt_pos]
  dt[, sd:=sqrt(sd_disch_pre + sd_disch_pos)]
  dt[, mean:=mean_disch_pos - mean_disch_pre]
  
  dt_ESE=data_pre[SUBCUENC_2==River_tr & month < 5]
  
  # l=ggplot(data = dt, aes(x=factor(month), y=mean))+
  #   geom_errorbar(aes(ymin=mean - 1.96 * sd , ymax=mean + 1.96 * sd), width=.1) +
  #   geom_line(aes(group=1)) +
  #   geom_point() +
  #   geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge), 
  #                 width=.2, col="red") +
  #   xlab("month")+
  #   ylab("Diff in Diff Coefficient")+
  #   ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
  l=ggplot(data = dt_ESE, aes(x=factor(month)))+
    geom_errorbar(aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge), 
                  width=.2, col="red") +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
  
    return(l)
}


# dt=data_reg_copy_daily[River=="Trinidad" & month < 5]
# dt = dt[year < 2008 | year > 2012][ ,year_treat:= ifelse(year>2007, 1, 0)]


## West Side
a1=FN_before_after_by_river(2008, 2012, "Trinidad")
a3=FN_before_after_by_river(2008, 2012, "Ciri grande")
a2=FN_before_after_by_river(2008, 2012, "Cano Quebrado")

multiplot(a1, a2, a3, cols=2)


## East Side
FN_before_after_by_river(2008, 2012, "Boqueron")
FN_before_after_by_river(2008, 2012, "Chagres")
FN_before_after_by_river(2008, 2012, "Gatun")
FN_before_after_by_river(2008, 2012, "Pequeni")


FN_coefficients_DiD_west(2008, 2012)
FN_coefficients_DiD_west(2008, 2013)
FN_coefficients_DiD_west(2007, 2013)
FN_coefficients_DiD_west(2007, 2012)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   0. diff-in-diff with daily data for the entire year ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############

FN_coefficients_DiD_monthly_nobs     <- function(year_pre, year_post) {
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=copy(data_reg_copy_daily)
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                     + factor(year), data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),1]
    std_error= summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  coefs_DiD[, CI_lower:=coefficient_DiD - 2.866 * standard_error_DiD]
  coefs_DiD[, CI_upper:=coefficient_DiD + 2.866 * standard_error_DiD]
  nnss=nobs(DiD_prepost)
  
  coefs_DiD[, m:= ceiling(((2.866 * standard_error_DiD) * sqrt(nnss-2)/coefficient_DiD)^2) + 2]
  coefs_DiD[, CI_lower_corrected:=coefficient_DiD - 2.866 * standard_error_DiD  * sqrt(nnss-2)/sqrt(m-2)]
  coefs_DiD[, CI_upper_corrected:=coefficient_DiD + 2.866 * standard_error_DiD  * sqrt(nnss-2)/sqrt(m-2)]
  coefs_DiD[month > 4, c("CI_lower_corrected", "CI_upper_corrected"):=0]
  
  l=ggplot(data = coefs_DiD[month<5], aes(x=factor(month), y=coefficient_DiD, group=1))+
    geom_errorbar(aes(ymin=CI_lower_corrected , ymax=CI_upper_corrected), width=.1) +
    geom_line() +
    geom_point() +
    geom_text(aes(label=m),vjust=-20)+
    theme_bw()+
    # geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge),
    #               width=.2, col="red") +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
  return(l)
}

FN_coefficients_DiD_monthly_nobs(2008, 2012)




dt=copy(data_reg_copy_daily)
dt = dt[year < 2008 | year > 2012][ ,year_treat:= ifelse(year>2007, 1, 0)]
dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
for(i in 1:12){
  DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                   + factor(year), data = dt[ month == i])
  beta_DiD = summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),1]
  std_error= summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),2]
  coefs_DiD[i, 2] = beta_DiD
  coefs_DiD[i, 3] = std_error
}
coefs_DiD[, CI_lower:=coefficient_DiD - 2.866 * standard_error_DiD]
coefs_DiD[, CI_upper:=coefficient_DiD + 2.866 * standard_error_DiD]
nnss=nobs(DiD_prepost)

coefs_DiD[, m:= ceiling(((2.866 * standard_error_DiD) * sqrt(nnss-2)/coefficient_DiD)^2) + 2]
coefs_DiD[, CI_lower_corrected:=coefficient_DiD - 2.866 * standard_error_DiD  * sqrt(nnss-2)/sqrt(m-2)]
coefs_DiD[, CI_upper_corrected:=coefficient_DiD + 2.866 * standard_error_DiD  * sqrt(nnss-2)/sqrt(m-2)]
coefs_DiD[month > 4, c("CI_lower_corrected", "CI_upper_corrected"):=0]

l=ggplot(data = coefs_DiD[month<5], aes(x=factor(month), y=coefficient_DiD, group=1))+
  geom_errorbar(aes(ymin=CI_lower_corrected , ymax=CI_upper_corrected), width=.1) +
  geom_line() +
  geom_point() +
  geom_text(aes(label=m),vjust=-20)+
  theme_bw()+
  # geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge),
  #               width=.2, col="red") +
  xlab("month")+
  ylab("Diff in Diff Coefficient")


+
  ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))












FN_coefficients_DiD_annual     <- function(year_pre, year_post) {
  # dt=data_reg_copy_daily
  dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                   + factor(year) + factor(dry_dummy), 
                   data = dt)
  a=summary(DiD_prepost)
  return(a)
}

FN_coefficients_DiD_annual_dry <- function(year_pre, year_post) {
  # dt=data_reg_daily_dry
  dt=data_reg_daily_dry[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                   + factor(year), 
                   data = dt)
  a=summary(DiD_prepost)
  return(a)
}

FN_coefficients_DiD_annual_wet <- function(year_pre, year_post) {
  # dt=data_reg_daily_dry
  dt=data_reg_daily_wet[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                   + factor(year), 
                   data = dt)
  a=summary(DiD_prepost)
  return(a)
}




FN_coefficients_DiD_annual(2008, 2012) 
FN_coefficients_DiD_annual_dry(2008, 2012)
FN_coefficients_DiD_annual_wet(2008, 2012)

FN_coefficients_DiD_annual(2008,     2013) 
FN_coefficients_DiD_annual_dry(2008, 2013)
FN_coefficients_DiD_annual_wet(2008, 2013)

# Total discharge decreases as a result of treatment. This is consistent with the SE hypothesis. However, both rainy season and dry season discharge are
# also decreasing.



### year by season fixed effects

data_reg_copy_daily[, year_by_season := .GRP, by=c("year", "dry_dummy")]

FN_coefficients_DiD_annual     <- function(year_pre, year_post) {
  # dt=data_reg_copy_daily
  dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                   + factor(year_by_season), 
                   data = dt)
  a=summary(DiD_prepost)
  return(a)
}

FN_coefficients_DiD_annual(2008, 2012) 
FN_coefficients_DiD_annual(2008, 2013) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   I. diff-in-diff with daily data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   a. with lagged discharge ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####


FN_coefficients_DiD_monthly     <- function(year_pre, year_post) {
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=copy(data_reg_copy_daily)
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    # DiD_prepost = lm(tot_discharge ~ year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge,
    #                  data = dt[ month == i])
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                     + factor(year)
                     , data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[12,1]
    std_error= summary(DiD_prepost)$coefficients[12,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  # dt_ESE=data_pre[SUBCUENC_2==River_tr]
  
  l=ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD, group=1))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    # geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge),
    #               width=.2, col="red") +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
  return(l)
}

FN_coefficients_DiD_monthly(2008, 2012)


FN_coefficients_DiD_west     <- function(year_pre, year_post) {
  # dt=data_reg_daily_dry[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    # geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, "with lagged discharge", sep = " "))
}

FN_coefficients_DiD_west_dry <- function(year_pre, year_post) {
  dt=data_reg_daily_dry[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:4, coefficient_DiD = rep(0 , 4), standard_error_DiD= rep(0, 4))
  for(i in 1:4){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Dry season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_west_wet <- function(year_pre, year_post) {
  dt=data_reg_daily_wet[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:8, coefficient_DiD = rep(0 , 8), standard_error_DiD= rep(0, 8))
  for(i in 1:8){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == 4+i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Rainy season - Before", year_pre, "and after", year_post, sep = " "))
}


FN_coefficients_DiD_west(2008, 2012)
FN_coefficients_DiD_west(2008, 2013)
FN_coefficients_DiD_west(2007, 2013)
FN_coefficients_DiD_west(2007, 2012)

# FN_coefficients_DiD_west_dry(2008, 2012)
# FN_coefficients_DiD_west_dry(2008, 2013)
# FN_coefficients_DiD_west_dry(2007, 2013)
# FN_coefficients_DiD_west_dry(2007, 2012)
# 
# FN_coefficients_DiD_west_wet(2008, 2012)
# FN_coefficients_DiD_west_wet(2008, 2013)
# FN_coefficients_DiD_west_wet(2007, 2013)
# FN_coefficients_DiD_west_wet(2007, 2012)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#                   b. without lagged discharge ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

FN_coefficients_DiD_west     <- function(year_pre, year_post) {
  dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[7,1]
    std_error= summary(DiD_prepost)$coefficients[7,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, "without lagged discharge", sep = " "))
}

FN_coefficients_DiD_west_dry <- function(year_pre, year_post) {
  dt=data_reg_daily_dry[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:4, coefficient_DiD = rep(0 , 4), standard_error_DiD= rep(0, 4))
  for(i in 1:4){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[7,1]
    std_error= summary(DiD_prepost)$coefficients[7,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Dry season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_west_wet <- function(year_pre, year_post) {
  dt=data_reg_daily_wet[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:8, coefficient_DiD = rep(0 , 8), standard_error_DiD= rep(0, 8))
  for(i in 1:8){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == 4+i])
    beta_DiD = summary(DiD_prepost)$coefficients[7,1]
    std_error= summary(DiD_prepost)$coefficients[7,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Rainy season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_west(2008, 2012)
FN_coefficients_DiD_west(2008, 2013)
FN_coefficients_DiD_west(2007, 2013)
FN_coefficients_DiD_west(2007, 2012)

# FN_coefficients_DiD_west_dry(2008, 2012)
# FN_coefficients_DiD_west_dry(2008, 2013)
# FN_coefficients_DiD_west_dry(2007, 2013)
# FN_coefficients_DiD_west_dry(2007, 2012)
# 
# FN_coefficients_DiD_west_wet(2008, 2012)
# FN_coefficients_DiD_west_wet(2008, 2013)
# FN_coefficients_DiD_west_wet(2007, 2013)
# FN_coefficients_DiD_west_wet(2007, 2012)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#                   c. without lagged discharge + squared day of season ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

FN_coefficients_DiD_west     <- function(year_pre, year_post) {
  dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + I(day_of_season^2), 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, "with squared dry season", sep = " "))
}

FN_coefficients_DiD_west_dry <- function(year_pre, year_post) {
  dt=data_reg_daily_dry[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:4, coefficient_DiD = rep(0 , 4), standard_error_DiD= rep(0, 4))
  for(i in 1:4){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + I(day_of_season^2), 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Dry season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_west_wet <- function(year_pre, year_post) {
  dt=data_reg_daily_wet[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:8, coefficient_DiD = rep(0 , 8), standard_error_DiD= rep(0, 8))
  for(i in 1:8){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + I(day_of_season^2), 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Rainy season - Before", year_pre, "and after", year_post, sep = " "))
}


FN_coefficients_DiD_west(2008, 2012)
FN_coefficients_DiD_west(2008, 2013)
FN_coefficients_DiD_west(2007, 2013)
FN_coefficients_DiD_west(2007, 2012)

# FN_coefficients_DiD_west_dry(2008, 2012)
# FN_coefficients_DiD_west_dry(2008, 2013)
# FN_coefficients_DiD_west_dry(2007, 2013)
# FN_coefficients_DiD_west_dry(2007, 2012)
# 
# FN_coefficients_DiD_west_wet(2008, 2012)
# FN_coefficients_DiD_west_wet(2008, 2013)
# FN_coefficients_DiD_west_wet(2007, 2013)
# FN_coefficients_DiD_west_wet(2007, 2012)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#                   d. Change treatment to 2008
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

FN_coefficients_DiD_west     <- function(year_pre, year_post) {
  dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2008, 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, "with treatment at 2008", sep = " "))
}

FN_coefficients_DiD_west_dry <- function(year_pre, year_post) {
  dt=data_reg_daily_dry[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2008, 1, 0)]
  coefs_DiD = data.table(month=1:4, coefficient_DiD = rep(0 , 4), standard_error_DiD= rep(0, 4))
  for(i in 1:4){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Dry season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_west_wet <- function(year_pre, year_post) {
  dt=data_reg_daily_wet[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2008, 1, 0)]
  coefs_DiD = data.table(month=1:8, coefficient_DiD = rep(0 , 8), standard_error_DiD= rep(0, 8))
  for(i in 1:8){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == 4+i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Rainy season - Before", year_pre, "and after", year_post, sep = " "))
}


FN_coefficients_DiD_west(2008, 2012)
FN_coefficients_DiD_west(2008, 2013)
FN_coefficients_DiD_west(2009, 2013)
FN_coefficients_DiD_west(2009, 2012)

# FN_coefficients_DiD_west_dry(2008, 2012)
# FN_coefficients_DiD_west_dry(2008, 2013)
# FN_coefficients_DiD_west_dry(2007, 2013)
# FN_coefficients_DiD_west_dry(2007, 2012)
# 
# FN_coefficients_DiD_west_wet(2008, 2012)
# FN_coefficients_DiD_west_wet(2008, 2013)
# FN_coefficients_DiD_west_wet(2007, 2013)
# FN_coefficients_DiD_west_wet(2007, 2012)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   e. Entire Watershed ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

FN_coefficients_DiD_all     <- function(year_pre, year_post) {
  dt=data_reg_copy_daily
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[7,1]
    std_error= summary(DiD_prepost)$coefficients[7,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post,"for entire watershed" ,sep = " "))
}

FN_coefficients_DiD_all_dry <- function(year_pre, year_post) {
  dt=data_reg_daily_dry
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:4, coefficient_DiD = rep(0 , 4), standard_error_DiD= rep(0, 4))
  for(i in 1:4){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[7,1]
    std_error= summary(DiD_prepost)$coefficients[,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Dry season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_all_wet <- function(year_pre, year_post) {
  dt=data_reg_daily_wet
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=data_reg_copy_daily
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:8, coefficient_DiD = rep(0 , 8), standard_error_DiD= rep(0, 8))
  for(i in 1:8){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge, 
                     data = dt[ month == 4+i])
    beta_DiD = summary(DiD_prepost)$coefficients[8,1]
    std_error= summary(DiD_prepost)$coefficients[8,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Rainy season - Before", year_pre, "and after", year_post, sep = " "))
}


FN_coefficients_DiD_all(2008, 2012)
FN_coefficients_DiD_all(2008, 2013)
FN_coefficients_DiD_all(2007, 2013)
FN_coefficients_DiD_all(2007, 2012)

# FN_coefficients_DiD_west_dry(2008, 2012)
# FN_coefficients_DiD_west_dry(2008, 2013)
# FN_coefficients_DiD_west_dry(2007, 2013)
# FN_coefficients_DiD_west_dry(2007, 2012)
# 
# FN_coefficients_DiD_west_wet(2008, 2012)
# FN_coefficients_DiD_west_wet(2008, 2013)
# FN_coefficients_DiD_west_wet(2007, 2013)
# FN_coefficients_DiD_west_wet(2007, 2012)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   f. Entire Watershed with year fixed effects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

FN_coefficients_DiD_all     <- function(year_pre, year_post) {
  dt=data_reg_copy_daily
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + factor(year) + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[11,1]
    std_error= summary(DiD_prepost)$coefficients[11,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, "with year FE", sep = " "))
}

FN_coefficients_DiD_all_dry <- function(year_pre, year_post) {
  dt=data_reg_daily_dry
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:4, coefficient_DiD = rep(0 , 4), standard_error_DiD= rep(0, 4))
  for(i in 1:4){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + factor(year) + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[11,1]
    std_error= summary(DiD_prepost)$coefficients[11,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Dry season - Before", year_pre, "and after", year_post, sep = " "))
}

FN_coefficients_DiD_all_wet <- function(year_pre, year_post) {
  dt=data_reg_daily_wet
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=data_reg_copy_daily
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  coefs_DiD = data.table(month=1:8, coefficient_DiD = rep(0 , 8), standard_error_DiD= rep(0, 8))
  for(i in 1:8){
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + factor(year) + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season, 
                     data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[11,1]
    std_error= summary(DiD_prepost)$coefficients[11,2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  
  ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Rainy season - Before", year_pre, "and after", year_post, sep = " "))
}


FN_coefficients_DiD_all(2008, 2012)

FN_coefficients_DiD_all_dry(2008, 2012)
FN_coefficients_DiD_all_wet(2008, 2012)











FN_coefficients_DiD_monthly     <- function(year_pre, year_post) {
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=copy(data_reg_copy_daily)
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    # DiD_prepost = lm(tot_discharge ~ year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge,
    #                  data = dt[ month == i])
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                     + factor(year)  + lag.discharge
                     , data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),1]
    std_error= summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  # dt_ESE=data_pre[SUBCUENC_2==River_tr]
  
  l=ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD, group=1))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    # geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge),
    #               width=.2, col="red") +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
  return(l)
}

FN_coefficients_DiD_monthly(2008, 2012)
FN_coefficients_DiD_monthly(2007, 2012)
FN_coefficients_DiD_monthly(2008, 2013)
FN_coefficients_DiD_monthly(2007, 2013)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
#                   f. Entire Watershed for annual effects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

FN_coefficients_DiD_monthly     <- function(year_pre, year_post) {
  # dt=data_reg_copy_daily[(River=="Trinidad" | River=="Ciri grande" | River=="Cano Quebrado")]
  dt=copy(data_reg_copy_daily)
  dt = dt[year < year_pre | year > year_post][ ,year_treat:= ifelse(year>2007, 1, 0)]
  dt=dt[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
  coefs_DiD = data.table(month=1:12, coefficient_DiD = rep(0 , 12), standard_error_DiD= rep(0, 12))
  for(i in 1:12){
    # DiD_prepost = lm(tot_discharge ~ year_treat + total_rain_day + lag.rain_day + day_of_season + lag.discharge,
    #                  data = dt[ month == i])
    DiD_prepost = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                     + factor(year)  + lag.discharge
                     , data = dt[ month == i])
    beta_DiD = summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),1]
    std_error= summary(DiD_prepost)$coefficients[nrow(summary(DiD_prepost)$coefficients),2]
    coefs_DiD[i, 2] = beta_DiD
    coefs_DiD[i, 3] = std_error
  }
  # dt_ESE=data_pre[SUBCUENC_2==River_tr]
  
  l=ggplot(data = coefs_DiD, aes(x=factor(month), y=coefficient_DiD, group=1))+
    geom_errorbar(aes(ymin=coefficient_DiD - 2.866 * standard_error_DiD , ymax=coefficient_DiD + 2.866 * standard_error_DiD), width=.1) +
    geom_line() +
    geom_point() +
    # geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge),
    #               width=.2, col="red") +
    xlab("month")+
    ylab("Diff in Diff Coefficient")+
    ggtitle(paste("Before", year_pre, "and after", year_post, sep = " "))
  return(l)
}

FN_coefficients_DiD_monthly(2008, 2012)
FN_coefficients_DiD_monthly(2007, 2012)
FN_coefficients_DiD_monthly(2008, 2013)
FN_coefficients_DiD_monthly(2007, 2013)










FN_coefficients_DiD_annual     <- function(year_pre) {
  dt=copy(data_reg_copy_daily)
  dt=dt[dry_dummy==1]
  coefs_DiD=data.table()
  for (i in year_pre:2015) {
    dt_1=rbind(dt[year< year_pre], dt[year== i])
    dt_1[, year_treat:=ifelse(year==i, 1, 0)]
    dt_1=dt_1[, treatment:=ifelse(River=="Trinidad" | River=="Ciri grande", 1, 0)]
    
    DiD_annual = lm(tot_discharge ~ treatment + year_treat + treatment:year_treat + total_rain_day + lag.rain_day + day_of_season
                    # + factor(year) + lag.discharge
                    , data = dt_1)
    beta_DiD = summary(DiD_annual)$coefficients[nrow(summary(DiD_annual)$coefficients),1]
    std_DiD  = summary(DiD_annual)$coefficients[nrow(summary(DiD_annual)$coefficients),2]
    dt_2=data.table(year=i, beta_DiD, std_DiD)
    coefs_DiD=rbind(coefs, dt_2)
  }
  l=ggplot(data = coefs_DiD, aes(x=factor(year), y=beta_DiD, group=1))+
    geom_errorbar(aes(ymin=beta_DiD - 1.96 * std_DiD , ymax=beta_DiD + 1.96 * std_DiD), width=.1) +
    geom_line() +
    geom_point() +
    # geom_errorbar(data = dt_ESE, aes(y=mean_extra_discharge, ymin=mean_extra_discharge - 1.96 * sd_extra_discharge , ymax=mean_extra_discharge + 1.96 * sd_extra_discharge),
    #               width=.2, col="red") +
    xlab("month")+
    ylab("Diff in Diff Coefficient")
  return(l)
}

FN_coefficients_DiD_annual(2008)







