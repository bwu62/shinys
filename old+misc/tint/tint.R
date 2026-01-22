# remotes::install_version("ggplot2",version="3.4.3")
# remotes::install_version("broom",version="1.0.4")
pacman::p_load(tidyverse,rvest,rgdal,maptools,rgeos,mapproj,magrittr)

### read data

tint = "https://www.raynofilm.com/blog/automotive-window-tint-laws-by-state" %>% 
  read_html %>% 
  html_element(xpath="//table[contains(@class,'table-fixed')]") %>% 
  html_table(header=T) %>% 
  set_names(c("state","front_side","back_side","rear",
              "windshield","reflectivity","other")) %>% 
  mutate_at(c("front_side","back_side","rear"),
            . %>% str_replace("%","") %>% str_replace("any","0") %>% 
              str_replace("no tinting allowed","100")) %>% 
  mutate(front_side_notes = "", windshield_notes = "") %>% 
  mutate(windshield = case_when(
    str_detect(windshield,"(4|5|6) inches|as-1 line") ~ "4-6 inches OR AS-1",
    str_detect(windshield,"70%") ~ "No more than 70%",
    str_detect(windshield,"(no tint|none).*allowed") ~ "No tint allowed"),
    windshield = factor(windshield,levels=c(
      "No tint allowed","4-6 inches OR AS-1","No more than 70%"))) %>% 
  column_to_rownames("state")

tint["Michigan","front_side"] = 100
tint["Michigan","front_side_notes"] = 
  "any percent, but only 4 inches from the top of window"

tint_car = tint_mpv = tint ; rm(tint)
tint_car["Washington D.C.","back_side"] = tint_car["Washington D.C.","rear"] = "50"
tint_mpv["Washington D.C.","back_side"] = tint_mpv["Washington D.C.","rear"] = "35"

tint_car %<>% mutate_at(c("front_side","back_side","rear"),~{as.numeric(.x)/100}) %>% 
  rownames_to_column("state")

tint_mpv %<>% mutate_at(c("front_side","back_side","rear"),~{as.numeric(.x)/100}) %>% 
  rownames_to_column("state")

### hex plot

spdf = readOGR("us_states_hexgrid.geojson.json",verbose=F) %>% suppressWarnings
spdf@data %<>% mutate(state = gsub(" \\(.*\\)", "", google_name),id=iso3166_2)
spdf@data$state[spdf@data$state=="District of Columbia"] = "Washington D.C."

# spdf_fortified %>% filter(id%in%c("ME","NH","RI","AK","HI","OR","NV")) %>%
#   select(id,state,long,lat) %>% group_by(id) %>% summarize(
#     long=paste(sort(unique(long)),collapse=","),
#     lat=paste(sort(unique(lat)),collapse=",")) %>%
#   separate(long,sep=",",into=str_c("x.",1:3)) %>%
#   separate(lat,sep=",",into=str_c("y.",1:4))

shifts = tibble(STATE=c("AK","HI","ME"),X=c(5.446,8.169,2.723),Y=c(-1.1,1,-2.838))
for(i in 1:nrow(shifts)){
  j = which(spdf@data$iso3166_2==shifts$STATE[i])
  spdf@polygons[[j]]@Polygons[[1]]@coords[,1] = 
    spdf@polygons[[j]]@Polygons[[1]]@coords[,1]+shifts$X[i]
  spdf@polygons[[j]]@Polygons[[1]]@coords[,2] = 
    spdf@polygons[[j]]@Polygons[[1]]@coords[,2]+shifts$Y[i]}

spdf_fortified_car = fortify(spdf, region="state") %>% rename(state=id) %>% 
  full_join(select(tint_car,state:windshield)) %>% 
  full_join(select(spdf@data,state,id)) %>% 
  relocate(id,state,long,lat,front_side:last_col(),everything())
centers_car = cbind.data.frame(data.frame(
  gCentroid(spdf, byid=TRUE),id=spdf@data$id,state=spdf@data$state)) %>% 
  full_join(select(tint_car,state:windshield))

spdf_fortified_mpv = fortify(spdf, region="state") %>% rename(state=id) %>% 
  full_join(select(tint_mpv,state:windshield)) %>% 
  full_join(select(spdf@data,state,id)) %>% 
  relocate(id,state,long,lat,front_side:last_col(),everything())
centers_mpv = cbind.data.frame(data.frame(
  gCentroid(spdf, byid=TRUE),id=spdf@data$id,state=spdf@data$state)) %>% 
  full_join(select(tint_mpv,state:windshield))

spdf_fortified_car[spdf_fortified_car$id=="ME",][4,c("long","lat")] = 
  spdf_fortified_car[spdf_fortified_car$id=="RI",][5,c("long","lat")]
spdf_fortified_mpv[spdf_fortified_mpv$id=="ME",][4,c("long","lat")] = 
  spdf_fortified_mpv[spdf_fortified_mpv$id=="RI",][5,c("long","lat")]

spdf_fortified = list(car=spdf_fortified_car,mpv=spdf_fortified_mpv)
centers = list(car=centers_car,mpv=centers_mpv)

type = "car"

(p1 = ggplot(spdf_fortified[[type]]) + 
  geom_polygon(aes(x=long,y=lat,group=state,fill=front_side), color="gray50") +
  geom_text(data=centers[[type]],aes(x,y,label=id,color=ifelse(front_side>.5,"a","b")),size=4) +
  scale_fill_gradientn(limits=c(0,1),colors=c("black","white"),
                       name="VLT",labels = scales::percent(0.25*0:4)) + 
  scale_color_manual(values=c("black","white"),guide=NULL) + ggtitle("Front side") +
  theme_void() + coord_map())

(p2 = ggplot(spdf_fortified[[type]]) + 
  geom_polygon(aes(x=long,y=lat,group=state,fill=back_side), color="gray50") +
  geom_text(data=centers[[type]],aes(x,y,label=id,color=ifelse(back_side>.5,"a","b")),size=4) +
  scale_fill_gradientn(limits=c(0,1),colors=c("black","white"),
                       name="VLT",labels = scales::percent(0.25*0:4)) + 
  scale_color_manual(values=c("black","white"),guide=NULL) + ggtitle("Back side") +
  theme_void() + coord_map())

(p3 = ggplot(spdf_fortified[[type]]) + 
  geom_polygon(aes(x=long,y=lat,group=state,fill=rear), color="gray50") +
  geom_text(data=centers[[type]],aes(x,y,label=id,color=ifelse(rear>.5,"a","b")),size=4) +
  scale_fill_gradientn(limits=c(0,1),colors=c("black","white"),
                       name="VLT",labels = scales::percent(0.25*0:4)) + 
  scale_color_manual(values=c("black","white"),guide=NULL) + ggtitle("Rear") + 
  theme_void() + coord_map())

(p4 = ggplot(spdf_fortified[[type]]) + 
  geom_polygon(aes(x=long,y=lat,group=state,fill=windshield), color="black") +
  scale_fill_manual(values=c("white","#7fcdbb","gray30"),name="Windshield rules") + 
  geom_text(data=centers[[type]],aes(x,y,label=id),size=4,color="black") + ggtitle("Windshield") + 
  theme_void() + coord_map())

library(gridExtra)
grid.arrange(p1,p2,p3,p4,ncol=2)

save(spdf_fortified,centers,file="tint_data.Rdata")
