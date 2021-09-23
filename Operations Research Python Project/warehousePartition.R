library(leaflet)
library(mapview)
warehouseLocations = read.csv('WarehouseLocations.csv')
map = leaflet(warehouseLocations) %>% 
  addTiles() %>% 
  addCircleMarkers(stroke=FALSE,radius=5, fillOpacity = 1, color = ifelse(warehouseLocations['Type'] == 'Distribution','red',ifelse(warehouseLocations['Type'] == 'The Warehouse','blue','green'))) %>% 
  addLegend("bottomright", colors = c('red','blue','green'), labels = c('Distribution','The Warehouse','Noel Leeming')) %>%
  addPolylines(lng = c(174.7301,174.8312), lat = c(-36.9388,-36.8478), color = 'purple')
mapshot(map, file = "warehousePartition.png")

