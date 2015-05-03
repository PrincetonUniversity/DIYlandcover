# Title      : polyFromWkt.R
# Purpose    : Subset of functions used in evaluating map accuracy
# Author     : Lyndon Estes

polyFromWkt <- function(geom.tab, crs) {
  # Creates spatialPolygonsDataFrame from geometry table read in from Postgis
  # Args: 
  #   geom.tab: Dataframe with geometry and identifiers in it.
  #   crs: Coordinate reference system
  # Returns: 
  #   A SpatialPolygonsDataFrame
  # Notes:  Identifier must be 1st column, geometries 2nd col  
  polys <- tst <- sapply(1:nrow(geom.tab), function(x) {
    poly <- as(readWKT(geom.tab[x, 2], p4s = crs), "SpatialPolygonsDataFrame")
    poly@data$ID <- geom.tab[x, 1]
    newid <- paste(x)
    poly <- spChFIDs(poly, newid)
    return(poly)
  })
  polyspdf <- do.call("rbind", polys)
  return(polyspdf)
}

