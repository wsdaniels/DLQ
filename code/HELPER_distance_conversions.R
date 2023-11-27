latlon.to.zone.number <- function(latitude, longitude){
  if (56 <= latitude & latitude < 64 & 3 <= longitude & longitude < 12){
    return(32)
  }
  if (72 <= latitude & latitude <= 84 & longitude >= 0){
    if (longitude <= 9){
      return(31)
    } else if (longitude <= 21){
      return(33)
    } else if (longitude <= 33){
      return(35)
    } else if (longitude <= 42){
      return(37)
    }
  }
  return(as.integer((longitude + 180) / 6) + 1)
}

zone.number.to.central.longitude <- function(zone_number){
  return((zone_number - 1) * 6 - 180 + 3)
}

latitude.to.zone.letter <- function(latitude){
  ZONE_LETTERS = c("C", "D", "E", "F", "G", "H", "J",
                   "K", "L", "M", "N", "P", "Q", "R", 
                   "S", "T", "U", "V", "W", "X", "X")
  if (-80 <= latitude & latitude <= 84){
    return(ZONE_LETTERS[bitwShiftR(as.integer(latitude + 80),3)+1])
  } else {
    return(NULL)
  }
}

latlon.to.utm <- function(latitude, longitude, force_zone_number=NULL, R=6378137, E=0.00669438){
  #This function convert Latitude and Longitude to UTM coordinate
  #      Parameters
  #      ----------
  #      latitude: float
  #          Latitude between 80 deg S and 84 deg N, e.g. (-80.0 to 84.0)
  #      longitude: float
  #          Longitude between 180 deg W and 180 deg E, e.g. (-180.0 to 180.0).
  #      force_zone number: int
  #          Zone Number is represented with global map numbers of an UTM Zone
  #          Numbers Map. You may force conversion including one UTM Zone Number.
  #          More information see http://www.jaworski.ca/utmzones.htm
  #      R: float
  #          semimajor axis of the planet, default Earth
  #      E: float
  #          eccentricity of the planet, default Earth
  
  K0 = 0.9996 # UTM scale on the central meridian
  E2 = E * E
  E3 = E2 * E
  E_P2 = E / (1.0 - E)
  M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
  M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
  M3 = (15 * E2 / 256 + 45 * E3 / 1024)
  M4 = (35 * E3 / 3072)
  lat_rad = pi*(latitude)/180
  lat_sin = sin(lat_rad)
  lat_cos = cos(lat_rad)
  lat_tan = lat_sin / lat_cos
  lat_tan2 = lat_tan * lat_tan
  lat_tan4 = lat_tan2 * lat_tan2
  
  if (is.null(force_zone_number)){
    zone_number = latlon.to.zone.number(latitude, longitude)
  } else {
    zone_number = force_zone_number
  }
  
  zone_letter = latitude.to.zone.letter(latitude)
  lon_rad = pi*(longitude)/180
  central_lon = zone.number.to.central.longitude(zone_number)
  central_lon_rad = pi*(central_lon)/180
  
  n = R / sqrt(1 - E * lat_sin**2)
  c = E_P2 * lat_cos**2
  a = lat_cos * (lon_rad - central_lon_rad)
  a2 = a * a
  a3 = a2 * a
  a4 = a3 * a
  a5 = a4 * a
  a6 = a5 * a
  m = R * (M1 * lat_rad -
             M2 * sin(2 * lat_rad) +
             M3 * sin(4 * lat_rad) -
             M4 * sin(6 * lat_rad))
  easting = K0 * n * (a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000
  northing = K0 * (m + n * lat_tan * (a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c**2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))
  #if latitude < 0:
  #    northing += 10000000
  return(c(easting, northing, zone_number, zone_letter))
}