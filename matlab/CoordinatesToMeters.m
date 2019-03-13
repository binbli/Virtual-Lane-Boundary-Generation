function dist = CoordinatesToMeters(latitude1,longitude1,latitude2,longitude2) 

m_earthDiameterMeters = 6365.998 * 2 * 1000;
latitude1  = deg2rad(latitude1);
longitude1 = deg2rad(longitude1);
latitude2  = deg2rad(latitude2);
longitude2 = deg2rad(longitude2);

x = sin((latitude2 - latitude1)/2);
y = sin((longitude2 -longitude1)/2);

dist =  m_earthDiameterMeters * asin(sqrt((x*x) + (cos(latitude1)*cos(latitude2)*y*y)));


