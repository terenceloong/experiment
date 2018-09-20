kmlwriteline('./kml/a.kml', pos(:,2),pos(:,3), 'Color','g', 'Width',2);

kmlwriteline('./kml/b.kml', nav_gps(:,2),nav_gps(:,3), 'color','b', 'Width',2);

kmlwriteline('./kml/c.kml', nav(:,1),nav(:,2), 'color','r', 'Width',2);