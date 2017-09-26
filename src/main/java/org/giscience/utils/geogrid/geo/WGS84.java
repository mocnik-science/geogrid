package org.giscience.utils.geogrid.geo;

/**
 * Data about the World Geodetic System (WGS) 84 reference ellipsoid
 *
 * @author Franz-Benjamin Mocnik
 */
public class WGS84 {
    public static final double semiMajorAxis = 6378.137; // [km]
    public static final double inverseFlattening = 1 / 298.257223563; // []
    public static final double semiMinorAxis = (1 - inverseFlattening) * WGS84.semiMajorAxis; // [km]
    public static final double radiusAuthalic = 6371.0071809184728409; // [km]
    public static final double areaOfEarth = 510065621.72408837080; // [km^2]
}
