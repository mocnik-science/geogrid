package org.giscience.utils.geogrid.geometry;

import org.giscience.utils.geogrid.generic.Trigonometric;
import org.giscience.utils.geogrid.geo.WGS84;

/**
 * Geographic coordinates of a location on Earth.
 *
 * @author Franz-Benjamin Mocnik
 */
public class GeoCoordinates implements Comparable<GeoCoordinates> {
    private final Double _lat;
    private final Double _lon;

    public GeoCoordinates(Double lat, Double lon) throws Exception {
        if (lat < -90 || lat > 90) throw new Exception("invalid latitude");
        lon %= 360;
        if (lon > 180) lon -= 360;
        else if (lon < -180) lon += 360;
        this._lat = lat;
        this._lon = lon;
    }

    public Double getLat() {
        return this._lat;
    }

    public Double getLon() {
        return this._lon;
    }

    public Double distanceTo(GeoCoordinates other) {
        return WGS84.radiusAuthalic * 2 * Math.asin(Math.sqrt(Math.pow(Trigonometric.sin((this.getLat() - other.getLat()) / 2), 2) + Math.pow(Trigonometric.sin((this.getLon() - other.getLon()) / 2), 2) * Trigonometric.cos(this.getLat()) * Trigonometric.cos(other.getLat())));
    }

    @Override
    public String toString() {
        return String.format("lat %f lon %f", this._lat, this._lon);
    }

    @Override
    public int compareTo(GeoCoordinates o) {
        int d = Double.compare(this._lat, o._lat);
        if (d != 0) return d;
        return Double.compare(this._lon, o._lon);
    }
}
