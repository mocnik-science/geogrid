/*
 * Copyright (C) 2017 Franz-Benjamin Mocnik, Heidelberg University
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.giscience.utils.geogrid.geometry;

/**
 * Geographic coordinates of a location on Earth.
 *
 * @author Franz-Benjamin Mocnik
 */
public class GridCell implements Comparable<GridCell> {
    private final double _precision = 1e-9;
    private final Integer _resolution;
    private final double _lat;
    private final double _lon;
    private final boolean _isPentagon;
    private Long _id = null;

    public GridCell(int resolution, double lat, double lon, boolean isPentagon) throws Exception {
        if (resolution < 0 || resolution > 18) throw new Exception("resolution must be between 0 and 18");
        this._resolution = resolution;
        if (lat < -90 || lat > 90) throw new Exception("invalid latitude");
        lon = lon % 360;
        if (lon > 180) lon -= 360;
        else if (lon < -180) lon += 360;
        this._lat = lat;
        this._lon = lon;
        this._isPentagon = isPentagon;
    }

    public GridCell(int resolution, GeoCoordinates c, boolean isPentagon) throws Exception {
        this(resolution, c.getLat(), c.getLon(), isPentagon);
    }

    public Integer getResolution() {
        return this._resolution;
    }

    public Double getLat() {
        return this._lat;
    }

    public Double getLon() {
        return this._lon;
    }

    public Boolean isPentagon() {
        return this._isPentagon;
    }

    /**
     * Returns an id of the cell. The id consists of the following elements:
     * <ul>
     *     <li>The leading sign is positive in case of a hexagon, and negative in case of a pentagon.</li>
     *     <li>This first two digits consist of the resolution incremented by 20 in case of negative latitude, by 40 in
     *         case of negative longitude, and by 60 in case of negative latitude and longitude.</li>
     *     <li>The consecutive digits consist of the latitude, with two pre-decimal and six decimal places.</li>
     *     <li>The consecutive digits consist of the longitude, with three pre-decimal and six decimal places.</li>
     * </ul>
     *
     * The id is only valid for resolution smaller less or equal 18.
     *
     * @return id of the cell
     */
    public Long getId() {
        if (this._id == null) {
            long sgnLat = (this._lat < 0) ? 20 : 0;
            long sgnLon = (this._lon < 0) ? 40 : 0;
            this._id = (this._isPentagon ? -1 : 1) * ((this._resolution.longValue() + sgnLat + sgnLon) * (long) 1e17 + Math.abs(Math.round(this._lat * 1e6)) * (long) 1e9 + Math.abs(Math.round(this._lon * 1e6)));
        }
        return this._id;
    }

    public boolean equals(Object o) {
        return (o instanceof GridCell) && ((GridCell) o).getId().equals(this.getId());
    }

    public int hashCode() {
        return Long.hashCode(this.getId());
    }

    @Override
    public String toString() {
        return String.format("resolution %d lat %f lon %f - %d", this._resolution, this._lat, this._lon, this.getId());
    }

    @Override
    public int compareTo(GridCell o) {
        int d = Integer.compare(this._resolution, o._resolution);
        if (d != 0) return d;
        d = (Math.abs(this._lat - o._lat) < this._precision) ? 0 : Double.compare(this._lat, o._lat);
        if (d != 0) return d;
        return (Math.abs(this._lon - o._lon) < this._precision) ? 0 : Double.compare(this._lon, o._lon);
    }
}
