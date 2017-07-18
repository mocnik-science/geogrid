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
package com.giscience.utils.geogrid.geo;

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
