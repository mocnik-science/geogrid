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
package com.giscience.utils.geogrid.examples;

import com.giscience.utils.geogrid.coordinates.FaceCoordinates;
import com.giscience.utils.geogrid.coordinates.GeoCoordinates;
import com.giscience.utils.geogrid.projections.ISEAProjection;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class example {
    public static void main(String[] args) throws Exception {
        ISEAProjection p = new ISEAProjection();
//        p.setOrientationSymmetricEquator();
        int iMax = 1000000;
        long start = System.currentTimeMillis();
        double x = 0;
        for (int i = 0; i < iMax; i++) {
            GeoCoordinates c = new GeoCoordinates(Math.random() * 170 - 85, Math.random() * 360);
            FaceCoordinates c2 = p.sphereToIcosahedron(c);
            GeoCoordinates c3 = p.icosahedronToSphere(c2);
            x += c3.getLat() - c.getLat();
//            System.out.format("lat %.20f lon %.20f%n", c3.getLat() - c.getLat(), c3.getLon() - c.getLon());
/*            if (c3.getLon() - c.getLon() > 0.001) {
                System.out.println("======");
                System.out.println(c);
                System.out.println(c2);
                System.out.format("lat %f lon %f%n", c3.getLat() - c.getLat(), c3.getLon() - c.getLon());
            }*/
        }
        long end = System.currentTimeMillis();
        System.out.format("%f sec%n", (end - start) / 1000.);
        System.out.println(x);
    }
}


/*

icosahedron
north and south points on the edge midpoints (symmetric around the equator)
hexagons, aperture 3 as the smallest possible aligned hexagon aperture
inverse ISEA projection
DGGS point at the center of the corresponding planar cell region
-> ISEA Aperture 3 Hexagon DGGS (ISEA3H)

*/