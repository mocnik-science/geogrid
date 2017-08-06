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
package org.giscience.utils.geogrid.examples;

import org.giscience.utils.geogrid.geometry.FaceCoordinates;
import org.giscience.utils.geogrid.geometry.GeoCoordinates;
import org.giscience.utils.geogrid.geometry.GridCell;
import org.giscience.utils.geogrid.grids.ISEA3H;
import org.giscience.utils.geogrid.projections.ISEAProjection;

import java.util.Collection;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class Example {
    public static void main(String[] args) throws Exception {

        // PROJECTION
        ISEAProjection p = new ISEAProjection();

        // project coordinates from the sphere to the icosahedron and back
        for (int i = 0; i < 10; i++) {
            GeoCoordinates c = new GeoCoordinates(Math.random() * 180 - 90, Math.random() * 360);
            FaceCoordinates c2 = p.sphereToIcosahedron(c);
            GeoCoordinates c3 = p.icosahedronToSphere(c2);
            System.out.println(c);
            System.out.println(c2);
            System.out.println(c3);
            System.out.println("------");
        }

        // GRID
        ISEA3H g = new ISEA3H(14);

        // print properties of the grid
        System.out.format("number of hexagon cells: %d%n", g.numberOfHexagonCells());
        System.out.format("number of pentagon cells: %d%n", g.numberOfPentagonCells());
        System.out.format("diameter of a hexagon cell: %f%n", g.diameterOfHexagonCellOnIcosahedron());
        System.out.format("area of a hexagon cell: %f%n", g.areaOfHexagonCell());
        System.out.println("------");

        // get cells in given bounds
        Collection<GridCell> cells = g.cellsForBound(41, 42, 6, 7);
        System.out.println(cells.size());
        System.out.println("------");

        // determine cell for geographic coordinates
        for (int i = 0; i < 10; i++) {
            double lat = Math.random() * 180 - 90;
            double lon = Math.random() * 360;
            GridCell c = g.cellForLocation(lat, lon);
            System.out.format("lat %f - lon %f%n", lat, lon);
            System.out.println(c);
            System.out.println("------");
        }
        // alternatively, cellForLocation accepts a GeoCoordinates object
        // the method cellForCentroid can be used to compute the grid cell for a jts Geometry
    }
}


/*

icosahedron
north and south points on the edge midpoints (symmetric around the equator)
hexagons, aperture 3 as the smallest possible aligned hexagon aperture
inverse ISEA projection
DGGS point at the center of the corresponding planar cell region
-> ISEA Aperture 3 Hexagon DGGS (ISEA3H)


        // -9 223 372 036 854 775 808
        // -92 23372036 854775808
        // -92 12345678 123456789
        // 90.372,036 lat
        // 180.775,808 lon
        // resolution: 111 km = 1 <=> 111 m = 0.001 <=> 0,111 m = 0.000,001

40,075 km / 5 faces = 8,015 km / face

8,015 km * 2/3 * (1/sqrt(3))^15 = 1.4106 m
8,015 km * 2/3 * (1/sqrt(3))^16 = 0.8144 m
8,015 km * 2/3 * (1/sqrt(3))^17 = 0.4702 m
8,015 km * 2/3 * (1/sqrt(3))^18 = 0.2715 m
8,015 km * 2/3 * (1/sqrt(3))^19 = 0.1567 m
8,015 km * 2/3 * (1/sqrt(3))^20 = 0.0905 m

resolution between 0 and 18

erdoberfläche: 510,100,000 km^2
-> radius: 6371.2218793249 km
-> umfang: 40031.5677009551 km

Documentation

* result: only between -180 and 180 lon

*/
