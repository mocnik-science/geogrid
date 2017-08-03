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
import org.giscience.utils.geogrid.geometry.GridCell;
import org.giscience.utils.geogrid.grids.ISEA3H;

import java.util.Collection;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class example {
    public static void main(String[] args) throws Exception {

        ISEA3H grid3 = new ISEA3H(11);
        Collection<GridCell> cells = grid3.cellsForBound(42, 43, 7, 8);
        System.out.println(cells.size());
        for (GridCell c : cells) System.out.println(c);
        System.out.println("==.==.==");
        for (GridCell c : cells) System.out.format("{lat: %f, lon: %f},", c.getLat(), c.getLon());




        /*

        ISEA3H grid2 = new ISEA3H(1);
        System.out.println(grid2.numberOfHexagonCells() + grid2.numberOfPentagonCells());
        System.out.println(grid2.areaOfHexagonCell());
        System.out.println(grid2.numberOfHexagonCells() * grid2.areaOfHexagonCell() + grid2.numberOfPentagonCells() * grid2.areaOfPentagonCell() - WGS84.areaOfEarth);




        ISEAProjection p = new ISEAProjection();
        p.setOrientationSymmetricEquator();
        int iMax = 1; //000000;
        long start = System.currentTimeMillis();
        for (int i = 0; i < iMax; i++) {
            GeoCoordinates c = new GeoCoordinates(Math.random() * 170 - 85, Math.random() * 360);
            FaceCoordinates c2 = p.sphereToIcosahedron(c);
            GeoCoordinates c3 = p.icosahedronToSphere(c2);
        }
        long end = System.currentTimeMillis();
        System.out.format("%f sec%n", (end - start) / 1000.);


        double lat = Math.random() * 180 - 90;
        double lon = Math.random() * 360;
        System.out.println(new GridCell(12, lat, lon));

        double l = 1;
        double inverseSqrt3 = 1 / Math.sqrt(3);

        ISEA3H grid = new ISEA3H(1);


//        System.out.println(grid.cellForLocation(new FaceCoordinates(1, 0., 0.)));
        example.testGrid(grid, new FaceCoordinates(1, 0., inverseSqrt3 / 2 - .01));
        example.testGrid(grid, new FaceCoordinates(1, 1.802799, 35.548452));
//        System.out.println(grid.cellForLocation(new FaceCoordinates(1, 1 / 2.9, 0.)));

        System.out.println(grid.diameterOfCellOnIcosahedron());

        for (int i = 0; i < 0; i++) {
            FaceCoordinates c = new FaceCoordinates(1, Math.random() * 100 - 50, Math.random() * 100 - 50);
            double d = c.distanceTo(grid.cellForLocation(c));
            System.out.println(d <= grid.diameterOfCellOnIcosahedron() / 2.);
        }
*/
    }

    public static void testGrid(ISEA3H grid, FaceCoordinates c) throws Exception {
        System.out.println("======");
        System.out.println(c);
        System.out.println(grid.cellForLocation(c));
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