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
package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.geometry.FaceCoordinates;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class ISEA3HTest {
    private final double _precision = 1e-9;
    private final double _precision2 = 1e-15;
    private final int _iterations = 1000000;

    @Test
    public void numberOfCellsByNumber() {
        this._numberOfCellsByNumber(1, 20);
        this._numberOfCellsByNumber(2, 80);
        this._numberOfCellsByNumber(3, 260);
        this._numberOfCellsByNumber(4, 800);
        this._numberOfCellsByNumber(15, 143489060);
        this._numberOfCellsByNumber(16, 430467200);
    }
    public void _numberOfCellsByNumber(int resolution, int numberOfHexagonCells) {
        ISEA3H grid = new ISEA3H(resolution);
        assertEquals(grid.numberOfHexagonCells(), numberOfHexagonCells);
        assertEquals(grid.numberOfPentagonCells(), 12);
    }

    @Test
    public void numberOfCellsByArea() {
        for (int r = 1; r < 19; r++) this._numberOfCellsByArea(r);
    }
    public void _numberOfCellsByArea(int resolution) {
        ISEA3H grid = new ISEA3H(resolution);
        assertTrue(grid.numberOfHexagonCells() * grid.areaOfHexagonCell() + grid.numberOfPentagonCells() * grid.areaOfPentagonCell() - WGS84.areaOfEarth < WGS84.areaOfEarth * this._precision2);
    }

    @Test
    public void numberOfCellsByGrid() throws Exception {
        for (int r = 1; r < 13; r++) this._numberOfCellsByGrid(r);
    }
    public void _numberOfCellsByGrid(int resolution) throws Exception {
        System.out.println(resolution);
        ISEA3H grid = new ISEA3H(resolution);
        System.out.println(grid.cellsForBound(-90, 90, -180, 180).size() + " " + (grid.numberOfHexagonCells() + grid.numberOfPentagonCells()));
        assertEquals(grid.cellsForBound(-90, 90, -180, 180).size(), grid.numberOfHexagonCells() + grid.numberOfPentagonCells());
    }

    @Test
    public void pointsInGridCells1() throws Exception {
        this._pointsInGridCells(1);
    }
    @Test
    public void pointsInGridCells12() throws Exception {
        this._pointsInGridCells(12);
    }
    @Test
    public void pointsInGridCells13() throws Exception {
        this._pointsInGridCells(13);
    }
    @Test
    public void pointsInGridCells16() throws Exception {
        this._pointsInGridCells(16);
    }
    private void _pointsInGridCells(int resolution) throws Exception {
        ISEA3H grid = new ISEA3H(resolution);
        for (int i = 0; i < this._iterations; i++) {
            FaceCoordinates c = new FaceCoordinates(1, Math.random() * 100 - 50, Math.random() * 100 - 50);
            assertTrue(c.distanceTo(grid.cellForLocation(c)) <= grid.diameterOfCellOnIcosahedron() / 2. + this._precision);
        }
    }
}
