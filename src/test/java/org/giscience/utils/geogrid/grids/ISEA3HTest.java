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
    public void _numberOfCellsByNumber(int resolution, int numberOfHexagonalCells) {
        ISEA3H grid = new ISEA3H(resolution);
        assertEquals(grid.numberOfHexagonalCells(), numberOfHexagonalCells);
        assertEquals(grid.numberOfPentagonalCells(), 12);
    }

    @Test
    public void numberOfCellsByArea() {
        for (int r = 1; r < 19; r++) this._numberOfCellsByArea(r);
    }
    public void _numberOfCellsByArea(int resolution) {
        ISEA3H grid = new ISEA3H(resolution);
        assertTrue(grid.numberOfHexagonalCells() * grid.areaOfAHexagonalCell() + grid.numberOfPentagonalCells() * grid.areaOfAPentagonalCell() - WGS84.areaOfEarth < WGS84.areaOfEarth * this._precision2);
    }

    @Test
    public void numberOfCellsByGrid() throws Exception {
        for (int r = 1; r <= 12; r++) this._numberOfCellsByGrid(r);
    }
    public void _numberOfCellsByGrid(int resolution) throws Exception {
        ISEA3H grid = new ISEA3H(resolution);
        long cells = grid.cells().size();
        long cellsExpected = grid.numberOfHexagonalCells() + grid.numberOfPentagonalCells();
        System.out.format("resolution %d - %d - %d\n", resolution, cells, cellsExpected);
        assertEquals(cells, cellsExpected);
    }

    @Test
    public void pointsInGridCells() throws Exception {
        this._pointsInGridCells(1);
        this._pointsInGridCells(12);
        this._pointsInGridCells(13);
        this._pointsInGridCells(16);
    }
    private void _pointsInGridCells(int resolution) throws Exception {
        ISEA3H grid = new ISEA3H(resolution);
        for (int i = 0; i < this._iterations; i++) {
            FaceCoordinates c = new FaceCoordinates(1, Math.random() * 100 - 50, Math.random() * 100 - 50);
            assertTrue(c.distanceTo(grid.cellForLocation(c)) <= grid.diameterOfHexagonalCellOnIcosahedron() / 2. + this._precision);
        }
    }
}
