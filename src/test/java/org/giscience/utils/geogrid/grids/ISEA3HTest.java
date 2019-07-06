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
        this.numberOfCellsByNumber(false);
    }
    @Test
    public void numberOfCellsByNumberRotated() {
        this.numberOfCellsByNumber(true);
    }
    public void numberOfCellsByNumber(boolean rotatedProjection) {
        this._numberOfCellsByNumber(1, 20, rotatedProjection);
        this._numberOfCellsByNumber(2, 80, rotatedProjection);
        this._numberOfCellsByNumber(3, 260, rotatedProjection);
        this._numberOfCellsByNumber(4, 800, rotatedProjection);
        this._numberOfCellsByNumber(15, 143489060, rotatedProjection);
        this._numberOfCellsByNumber(16, 430467200, rotatedProjection);
    }
    public void _numberOfCellsByNumber(int resolution, int numberOfHexagonalCells, boolean rotatedProjection) {
        ISEA3H grid = new ISEA3H(resolution, rotatedProjection);
        assertEquals(grid.numberOfHexagonalCells(), numberOfHexagonalCells);
        assertEquals(grid.numberOfPentagonalCells(), 12);
    }

    @Test
    public void numberOfCellsByArea() {
        this.numberOfCellsByArea(false);
    }
    @Test
    public void numberOfCellsByAreaRotated() {
        this.numberOfCellsByArea(true);
    }
    public void numberOfCellsByArea(boolean rotatedProjection) {
        for (int r = 1; r < 19; r++) this._numberOfCellsByArea(r, rotatedProjection);
    }
    public void _numberOfCellsByArea(int resolution, boolean rotatedProjection) {
        ISEA3H grid = new ISEA3H(resolution, rotatedProjection);
        assertTrue(grid.numberOfHexagonalCells() * grid.areaOfAHexagonalCell() + grid.numberOfPentagonalCells() * grid.areaOfAPentagonalCell() - WGS84.areaOfEarth < WGS84.areaOfEarth * this._precision2);
    }

    @Test
    public void numberOfCellsByGrid() throws Exception {
        this.numberOfCellsByGrid(false);
    }
    @Test
    public void numberOfCellsByGridRotated() throws Exception {
        this.numberOfCellsByGrid(true);
    }
    public void numberOfCellsByGrid(boolean rotatedProjection) throws Exception {
        for (int r = 1; r <= 12; r++) this._numberOfCellsByGrid(r, rotatedProjection);
    }
    public void _numberOfCellsByGrid(int resolution, boolean rotatedProjection) throws Exception {
        ISEA3H grid = new ISEA3H(resolution, rotatedProjection);
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
