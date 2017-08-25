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

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import org.giscience.utils.geogrid.generic.Tuple;
import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.geometry.FaceCoordinates;
import org.giscience.utils.geogrid.geometry.GeoCoordinates;
import org.giscience.utils.geogrid.geometry.GridCell;
import org.giscience.utils.geogrid.projections.ISEAProjection;

import java.util.*;

/**
 * ISEA Aperture 3 Hexagon (ISEA3H) Discrete Global Grid System (DGGS)
 *
 * The ISEA3H grid is constructed by using the icosahedron Snyder equal-area (ISEA) projection to map the surface of the
 * globe to the icosahedron. Thereby, the orientation of the icosahedron is chosen such that the north and the south
 * poles are mapped to the edge midpoints of the icosahedron. The equator is thus mapped symmetrically. A grid
 * (aperture 3) is constructed on the icosahedron, and this grid is mapped back by the inverse projection to the globe.
 *
 * The cells of the grid are identified by the resolution and their center points.
 *
 * The ISEA3H has been proposed by:
 *
 * Kevin Sahr, Denis White, and A. Jon Kimerling: Geodesic Discrete Global Grid Systems. Cartography and Geographic
 * Information Science, 30(2), 121–134, 2003.
 *
 * @author Franz-Benjamin Mocnik
 */
public class ISEA3H {
    private final double _precision = 1e-9;
    private final ISEAProjection _projection = new ISEAProjection();
    private final int _resolution; // resolution - 1
    private final long _numberOfHexagonCells;
    private final int _numberOfPentagonCells = 12;
    private final double _l0; // length of the triangle base at resolution 0
    private final double _inverseSqrt3l0; // 1 / \sqrt{3} * l_0
    private final double _l; // length of the triangle base at the given resolution
    private final double _l2; // l / 2
    private final double _l6; // l / 6
    private final double _l23; // l * 2 / 3
    private final double _inverseSqrt3 = 1 / Math.sqrt(3); // 1 / \sqrt{3}
    private final double _inverseSqrt3l; // 1 / \sqrt{3} * l
    private final double _inverseSqrt3l2; // 1 / (2 \sqrt{3}) * l
    private final double _triangleA; // l_0 / 2 // half base
    private final double _triangleB; // 1/\sqrt{3} * l_0 // distance center point to tip
    private final double _triangleC; // 1/(2 \sqrt{3}) * l_0 // distance base to center point
    private final double _triangleBCA; // (_triangleB + _triangleC) / _triangleA

    public ISEA3H(int resolution) {
        this(resolution, true);
    }

    public ISEA3H(int resolution, boolean rotatedProjection) {
        if (rotatedProjection) this._projection.setOrientationSymmetricEquator();
        this._resolution = resolution - 1;
        long numberOfHexagonCells = 1;
        for (int i = 0; i < this._resolution; i++) numberOfHexagonCells = 3 * numberOfHexagonCells + 1;
        this._numberOfHexagonCells = 20 * numberOfHexagonCells;
        this._l0 = this._projection.lengthOfTriangleBase();
        this._inverseSqrt3l0 = this._inverseSqrt3 * this._l0;
        this._l = Math.pow(this._inverseSqrt3, this._resolution) * this._l0;
        this._l2 = this._l / 2.;
        this._l6 = this._l / 6.;
        this._l23 = this._l * 2 / 3.;
        this._inverseSqrt3l = this._inverseSqrt3 * this._l;
        this._inverseSqrt3l2 = this._inverseSqrt3l / 2.;
        this._triangleA = this._l0 / 2.;
        this._triangleB = this._inverseSqrt3l0;
        this._triangleC = this._inverseSqrt3l0 / 2.;
        this._triangleBCA = (this._triangleB + this._triangleC) / this._triangleA;
    }

    /**
     * @return diameter of a hexagonal cell on the icosahedron
     */
    public double diameterOfHexagonCellOnIcosahedron() {
        return this._l23;
    }

    /**
     * @return length of a side of a hexagonal cell on the icosahedron
     */
    public double lengthOfASideOfHexagonCellOnIcosahedron() {
        return this._l / 3;
    }

    /**
     * @return lower bound for the length of a side of a hexagonal cell on the sphere
     */
    public double lowerBoundForLengthOfASideOfHexagonCellOnSphere() {
        return this._projection.sphericalDistanceFromCenterToVerticesOnSphere() * WGS84.radiusAuthalic * 2 * Math.PI / (360 * Math.sqrt(Math.pow(3, this._resolution) * 5));
    }

    /**
     * Returns the area of a hexagonal cell. The cells should all have the same area by construction, because the ISEA
     * projection is equal-area.
     *
     * @return area of a hexagonal cell
     */
    public double areaOfAHexagonCell() {
        return WGS84.areaOfEarth / (this._numberOfHexagonCells + 5 / 6. * this._numberOfPentagonCells);
    }

    /**
     * Returns the area of a pentagonal cell. The cells should all have the same area by construction, because the ISEA
     * projection is equal-area.
     *
     * @return area of a pentagonal cell
     */
    public double areaOfAPentagonCell() {
        return 5 / 6. * this.areaOfAHexagonCell();
    }

    /**
     * @return number of hexagonal cells
     */
    public long numberOfHexagonCells() {
        return this._numberOfHexagonCells;
    }

    /**
     * @return number of pentagonal cells
     */
    public int numberOfPentagonCells() {
        return this._numberOfPentagonCells;
    }

    /**
     * Returns the grid cell for a given location
     *
     * @param lat latitude
     * @param lon longitude
     * @return corresponding grid cell
     * @throws Exception
     */
    public GridCell cellForLocation(double lat, double lon) throws Exception {
        return this.cellForLocation(new GeoCoordinates(lat, lon));
    }

    /**
     * Returns the grid cell for a given location
     *
     * @param c geographic coordinates
     * @return corresponding grid cell
     * @throws Exception
     */
    public GridCell cellForLocation(GeoCoordinates c) throws Exception {
        FaceCoordinates fc = this.cellForLocation(this._projection.sphereToIcosahedron(c));
        return this._newGridCell(this._projection.icosahedronToSphere(fc), fc);
    }

    /**
     * Returns the coordinates of the center of the corresponding grid cell for given coordinates in the face
     *
     * @param c face coordinates
     * @return face coordinates of the center of the corresponding grid cell
     * @throws Exception
     */
    public FaceCoordinates cellForLocation(FaceCoordinates c) throws Exception {
        double x = (this._coordinatesNotSwapped()) ? c.getX() : c.getY();
        double y = (this._coordinatesNotSwapped()) ? c.getY() : c.getX();
        int nxCenter = (int) Math.round(x / this._l2);
        double xCenter = nxCenter * this._l2;
        if (Math.abs(x - xCenter) <= this._l6) {
            int nyCenter12 = (int) Math.round((((nxCenter % 2 == 0) ? y : y - this._inverseSqrt3l2)) / this._inverseSqrt3l);
            double yCenter12 = nyCenter12 * this._inverseSqrt3l + ((nxCenter % 2 == 0) ? 0 : this._inverseSqrt3l2);
            return this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, yCenter12);
        }
        else {
            int nyCenter1 = (int) Math.round(y / this._inverseSqrt3l);
            double yCenter1 = nyCenter1 * this._inverseSqrt3l;
            int nyCenter2 = (int) Math.round((y - this._inverseSqrt3l2) / this._inverseSqrt3l);
            double yCenter2 = nyCenter2 * this._inverseSqrt3l + this._inverseSqrt3l2;
            FaceCoordinates cCandidate1 = this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, (nxCenter % 2 == 0) ? yCenter1 : yCenter2);
            FaceCoordinates cCandidate2 = this._faceCoordinatesSwapByResolution(c.getFace(), (x > xCenter) ? xCenter + this._l2 : xCenter - this._l2, (nxCenter % 2 != 0) ? yCenter1 : yCenter2);
            return (c.distanceTo(cCandidate1) < c.distanceTo(cCandidate2)) ? cCandidate1 : cCandidate2;
        }
    }

    /**
     * The following equality holds: cellForLocation = _getCoordinatesOfCenter . _integerForFaceCoordinates
     *
     * However, the method cellForLocation is defined separately to speed up computations.
     *
     * @param face face to compute the coordinates on
     * @param nx steps into the direction of the vertex of the hexagon
     * @param ny steps into the direction of the edge of the hexagon
     * @return coordinates on the face
     */
    private FaceCoordinates _getCoordinatesOfCenter(int face, int nx, int ny) {
        double x = nx * this._l2;
        double y = (ny + ((nx % 2 == 0) ? 0 : .5)) * this._inverseSqrt3l;
        return this._faceCoordinatesSwapByResolution(face, x, y);
    }

    private Tuple<Integer, Integer> _integerForFaceCoordinates(FaceCoordinates c) {
        double x = (this._coordinatesNotSwapped()) ? c.getX() : c.getY();
        double y = (this._coordinatesNotSwapped()) ? c.getY() : c.getX();
        int nxCenter = (int) Math.round(x / this._l2);
        double xCenter = nxCenter * this._l2;
        if (Math.abs(x - xCenter) <= this._l6) return new Tuple(nxCenter, (int) Math.round((((nxCenter % 2 == 0) ? y : y - this._inverseSqrt3l2)) / this._inverseSqrt3l));
        else {
            int nyCenter1 = (int) Math.round(y / this._inverseSqrt3l);
            double yCenter1 = nyCenter1 * this._inverseSqrt3l;
            int nyCenter2 = (int) Math.round((y - this._inverseSqrt3l2) / this._inverseSqrt3l);
            double yCenter2 = nyCenter2 * this._inverseSqrt3l + this._inverseSqrt3l2;
            FaceCoordinates cCandidate1 = this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, (nxCenter % 2 == 0) ? yCenter1 : yCenter2);
            FaceCoordinates cCandidate2 = this._faceCoordinatesSwapByResolution(c.getFace(), (x > xCenter) ? xCenter + this._l2 : xCenter - this._l2, (nxCenter % 2 != 0) ? yCenter1 : yCenter2);
            return (c.distanceTo(cCandidate1) < c.distanceTo(cCandidate2)) ? new Tuple(nxCenter, (nxCenter % 2 == 0) ? nyCenter1 : nyCenter2) : new Tuple((x > xCenter) ? nxCenter + 1 : nxCenter - 1, (nxCenter % 2 != 0) ? nyCenter1 : nyCenter2);
        }
    }

    /**
     * Returns the grid cell for the centroid of a given geometry
     *
     * @param g geometry
     * @return corresponding grid cell
     * @throws Exception
     */
    public GridCell cellForCentroid(Geometry g) throws Exception {
        Coordinate c = g.getCentroid().getCoordinate();
        return this.cellForLocation(c.y, c.x);
    }

    /**
     * Returns cells that are inside the bounds, or at least very near.
     *
     * Note that the result should, in fact, include all cells whose center points are inside the given bounds, but also
     * cells nearby. Observe that it can, however, not be guaranteed that all such cells are returned if the longitude
     * range is less than 360 degrees or the latitude range is less than 180 degrees.
     *
     * @param lat0
     * @param lat1
     * @param lon0
     * @param lon1
     * @return cells inside the bounds
     */
    public Collection<GridCell> cellsForBound(double lat0, double lat1, double lon0, double lon1) throws Exception {
        return _cellsForBound(new CellAggregatorByCells(), lat0, lat1, lon0, lon1).cellAggregator.getCells();
    }

    /**
     * Returns cell ids that are inside the bounds, or at least very near.
     *
     * Same as cellsForBound, apart that this method returns only ids and is much more memory efficient.
     *
     * @param lat0
     * @param lat1
     * @param lon0
     * @param lon1
     * @return ids of cells inside the bounds
     */
    public Collection<Long> cellIdsForBound(double lat0, double lat1, double lon0, double lon1) throws Exception {
        return _cellsForBound(new CellAggregatorByCellIds(), lat0, lat1, lon0, lon1).cellAggregator.getCellIds();
    }

    public <T extends CellAggregator> ResultCellForBound<T> _cellsForBound(T ca, double lat0, double lat1, double lon0, double lon1) throws Exception {
        // compute center of bounding box
        double lat = (lat0 + lat1) / 2.;
        double lon = (lon0 + lon1 + ((lon0 <= lon1) ? 0 : 360)) / 2.;
        if (lon > 360) lon -= 360;
        FaceCoordinates fc = this._projection.sphereToIcosahedron(new GeoCoordinates(lat, lon));
        // compute
        return _cellsForBound(new ResultCellForBound<>(ca), fc, lat0, lat1, lon0, lon1);
    }

    private <T extends CellAggregator> ResultCellForBound<T> _cellsForBound(ResultCellForBound<T> result, FaceCoordinates fcStart, double lat0, double lat1, double lon0, double lon1) throws Exception {
        // if fcStart is already in result, skip the computation
        if (result.visitedCells.contains(fcStart.getFace())) return result;
        result.visitedCells.add(fcStart.getFace());
        // prepare dNs
        List<Tuple<Integer, Integer>> dNs = new ArrayList<>();
        dNs.add(new Tuple<>(1, 1));
        dNs.add(new Tuple<>(-1, 1));
        dNs.add(new Tuple<>(1, -1));
        dNs.add(new Tuple<>(-1, -1));
        // compute
        Tuple<Integer, Integer> fcn = this._integerForFaceCoordinates(fcStart);
        for (Tuple<Integer, Integer> dN : dNs) {
            FaceCoordinates fc = null;
            GeoCoordinates gc = null;
            Set<Integer> success = new HashSet();
            Set<Integer> successLast;
            boolean hasFoundInside;
            Map<Integer, FaceCoordinates> faceTodo = new HashMap();
            int nx = fcn._1 - dN._1;
            while (true) {
                nx += dN._1;
                int ny = fcn._2 - dN._2;
                successLast = success;
                Integer maxMinValue = (!successLast.isEmpty()) ? ((dN._2 >= 0) ? Collections.max(successLast) : Collections.min(successLast)) : null;
                success = new HashSet<>();
                hasFoundInside = false;
                while (true) {
                    ny += dN._2;
                    fc = this._getCoordinatesOfCenter(fcStart.getFace(), nx, ny);
                    gc = this._projection.icosahedronToSphere(fc);
                    if (this._isInside(gc, lat0, lat1, lon0, lon1)) hasFoundInside = true;
                    else if (hasFoundInside || maxMinValue == null || ((dN._2 >= 0) ? ny > maxMinValue : ny < maxMinValue)) break;
                    else continue;
                    if (this._isCoordinatesInFace(fc)) {
                        success.add(ny);
                        result.cellAggregator.add(this._newGridCell(gc, fc));
                    } else {
                        FaceCoordinates fc2 = this._projection.sphereToIcosahedron(gc);
                        if (faceTodo.containsKey(fc2.getFace())) faceTodo.put(fc2.getFace(), fc2);
                        else {
                            int sizeBefore = result.cellAggregator.size();
                            result = this._cellsForBound(result, fc2, lat0, lat1, lon0, lon1);
                            if (result.cellAggregator.size() != sizeBefore) faceTodo.put(fc2.getFace(), null);
                        }
                    }
                }
                if (success.isEmpty()) break;
            }
            for (Map.Entry<Integer, FaceCoordinates> e : faceTodo.entrySet()) if (e.getValue() != null && !result.cellAggregator.contains(this.cellForLocation(this._projection.icosahedronToSphere(e.getValue())))) result = this._cellsForBound(result, e.getValue(), lat0, lat1, lon0, lon1);
        }
        return result;
    }

    private class ResultCellForBound<T extends CellAggregator> {
        public T cellAggregator;
        public Set<Integer> visitedCells = new HashSet();

        public ResultCellForBound(T cellAggregator) {
            this.cellAggregator = cellAggregator;
        }
    }

    private abstract class CellAggregator {
        public abstract void add(GridCell gc);
        public abstract int size();
        public abstract boolean contains(GridCell gc);
    }

    private class CellAggregatorByCells extends CellAggregator {
        private Set<GridCell> _cells = new HashSet();

        @Override
        public void add(GridCell gc) {
            this._cells.add(gc);
        }

        @Override
        public int size() {
            return this._cells.size();
        }

        @Override
        public boolean contains(GridCell gc) {
            return this._cells.contains(gc);
        }

        public Set<GridCell> getCells() {
            return this._cells;
        }
    }

    private class CellAggregatorByCellIds extends CellAggregator {
        private Set<Long> _cells = new HashSet();

        @Override
        public void add(GridCell gc) {
            this._cells.add(gc.getId());
        }

        @Override
        public int size() {
            return this._cells.size();
        }

        @Override
        public boolean contains(GridCell gc) {
            return this._cells.contains(gc);
        }

        public Set<Long> getCellIds() {
            return this._cells;
        }
    }

    private boolean _isInside(GeoCoordinates c, double lat0, double lat1, double lon0, double lon1) {
        boolean bLat = lat0 <= c.getLat() && c.getLat() <= lat1;
        boolean bLon = (lon0 <= lon1) ? lon0 <= c.getLon() && c.getLon() <= lon1 : lon0 < c.getLon() || c.getLon() < lon1;
        return bLat && bLon;
    }

    private GridCell _newGridCell(GeoCoordinates gc, FaceCoordinates fc) throws Exception {
        int d = this._projection.faceOrientation(fc);
        boolean isPentagon = (Math.abs(Math.abs(fc.getX()) - this._triangleA) < this._precision && Math.abs(fc.getY() + d * this._triangleC) < this._precision) || (Math.abs(fc.getX()) < this._precision && Math.abs(fc.getY() - d * this._triangleB) < this._precision);
        return new GridCell(this._resolution, gc, isPentagon);
    }

    /**
     * Returns a buffer for a bounding box of geographic coordinates that needs to be considered in order to ensure
     * that all grid cells, whose center point is within the given bounding box, are contained in the bounding box or
     * the buffer
     *
     * @param lat0
     * @param lat1
     * @param lon0
     * @param lon1
     * @return buffer
     * @throws Exception
     */
    public Double bufferEstimator(double lat0, double lat1, double lon0, double lon1) throws Exception {
        FaceCoordinates southWest = this._projection.sphereToIcosahedron(new GeoCoordinates(lat0, lon0));
        FaceCoordinates northEast = this._projection.sphereToIcosahedron(new GeoCoordinates(lat0, lon0));
        GeoCoordinates southWest2 = this._projection.icosahedronToSphere(new FaceCoordinates(southWest.getFace(), southWest.getX() - this.diameterOfHexagonCellOnIcosahedron() / 2, southWest.getY() - this.diameterOfHexagonCellOnIcosahedron() / 2));
        GeoCoordinates northEast2 = this._projection.icosahedronToSphere(new FaceCoordinates(northEast.getFace(), northEast.getX() - this.diameterOfHexagonCellOnIcosahedron() / 2, northEast.getY() - this.diameterOfHexagonCellOnIcosahedron() / 2));
        List<Double> l = new ArrayList<>();
        l.add(Math.abs(southWest2.getLat() - lat0));
        l.add(Math.abs(southWest2.getLon() - lon0));
        l.add(Math.abs(northEast2.getLat() - lat1));
        l.add(Math.abs(northEast2.getLon() - lon1));
        return l.stream().max(Double::compareTo).get();
    }

    private boolean _isCoordinatesInFace(FaceCoordinates c) {
        int d = this._projection.faceOrientation(c);

        // test whether coordinate is left of the triangle, right of the triangle, or below the triangle
        if (d * c.getY() > c.getX() * this._triangleBCA + this._triangleB + this._precision) return false;
        if (d * c.getY() > -c.getX() * this._triangleBCA + this._triangleB + this._precision) return false;
        if (d * c.getY() < -this._triangleC - this._precision) return false;

        return true;
    }

    private FaceCoordinates _faceCoordinatesSwapByResolution(int face, double x, double y) {
        return new FaceCoordinates(face, this._coordinatesNotSwapped() ? x : y, this._coordinatesNotSwapped() ? y : x);
    }

    private boolean _coordinatesNotSwapped() {
        return this._resolution % 2 == 0;
    }
}
