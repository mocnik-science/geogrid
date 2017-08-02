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
 * Earth to the icosahedron. Thereby, the orientation of the icosahedron is chosen such that the north and the south
 * poles are mapped to the edge midpoints of the icosahedron. The equator is thus mapped symmetrically. A grid
 * (aperture 3) is constructed on the icosahedron, and this grid is mapped back by the inverse projection to the Earth.
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
    private final int _numberOfHexagonCells;
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
        this._projection.setOrientationSymmetricEquator();
        this._resolution = resolution - 1;
        int numberOfHexagonCells = 1;
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
     * @return diameter of a cell
     */
    public double diameterOfCellOnIcosahedron() {
        return this._l23;
    }

    /**
     * Returns the area of a hexagon cell. The cells should all have the same area by construction, because the ISEA
     * projection is equal-area.
     *
     * @return area of a hexagon cell
     */
    public double areaOfHexagonCell() {
        return WGS84.areaOfEarth / (this._numberOfHexagonCells + 5 / 6. * this._numberOfPentagonCells);
    }

    /**
     * Returns the area of a pentagon cell. The cells should all have the same area by construction, because the ISEA
     * projection is equal-area.
     *
     * @return area of a pentagoncell
     */
    public double areaOfPentagonCell() {
        return 5 / 6. * this.areaOfHexagonCell();
    }

    /**
     * @return number of hexagon cells
     */
    public int numberOfHexagonCells() {
        return this._numberOfHexagonCells;
    }

    /**
     * @return number of pentagon cells
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
        double nxCenter = Math.round(x / (this._l2));
        double xCenter = nxCenter * this._l2;
        if (Math.abs(x - xCenter) <= this._l6) {
            double nyCenter12;
            double yCenter12;
            if (nxCenter % 2 == 0) {
                nyCenter12 = Math.round(y / this._inverseSqrt3l);
                yCenter12 = nyCenter12 * this._inverseSqrt3l;
            } else {
                nyCenter12 = Math.round((y - this._inverseSqrt3l2) / this._inverseSqrt3l);
                yCenter12 = nyCenter12 * this._inverseSqrt3l + this._inverseSqrt3l2;
            }
            return this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, yCenter12);
        }
        else {
            double nyCenter1 = Math.round(y / this._inverseSqrt3l);
            double yCenter1 = nyCenter1 * this._inverseSqrt3l;
            double nyCenter2 = Math.round((y - this._inverseSqrt3l2) / this._inverseSqrt3l);
            double yCenter2 = nyCenter2 * this._inverseSqrt3l + this._inverseSqrt3l2;
            FaceCoordinates cCandidate1 = this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, (nxCenter % 2 == 0) ? yCenter1 : yCenter2);
            FaceCoordinates cCandidate2 = this._faceCoordinatesSwapByResolution(c.getFace(), (x > xCenter) ? xCenter + this._l2 : xCenter - this._l2, (nxCenter % 2 != 0) ? yCenter1 : yCenter2);
            return (c.distanceTo(cCandidate1) < c.distanceTo(cCandidate2)) ? cCandidate1 : cCandidate2;
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
     * Returns cells that are inside the bounds, or at least very near. Note that, in fact, all cells are included,
     * whose center points are less than
     *
     * @param lat0
     * @param lat1
     * @param lon0
     * @param lon1
     * @return cells inside the bounds
     */
    public Collection<GridCell> cellsForBound(double lat0, double lat1, double lon0, double lon1) throws Exception {
        Set<GridCell> cells = new HashSet<>();
        for (int f = 0; f < this._projection.numberOfFaces(); f++) cells.addAll(this._cellsForBound(f, lat0, lat1, lon0, lon1));
        return cells;
    }

    private Tuple<Integer, Integer> _integerForFaceCoordinates(FaceCoordinates c) {
        double x = (this._coordinatesNotSwapped()) ? c.getX() : c.getY();
        double y = (this._coordinatesNotSwapped()) ? c.getY() : c.getX();
        int nx = (int)Math.round(x / this._l2);
        int ny = (int)Math.round(y / this._inverseSqrt3l - ((nx % 2 == 0) ? 0 : .5));
        return new Tuple(nx, ny);
    }

    private FaceCoordinates _cellCoordinatesForLocationAndFace(int face, GeoCoordinates c) throws Exception {
        return this.cellForLocation(this._projection.sphereToPlanesOfTheFacesOfTheIcosahedron(face, c));
    }

    private Double _bufferEstimator(int face, int nx, int ny) throws Exception {
        GeoCoordinates x = this._projection.icosahedronToSphere(this._getCoordinatesOfCenter(face, nx, ny));
        GeoCoordinates y0 = this._projection.icosahedronToSphere(this._getCoordinatesOfCenter(face, nx + 1, ny));
        GeoCoordinates y1 = this._projection.icosahedronToSphere(this._getCoordinatesOfCenter(face, nx, ny + 1));
        double d0 = Math.pow(x.getLat() - y0.getLat(), 2) + Math.pow(x.getLon() - y0.getLon(), 2);
        double d1 = Math.pow(x.getLat() - y1.getLat(), 2) + Math.pow(x.getLon() - y1.getLon(), 2);
        double d = (d0 > d1) ? Math.sqrt(d0) : Math.sqrt(d1);
        return 2 * d;
    }

    private Collection<GridCell> _cellsForBound(int face, double lat0, double lat1, double lon0, double lon1) throws Exception {
        Set<GridCell> cells = new HashSet<>();

        // bounding box of face (triangle)
        double size1 = this._l0 / 3.;
        double size2 = this._l0 * this._inverseSqrt3 / 2.;
        double sizeX = this._coordinatesNotSwapped() ? size1 : size2;
        double sizeY = this._coordinatesNotSwapped() ? size2 : size1;

        // coordinates for vertices of bbox
        FaceCoordinates fc1 = this._cellCoordinatesForLocationAndFace(face, new GeoCoordinates(lat0, lon0));
        FaceCoordinates fc2 = this._cellCoordinatesForLocationAndFace(face, new GeoCoordinates(lat0, lon1));
        FaceCoordinates fc3 = this._cellCoordinatesForLocationAndFace(face, new GeoCoordinates(lat1, lon0));
        FaceCoordinates fc4 = this._cellCoordinatesForLocationAndFace(face, new GeoCoordinates(lat1, lon1));

        // find minimum and maximum values
        double xMin = Math.min(fc1.getX(), Math.min(fc2.getX(), Math.min(fc3.getX(), fc4.getX())));
        double yMin = Math.min(fc1.getY(), Math.min(fc2.getY(), Math.min(fc3.getY(), fc4.getY())));
        double xMax = Math.max(fc1.getX(), Math.max(fc2.getX(), Math.max(fc3.getX(), fc4.getX())));
        double yMax = Math.max(fc1.getY(), Math.max(fc2.getY(), Math.max(fc3.getY(), fc4.getY())));

        // check whether bbox intersects face
        double buffer = this._l;
        if (xMin - buffer > sizeX || xMax + buffer < -sizeX || yMin - buffer > sizeY || yMax + buffer < -sizeY) return cells;

        // compute cells
        Tuple<Integer, Integer> fcMinN = this._integerForFaceCoordinates(this.cellForLocation(new FaceCoordinates(face, xMin, yMin)));
        Tuple<Integer, Integer> fcMaxN = this._integerForFaceCoordinates(this.cellForLocation(new FaceCoordinates(face, xMax, yMax)));
        Double buffer2 = this._bufferEstimator(face, Math.round((fcMaxN._1 - fcMinN._1 + 1) / 2), Math.round((fcMaxN._2 - fcMinN._2 + 1) / 2));
        for (int nx = fcMinN._1 - 1; nx <= fcMaxN._1 + 1; nx++) {
            for (int ny = fcMinN._2 - 1; ny <= fcMaxN._2 + 1; ny++) {
                FaceCoordinates fc = this._getCoordinatesOfCenter(face, nx, ny);
                if (this._isCoordinatesInFace(fc)) {
                    GeoCoordinates gc = this._projection.icosahedronToSphere(fc);
                    if (lat0 - buffer2 <= gc.getLat() && gc.getLat() <= lat1 + buffer2 && lon0 - buffer2 <= gc.getLon() && gc.getLon() <= lon1 + buffer2) cells.add(this._newGridCell(gc, fc));
                }
            }
        }

        return cells;
    }

    private GridCell _newGridCell(GeoCoordinates gc, FaceCoordinates fc) throws Exception {
        int d = this._projection.faceOrientation(fc);
        boolean isPentagon = (Math.abs(Math.abs(fc.getX()) - this._triangleA) < this._precision && Math.abs(fc.getY() + d * this._triangleC) < this._precision) || (Math.abs(fc.getX()) < this._precision && Math.abs(fc.getY() - d * this._triangleB) < this._precision);
        return new GridCell(this._resolution, gc, isPentagon, fc.getFace());
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
        GeoCoordinates southWest2 = this._projection.icosahedronToSphere(new FaceCoordinates(southWest.getFace(), southWest.getX() - this.diameterOfCellOnIcosahedron() / 2, southWest.getY() - this.diameterOfCellOnIcosahedron() / 2));
        GeoCoordinates northEast2 = this._projection.icosahedronToSphere(new FaceCoordinates(northEast.getFace(), northEast.getX() - this.diameterOfCellOnIcosahedron() / 2, northEast.getY() - this.diameterOfCellOnIcosahedron() / 2));
        List<Double> l = new ArrayList<>();
        l.add(Math.abs(southWest2.getLat() - lat0));
        l.add(Math.abs(southWest2.getLon() - lon0));
        l.add(Math.abs(northEast2.getLat() - lat1));
        l.add(Math.abs(northEast2.getLon() - lon1));
        return l.stream().max(Double::compareTo).get();
    }

    /**
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

    private boolean _isCoordinatesInFace(FaceCoordinates fc) {
        double x = this._coordinatesNotSwapped() ? fc.getX() : fc.getY();
        double y = this._coordinatesNotSwapped() ? fc.getY() : fc.getX();

        // cell orientation
        int d = this._projection.faceOrientation(fc);

        // test whether coordinate is left of the triangle, right of the triangle, or below the triangle
        if (d * y > x * this._triangleBCA + this._triangleB) return false;
        if (d * y > -x * this._triangleBCA + this._triangleB) return false;
        if (d * y < -this._triangleC) return false;

        return true;
    }

    private FaceCoordinates _faceCoordinatesSwapByResolution(int face, double x, double y) {
        return new FaceCoordinates(face, this._coordinatesNotSwapped() ? x : y, this._coordinatesNotSwapped() ? y : x);
    }

    private boolean _coordinatesNotSwapped() {
        return this._resolution % 2 == 0;
    }
}
