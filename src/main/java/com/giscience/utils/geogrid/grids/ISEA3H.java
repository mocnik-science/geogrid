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
package com.giscience.utils.geogrid.grids;

import com.giscience.utils.geogrid.geometry.FaceCoordinates;
import com.giscience.utils.geogrid.geometry.GeoCoordinates;
import com.giscience.utils.geogrid.geometry.GridCell;
import com.giscience.utils.geogrid.projections.ISEAProjection;

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
 * @author Franz-Benjamin Mocnik
 */
public class ISEA3H {
    private final int _resolution;
    private final double _l0; // length of the triangle base at resolution 0
    private final double _l; // length of the triangle base at the given resolution
    private final double _l2; // l / 2
    private final double _l3; // l / 3
    private final double _l6; // l / 6
    private final double _l23; // l * 2 / 3
    private final ISEAProjection _projection = new ISEAProjection();
    private final double _inverseSqrt3 = 1 / Math.sqrt(3);
    private final double _inverseSqrt3l;
    private final double _inverseSqrt3l2;

    public ISEA3H(int resolution) {
        this._resolution = resolution;
        this._projection.setOrientationSymmetricEquator();
        this._l0 = this._projection.getLengthOfTriangleBase();
        this._l = Math.pow(this._inverseSqrt3, this._resolution) * this._l0;
        this._l2 = this._l / 2.;
        this._l3 = this._l / 3.;
        this._l6 = this._l / 6.;
        this._l23 = this._l * 2 / 3.;
        this._inverseSqrt3l = this._inverseSqrt3 * this._l;
        this._inverseSqrt3l2 = this._inverseSqrt3l / 2.;
    }

    /**
     * Returns the diameter of a cell
     * 
     * @return diameter of a cell
     */
    public double getDiameterOfCell() {
        return this._l23;
    }

    /**
     * Returns the area of a cell. The cells should all have the same area by construction, because the ISEA projection
     * is equal-area.
     * 
     * @return area of a cell
     */
    public double getAreaOfCell() {
        return this._l3;
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
        return new GridCell(this._resolution, this._projection.icosahedronToSphere(this.cellForLocation(this._projection.sphereToIcosahedron(c))));
    }

    /**
     * Returns the coordinates of the center of the corresponding grid cell for given coordinates in the face
     * 
     * @param c face coordinates
     * @return face coordinates of the center of the corresponding grid cell
     * @throws Exception
     */
    public FaceCoordinates cellForLocation(FaceCoordinates c) throws Exception {
        double nxCenter = Math.round(c.getX() / (this._l2));
        double xCenter = nxCenter * this._l2;
        double nyCenter = Math.round(c.getY() / this._inverseSqrt3l);
        double yCenter = nyCenter * this._inverseSqrt3l;
        if (Math.abs(c.getX() - xCenter) <= this._l6) return this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, yCenter);
        if (Math.abs(c.getX() - xCenter) > this._l3) return this._faceCoordinatesSwapByResolution(c.getFace(), (c.getX() > xCenter) ? xCenter + this._l2 : xCenter - this._l2, (c.getY() > yCenter) ? yCenter + this._inverseSqrt3l2 : yCenter - this._inverseSqrt3l2);
        else {
            FaceCoordinates cCandidate1 = this._faceCoordinatesSwapByResolution(c.getFace(), xCenter, yCenter);
            FaceCoordinates cCandidate2 = this._faceCoordinatesSwapByResolution(c.getFace(), (c.getX() > xCenter) ? xCenter + this._l2 : xCenter - this._l2, (c.getY() > yCenter) ? yCenter + this._inverseSqrt3l2 : yCenter - this._inverseSqrt3l2);
            return (c.distanceTo(cCandidate1) < c.distanceTo(cCandidate2)) ? cCandidate1 : cCandidate2;
        }
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

        // check for out of scope

    }

    private FaceCoordinates _faceCoordinatesSwapByResolution(int face, double x, double y) {
        return new FaceCoordinates(face, (this._resolution % 2 == 0) ? x : y, (this._resolution % 2 == 0) ? y : x);
    }
}
