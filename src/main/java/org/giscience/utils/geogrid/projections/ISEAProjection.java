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
package org.giscience.utils.geogrid.projections;

import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.geometry.GeoCoordinates;
import org.giscience.utils.geogrid.geometry.FaceCoordinates;
import org.giscience.utils.geogrid.generic.Trigonometric;

/**
 * Icosahedron Snyder equal-area (ISEA) projection
 *
 * The ISEA projection is a projects a sphere on the icosahedron. Thereby the size of areas mapped to the icosahedron
 * are preserved. Angles and distances are however slightly distorted. The angular distortion is below 17.27 degree, and
 * the scale variation is less than 16.3 per cent.
 *
 * The projection has been proposed and has been described in detail by:
 *
 * John P. Snyder: An equal-area map projection for polyhedral globes. Cartographica, 29(1), 10–21, 1992.
 * doi:10.3138/27H7-8K88-4882-1752
 *
 * Another description and improvements can be found in:
 *
 * Erika Harrison, Ali Mahdavi-Amiri, and Faramarz Samavati: Optimization of inverse Snyder polyhedral projection.
 * International Conference on Cyberworlds 2011. doi:10.1109/CW.2011.36
 *
 * Erika Harrison, Ali Mahdavi-Amiri, and Faramarz Samavati: Analysis of inverse Snyder optimizations.
 * In: Marina L. Gavrilova, and C. J. Kenneth Tan (Eds): Transactions on Computational Science XVI. Heidelberg,
 * Springer, 2012. pp. 134–148. doi:10.1007/978-3-642-32663-9_8
 *
 * @author Franz-Benjamin Mocnik
 */
public class ISEAProjection {
    // radius
    private final double _RR_earth = Math.sqrt((this._G - this._theta) * Math.PI / (90 * Trigonometric.sin(this._theta) * Trigonometric.cos(this._theta))) / Trigonometric.tan(this._g); // R' / R
    private final double _R = this._RR_earth * this._R_earth; // R'
    private final double _R_earth = WGS84.radiusAuthalic; // R // authalic sphere radius for WGS84 [km]
    // faces
    private final int _numberOfFaces = 20;
    // orientation
    private double _orientationLat = 0;
    private double _orientationLon = 0;
    // spherical constants
    private final double _g = 37.37736814; // g
    private final int _G = 36; // G
    private final int _theta = 30; // \theta
    // distortion
    private final double _omega = 17.27; // \omega
    private final double _a = 1.163; // a
    private final double _b = .860; // b
    // face constants
    private final double __E = 90 - this._g; // E
    private final double __F = 10.81231696; // F
    private final double __G = this._R * Trigonometric.tan(this._g) * Math.sqrt(3) / 2.; // G // this value incorporates R', and not R, as is stated wrongly in the paper by Snyder
    private final int __X = 36; // half the difference in latitude between two horizontally adjacent faces
    private final double[] __lats = new double[20];
    private final int[] __lons = new int[20];
    // precision
    private final double _precision = 1e-9;
    // computed values
    private final double _2R = 2 * this._R; // 2 R'
    private final double __EF = this.__E - this.__F; // E - F
    private final int _AzMax = 2 * (90 - this._theta); // 2 (90 - \theta)
    private final double _tan_g = Trigonometric.tan(this._g); // \tan g
    private final double _cosG = Trigonometric.cos(this._G); // \cos G
    private final double _cotTheta = Trigonometric.cot(this._theta); // \cot \theta
    private final double _2cotTheta = 2 * this._cotTheta; // 2 \cot \theta
    private final double _pi_R_earth2_180 = Math.PI * Math.pow(this._R_earth, 2) / 180; // \pi R^2 / 180
    private final double _R_tan_g = this._R * this._tan_g; // R' \tan g
    private final double _R_tan_g_2 = Math.pow(this._R_tan_g, 2); // R'^2 * \tan^2 g
    private final double _sinG_cos_g = Trigonometric.sin(this._G) * Trigonometric.cos(this._g); // \sin G \cos g
    private final double _G_180 = this._G - 180; // G - 180

    public ISEAProjection() {
        this.__lats[0] = this.__E;
        this.__lats[1] = this.__E;
        this.__lats[2] = this.__E;
        this.__lats[3] = this.__E;
        this.__lats[4] = this.__E;
        this.__lats[5] = this.__F;
        this.__lats[6] = this.__F;
        this.__lats[7] = this.__F;
        this.__lats[8] = this.__F;
        this.__lats[9] = this.__F;
        this.__lats[10] = -this.__F;
        this.__lats[11] = -this.__F;
        this.__lats[12] = -this.__F;
        this.__lats[13] = -this.__F;
        this.__lats[14] = -this.__F;
        this.__lats[15] = -this.__E;
        this.__lats[16] = -this.__E;
        this.__lats[17] = -this.__E;
        this.__lats[18] = -this.__E;
        this.__lats[19] = -this.__E;
        this.__lons[0] = -4 * this.__X;
        this.__lons[1] = -2 * this.__X;
        this.__lons[2] = 0;
        this.__lons[3] = 2 * this.__X;
        this.__lons[4] = 4 * this.__X;
        this.__lons[5] = -4 * this.__X;
        this.__lons[6] = -2 * this.__X;
        this.__lons[7] = 0;
        this.__lons[8] = 2 * this.__X;
        this.__lons[9] = 4 * this.__X;
        this.__lons[10] = -3 * this.__X;
        this.__lons[11] = -this.__X;
        this.__lons[12] = this.__X;
        this.__lons[13] = 3 * this.__X;
        this.__lons[14] = 5 * this.__X;
        this.__lons[15] = -3 * this.__X;
        this.__lons[16] = -this.__X;
        this.__lons[17] = this.__X;
        this.__lons[18] = 3 * this.__X;
        this.__lons[19] = 5 * this.__X;
    }

    /**
     * @return number of faces of the icosahedron
     */
    public int numberOfFaces() {
        return this._numberOfFaces;
    }

    /**
     * Sets the orientation of the icosahedron.
     *
     * One corner of the icosahedron is, by default, facing to the north pole, and one to the south pole. The provided
     * orientation is relative to the default orientation.
     *
     * The orientation shifts every location by the angle <code>orientationLon</code> in direction of positive
     * longitude, and thereafter by the angle <code>orientationLat</code> in direction of positive latitude.
     *
     * @param orientationLat
     * @param orientationLon
     */
    public void setOrientation(double orientationLat, double orientationLon) {
        this._orientationLat = orientationLat;
        this._orientationLon = orientationLon;
    }

    /**
     * Sets the orientation of the icosahedron such that the north and the south poles are mapped to the edge midpoints
     * of the icosahedron. The equator is thus mapped symmetrically.
     */
    public void setOrientationSymmetricEquator() {
        this.setOrientation((90 + this.__E) / 2., 36.);
    }

    private GeoCoordinates _changeOrientation(GeoCoordinates c) throws Exception {
        double sinOrientationLat = Trigonometric.sin(this._orientationLat);
        double cosOrientationLat = Trigonometric.cos(this._orientationLat);
        double sinLat1 = Trigonometric.sin(c.getLat());
        double cosLat1 = Trigonometric.cos(c.getLat());
        double lon1 = c.getLon() + this._orientationLon;
        double sinLon1 = Trigonometric.sin(lon1);
        double cosLon1 = Trigonometric.cos(lon1);
        double lat2 = Trigonometric.asin(sinLat1 * cosOrientationLat + cosLon1 * cosLat1 * sinOrientationLat);
        double lon2 = Trigonometric.atan2(sinLon1 * cosLat1, cosLon1 * cosLat1 * cosOrientationLat - sinLat1 * sinOrientationLat);
        return new GeoCoordinates(lat2, lon2);
    }

    private GeoCoordinates _revertOrientation(GeoCoordinates c) throws Exception {
        double sinOrientationLat = Trigonometric.sin(-this._orientationLat);
        double cosOrientationLat = Trigonometric.cos(this._orientationLat);
        double sinLat1 = Trigonometric.sin(c.getLat());
        double cosLat1 = Trigonometric.cos(c.getLat());
        double sinLon1 = Trigonometric.sin(c.getLon());
        double cosLon1 = Trigonometric.cos(c.getLon());
        double lat2 = Trigonometric.asin(sinLat1 * cosOrientationLat + cosLon1 * cosLat1 * sinOrientationLat);
        double lon2 = Trigonometric.atan2(sinLon1 * cosLat1, cosLon1 * cosLat1 * cosOrientationLat - sinLat1 * sinOrientationLat);
        return new GeoCoordinates(lat2, lon2 - this._orientationLon);
    }

    /**
     * The projection distorts angles. This method returns the maximum angular distortion.
     *
     * @return maximum angular distortion
     */
    public double maximumAngularDistortion() {
        return this._omega;
    }

    /**
     * The projection distorts scale. This method returns the maximum scale variation.
     *
     * @return maximum scale variation
     */
    public double maximumScaleVariation() {
        return this._a;
    }

    /**
     * The projection distorts scale. This method returns the minimum scale variation.
     *
     * @return minimum scale variation
     */
    public double miniumScaleVariation() {
        return this._b;
    }

    /**
     * Returns the length of the bases of the triangles of the icosahedron
     *
     * @return length of the bases of the triangles
     */
    public double lengthOfTriangleBase() {
        return 2 * this.__G;
    }

    /**
     * Converts geographic coordinates to coordinates on the icosahedron.
     *
     * @param c geographic coordinates
     * @return coordinates on the icosahedron
     * @throws Exception
     */
    public FaceCoordinates sphereToIcosahedron(GeoCoordinates c) throws Exception {
        c = this._changeOrientation(c);
        Face face = this._sphereToFace(c);
        return this.sphereToFace(c, face);
    }

    /**
     * Converts geographic coordinates to the coordinate system of a plane which contains a face of the icosahedron.
     *
     * <b>CAUTION</b>: This function should only be used very rarely.  It generates coordinates which may not lay on the
     * faces themselves but only on the planes of the faces.  In case the given face equals the one that the given
     * coordinates would be projected to, the coordinates coincide.
     *
     * This function is useful for computing the overlap of bounding boxes or areas and a given face.
     *
     * @param c geographic coordinates
     * @return coordinates on the plane of the faces of the icosahedron
     * @throws Exception
     */
    public FaceCoordinates sphereToPlanesOfTheFacesOfTheIcosahedron(int face, GeoCoordinates c) throws Exception {
        c = this._changeOrientation(c);
        return this.sphereToFace(c, new Face(face, this, c));
    }

    /**
     * Converts geographic coordinates to coordinates on the icosahedron.
     *
     * @param c coordinates on the icosahedron
     * @return geographic coordinates
     * @throws Exception
     */
    public GeoCoordinates icosahedronToSphere(FaceCoordinates c) throws Exception {
        GeoCoordinates c2 = this._faceToSphere(c);
        return this._revertOrientation(c2);
    }

    private FaceCoordinates sphereToFace(GeoCoordinates c, Face face) {
        double Az_earth = Trigonometric.atan2(face.cosLat() * face.sinLonLon0(), face.cosLat0() * face.sinLat() - face.sinLat0() * face.cosLat() * face.cosLonLon0()); // Az
        double AzAdjustment = (this.faceOrientation(face.getFace()) > 0) ? 0 : 180;
        Az_earth += AzAdjustment;
        while (Az_earth < 0) {
            AzAdjustment += this._AzMax;
            Az_earth += this._AzMax;
        }
        while (Az_earth > this._AzMax) {
            AzAdjustment -= this._AzMax;
            Az_earth -= this._AzMax;
        }
        double sinAz_earth = Trigonometric.sin(Az_earth); // \sin Az
        double cosAz_earth = Trigonometric.cos(Az_earth); // \cos Az
        double H = this._compute_H(sinAz_earth, cosAz_earth); // H
        double area = (Az_earth + this._G_180 + H) * this._pi_R_earth2_180; // A_G and A_{ABD}
        double Az = Trigonometric.atan2(2 * area, this._R_tan_g_2 - area * this._2cotTheta); // Az'
        double sinAz = Trigonometric.sin(Az); // \sin Az'
        double cosAz = Trigonometric.cos(Az); // \cos Az'
        double f = this._compute_f(sinAz, cosAz, sinAz_earth, cosAz_earth); // f
        double rho = this._2R * f * Trigonometric.sin(face.z() / 2.); // \rho
        Az -= AzAdjustment;
        double x = rho * Trigonometric.sin(Az); // x
        double y = rho * Trigonometric.cos(Az); // y
        return new FaceCoordinates(face.getFace(), x, y);
    }

    /**
     * Returns the corresponding face for given coordinates
     *
     * @param c coordinates
     * @return face that the coordinates belongs to
     */
    public int sphereToFace(GeoCoordinates c) {
        return this._sphereToFace(c).getFace();
    }

    private Face _sphereToFace(GeoCoordinates c) {
        if (c.getLat() > this.__EF) {
            if (c.getLon() < -108) return new Face(0, this, c);
            else if (c.getLon() < -36) return new Face(1, this, c);
            else if (c.getLon() < 36) return new Face(2, this, c);
            else if (c.getLon() < 108) return new Face(3, this, c);
            else return new Face(4, this, c);
        } else if (c.getLat() < -this.__EF) {
            if (c.getLon() < -144) return new Face(19, this, c);
            else if (c.getLon() < -72) return new Face(15, this, c);
            else if (c.getLon() < 0) return new Face(16, this, c);
            else if (c.getLon() < 72) return new Face(17, this, c);
            else if (c.getLon() < 144) return new Face(18, this, c);
            else return new Face(19, this, c);
        } else {
            if (c.getLon() < -144) return this._sphereToFaceTestFaces(c, 5, 14);
            else if (c.getLon() < -108) return this._sphereToFaceTestFaces(c, 5, 10);
            else if (c.getLon() < -72) return this._sphereToFaceTestFaces(c, 6, 10);
            else if (c.getLon() < -36) return this._sphereToFaceTestFaces(c, 6, 11);
            else if (c.getLon() < 0) return this._sphereToFaceTestFaces(c, 7, 11);
            else if (c.getLon() < 36) return this._sphereToFaceTestFaces(c, 7, 12);
            else if (c.getLon() < 72) return this._sphereToFaceTestFaces(c, 8, 12);
            else if (c.getLon() < 108) return this._sphereToFaceTestFaces(c, 8, 13);
            else if (c.getLon() < 144) return this._sphereToFaceTestFaces(c, 9, 13);
            return this._sphereToFaceTestFaces(c, 9, 14);
        }
    }

    private GeoCoordinates _faceToSphere(FaceCoordinates c) throws Exception {
        double Az = Trigonometric.atan2(c.getX(), c.getY()); // Az'
        double rho = Math.sqrt(Math.pow(c.getX(), 2) + Math.pow(c.getY(), 2)); // \rho
        double AzAdjustment = (this.faceOrientation(c) > 0) ? 0 : 180;
        Az += AzAdjustment;
        while (Az < 0) {
            AzAdjustment += this._AzMax;
            Az += this._AzMax;
        }
        while (Az > this._AzMax) {
            AzAdjustment -= this._AzMax;
            Az -= this._AzMax;
        }
        double sinAz = Trigonometric.sin(Az); // \sin Az'
        double cosAz = Trigonometric.cos(Az); // \cos Az'
        double cotAz = cosAz / sinAz; // \cot Az'
        double area = this._R_tan_g_2 / (2 * (cotAz + this._cotTheta)); // A_G or A_{ABD}
        double deltaAz = 10 * this._precision;
        double area_pi_R_earth2_180_G_180 = area / this._pi_R_earth2_180 - this._G_180;
        double Az_earth = Az;
        while (Math.abs(deltaAz) > this._precision) {
            double H = this._compute_H(Trigonometric.sin(Az_earth), Trigonometric.cos(Az_earth)); // H
            double FAz_earth = area_pi_R_earth2_180_G_180 - H - Az_earth; // F(Az) or g(Az)
            double F2Az_earth = (Trigonometric.cos(Az_earth) * this._sinG_cos_g + Trigonometric.sin(Az_earth) * this._cosG) / Trigonometric.sin(H) - 1; // F'(Az) or g'(Az)
            deltaAz = - FAz_earth / F2Az_earth; // \Delta Az^0 or \Delta Az
            Az_earth += deltaAz;
        }
        double sinAz_earth = Trigonometric.sin(Az_earth); // \sin Az
        double cosAz_earth = Trigonometric.cos(Az_earth); // \cos Az
        double f = this._compute_f(sinAz, cosAz, sinAz_earth, cosAz_earth); // f
        double z = 2 * Trigonometric.asin(rho / (this._2R * f)); // z
        Az_earth -= AzAdjustment;
        double sinLat0 = Trigonometric.sin(this._getLat(c)); // \sin \phi_0
        double cosLat0 = Trigonometric.cos(this._getLat(c)); // \cos \phi_0
        double sinZ = Trigonometric.sin(z); // \sin z
        double cosZ = Trigonometric.cos(z); // \cos z
        double lat = Trigonometric.asin(sinLat0 * cosZ + cosLat0 * sinZ * Trigonometric.cos(Az_earth)); // \phi
        double lon = this._getLon(c) + Trigonometric.atan2(Trigonometric.sin(Az_earth) * sinZ * cosLat0, cosZ - sinLat0 * Trigonometric.sin(lat)); // \lambda
        return new GeoCoordinates(lat, lon);
    }

    /**
     * Tests whether coordinates c belong to face1 or face2. Assumes that the coordinates c belongs to one of these
     * faces.
     *
     * @param c
     * @param face1
     * @param face2
     * @return
     */
    private Face _sphereToFaceTestFaces(GeoCoordinates c, int face1, int face2) {
        Face f1 = new Face(face1, this, c);
        Face f2 = new Face(face2, this, c);
        return (f1.z() <= f2.z()) ? f1 : f2;
    }

    /**
     * Returns orientation of a face.
     *
     * @param fc
     * @return 1 for upright, and -1 for upside down
     */
    public int faceOrientation(FaceCoordinates fc) {
        return this.faceOrientation(fc.getFace());
    }

    /**
     * Returns orientation of a face.
     *
     * @param face
     * @return 1 for upright, and -1 for upside down
     */
    public int faceOrientation(int face) {
        return (face <= 4 || (10 <= face && face <= 14)) ? 1 : -1;
    }

    private double _compute_H(double sinAz_earth, double cosAz_earth) {
        return Trigonometric.acos(sinAz_earth * this._sinG_cos_g - cosAz_earth * this._cosG); // H
    }
    private double _compute_f(double sinAz, double cosAz, double sinAz_earth, double cosAz_earth) {
        return this._compute_d(sinAz, cosAz) / (2 * this._R * Trigonometric.sin(this._compute_q(sinAz_earth, cosAz_earth) / 2)); // f
    }
    private double _compute_d(double sinAz, double cosAz) {
        return this._R_tan_g / (cosAz + sinAz * this._cotTheta); // d'
    }
    private double _compute_q(double sinAz_earth, double cosAz_earth) {
        return Trigonometric.atan2(this._tan_g, (cosAz_earth + sinAz_earth * this._cotTheta)); // q
    }

    private double _getLat(Face f) {
        return this.getLat(f.getFace());
    }
    private double _getLon(Face f) {
        return this.getLon(f.getFace());
    }
    private double _getLat(FaceCoordinates c) {
        return this.getLat(c.getFace());
    }
    private double _getLon(FaceCoordinates c) {
        return this.getLon(c.getFace());
    }

    /**
     * Returns latitude for face
     *
     * @param face
     * @return latitude
     */
    public double getLat(int face) {
        return this.__lats[face];
    }

    /**
     * Returns longitude for face
     *
     * @param face
     * @return longitude
     */
    public double getLon(int face) {
        return this.__lons[face];
    }

    /**
     * Returns minimium latitude of the given face
     *
     * @param face
     * @return
     */
    public double getLonMin(int face) {
        double lon = this.getLon(face) - this.__X;
        if (lon < -180) lon += 360;
        return lon;
    }

    /**
     * Returns maximum latitude of the given face
     *
     * @param face
     * @return
     */
    public double getLonMax(int face) {
        double lon = this.getLon(face) + this.__X;
        if (lon > 180) lon -= 360;
        return lon;
    }

    private class Face {
        private final int _face;
        private final ISEAProjection _iseaProjection;
        private final GeoCoordinates _c;
        private Double _z = null;
        private Double _sinLat = null;
        private Double _cosLat = null;
        private Double _sinLat0 = null;
        private Double _cosLat0 = null;
        private Double _sinLonLon0 = null;
        private Double _cosLonLon0 = null;

        public Face(int face, ISEAProjection iseaProjection, GeoCoordinates c) {
            this._face = face;
            this._iseaProjection = iseaProjection;
            this._c = c;
        }

        public int getFace() {
            return this._face;
        }

        protected double z() {
            if (this._z == null) this._z = Trigonometric.acos(this.sinLat0() * this.sinLat() + this.cosLat0() * this.cosLat() * this.cosLonLon0());
            return this._z;
        }
        protected double sinLat() {
            if (this._sinLat == null) this._sinLat = Trigonometric.sin(this._c.getLat());
            return this._sinLat;
        }
        protected double cosLat() {
            if (this._cosLat == null) this._cosLat = Trigonometric.cos(this._c.getLat());
            return this._cosLat;
        }
        protected double sinLat0() {
            if (this._sinLat0 == null) this._sinLat0 = Trigonometric.sin(this._iseaProjection._getLat(this));
            return this._sinLat0;
        }
        protected double cosLat0() {
            if (this._cosLat0 == null) this._cosLat0 = Trigonometric.cos(this._iseaProjection._getLat(this));
            return this._cosLat0;
        }
        protected double sinLonLon0() {
            if (this._sinLonLon0 == null) this._sinLonLon0 = Trigonometric.sin(this._c.getLon() - this._iseaProjection._getLon(this));
            return this._sinLonLon0;
        }
        protected double cosLonLon0() {
            if (this._cosLonLon0 == null) this._cosLonLon0 = Trigonometric.cos(this._c.getLon() - this._iseaProjection._getLon(this));
            return this._cosLonLon0;
        }
    }
}
