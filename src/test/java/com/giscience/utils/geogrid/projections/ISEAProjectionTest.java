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
package com.giscience.utils.geogrid.projections;

import com.giscience.utils.geogrid.coordinates.FaceCoordinates;
import com.giscience.utils.geogrid.coordinates.GeoCoordinates;
import static org.junit.Assert.assertTrue;

import com.giscience.utils.geogrid.generic.Trigonometric;
import org.junit.Test;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class ISEAProjectionTest {
    private final double _precision = 10E-10;
    private final int _iterations = 1000000;
    
    @Test
    public void projectionApartFromPoles() throws Exception {
        ISEAProjection p = new ISEAProjection();
        p.setOrientationSymmetricEquator();
        for (int i = 0; i < this._iterations; i++) {
            GeoCoordinates c = new GeoCoordinates(Math.random() * 179.99 - 89.995, Math.random() * 360);
            FaceCoordinates c2 = p.sphereToIcosahedron(c);
            GeoCoordinates c3 = p.icosahedronToSphere(c2);
            assertTrue(Math.abs(c3.getLat() - c.getLat()) < this._precision);
            assertTrue((Math.abs(c3.getLon() - c.getLon()) % 360) * Trigonometric.cos(c.getLat()) < this._precision);
        }
    }
}
