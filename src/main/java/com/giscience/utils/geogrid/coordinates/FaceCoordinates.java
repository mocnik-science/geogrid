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
package com.giscience.utils.geogrid.coordinates;

/**
 * Cartesian Coordinates of a location on a face of a platonic solid.
 * 
 * @author Franz-Benjamin Mocnik
 */
public class FaceCoordinates {
    private final Integer _face;
    private final Double _x;
    private final Double _y;
    
    public FaceCoordinates(Integer face, Double x, Double y) {
        this._face = face;
        this._x = x;
        this._y = y;
    }
    
    public Integer getFace() {
        return this._face;
    }
    
    public Double getX() {
        return this._x;
    }
    
    public Double getY() {
        return this._y;
    }
    
    @Override
    public String toString() {
        return String.format("face %d x %f y %f", this._face, this._x, this._y);
    }
}
