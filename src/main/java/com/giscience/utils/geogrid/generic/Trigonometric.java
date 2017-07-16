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
package com.giscience.utils.geogrid.generic;

/**
 * Trigonometric functions with degree instead of radians
 * 
 * @author Franz-Benjamin Mocnik
 */
public class Trigonometric {
    public static double sin(double x) {
        return Math.sin(Math.toRadians(x));
    }
    public static double cos(double x) {
        return Math.cos(Math.toRadians(x));
    }
    public static double tan(double x) {
        return Math.tan(Math.toRadians(x));
    }
    public static double cot(double x) {
        return 1 / Trigonometric.tan(x);
    }
    public static double asin(double x) {
        return Math.toDegrees(Math.asin(x));
    }
    public static double acos(double x) {
        return Math.toDegrees(Math.acos(x));
    }
    public static double atan(double x) {
        return Math.toDegrees(Math.atan(x));
    }
    public static double atan2(double x, double y) {
        return Math.toDegrees(Math.atan2(x, y));
    }
}
