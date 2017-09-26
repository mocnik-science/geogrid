package org.giscience.utils.geogrid.generic;

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
    public static double atan2(double y, double x) {
        return Math.toDegrees(Math.atan2(y, x));
    }
}
