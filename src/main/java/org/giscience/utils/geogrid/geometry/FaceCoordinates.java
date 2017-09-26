package org.giscience.utils.geogrid.geometry;

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

    public Double distanceTo(FaceCoordinates c) {
        return Math.sqrt(Math.pow(this._x - c.getX(), 2) + Math.pow(this._y - c.getY(), 2));
    }
    @Override
    public String toString() {
        return String.format("face %d x %f y %f", this._face, this._x, this._y);
    }
}
