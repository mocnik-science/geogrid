package org.giscience.utils.geogrid.generic;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class Tuple<T1, T2> {
    public final T1 _1;
    public final T2 _2;

    public Tuple(T1 arg1, T2 arg2){
        super();
        this._1 = arg1;
        this._2 = arg2;
    }

    @Override
    public String toString() {
        return String.format("(%s, %s)", this._1, this._2);
    }
}
