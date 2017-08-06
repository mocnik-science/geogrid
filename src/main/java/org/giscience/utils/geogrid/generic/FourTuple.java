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
package org.giscience.utils.geogrid.generic;

/**
 *
 * @author Franz-Benjamin Mocnik
 */
public class FourTuple<T1, T2, T3, T4> {
    public final T1 _1;
    public final T2 _2;
    public final T3 _3;
    public final T4 _4;

    public FourTuple(T1 arg1, T2 arg2, T3 arg3, T4 arg4){
        super();
        this._1 = arg1;
        this._2 = arg2;
        this._3 = arg3;
        this._4 = arg4;
    }

    @Override
    public String toString() {
        return String.format("(%s, %s, %s, %s)", this._1, this._2, this._3, this._4);
    }
}
