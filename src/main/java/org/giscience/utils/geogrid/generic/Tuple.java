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
