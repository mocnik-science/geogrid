package org.giscience.utils.geogrid.cells;

import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.grids.ISEA3H;

import java.util.Map;
import java.util.TreeMap;

/**
 * Meta data for a grid cell
 *
 * @author Franz-Benjamin Mocnik
 */public class GridCellMetaData {
    private static GridCellMetaData _gridCellMetaData = new GridCellMetaData();
    private static final int _maxNumberOfConsecutiveDigits = 6;
    private Map<Integer, Integer> _numberOfConsecutiveDigits = new TreeMap();

    private GridCellMetaData() {}

    /**
     * Get an instance of this (singleton) class
     *
     * @return singleton class
     */
    public static GridCellMetaData getInstance() {
        return GridCellMetaData._gridCellMetaData;
    }

    /**
     * Returns the number of consecutive digits to be used for the ID for a given resolution and a given ID type
     *
     * @param resolution
     * @param gridCellIDType
     * @return number of consecutive digits
     */
    public int numberOfConsecutiveDigits(int resolution, GridCellIDType gridCellIDType) {
        if (gridCellIDType == GridCellIDType.NON_ADAPTIVE) return this._maxNumberOfConsecutiveDigits;
        int nocd;
        if (this._numberOfConsecutiveDigits.containsKey(resolution)) nocd = this._numberOfConsecutiveDigits.get(resolution);
        else {
            double distBetweenCells = 2 * (new ISEA3H(resolution)).lowerBoundForLengthOfASideOfHexagonalCellOnSphere();
            nocd = (int)Math.ceil(-Math.log10(distBetweenCells / (2 * Math.PI * WGS84.radiusAuthalic / 360)));
            this._numberOfConsecutiveDigits.put(resolution, nocd);
        }
        switch (gridCellIDType) {
            case ADAPTIVE_UNIQUE:
                return Math.max(Math.min(nocd, this._maxNumberOfConsecutiveDigits), 0);
            case ADAPTIVE_1_PERCENT:
                return Math.max(Math.min(nocd + 2, this._maxNumberOfConsecutiveDigits), 0);
            default:
                return this._maxNumberOfConsecutiveDigits;
        }
    }
}
