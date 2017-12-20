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
    private static final int _maxNumberOfDecimalPlaces = 6;
    private Map<Integer, Integer> _numberOfDecimalPlaces = new TreeMap();

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
     * Returns the number of decimal places to be used for the ID for a given resolution and a given ID type
     *
     * @param resolution
     * @param gridCellIDType
     * @return number of decimal places
     */
    public int numberOfDecimalPlaces(int resolution, GridCellIDType gridCellIDType) {
        if (gridCellIDType == GridCellIDType.NON_ADAPTIVE) return this._maxNumberOfDecimalPlaces;
        int nodp;
        if (this._numberOfDecimalPlaces.containsKey(resolution)) nodp = this._numberOfDecimalPlaces.get(resolution);
        else {
            double distBetweenCells = 2 * (new ISEA3H(resolution)).lowerBoundForLengthOfASideOfHexagonalCellOnSphere();
            nodp = (int)Math.ceil(-Math.log10(distBetweenCells / (2 * Math.PI * WGS84.radiusAuthalic / 360)));
            this._numberOfDecimalPlaces.put(resolution, nodp);
        }
        switch (gridCellIDType) {
            case ADAPTIVE_UNIQUE:
                return Math.max(Math.min(nodp, this._maxNumberOfDecimalPlaces), 0);
            case ADAPTIVE_1_PERCENT:
                return Math.max(Math.min(nodp + 2, this._maxNumberOfDecimalPlaces), 0);
            default:
                return this._maxNumberOfDecimalPlaces;
        }
    }
}
