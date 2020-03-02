package org.giscience.utils.geogrid.examples;

import org.giscience.utils.geogrid.cells.GridCellIDType;
import org.giscience.utils.geogrid.cells.GridCellMetaData;

/**
 * Computes the number of decimal places needed at a given resolution, depending on the type of cell ID.
 *
 * @author Franz-Benjamin Mocnik
 */
public class AdaptiveIDs {
    public static void main(String[] args) throws Exception {
        for (int resolution = 1; resolution <= 22; resolution++) {
            System.out.println("====== " + resolution);
            System.out.println(GridCellMetaData.getInstance().numberOfDecimalPlaces(resolution, GridCellIDType.NON_ADAPTIVE));
            System.out.println(GridCellMetaData.getInstance().numberOfDecimalPlaces(resolution, GridCellIDType.ADAPTIVE_UNIQUE));
            System.out.println(GridCellMetaData.getInstance().numberOfDecimalPlaces(resolution, GridCellIDType.ADAPTIVE_1_PERCENT));
        }
    }
}
