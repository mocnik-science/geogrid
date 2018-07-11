# geogrid

The library `geogrid` provides methods to generate and handle the ISEA Aperture 3 Hexagon Discrete Global Grid System (ISEA3H DGGS) using the Inverse Snyder Equal-Area Projection (ISEA).

## Scientific Publications

An overview of the ISEA3H DGGS and the identifier used in this library can be found in:

* F-B Mocnik: [**A Novel Identifier Scheme for the ISEA Aperture 3 Hexagon Discrete Global Grid System.**](http://doi.org/10.1080/15230406.2018.1455157) Cartography and Geographic Information science, 2018

The ISEA3H has originally been proposed by:

* K Sahr, D White, and A J Kimerling: [**Geodesic Discrete Global Grid Systems**](http://doi.org/10.1559/152304003100011090). Cartography and Geographic Information Science, 30(2), 2003, 121–134. doi:10.1559/152304003100011090

The projection has been proposed and has been described in detail by:

* J P Snyder: [**An equal-area map projection for polyhedral globes**](http://doi.org/http://doi.org/10.3138/27H7-8K88-4882-1752). Cartographica, 29(1), 1992, 10–21. doi:10.3138/27H7-8K88-4882-1752

## Related Software

Data aggregated by a DGGS can be presented in the browser by using the JavaScript library [**geogrid.js**](https://github.com/giscience/geogrid.js), which visualizes such data as a layer to [Leaflet](http://leafletjs.com).  Data aggregation can even be performed automatically by the framework [**Measures REST**](https://github.com/giscience/measures-rest), which even provides the aggregated data by a REST interface.

## Short Overview

### ISEA Projection

The ISEA projection maps the regular icosahedron to the sphere.  A hexagonal grid introduced on the icosahedron can thus be mapped to the sphere as well. As the ISEA projection and the ISEA3H grid are thus strongly interlinked, this library offers an implementation of both, the ISEA projection and the ISEA3H grid.

A location on the icosahedron can be mapped to the sphere, as in the following example of the coordinate `(20.3, 12.5)` on face `2`:

```java
ISEAProjection projection = new ISEAProjection();
GeoCoordinates c = projection.icosahedronToSphere(new FaceCoordinates(2, 20.3, 12.5));
c.getLat(); // 52.739063
c.getLon(); // 0.315051
```

The resulting geographic coordinates can be projected back to the icosahedron:

```java
FaceCoordinates c2 = projection.sphereToIcosahedron(c);
c2.getFace(); // 2
c2.getX(); // 20.3
c2.getY(); // 12.5
```

### ISEA3H DGGS: Determining the Cell for a Given Location

When the grid cell containing a given location should be determined, the method `cellForLocation` can be used, either for a pair of coordinates, or an instance of `GeoCoordinates`:

```java
GridCell gc1 = grid.cellForLocation(41.5, 6.5);
GridCell gc2 = grid.cellForLocation(new GeoCoordinates(41.5, 6.5));
```

The cell containing the centroid of a `Geometry` (from the [JTS Topology Suite](https://github.com/locationtech/jts)) can be determined accordingly:

```java
grid.cellForCentroid(g);
```

### ISEA3H DGGS: Computing the Cells

The grid cells of the ISEA3H for a given boundary box can easily be computed as follows:

```java
ISEA3H grid = new ISEA3H(10);
Collection<GridCell> cells = grid.cellsForBound(41, 42, 6, 7);
```

Observe that the constructor of the class `ISEA3H` accepts the grid resolution as a parameter.  The method `cellsForBound` computes, in this case, all cells that are (about) contained in the area between `41` and `42`  degree latitude and `6` and `7` degree longitude.  If all cells should be computed, the method `cells` can be used instead:

```java
Collection<GridCell> cellsAll = grid.cells();
```

The resulting collection of cells contains instances of `GridCell`.  Such an instance contains information about the centroid of the grid cell, the grid resolution, whether it is a hexagon (or a pentagon) cell, and an identifier (ID):

```java
for (GridCell gc : cells) {
    gc.getLat();
    gc.getLon();
    gc.getResolution();
    gc.isPentagon();
    gc.getID();
}
```

Different types of cell IDs can be used:

| Type of ID | Description | Typical use |
| ---------- | ----------- | ----------- |
| `NON_ADAPTIVE` | *Non-adaptive IDs* include the geometry of the cell in a heigh precision. | Perfect rendering of the grid cells |
| `ADAPTIVE_1_PERCENT` | *Adaptive 1% IDs* include the geometry in a scale dependent but yet reasonable precision. | Rendering of the grid cells on a website |
| `ADAPTIVE_1_PERCENT` | *Adaptive unique IDs* aim at being as short as possible while preserving uniqueness. | Comparison of whether two cells are identical |

Further information can be found in:

tba

By default, non-adaptive IDs are returned.  The type of ID can though be inserted explicitly as follows:

```java
gc.getID(NON_ADAPTIVE);
gc.getID(ADAPTIVE_1_PERCENT);
gc.getID(ADAPTIVE_UNIQUE);
```

### ISEA3H DGGS: Computing Cell IDs only

In case that just the identifiers (IDs) of the cells should be computed, two more efficient methods can be used:

```java
ISEA3H grid = new ISEA3H(10);
Collection<Long> cells = grid.cellIDsForBound(41, 42, 6, 7);
Collection<Long> cells2 = grid.cellIDsForBound(41, 42, 6, 7, ADAPTIVE_UNIQUE);
Collection<Long> cellsAll = grid.cellIDs();
Collection<Long> cellsAll2 = grid.cellIDs(ADAPTIVE_UNIQUE);
```

The IDs can even be saved to the hard disc:

```java
ISEA3H grid = new ISEA3H(10);
grid.cellIDs("path/filename");
grid.cellIDs("path/filename2", ADAPTIVE_UNIQUE);
```

Each of the created files contains only IDs that refer to one face of the icosahedron, and the data is divided into a number of chunks.  The same ID can occur multiple times (for different faces) when it refers to a cell that is partially contained in different faces.  The files are named as `filename.face.chunk_number`.  The first part of this name refers to the filename provided in the method `cellsIDs`; the second part, to the face to which the IDs contained in the file belong to; and the third part consists of consecutive numbers to distinguish the chunks.

## ISEA Projection

### Properties

Many properties of the ISEA projection are provided by an instance of the `ISEAProjection`. The following table provides an overview over some of them:

| Function | Description |
| -------- | ----------- |
| `numberOfFaces()` | number of faces of the icosahedron |
| `maximumAngularDistortion()` | maximum angular distortion |
| `maximumScaleVariation()` | maximum scale variation |
| `miniumScaleVariation()` | minimum scale variation |
| `lengthOfTriangleBase()` | length of the bases of the triangles of the icosahedron |
| `sphericalDistanceFromCenterToVerticesOnSphere()` | spherical distance from center of a face to any of its vertices on the sphere; in degrees |

### Orientation

The ISEA projection maps from the icosahedron to the sphere.  By default, two of the vertices of the faces of the icosahedron are mapped to the poles.  Another orientation might, however, be useful in many cases, as, for example, in the case of the ISEA3H DGGS: distortions are especially large at the vertices of the faces, and one aims at mapping the vertices to places that are within the oceans.  Such an ‘ideal’ choice of an orientation can be reached by:

```java
ISEAProjection projection = new ISEAProjection();
projection.setOrientationSymmetricEquator();
```

If the orientation is to be set manually by shifting the orientation relative to the standard orientation (two of the vertices are mapped to the poles).  The following method shifts every location by the angle `orientationLon` in the direction of positive longitude, and thereafter by the angle `orientationLat` in direction of positive latitude:

```java
projection.setOrientation(orientationLat, orientationLon);
```

## ISEA3H DGGS

### Properties

Many properties of the ISEA3H DGGS are provided by an instance of the `ISEA3H`. The following table provides an overview over some of them:

| Function | Description |
| -------- | ----------- |
| `diameterOfHexagonalCellOnIcosahedron()` | diameter of a hexagonal cell on the icosahedron, in kilometres |
| `lengthOfASideOfHexagonalCellOnIcosahedron()` | length of a side of a hexagonal cell on the icosahedron, in kilometres |
| `lowerBoundForLengthOfASideOfHexagonalCellOnSphere()` | lower bound for the length of a side of a hexagonal cell on the sphere, in kilometres |
| `areaOfAHexagonalCell()` | area of a hexagonal cell; the cells have all the same area by construction, because the ISEA projection is equal-area, in square kilometres |
| `areaOfAPentagonalCell()` | area of a pentagonal cell; the cells have all the same area by construction, because the ISEA projection is equal-area, in square kilometres |
| `numberOfHexagonalCells()` | number of hexagonal cells |
| `numberOfPentagonalCells()` | number of pentagonal cells |

### Threads

Some of the computations can be executed in parallel.  The number of threads can be set as follows:

```java
ISEA3H grid = new ISEA3H(10);
grid.setNumberOfThreads(8);
```

By default, 8 threads are used.

## Author

This software is written and maintained by Franz-Benjamin Mocnik, <mocnik@uni-heidelberg.de>, GIScience Research Group, Institute of Geography, Heidelberg University.

The development has been supported by the DFG project *A framework for measuring the fitness for purpose of OpenStreetMap data based on intrinsic quality indicators* (FA 1189/3-1).

(c) by Heidelberg University, 2017–2018.

## License

The code is licensed under the [MIT license](https://github.com/giscience/geogrid/blob/master/LICENSE).
