# ArcpyGlobeModel

Solution of exercise 2 of 2024 course of Mathematical Cartography. Creates a dodecahedral net globe model in gnomonic projection in oblique aspect. The net is created in an ArcGIS Pro project.

---
### Script requirements
- SW:
  - ArcGIS Pro must be installed, version 3.4.3 was used, corresponds to arcpy version 3.4. (runs on Python 3.10.10). Also compatible with ArcGIS Pro 3.5.1.
    - Used libraries present in arcpy: matplotlib 3.6.3, numpy: 1.24.3, scipy 1.9.3.
  - Libraries included in Python standard library (installed in Python 3.12.2): shutil 1.7.0, os (version unknown).
- Data:
  - All paths referenced below are relative to CWD the script was launched from (since 'paths relative to CWD' commit).
  - Continents data must be in `DataSmall`.
  - Input ArcGIS project must be present in `referenceProj`.

### Inputs
- Singlepart polygons (`DataSmall`) in WGS-84 used for _.SVG_ faces, 
- `referenceProj` contents: 
  - A map `face<1-12>` for each face. These maps don´t have the requested projection. They are only projected in North / South Pole Gnomonic projection.
  - An A3 portrait layout including: 
    - Map frames with names corresponding to map names, anchor is set to center point.
    - A transparent basemap citation object. It keeps the citations out of the faces. 
  - A symbology map with layers controling the look of boundary and graticule layers.

### Outputs 
- `faces` folder containing 12 _SVG_ face plots with graticule and continents. The faces must be rotated manually to asseble a dodecahedral net. 
- ArcGIS project `outputProj` including face maps and dodecahedral net layout.
- Exported layout `outputMap.pdf`.

### Globe modelling
Globe scale can be changed using `m` variable on row 7. Radius of a circle incribed to globe face is computed as `Earth radius / M`. Face edge length is based on this value. The script prints globe paramaters related to its size so correctness can be checked (see rows 229-234). 

Face map frame centers in the layout are scaled based on difference of current model face edge length and length of 92,68 mm. This edge length and face center locations are based on M = 100 000 000 (see rows 341-342, 243-245).

### Sources used
BAYER, T. (2024): Konstrukce glóbů na platónských tělesech, návod na cvičení. Přírodovědecká fakulta
Univerzity Karlovy (14. 4. 2024).

ESRI (2022): World Continents, Esri. https://hub.arcgis.com/datasets/esri::world-continents/about (14. 4. 2024).

ESRI (2025): OpenStreetMap Vector Basemap. https://www.arcgis.com/home/item.html?id=d6b18a2e774c4959ba855f6ac90952a2 (19. 6. 2025).

ESRI, TOMTOM, GARMIN, FAO, NOAA, USGS, OPENSTREETMAP CONTRIBUTORS, GIS USER COMMUNITY (2025): World Topographic map. https://esri.maps.arcgis.com/home/item.html?id=67372ff42cd145319639a99152b15bc3 (1. 6. 2025).

WIKIPEDIA (2025): Pětiúhelník. https://cs.wikipedia.org/wiki/P%C4%9Bti%C3%BAheln%C3%ADk (27. 5. 2025).