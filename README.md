[![](https://img.shields.io/badge/License-MIT-black.svg)](https://lbesson.mit-license.org/)
[![](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8-blue)](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8-blue)
![](https://img.shields.io/badge/version-1.0.0-blue.svg?cacheSeconds=2592000)

```python
>>> from geosample import QuadTree
>>> import geopandas as gpd
>>>
>>> samples = gpd.read_file('samples.gpkg')
>>>
>>> qt = QuadTree(samples)
>>>
>>> # Split until the quadrants are less than 5,000 meters
>>> qt.split_recursive(max_length=5000)
>>>
>>> # Get the actual quadrant length
>>> qt.qlen
>>>
>>> # Get the quadrants as a GeoDataFrame
>>> df = qt.to_frame()
>>>
>>> # Get 5 random points using the Generalized Random Tessellation Stratified (GRTS) method
>>> dfs = qt.sample(n=5)
```