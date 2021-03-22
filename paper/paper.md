---
title: 'GAMUT: A Geospatial Analysis of Multisector Urban Teleconnections'
tags:
  - R
  - Multisector
  - Dynamics
  - Urban
  - Watersheds
  - Geospatial
authors:
  - name: Kristian D. Nelson
    orcid: 0000-0002-6745-167X
    affiliation: 1
  - name: Sean W. Turner
    orcid: 0000-0003-4400-9800
    affiliation: 1
  - name: Chris Vernon
    orcid: 0000-0002-3406-6214
    affiliation: 1
  - name: Jennie Rice
    orcid: 0000-0002-7833-9456
    affiliation: 1
affiliations:
 - name: Pacific Northwest National Laboratory
   index: 1
date: 1 February 2020
bibliography: paper.bib
---

# Summary

The vulnerability of urban cities' source watersheds is a characteristic that fluctuates based on the changes of many land-water-energy systems. The connections between these sectors are highly complex, which leads to a vague understanding of how cities and their watersheds are affected through distal relationships of a variety of factors. These distal relationships between sectors are defined as "teleconnections" and help bridge the gap between human and natural systems. Some of these factors include, energy production through thermal or hydro generation, land cover, hydrological properties of the watersheds, and water demands of the surrounding population. The relationships between these systems can be quantified through the ``GAMUT`` package and can aid in classifying a city's vulnerability based on the complexity of its teleconnections.

# Statement of Need

``GAMUT`` is an geospatial R package that analyzes multisector urban teleconnections. The main input of this package is the Urban Water Blueprint dataset [@McDonald:2014] which supplies spatial polygons of source watersheds and intake points for 235 major US cities. The package uses other geospatial layers like high resolution land cover rasters, annual runoff simulations, population data, wastewater discharge, watershed pollution, and much more. These layers are used to analyze cities' vulnerability or other research interests like water quality. This is done through intensive overlay of layers and comprehensive analysis of features within urban water sources. ``GAMUT`` is able to quickly run through each city, conduct these overlays and analysis, and produce a multitude of results that describe the connections between human and natural systems.

This package is being utilized in a project that looks at potential surface water contamination in urban cities' drinking water supplies [@Turner:2020]. With the help of ``GAMUT``, this project is able to input layers like land cover or land use, wastewater treatment plants, runoff volumes, and river flows to conduct this analysis. By overlaying land cover and land use rasters with the watershed polygons, a fraction can be found of land that is developed or used for crops in a catchment. With this fraction, the amount of runoff flowing from those sources into the watershed is calculated to find the amount of polluted flow coming from non-natural lands. Another variable that is analyzed is the amount of waste water flowing back into the basin from upstream wastewater treatment plants. This is found by overlaying a treatment plant point file with the watershed polygons and calculating the total amount of recycled water flowing out of the plants. Using ``GAMUT`` as the primary tool, this project was successfully able to analyze urban city water quality for 116 cities and pinpoint places that might have drinking water pollution. This project can be found at the linked DOI: [@Turner:2020].

To download and use the ``GAMUT`` package, you can find it at this cited link: [] 


# Acknowledgements

This research was supported by the U.S. Department of Energy, Office of Science, as part of research in Multi-Sector Dynamics, Earth and Environmental System Modeling Program. 

# References
