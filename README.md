
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/IMMM-SFA/teleconnect.svg?branch=master)](https://travis-ci.org/IMMM-SFA/teleconnect)
[![DOI](https://zenodo.org/badge/203447802.svg)](https://zenodo.org/badge/latestdoi/203447802)

# teleconnect

#### An R package to identify multi-sector teleconnection complexity

## Description

The `teleconnect` package classifies US cities according to the
complexity of their distal relationships across water, energy, and land
sectors. The cities analyzed in this package are based off the Urban
Water Blueprint dataset.

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Get Started with `teleconnect`

`teleconnect` can be installed directly from its GitHub repository using
the R `devtools` package. From an R prompt, run the command:

``` r
devtools::install_github('IMMM-SFA/teleconnect')
```

##### Data Directory Set Up

In order for the package to function properly, the user needs to
correctly set up a data directory, and place all of the necessary
datasets in that directory. The steps below show how to set up the
directory so that the package can read in the input datasets correctly.

    Step 1: In a drive of your choice, create a folder called "data_dir".
    Step 2: Within the "data_dir" folder, create the 3 folders: "energy","land","water". 
            Follow this set up:
            
            Main Folder: data_dir
                  Sub Folder: energy
                  Sub Folder: land
                  Sub Folder: water
    
    Step 3: Follow the instructions in the "Data Downloading Instructions" section below to download the datasets. The             table below shows where to find each dataset, and what folder it needs to go in. 

| Teleconnect Name | Sub Folder Location | Full File Path                                                                      | Data Location                                                       |
| :--------------- | :------------------ | :---------------------------------------------------------------------------------- | :------------------------------------------------------------------ |
| waterseds        | water               | water/CWM\_v2\_2/World\_Watershed8.shp                                              | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>      |
| withdrawal       | water               | water/CWM\_v2\_2/Snapped\_Withdrawal\_Points.shp                                    | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>      |
| citypoint        | water               | water/CWM\_v2\_2/City\_Centroid.shp                                                 | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>      |
| powerplants      | water               | water/UCS-EW3-Energy-Water-Database.xlsx                                            | <https://www.ucsusa.org/resources/ucs-ew3-energy-water-database>    |
| crop             | land                | land/2016\_90m\_cdls/cdl\_lowres\_usa.img                                           | <https://www.nass.usda.gov/Research_and_Science/Cropland/Release/>  |
| crop\_attributes | land                | land/2016\_90m\_cdls/cdl\_lowres\_usa.img.vat.dbf                                   | NA                                                                  |
| irrigation       | land                | land/Version2\_USA\_Demeter.csv                                                     | NA                                                                  |
| nlud             | land                | land/usa\_nlud\_LR.tif                                                              | NA                                                                  |
| hydro            | energy              | energy/EHA\_Public\_PlantFY2019\_GIS\_6/ORNL\_EHAHydroPlant\_PublicFY2019final.xlsx | NA                                                                  |
| transfers        | water               | water/transfers/USIBTsHUC6\_Dickson.shp                                             | NA                                                                  |
| climate          | land                | land/kop\_climate\_classes.tif                                                      | <http://koeppen-geiger.vu-wien.ac.at/present.htm>                   |
| HUC4             | water               | water/USA\_HUC4/huc4\_to\_huc2.shp                                                  | NA                                                                  |
| population       | land                | land/pden2010\_block/pden2010\_60m.tif                                              | <https://www.sciencebase.gov/catalog/item/57753ebee4b07dd077c70868> |
| runoff           | water               | water/Historical\_Mean\_Runoff/USA\_Mean\_Runoff.tif                                | NA                                                                  |
| nhd\_flow        | water               | water/Watershed\_Flow\_Contributions/UWB\_Intake\_Flows.shp                         | NA                                                                  |
| contributions    | water               | water/Watershed\_Flow\_Contributions/Watershed\_Contributions.csv                   | NA                                                                  |

##### Data Downloading Instructions

Below are instructions on how to download each dataset above from the
online sources.

``` 
`watersheds`,`withdrawal`,`citypoint` datasets:
  1. Follow the link in the table above and click the Download button across from "CWM_v2_2.zip".
  2. This will download a data file named `knb.1002.1`. Right click and unzip this folder. Go inside the unzipped file, copy the folder named "CWM_v2_2.zip", and place it into the `water` folder.
  
`powerplants` dataset:
  1. Follow the link in the table above and download the Excel file at the bottom of the page.
  2. After it downloads, place it into the `water` folder.
  
`crop` and `crop_attributes` datasets:
  1. Follow the link in the table above and download the zipped 2019 CDL file.
  2. This layer requires some processing before it can function properly in the package. In order for the package to run this data, the file size needs to be decreased. We do this by changing the resolution from 30m cells to 90 meter cells.
  3. Open up ArcGIS and load in the CDL file. Run this file through the "Resample" function, and set the cell size to 90 meters by 90 meters.
  4. Once the output is created, run the output through the "Project Raster" function, and change the projection to `GCS_WGS_1984`.
  5. Export the output of this function to `land` folder under the file name of `cdl_lowres_usa`, and set the file type to an IMAGINE image file. This will also produce the `crop_attributes` file. 
  
```

##### Data information

Data Collected: 2019-2020

Geographic location: Continental United States

Projection: WGS 1984

## Usage

Once all the necessary data is organized in the correct format, the
package is ready to be used. The main function is
`count_watershed_teleconnections`. Within this function, set up your
data directory and select cities. The city names need to be in the
format of `City | State`, and for multiple cities, place them inside
`c()`.

Here is an example of what you would type into your console:

``` r
count_watershed_teleconnections(data_dir = "C:/data_dir/", cities = c("Portland | OR", "Charlotte | NC", "Austin | TX"))
```

The package will cycle through each city and their respective
watersheds, and produce a table with several columns of information. To
learn what each of these variables mean, scroll to the bottom of this
section and see the table of variables. The result of this function will
look something like this:

| city          | city\_population | n\_watersheds | n\_other\_cities | dependent\_city\_pop | watershed\_area\_sqkm | storage\_BCM | yield\_BCM | irr\_cons\_BCM | n\_climate\_zones | n\_transfers\_in | n\_transfers\_out | n\_transfers\_within | n\_hydro\_plants | n\_thermal\_plants | n\_fac\_agcrop | n\_fac\_aglivestock | n\_fac\_cnsmnf | n\_fac\_mining | n\_fac\_oilgas | n\_fac\_total | hydro\_gen\_MWh | thermal\_gen\_MWh | thermal\_cons\_BCM | thermal\_with\_BCM | n\_utilities | n\_ba | n\_crop\_classes | cropland\_fraction | developed\_fraction | ag\_runoff\_max | ag\_runoff\_av\_exgw | ag\_runoff\_av | dev\_runoff\_max | dev\_runoff\_av\_exgw | dev\_runoff\_av | np\_runoff\_max | np\_runoff\_av\_exgw | np\_runoff\_av\_exgw\_unweighted | np\_runoff\_av | n\_economic\_sectors | max\_withdr\_dist\_km | avg\_withdr\_dis\_km | n\_treatment\_plants | watershed\_pop | pop\_cons\_m3sec | av\_fl\_sur\_conc\_pct | av\_fl\_sur\_conc\_pct\_unweighted | av\_ro\_sur\_conc\_pct | av\_fl\_all\_conc\_pct | av\_ro\_all\_conc\_pct | av\_fl\_max\_conc\_pct | av\_ro\_max\_conc\_pct | surface\_contribution\_pct | importance\_of\_worst\_watershed\_pct |
| :------------ | ---------------: | ------------: | ---------------: | -------------------: | --------------------: | -----------: | ---------: | -------------: | ----------------: | ---------------: | ----------------: | -------------------: | ---------------: | -----------------: | -------------: | ------------------: | -------------: | -------------: | -------------: | ------------: | --------------: | ----------------: | -----------------: | -----------------: | -----------: | ----: | ---------------: | -----------------: | ------------------: | --------------: | -------------------: | -------------: | ---------------: | --------------------: | --------------: | --------------: | -------------------: | -------------------------------: | -------------: | -------------------: | --------------------: | -------------------: | -------------------: | -------------: | ---------------: | ---------------------: | ---------------------------------: | ---------------------: | ---------------------: | ---------------------: | ---------------------: | ---------------------: | -------------------------: | ------------------------------------: |
| Portland, OR  |           653115 |             1 |                0 |                    0 |              280.7526 |    0.0970009 |  0.3529533 |      0.0000000 |                 4 |               NA |                NA |                   NA |                1 |                  0 |              0 |                   0 |              0 |              0 |              0 |             0 |        55954.44 |                 0 |          0.0000000 |          0.0000000 |            0 |     0 |                0 |          0.0000000 |           0.0044529 |       0.0000000 |            0.0000000 |      0.0000000 |        0.0000172 |             0.0000172 |       0.0000172 |       0.0000172 |            0.0000172 |                        0.0000172 |      0.0000172 |                    4 |             40.923236 |            40.923236 |                    0 |       19.63243 |        0.0000454 |               0.000000 |                           0.000000 |               0.000000 |               0.000000 |               0.000000 |               0.000000 |               0.000000 |                        100 |                                   100 |
| Charlotte, NC |           872498 |             1 |                4 |              1084577 |             4875.1775 |    3.2929366 |  0.8235577 |      0.0173409 |                 3 |               NA |                NA |                   NA |                6 |                  3 |              0 |                   1 |            213 |             10 |              0 |           442 |       546288.01 |          23312996 |          0.0226557 |          3.4663151 |            1 |     1 |                8 |          0.0456691 |           0.1762985 |       0.0023421 |            0.0023421 |      0.0023421 |        0.0436056 |             0.0436056 |       0.0436056 |       0.0459477 |            0.0459477 |                        0.0459477 |      0.0459477 |                   12 |             17.065885 |            17.065885 |                    0 |   466823.11767 |        1.0806115 |               1.569033 |                           1.569033 |               1.688604 |               1.569033 |               1.688604 |               1.569033 |               1.688604 |                        100 |                                   100 |
| Austin, TX    |           964254 |             1 |                2 |               243567 |            97194.7485 |    6.4692813 |  0.8235577 |      0.6764835 |                 3 |               NA |                NA |                   NA |                6 |                 10 |              0 |                   6 |             10 |              0 |              1 |            72 |       227513.00 |           7425960 |          0.0098231 |          0.1109466 |            9 |     2 |                8 |          0.1451631 |           0.0520904 |       0.0103864 |            0.0103864 |      0.0103864 |        0.0060554 |             0.0060554 |       0.0060554 |       0.0164418 |            0.0164418 |                        0.0164418 |      0.0164418 |                   13 |              3.503236 |             3.503236 |                    0 |   939405.35947 |        2.1745543 |               1.158840 |                           1.158840 |               1.268298 |               1.158840 |               1.268298 |               1.158840 |               1.268298 |                        100 |                                   100 |

This table can be used to compare different variables between multiple
cities. Below is a graph comparing how much developed land are in
cities’ watersheds.

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

To get a closer look at the cities, their watersheds, and the type of
land within those watersheds, run the `plot_watershed` function.

Here is an example of what to type into the console:

``` r
plot_watershed(data_dir = "C:/data_dir/", city = "Los Angeles | CA")
```

Below are some examples of maps that can be produced:

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

The table below shows explanations for each of these variables that are
created through this function:

| Variable Name                         | Description                                                                             | Units                   |
| :------------------------------------ | :-------------------------------------------------------------------------------------- | :---------------------- |
| city\_population                      | The population of the city being analyzed                                               | people                  |
| n\_watershed                          | Number of watersheds that city uses to source drinking water                            | watersheds              |
| n\_other\_cities                      | Number of other cities pulling off the same watersheds                                  | cities                  |
| dependent\_city\_pop                  | Total population of people dependent on that city’s watersheds                          | people                  |
| watershed\_area\_sqkm                 | Combined area of all the source watersheds of a city                                    | square kilometres       |
| storage\_BCM                          | Combined storage capacity of all the city catchments                                    | billion cubic meters    |
| yield\_BCM                            | Combined yield capacity of all the city catchments                                      | billion cubic meters    |
| irr\_cons\_BCM                        | Combined water consumption that is used for irrigation with the watersheds              | billion cubic meters    |
| n\_climate\_zones                     | Number of climate zones that the source watersheds cover                                | zones                   |
| n\_transfers\_in                      | Number of interbasin transfers that flow into the source watersheds                     | transfers               |
| n\_transfers\_out                     | Number of interbasin transfers that flow out of the source watersheds                   | transfers               |
| n\_transfers\_within                  | Number of water transfers that occur within the source watersheds                       | transfers               |
| n\_hydro\_plants                      | Number of hydro electric power plants operating within the source watersheds            | plants                  |
| n\_thermal\_plants                    | Number of thermal power plants operating within the source watersheds                   | plants                  |
| n\_fac\_agcrop                        | Number of agricultural crop facilities within the source watersheds                     | facilities              |
| n\_fac\_aglivestock                   | Number of agicultural livestock facilities within the source watersheds                 | facilities              |
| n\_fac\_cnsmnf                        | NA                                                                                      | facilities              |
| n\_fac\_mining                        | Number of mining facilities within the source watersheds                                | facilities              |
| n\_fac\_oilgas                        | Number of oil and gas facilities within the source watersheds                           | facilities              |
| n\_fac\_total                         | Total number of facilities operating within the source watersheds                       | facilities              |
| hydro\_gen\_MWh                       | Combined hydro electric generation from all the facilities within the source watersheds | Megawatthours           |
| thermal\_gen\_MWh                     | Combined thermal generation from all the facilities within the source watersheds        | Megawatthours           |
| thermal\_cons\_BCM                    | Combined water consumption that is used for thermal generation                          | billion cubic meters    |
| thermal\_with\_BCM                    | Combined water withdrawal for thermal generation                                        | billion cubic meters    |
| n\_utilities                          | Number of electric utilities within the source watersheds                               | utilities               |
| n\_ba                                 | Number of balancing authorities within the source watersheds                            | balancing authorities   |
| n\_crop\_classes                      | Total number of different types of crops within the source watersheds                   | crops                   |
| cropland\_fraction                    | Fraction of land that is used for crops within the source watersheds                    | %                       |
| developed\_fraction                   | Fraction of land that is developed within the source watersheds                         | %                       |
| ag\_runoff\_max                       | Max amount of agricultural runoff within the source watersheds                          | NA                      |
| ag\_runoff\_av\_exgw                  | NA                                                                                      | NA                      |
| ag\_runoff\_av                        | Average runoff from agricultural lands                                                  | NA                      |
| dev\_runoff\_max                      | Max amount of agricultural runoff within the source watersheds                          | NA                      |
| dev\_runoff\_av\_exgw                 | NA                                                                                      | NA                      |
| dev\_runoff\_av                       | Average runoff from developed landsl.                                                   | NA                      |
| np\_runoff\_max                       | Max amount of non-point source runoff within the source watersheds                      | NA                      |
| np\_runoff\_av\_exgw                  | NA                                                                                      | NA                      |
| np\_runoff\_av\_exgw\_unweighted      | NA                                                                                      | NA                      |
| np\_runoff\_av                        | Average non-point source runoff.                                                        | NA                      |
| n\_economic\_sectors                  | Total number of different economic sectors within the source watersheds                 | sectors                 |
| max\_withdr\_dist\_km                 | Maximum distance between a city’s intake points                                         | kilometers              |
| avg\_withdr\_dis\_km                  | Average distance between a city’s intake points                                         | kilometers              |
| n\_treatment plants                   | Total number of waste water treatment plants operating within the source watersheds     | plants                  |
| watershed\_pop                        | Total number of people living within the source watershed boundaries                    | people                  |
| pop\_cons\_m3sec                      | Combined water consumption from the source watersheds that is used for people.          | cubic meters per second |
| av\_fl\_sur\_conc\_pct                | NA                                                                                      | %                       |
| av\_fl\_sur\_conc\_pct\_unweighted    | NA                                                                                      | %                       |
| av\_ro\_sur\_conc\_pct                | NA                                                                                      | %                       |
| av\_fl\_all\_conc\_pct                | NA                                                                                      | %                       |
| av\_ro\_all\_conc\_pct                | NA                                                                                      | %                       |
| av\_fl\_max\_conc\_pct                | NA                                                                                      | %                       |
| av\_ro\_max\_conc\_pct                | NA                                                                                      | %                       |
| surface\_contribution\_pct            | NA                                                                                      | %                       |
| importance\_of\_worst\_watershed\_pct | NA                                                                                      | %                       |

## Support

For any questions about the package, please contact any of the
contributors below:

Sean Turner: <sean.turner@pnnl.gov>

Kristian Nelson: <kristian.nelson@pnnl.gov>

## Authors and Acknowledgement

Authors: Sean Turner, Kristian Nelson, Chris Vernon, Jennie Rice, Casey
Burleyson, Ryan McManamay, Kerim Dickson

This research was supported by the US Department of Energy, Office of
Science, as part of research in the MultiSector Dynamics, Earth and
Environmental System Modeling Program.
