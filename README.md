
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

| Teleconnect Name | Sub Folder Location | Full File Path                                                                      | Data Location                                                                                            |
| :--------------- | :------------------ | :---------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------- |
| watersheds       | water               | water/CWM\_v2\_2/World\_Watershed8.shp                                              | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>                                           |
| withdrawal       | water               | water/CWM\_v2\_2/Snapped\_Withdrawal\_Points.shp                                    | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>                                           |
| citypoint        | water               | water/CWM\_v2\_2/City\_Centroid.shp                                                 | <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF1J67DWR>                                           |
| powerplants      | water               | water/UCS-EW3-Energy-Water-Database.xlsx                                            | <https://www.ucsusa.org/resources/ucs-ew3-energy-water-database>                                         |
| crop             | land                | land/2016\_90m\_cdls/cdl\_lowres\_usa.img                                           | <https://www.nass.usda.gov/Research_and_Science/Cropland/Release/>                                       |
| crop\_attributes | land                | land/2016\_90m\_cdls/cdl\_lowres\_usa.img.vat.dbf                                   | <https://www.nass.usda.gov/Research_and_Science/Cropland/Release/>                                       |
| irrigation       | land                | land/Version2\_USA\_Demeter.csv                                                     | NA                                                                                                       |
| nlud             | land                | land/usa\_nlud\_LR.tif                                                              | <https://drive.google.com/file/d/1vmNfwjcaLf0sZTYJ1wsB3liG37sN8gyC/view>                                 |
| hydro            | energy              | energy/EHA\_Public\_PlantFY2019\_GIS\_6/ORNL\_EHAHydroPlant\_PublicFY2019final.xlsx | <https://hydrosource.ornl.gov/node/250>                                                                  |
| transfers        | NA                  | NA                                                                                  | NA                                                                                                       |
| climate          | land                | land/kop\_climate\_classes.tif                                                      | <http://koeppen-geiger.vu-wien.ac.at/present.htm>                                                        |
| HUC4             | water               | water/USA\_HUC4/huc4\_to\_huc2.shp                                                  | <http://prd-tnm.s3-website-us-west-2.amazonaws.com/?prefix=StagedProducts/Hydrography/WBD/National/GDB/> |
| population       | land                | land/pden2010\_block/pden2010\_60m.tif                                              | <https://www.sciencebase.gov/catalog/item/57753ebee4b07dd077c70868>                                      |
| runoff           | water               | water/Historical\_Mean\_Runoff/USA\_Mean\_Runoff.tif                                | NA                                                                                                       |
| nhd\_flow        | water               | water/Watershed\_Flow\_Contributions/UWB\_Intake\_Flows.shp                         | NA                                                                                                       |
| contributions    | water               | water/Watershed\_Flow\_Contributions/Watershed\_Contributions.csv                   | NA                                                                                                       |

##### Data Downloading Instructions

Below are instructions on how to download each dataset and prepare them
for the package.

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
      
      `nlud` dataset:
      1: Download the zip file from the google drive listed in the table above. 
      2: Unzip the file, and the remove the TIFs listed below (these are empty):
          - NLUDv1_2010_30m-0000000000-0000000000.tif
          - NLUDv1_2010_30m-0000093184-0000232960.tif
          - NLUDv1_2010_30m-0000093184-0000000000.tif
          - NLUDv1_2010_30m-0000046592-0000232960.tif
          - NLUDv1_2010_30m-0000046592-0000000000.tif
          - NLUDv1_2010_30m-0000000000-0000232960.tif
      3: Change the resolution of these rasters to 90m by 90m, and then mosaic them together. Save the mosaic file to the file path in the table.
      4: (Efficient code to do this will be placed here.)
      
      `hyrdo` dataset:
      1. Follow the link in the table above and click the download button for the ORNL EHAHydroPlant FY2020revised (xlsx).
      2. After downloading place it into the `energy` folder.
      
      `climate` dataset:
      1. Follow the link in the table above and scroll to the High Resolution Map and Data section.
      2. Click the download link next to R code with raster file.
      3. Crop this raster file to the shape of the Continental United States, and save it as a new TIF file in the `land` folder under the name "koppen_climate_classes.tif".
      
      `population` dataset:
      1. Follow the link in the table above and download the attached file named "pden2010_block.zip". Place this file into the `land` folder.
      2. This raster will also need to be cropped to the shape of the Continental United States.

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

| city             | city\_population | n\_watersheds | n\_other\_cities | dependent\_city\_pop | watershed\_area\_sqkm | storage\_BCM | yield\_BCM | irr\_cons\_BCM | n\_climate\_zones | n\_transfers\_in | n\_transfers\_out | n\_transfers\_within | n\_hydro\_plants | n\_thermal\_plants | n\_fac\_agcrop | n\_fac\_aglivestock | n\_fac\_cnsmnf | n\_fac\_mining | n\_fac\_oilgas | n\_fac\_total | hydro\_gen\_MWh | thermal\_gen\_MWh | thermal\_cons\_BCM | thermal\_with\_BCM | n\_utilities | n\_ba | n\_crop\_classes | cropland\_fraction | developed\_fraction | ag\_runoff\_max | ag\_runoff\_av\_exgw | ag\_runoff\_av | dev\_runoff\_max | dev\_runoff\_av\_exgw | dev\_runoff\_av | np\_runoff\_max | np\_runoff\_av\_exgw | np\_runoff\_av\_exgw\_unweighted | np\_runoff\_av | n\_economic\_sectors | max\_withdr\_dist\_km | avg\_withdr\_dis\_km | n\_treatment\_plants | watershed\_pop | pop\_cons\_m3sec | av\_fl\_sur\_conc\_pct | av\_fl\_sur\_conc\_pct\_unweighted | av\_ro\_sur\_conc\_pct | av\_fl\_all\_conc\_pct | av\_ro\_all\_conc\_pct | av\_fl\_max\_conc\_pct | av\_ro\_max\_conc\_pct | surface\_contribution\_pct | importance\_of\_worst\_watershed\_pct |
| :--------------- | ---------------: | ------------: | ---------------: | -------------------: | --------------------: | -----------: | ---------: | -------------: | ----------------: | ---------------: | ----------------: | -------------------: | ---------------: | -----------------: | -------------: | ------------------: | -------------: | -------------: | -------------: | ------------: | --------------: | ----------------: | -----------------: | -----------------: | -----------: | ----: | ---------------: | -----------------: | ------------------: | --------------: | -------------------: | -------------: | ---------------: | --------------------: | --------------: | --------------: | -------------------: | -------------------------------: | -------------: | -------------------: | --------------------: | -------------------: | -------------------: | -------------: | ---------------: | ---------------------: | ---------------------------------: | ---------------------: | ---------------------: | ---------------------: | ---------------------: | ---------------------: | -------------------------: | ------------------------------------: |
| Atlanta, GA      |           498044 |             1 |                0 |                    0 |             4109.2426 |    2.4119697 |  0.0710615 |      0.0054521 |                 2 |               NA |                NA |                   NA |                2 |                  0 |              0 |                   2 |            361 |             14 |              0 |          9402 |       185763.17 |                 0 |          0.0000000 |          0.0000000 |            0 |     0 |                7 |          0.0094468 |           0.3660531 |       0.0001360 |            0.0001360 |      0.0001360 |        0.2660990 |             0.2660990 |       0.2660990 |       0.2662349 |            0.2662349 |                        0.2662349 |      0.2662349 |                   13 |             10.798758 |            10.798758 |                    0 |   1.530574e+06 |        3.5430035 |              4.6432565 |                          4.6432565 |              5.8426343 |              4.6432565 |              5.8426343 |               4.643256 |               5.842634 |                      100.0 |                                 100.0 |
| Portland, OR     |           653115 |             1 |                1 |               653115 |              280.7526 |    0.0970009 |  0.0710615 |      0.0000000 |                 4 |               NA |                NA |                   NA |                1 |                  0 |              0 |                   0 |              0 |              0 |              0 |             0 |        55263.24 |                 0 |          0.0000000 |          0.0000000 |            0 |     0 |                0 |          0.0000000 |           0.0044529 |       0.0000000 |            0.0000000 |      0.0000000 |        0.0000172 |             0.0000172 |       0.0000172 |       0.0000172 |            0.0000172 |                        0.0000172 |      0.0000172 |                    4 |             40.923236 |            40.923236 |                    0 |   1.963243e+01 |        0.0000454 |              0.0000000 |                          0.0000000 |              0.0000000 |              0.0000000 |              0.0000000 |               0.000000 |               0.000000 |                      100.0 |                                 100.0 |
| Charlotte, NC    |           872498 |             1 |                3 |               212079 |             4875.1775 |    3.2929366 |  0.0710615 |      0.0173409 |                 3 |               NA |                NA |                   NA |                6 |                  3 |              0 |                   1 |            213 |             10 |              0 |           442 |       574511.46 |          23312996 |          0.0226557 |          3.4663151 |            1 |     1 |                8 |          0.0456691 |           0.1762985 |       0.0023421 |            0.0023421 |      0.0023421 |        0.0436056 |             0.0436056 |       0.0436056 |       0.0459477 |            0.0459477 |                        0.0459477 |      0.0459477 |                   12 |             17.065885 |            17.065885 |                    0 |   4.668231e+05 |        1.0806115 |              1.5690333 |                          1.5690333 |              1.6886042 |              1.5690333 |              1.6886042 |               1.569033 |               1.688604 |                      100.0 |                                 100.0 |
| Austin, TX       |           964254 |             1 |                3 |              1207821 |            97194.7485 |    6.4692813 |  0.0710615 |      0.6764835 |                 3 |               NA |                NA |                   NA |                6 |                 10 |              0 |                   6 |             10 |              0 |              1 |            72 |       231050.88 |           7425960 |          0.0098231 |          0.1109466 |            9 |     2 |                8 |          0.1451631 |           0.0520904 |       0.0103864 |            0.0103864 |      0.0103864 |        0.0060554 |             0.0060554 |       0.0060554 |       0.0164418 |            0.0164418 |                        0.0164418 |      0.0164418 |                   13 |              3.503236 |             3.503236 |                    0 |   9.394054e+05 |        2.1745543 |              1.1588398 |                          1.1588398 |              1.2682985 |              1.1588398 |              1.2682985 |               1.158840 |               1.268298 |                      100.0 |                                 100.0 |
| Macon, GA        |           153095 |             1 |                1 |               153095 |             5805.3543 |    0.5672585 |  0.0710615 |      0.0052724 |                 1 |               NA |                NA |                   NA |                3 |                  1 |              0 |                   2 |            358 |             33 |              0 |         12440 |        68416.13 |          24348774 |          0.0633148 |          0.0926290 |            1 |     1 |                8 |          0.0179467 |           0.3019529 |       0.0004749 |            0.0004749 |      0.0004749 |        0.2396321 |             0.2396321 |       0.2396321 |       0.2401070 |            0.2401070 |                        0.2401070 |      0.2401070 |                   13 |              3.430790 |             3.430790 |                    0 |   1.674266e+06 |        3.8756251 |              5.7257250 |                          5.7257250 |              5.7789911 |              5.7257250 |              5.7789911 |               5.725725 |               5.778991 |                      100.0 |                                 100.0 |
| Philadelphia, PA |          1584138 |             2 |                4 |             10155355 |            25000.1122 |    2.4580071 |  0.1421230 |      0.0413309 |                 2 |               NA |                NA |                   NA |               10 |                 16 |              1 |                  17 |            301 |             29 |              1 |          1764 |       249723.16 |          23756210 |          0.0431045 |          0.7471718 |           14 |     1 |                8 |          0.1128719 |           0.1799514 |       0.1019802 |            0.0494111 |      0.0494111 |        0.2038389 |             0.1076844 |       0.1076844 |       0.3058191 |            0.1570955 |                        0.1818828 |      0.1570955 |                   13 |             15.865748 |             9.745459 |                    0 |   4.083398e+06 |        9.4523306 |              3.7301834 |                          4.1224094 |              4.0931675 |              3.7301834 |              4.0931675 |               6.083539 |               6.462299 |                      100.0 |                                  40.0 |
| Houston, TX      |          2325502 |             2 |                3 |              4565557 |            50273.6660 |   10.0813681 |  0.1421230 |      0.0057792 |                 1 |               NA |                NA |                   NA |                2 |                 18 |              0 |                   0 |              2 |              0 |              0 |             8 |         5079.73 |          25375485 |          0.0307072 |          0.5044097 |           11 |     2 |               10 |          0.0732626 |           0.1768837 |       0.0122991 |            0.0108832 |      0.0093595 |        0.1147764 |             0.0657315 |       0.0565291 |       0.1148982 |            0.0766146 |                        0.0932378 |      0.0658886 |                   13 |            102.672878 |            65.669248 |                    0 |   7.679998e+06 |       17.7778126 |             16.7718115 |                         10.0307609 |              9.9591510 |             14.4237579 |              8.5648699 |              18.814554 |              11.132603 |                       86.0 |                                  76.0 |
| Raleigh, NC      |           469298 |             1 |                2 |               743589 |             1987.8516 |    0.1921883 |  0.0710615 |      0.0041787 |                 1 |               NA |                NA |                   NA |                0 |                  0 |              0 |                   0 |             38 |              6 |              0 |           657 |            0.00 |                 0 |          0.0000000 |          0.0000000 |            0 |     0 |                8 |          0.0939267 |           0.1669958 |       0.0161245 |            0.0161245 |      0.0161245 |        0.0529000 |             0.0529000 |       0.0529000 |       0.0690245 |            0.0690245 |                        0.0690245 |      0.0690245 |                   12 |             18.713965 |            18.713965 |                    0 |   2.537039e+05 |        0.5872788 |              2.0842711 |                          2.0842711 |              3.3143195 |              2.0842711 |              3.3143195 |               2.084271 |               3.314319 |                      100.0 |                                 100.0 |
| Knoxville, TN    |           187500 |             1 |                3 |               346730 |            23196.0276 |    5.1667173 |  0.0710615 |      0.0132724 |                 3 |               NA |                NA |                   NA |               13 |                  4 |              0 |                   7 |            149 |             20 |              0 |           569 |      1633153.23 |           8595014 |          0.0105083 |          0.9841422 |            4 |     4 |                7 |          0.0392887 |           0.1217025 |       0.0017350 |            0.0017350 |      0.0017350 |        0.0187317 |             0.0187317 |       0.0187317 |       0.0204667 |            0.0204667 |                        0.0204667 |      0.0204667 |                   13 |              5.942854 |             5.942854 |                    0 |   1.559519e+06 |        3.6100047 |              1.4421179 |                          1.4421179 |              1.7276066 |              1.4421179 |              1.7276066 |               1.442118 |               1.727607 |                      100.0 |                                 100.0 |
| NewYork, NY      |          8398748 |             8 |                1 |              8398748 |             5203.8799 |    2.3480124 |  0.5684919 |      0.0046190 |                 2 |               NA |                NA |                   NA |                5 |                  0 |              0 |                   3 |             45 |              7 |              0 |           377 |       155873.63 |                 0 |          0.0000000 |          0.0000000 |            0 |     0 |                7 |          0.0473016 |           0.0822482 |       0.0207107 |            0.0034037 |      0.0034037 |        0.3966929 |             0.0263228 |       0.0263228 |       0.3966933 |            0.0297265 |                        0.0667699 |      0.0297265 |                   13 |            191.509637 |           130.711473 |                    0 |   2.566443e+05 |        0.5940853 |              0.1757086 |                          0.3212001 |              0.1946581 |              0.1757086 |              0.1946581 |               1.363024 |               1.150314 |                      100.0 |                                   5.0 |
| LosAngeles, CA   |          3990456 |            30 |               29 |             14790722 |           499549.8199 |  129.2864648 |  2.4160906 |      0.9946673 |                15 |               NA |                NA |                   NA |              196 |                 38 |              7 |                  12 |           3312 |            445 |            211 |          8436 |     26779746.24 |         132000904 |          0.3249964 |          1.2873801 |           23 |    11 |               10 |          0.0125410 |           0.0110604 |       0.0594559 |            0.0010065 |      0.0009210 |        0.3849221 |             0.0004974 |       0.0004551 |       0.3849221 |            0.0015039 |                        0.0211506 |      0.0013761 |                   13 |            843.137141 |           331.178945 |                    0 |   3.809702e+06 |        8.8187738 |              0.5737800 |                          0.3949688 |              0.7595693 |              0.5250087 |              0.6950059 |              10.872884 |              10.872884 |                       91.5 |                                   1.5 |

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

Below is the map that is produced:

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

The table below shows explanations for each of these variables that are
created through this function:

| Variable Name                         | Description                                                                             | Units                 |
| :------------------------------------ | :-------------------------------------------------------------------------------------- | :-------------------- |
| city\_population                      | The population of the city being analyzed                                               | people                |
| n\_watersheds                         | Number of watersheds that city uses to source drinking water                            | watersheds            |
| n\_other\_cities                      | Number of other cities pulling off the same watersheds                                  | cities                |
| dependent\_city\_pop                  | Total population of people dependent on that city’s watersheds                          | people                |
| watershed\_area\_sqkm                 | Combined area of all the source watersheds of a city                                    | square kilometers     |
| storage\_BCM                          | Combined storage capacity of all the city catchments                                    | billion cubic meters  |
| yield\_BCM                            | Combined yield capacity of all the city catchments                                      | billion cubic meters  |
| irr\_cons\_BCM                        | Combined water consumption that is used for irrigation with the watersheds              | billion cubic meters  |
| n\_climate\_zones                     | Number of climate zones that the source watersheds cover                                | zones                 |
| n\_transfers\_in                      | Number of interbasin transfers that flow into the source watersheds                     | transfers             |
| n\_transfers\_out                     | Number of interbasin transfers that flow out of the source watersheds                   | transfers             |
| n\_transfers\_within                  | Number of water transfers that occur within the source watersheds                       | transfers             |
| n\_hydro\_plants                      | Number of hydro electric power plants operating within the source watersheds            | plants                |
| n\_thermal\_plants                    | Number of thermal power plants operating within the source watersheds                   | plants                |
| n\_fac\_agcrop                        | Number of agricultural crop facilities within the source watersheds                     | facilities            |
| n\_fac\_aglivestock                   | Number of agicultural livestock facilities within the source watersheds                 | facilities            |
| n\_fac\_cnsmnf                        |                                                                                         | facilities            |
| n\_fac\_mining                        | Number of mining facilities within the source watersheds                                | facilities            |
| n\_fac\_oilgas                        | Number of oil and gas facilities within the source watersheds                           | facilities            |
| n\_fac\_total                         | Total number of facilities operating within the source watersheds                       | facilities            |
| hydro\_gen\_MWh                       | Combined hydro electric generation from all the facilities within the source watersheds | Megawatthours         |
| thermal\_gen\_MWh                     | Combined thermal generation from all the facilities within the source watersheds        | Megawatthours         |
| thermal\_cons\_BCM                    | Combined water consumption that is used for thermal generation                          | billion cubic meters  |
| thermal\_with\_BCM                    | Combined water withdrawal for thermal generation                                        | billion cubic meters  |
| n\_utilities                          | Number of electric utilities within the source watersheds                               | utilities             |
| n\_ba                                 | Number of balancing authorities within the source watersheds                            | balancing authorities |
| n\_crop\_classes                      | Total number of different types of crops within the source watersheds                   | crops                 |
| cropland\_fraction                    | Fraction of land that is used for crops within the source watersheds                    | NA                    |
| developed\_fraction                   | Fraction of land that is developed within the source watersheds                         | NA                    |
| ag\_runoff\_max                       | Max amount of agricultural runoff within the source watersheds                          |                       |
| ag\_runoff\_av\_exgw                  |                                                                                         |                       |
| ag\_runoff\_av                        | Average runoff from agricultural lands                                                  |                       |
| dev\_runof\_max                       | Max amount of agricultural runoff within the source watersheds                          |                       |
| dev\_runof\_av\_exgw                  |                                                                                         |                       |
| dev\_runof\_av                        | Average runoff from developed lands                                                     |                       |
| np\_runoff\_max                       | Max amount of non-point source runoff within the source watersheds                      |                       |
| np\_runoff\_av\_exgw                  |                                                                                         |                       |
| np\_runoff\_av\_exgw\_unweighted      |                                                                                         |                       |
| np\_runoff\_av                        | Average non-point source runoff.                                                        |                       |
| n\_economic\_sectors                  | Total number of different economic sectors within the source watersheds                 | sectors               |
| max\_withdr\_dis\_km                  | Maximum distance between a city’s intake points                                         | kilometers            |
| avg\_withdr\_dis\_km                  | Average distance between a city’s intake points                                         | kilometers            |
| n\_treatment\_plants                  | Total number of waste water treatment plants operating within the source watersheds     | plants                |
| watershed\_pop                        | Total number of people living within the source watershed boundaries                    | people                |
| pop\_cons\_m3sec                      | Combined water consumption from the source watersheds that is used for people           | m3/sec                |
| av\_fl\_sur\_conc\_pct                |                                                                                         | %                     |
| av\_fl\_sur\_conc\_pct\_unweighted    |                                                                                         | %                     |
| av\_ro\_sur\_conc\_pct                |                                                                                         | %                     |
| av\_fl\_all\_conc\_pct                |                                                                                         | %                     |
| av\_ro\_all\_conc\_pct                |                                                                                         | %                     |
| av\_fl\_max\_conc\_pct                |                                                                                         | %                     |
| av\_ro\_max\_conc\_pct                |                                                                                         | %                     |
| surface\_contribution\_pct            |                                                                                         | %                     |
| importance\_of\_worst\_watershed\_pct |                                                                                         | %                     |

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
