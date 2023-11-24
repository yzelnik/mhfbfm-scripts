This folder demonstrates simulations used for the manuscript "Managing for heterogeneity reduces fire risk in boreal forest landscapes – a model analysis".

The main script to run is SomeExamples.m, where it the script is partioned to several parts that can (and should) be run separately. This file directly runs 3 functions:
SetupLandscape.m   -- define landscape structure
GetPreFireState.m  -- define conditions before fire, based on landscape and other non-spatial conditions
RunFire.m          -- run a fire simulation, based on the output of the two previous functions

These 3 functions in turn use other functions in the folder.
