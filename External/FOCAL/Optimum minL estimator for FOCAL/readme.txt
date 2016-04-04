readme.txt 2015/8/17

To run the FOCAL analysis

1. Run FOCALmain.m in matlab. 
2. In the dialogue box select the sample localization table. The table is produced by the freely available software RapidSTORM and from a dSTORM acquisition of mouse cortex cell. Images were corrected for static background by a temporal median filter which is available in this webpage as well. After localization by RapidStorm  the drift is corrected and the anomalous localizations filtered out by a slicing filter as described in the supplementary material.
3. The program perfoms the uncertainty analysis and shows the result in a figure. This helps to decide which of the FOCAL generated localization table should be used to reconstruct super resolution image. 
