segTestAra folder

Contains python script files and results for testing the use of Phytotyping4D (set of python code) to create pipelines for leaf feature extraction using example images of arabidopsis.

List of python scripts:
* segTestAra.py: the script to read raw images and run segmentation function 
(depth images are fixed to 1 matrix of all element = 1). 

* P3D_segment.py: modified from P4D_segment.py

* P4D_help.py: set of useful functions for extracting relevant features and running batch jobs. (direct copy of P4D_help.py, minor changes in library module names)


List of folders:
Focus: for raw images. Currently only 9 png files mannually cut from one images containing 20 arabidopsis plants (description_for_plant_position_coding.png).  

Pics: for storing processed images after segmentation.

Nums: for storing plant leaf parameters extracted. 


