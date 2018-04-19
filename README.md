# first-project
This is the first project of my freshman year, the code here is messy and i didn't upload the finished code. The main function is to get disparity map of two rectified rgb images (normal focal length, wide and tele). If you want to execute stereo matching with color and monochrome cameras, maybe the first step is to grayscale (decolorization) the picture, and meanwhile preserve the contrast.

Prerequisites
Windows 7-10
OpenCV 3.1 (or above)


Pipline
cost computation->cost aggregation->post process

..CC.cpp are the methods for cost computation

For cost aggregation, non-local is performed well.
Yang Q. A non-local cost aggregation method for stereo matching[J]. 2012, 157(10):1402-1409.

Post process
Min D, Choi S, Lu J, et al. Fast Global Image Smoothing Based on Weighted Least Squares.[J]. IEEE Transactions on Image Processing, 2014, 23(12):5638-5653.


