
1. First thing first, you need to have reasonable amount of IDL knowledge and have IDL installed in order to use this package. 
2. Type "make" and you should be able to use the idl procedures in this folder
3. If you want to add these procedures to your idl library, you should modify your idl path to include this folder, which is not necessary.
4. To clean up the installation, type "make clean".

In this folder, there are three IDL procedures that are essential to this package:
(a) detect_source.pro. It detects sources from an image.
(b) purge_match.pro. It takes the output from detect_source.pro and trims those that are already in the user supplied coordinates list.
(c) der_residual.pro. It derives contamination subtracted residuals for sources.
For how these procedures work, read Zhang et al., 2014 in preparation.

There are also a few supporting procedures:
(a) flood_shed.pro, maximas.pro, fmm_inpaint_image.pro, purge_area.pro. IDL wrappers to the corresponding c++ routines. 
(b) arrange_stamps1.pro. It arranges the stamps of images obtained from der_residual.pro into arrays of images. Needed only for running wrapper_cat.pro and wrapper_noncat.pro.
(c) gradient.pro. This routine is not written by me. It derives gradient of an image. I use IDL code for this purpose because IDL is very fast with arrays.

There are three additional IDL procedures that illustrate how to use this package. You would need some IDL knowledge to understand what these procedures do. You can certainly modify these procedures as you wish.
(a) testime.pro. It runs the three procedures listed above and print the time used of each procedure. 
(b) wrapper_cat.pro. It re-processes data that have "blended" flags (flag=3) in a SExtractor catalog, and re-measure their photometry using SExtractor on the contamination subtracted images.
(c) wrapper_noncat.pro. It searches for un-detected sources on a image and measure their photometry using SExtractor on the contamination subtracted images.
For (b) and (c), you would also need a reasonable amount of SExtractor experience to understand them. I also included a minimum amount of SE set-up files in the example folder, but the setups are not optimized or heavily tested for various purposes.

There are also a few supporting procedures:
(a) flood_shed.pro, maximas.pro, fmm_inpaint_image.pro, purge_area.pro. IDL wrappers to the corresponding c++ routines.     
(b) arrange_stamps1.pro. It arranges the stamps of images obtained from der_residual.pro into arrays of images. Needed only for running wrapper_cat.pro and wrapper_noncat.pro.
(c) gradient.pro. This routine is not written by me. It derives gradient of an image. I use IDL code for this purpose because IDL is very fast with arrays.

If you have questions, suggestions, or are very frustrated with these procedures, contact ynzhang@umich.edu. If you find this package useful, please cite Zhang et al., 2014, in preparation.

If you are part of the DES collaboration, before getting into the hassles of running this package youself, I strongly recommend trying testing data from http://umdes1.physics.lsa.umich.edu/deblend/ first. Also, please send an email to ynzhang@umich.edu describing your purpose. I ask for this because we might have processed the data you need or have the plan to process them.
