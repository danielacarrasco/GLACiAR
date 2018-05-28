# GLACiAR

Overview 
----------

``GLACiAR`` is an open-source python tool for simulations of source recovery and completeness in galaxy surveys.
For a more detailed explanation, please revise the paper https://arxiv.org/abs/1805.08985 .

Requirements
----------

To install ``GLACiAR``, users should download the source from GitHub.

Python and SExtractor are required to run the program. We also recommend DS9 to visualise the images.

We suggest users download Anaconda https://www.anaconda.com/download/ which includes all the softwares needed to run ``GLACiAR``. 

``GLACiAR`` uses a module that may not be included in most python installations, which is pysysp. 

Running GLACiAR 
----------

1. Download the source code from GitHub
2. Modify the 'parameters.yaml' file accordingly.
3. Create a directory that contains the images for all different fields and bands required. In this directory, a subfolder 'Results' will be created by ``GLACiAR``. All the files created will be saved here. 
5. Modify 'dropouts.py' if needed.
6. Run 'completeness.py'

Parameters
----------

Modify the parameters file 'parameters.yaml'.

**Parameters files**

- *n_galaxies:* Number of galaxies to place in each iteration (default = 100). Type = int.
- *n_iterations:* Number of iterations, i.e., the number of times the simulation is going to be run on each image for galaxies with the same redshift and magnitude.(default = 100).
- *mag_bins:* The numbers of magnitude bins wanted. For a simulation run from m1 = 24.0 to m2 = 25.0 in steps of 0.2 magnitudes, there will be 6 bins (default = 20).
- *min_mag:* Brightest magnitude of the simulated galaxies (default = 24.1).
- *max_mag:* Faintest magnitude of the simulated galaxies (default = 27.9).
- *z_bins:* The numbers of redshift bins wanted. For a simulation run from z1 = 9.5 to m2 = 10.5 in steps of 0.2 magnitudes, there will be 6 bins (default = 15).
- *min_z:* Minimum redshift of the simulated galaxies (default = 9.0).
- *max_z:* Maximum redshift of the simulated galaxies (default = 11.9).
- *n_bands:* How many filters the images have been observed in. If not specified, it will raise an error.
- *detection_band:* This is the band in which the objects are identified. The images taken in this band will be where the simulated galaxies are first put in. If not specified, it will raise an error.
- *bands:* Name of the bands. The detection band has to go first. If not specified, it will raise an error.
- *zeropoints:* Zeropoint value corresponding to each band. The default value is 25 for each band.
- *gain_values:* Gain values for each band. If not specified, it will raise an error.
- *list_of_fields:* Text file containing the name of the fields where the simulation is going to be run. If not specified, it will raise an error.
- *R_eff:* Effective radius in pixels. It is the half light radius, i. e., the radius within half of the light emitted by the galaxy is enclosed (default = 2.5).
- *size_pix:* Pixel scale for the images (default = 0.8).
- *path_to_images:* Directory where the images are located. The program will create a folder inside it with the results. If not specified, it will raise an error.
- *image_name:* Name of the images. They all should have the same name with the band written at the end. For example: "image_name+band.fits"
- *sersic_indices:* Sérsic indices for the galaxies (default = [1,4]).
- *fraction_type_of_galaxies:* Fraction of galaxies corresponding the the Sérsic indices given (default = [0.5,0.5]). 
- *min_sn:* Minimum signal to noise ratio n the detection band for a galaxy to be considered detected by SExtractor. (default = 8.0)
- *dropouts:* Yes or no. Boolean that indicates whether the user wants to run a dropout selection (default = False).

Files
----------

Files required to run ``GLACiAR``.


- Science images: Files with the observed images of the survey including all the fields and filters. It typically includes several files of each type, although one of each is enough.

- List: Text file with the names of the fields from the survey. This list is given as an input parameter. Its minimum length is one, i.e. the name of one field.

- SExtractor parameters: As previously explained, one of the steps of the code involves running ``SExtractor`` on the images (original and with simulated galaxies). To run the software, a file defining the parameters is required. There is an example provided under ``SExtractor\_files``, but we recommend for the user to change it according to their data.

- RMS maps or weight maps: They describe the noise intensity at each pixel in relation to the science image. Although they are necessary only if required for the SExtractor parameters, it is strongly recommended that one of these is used in order to improve the source detection.

- PSF: Image of the PSF that we'll be used to convolve with the simulated galaxy. It depends on the instrument, telescope, and filter. Some PSFs are included with ``GLACiAR``, but users should add theirs if needed.

Modules
----------

The code is made of several modules which are called by a main module, ``completeness.py``. A description of them follows.

- ``write_conf_files.py``: This module is used in ``run_sextractor.py`` and it reads the SExtractor parameters file, which has the criteria to identify the sources. Subsequently, it writes new temporary configuration files for each band. These are used when SExtractor is run.

- ``run_sextractor.py``: This module is called by ``completeness.py``, and it calls ``write_conf_files.py`` to write the SExtractor parameters file for each band. It then distinguishes between the detection band and the other bands. It runs first on the detection band, identifies the sources, and runs then SExtractor in dual mode. This means that the photometry is performed on the location of the sources found in the detection band. This is crucial so then all the sources can be compared and have their information in all the bands.

- ``creation_of_galaxy.py``: This module performs most of the mathematics processes involved in the code. It calculates the flux for each of the pixels following the Sersic profile. Accordingly, it performs all the required operations for the profile. It also generates the random positions for the simulated sources, generates the mock spectrum and calculates the expected flux for that spectrum in a given band.

- ``blending.py``: It identifies the detection and blending status of a source. It does this by comparing segmentation maps, the ones from the original science image and the ones with the simulated sources. It also extracts information from the catalogues when needed in order to compare magnitudes for the blending status. It runs for each iteration, so it retrieves a list with all the galaxies in that iteration and their status.

- ``dropouts.py``: If the parameters dropouts is set to "Yes" by the user, this module is called by ``blending.py``. This is a more specific module as it is only useful in case of dropouts. The one included here has the selection criteria from BoRG (Bernard et al. 2016), but it can be modified. It receives information on the magnitudes and status of all sources and then, depending on their colours and signal-to-noise, classify the objects.

- ``plot_completeness.py``: This module produces the plots. It always produces the plot of the completeness *C(m)* as a function of redshift and magnitude. Depending on the requirements by the user, it can produce a plot of *S(s,m)* and *S(s,m)C(m)* as well.

- ``completeness.py``: Main module. It manages the files and calls ``creation_of_galaxy.py`` to perform the mathematical operations in order to calculate the flux for the simulated galaxies and expected magnitudes. It also calls ``run_sextractor.py`` to run the source identification software on the images as well as ``plot_completeness.py`` in order to produce plots. On its own, this module produces the features of the artificial galaxies according to the input parameters. It opens the images, create the stamps with the calculated fluxes, places the galaxies in the given positions, and adds them to the science images. Then calls the module to run SExtractor, and records the statistics regarding the recovery of sources and also the individual status and extracted properties of the simulated objects. It then produces the final tables and plots.



How it works
----------

In order to estimate the completeness of a survey, it is necessary to quantify the fraction of galaxies that are not being detected. 
``GLACiAR`` simulates artificial galaxies with particular features (inputted by the user through a parameters file) and adds them to the images of the survey. Afterwards, a source identification software is run over the original science images and the images with the simulated sources. Both catalogues are compared and the fraction of recovered artificial galaxies is measured. 

This process is repeated a given number of iterations for each magnitude bin, redshift bin, and number of fields.
For each image, there is also a number of galaxies placed in it.

An artificial galaxy stamp is created following the parameters given by the user, such as magnitude, redshift, and radius. This source is then added to the science image in random positions and the new image is saved as a new file. This file is created for each band for which the magnitude is calculated as well. To calculate the expected magnitude in each band, the program generates a synthetic spectrum, and then measures the expected flux in each band. All galaxies in one iteration have the same UV slope, which is done in order to save time.

Now, the new image information is ready to be extracted. In order to do so, \texttt{SExtractor} runs on the original science image first, creates a catalogue, and then runs over all the new images with the same SExtractor parameters. 

SExtractor runs in dual mode, this means that there is an identification band. All the sources are identified in this band first and then the information of the source (or absence of) in said position is calculated from the images on each band.

The way in which the catalogues are compared is by checking a grid centered in the input position of the galaxies. Depending on whether there is a source detected in that grid in either the segmentation map of the original science image or the segmentation map of the science image with the added artificial galaxies, a status will be assigned to the artificial galaxies. The possible statuses are:

- Detected and isolated.
- Detected and blended with a fainter object.
- Detected and blended with a brighter object.
- Detected below the S/N threshold.
- Not detected.

If required by the user, the program can run a redshift selection criteria which we set here as an example, and it is from Bernard et al (2016.)

Finally, plots are generated for C(m), S(z,m), and C(m)S(z,m).

Outputs
----------

The code produces a set of files, images and tables. Some of them are deleted for space reasons, while others are kept as a final result. We outline them below in order of appearance.

- **New images**.

The new images are the created images that contain the original observed data from the survey with the simulated galaxies added. Each iteration in the simulation produces one image. These are saved and then deleted immediately after SExtractor runs on them.

- **SExtractor Catalogues**

The source identification process generates two sets of files. One of them is the catalogues (described here), and the other one is segmentation maps.

The catalogues include a list with all the sources identified in the images and their parameters. Due to the structure of the code, the first file created for a field of the survey contains the information of the real sources as it is run on the original science image. Besides this, for each iteration of the simulation a catalogue with the information of the sources is created. The header of the catalogue is shown below:

1. NUMBER   Running object number   
2.   FLUX_ISO   Isophotal flux   [count]
3.   FLUXERR_ISO   RMS error for isophotal flux   [count] 
4.   MAG_ISO   Isophotal magnitude   [mag] 
5.   MAGERR_ISO   RMS error for isophotal magnitude   [mag] 
6.   FLUX_APER   Flux vector within fixed circular aperture(s)   [count] 
11.   FLUXERR_APER   RMS error vector for aperture flux(es)   [count] 
16.   MAG_APER   Fixed aperture magnitude vector   [mag] 
21.   MAGERR_APER   RMS error vector for fixed aperture mag.   [mag] 
26.   FLUX_AUTO    Flux within a Kron-like elliptical aperture   [count] 
27.   FLUXERR_AUTO   RMS error for AUTO flux   [count] 
28.   MAG_AUTO   Kron-like elliptical aperture magnitude   [mag] 
29.   MAGERR_AUTO   RMS error for AUTO magnitude   [mag] 
30.   KRON_RADIUS   Kron apertures in units of A or B                      
31.   BACKGROUND   Background at centroid position   [count] 
32.   X_IMAGE   Object position along x   [pixel] 
33.   Y_IMAGE    Object position along y   [pixel] 
34.   ALPHA_J2000   Right ascension of barycenter (J2000)   [deg] 
35.   DELTA_J2000   Declination of barycenter (J2000)   [deg] 
36.   A_IMAGE   Profile RMS along major axis   [pixel] 
37.   B_IMAGE   Profile RMS along minor axis   [pixel] 
38.   THETA_IMAGE   Position angle (CCW/x)   [deg] 
39.   FWHM_IMAGE   FWHM assuming a gaussian core   [pixel] 
40.   FWHM_WORLD   FWHM assuming a gaussian core   [deg] 
41.   FLAGS   Extraction flags 
42.   CLASS_STAR   S/G classifier output 
43.   FLUX_RADIUS   Fraction-of-light radii   [pixel]


- **Segmentation maps**

A segmentation map is a map with the definition of the location of a source and its borders found by a source identification software. In this case, they are produced by SExtractor. For ``GLACiAR``, they can be classified in two groups: the segmentation maps from the original science images of the survey, and the ones from the images that include the simulated galaxies.

Together with the old catalogues of the sources a segmentation map is created. This is a .fits file and it is kept.

The images with the simulated galaxies will produce new segmentation maps with the same characteristics. The only difference with the catalogue described above is the inclusion of the new detected sources.
This file is discarded in order to save space.

- **GLACiAR catalogues**

The main results produced by ``GLACiAR`` can be summarised in three tables. These contain the information of the inserted galaxies, including their input properties and their output status and properties as well. They are described below.

The first table contains information on the statistics of the results in terms of the recovered sources and the dropouts (if specified). It counts the amount of sources inserted per redshift bin and magnitude bin and it keeps track of the amount that corresponds to each detection status. It also calculates the amount of total recoveries over all the inserted simulated sources, C(m). If applicable, it counts the number of dropouts for each bin, and calculates the fraction of these over the total amount of inserted galaxies S(z,m), and over the recovered dropouts in the required redshift range over the number of recovered simulated galaxies C(m)S(z,m). Below there is an example of its structure with a brief description of the columns.

   
1. z: Input redshift of the simulated galaxy.
2. m: Median value of the magnitude bin.
3. N_Obj: Number of objects inserted for the corresponding redshift and magnitude bin in all the iterations.
4. S(0): Number of artificial sources recovered by \texttt{SExtractor} that were isolated.
5. S(2,1): Number of artificial sources recovered that were blended with a fainter object.
6. S(-1): Number of artificial sources recovered that were blended with a brighter object.
7. S(-2): Number of artificial sources that were detected by \texttt{SExtractor} with a $S/N$ under the required threshold.
8. S(-3): Number of artificial sources that were not detected by \texttt{SExtractor}.
9. N_Rec: Number of recovered artificial sources, S(0)+S(2,1).
10. N_Drop: Number of artificial sources that passed the dropout selection criteria.
11. Rec: Fraction of not recovered artificial sources : N_Rec/N_Obj.
12. Drop: Fraction of artificial sources that passed the selection criteria: N_Drop/N_Obj.


Complementary to the previous table, the algorithm produces a table with the each one of the inserted sources, their positions, the input magnitude, the blending status, and the detection status. Several tables (one for each redshift step) are produced with all the galaxies that were placed in the simulations at that redshift. It provides an insight to understand the characteristics or reasons to detect or miss an object. It also yields their detected magnitude in the detection band, and their size. In summary, it gives us information about how this source is detected instead of the input information about it. 

1. Initial Mag: Magnitude corresponding to the input flux for the star. This is not the same as Input Magnitude since the input magnitude changes depending on the beta value and size of the object.
2. Iteration: Iteration number.
3. ID Number: Identification number given by SExtractor after it runs on the image with the simulated galaxies. This number is unique for every iteration for a given magnitude and redshift.
4. Input Magnitude: Magnitude corresponding to the added flux inside all the pixels that the source includes.
5. Output Magnitude: Magnitude of the source found with SExtractor after it runs on the image with the simulated galaxies.
6. Id Status: Integer number that indicates whether a source has been recovered and/or is blended.


One last table, which is useful for redshift selection, is produced. Given that the number of bands is variable  this table is released in a Python-specific compact binary representation (using the pickle module). It does contain the ID of the object, input magnitude, status, magnitudes in all bands, and S/N for each band as well. This is an important file for redshift estimations/selection techniques. Similar to the previous tables, it contains information about the output/measured characteristics of the detected objects. It contains all the magnitudes in different bands, so it is especially useful for photometric estimations. Unlike the other two tables, this is not an ASCII file as that is not efficient and it requires too many resources.


- **Plots**

``GLACiAR`` also produces a plot of the completeness and two extra plots if the boolean dropouts parameters is set to True. The first plot corresponds to the completeness function C(m) as a function of the magnitude and the redshift. The second and third plot are the S(z,m) and S(z,m)C(m). This is only produced in case the dropout technique is applied, but given the tables produced by ``GLACiAR``, it is easily calculable with the final catalogues.
