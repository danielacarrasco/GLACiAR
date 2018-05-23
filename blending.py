import numpy as np
import pickle
import dropouts
import os
from astropy.io import fits


def main(path_to_cat, niter, detection_band, cat, m1, redshift, xpos, ypos,
         xpos_oc, ypos_oc, segm_science, m_oc, f_oc, id_oc, input_mag, zp,
         bands, min_sn, dp, fw, fw2):
    """
    Uses the information from the new and old catalogs and segmentation
    maps to find the new sources and label them according to their
    identification and blending statuses.
    Args:
        path_to_cat (string) = Path to the folder with the science images.
                               Given in the parameters file.
        niter (integer) = iteration number.
        detection_band (string) = Name of the detection band given in the
                                  parameters file.
        cat (string) = Name of the field for which the simulation is run.
        m1 (float) = Initial input magnitude for the simulated galaxy in the
                     detection band.
        redshift (float) = Redshift for the simulated galaxy.
        xpos (int array) = Position of the simulated galaxy in the x axis.
        ypos (int array) = Position of the simulated galaxy in the y axis.
        xpos_oc (float array) = Array with the position in the x axis of the
                                centre for all the sources identified in the
                                original science image.
        ypos_oc (float array) = Array with the position in the y axis of the
                                centre for all the sources identified in the
                                original science image.
        segm_science (array) = Segmentation map produced by SExtractor for the
                               identified sources from the science image.
        m_oc (float array) = Array with the AB AUTO magnitude for the sources
                             identified by SExtractor in the science image.
        f_oc (float array) = Array with the ISO flux for the sources
                             identified by SExtractor in the science image.
        id_oc (int array) = Array with the ID assigned by SExtractor for each
                            one of the sources identified in the science image.
        input_mag (float) = Expected magnitude of the artificial galaxies in
                            the detection band.
        zp (float array) = Zeropoint values for all bands as given in the
                           input parameters file.
        bands (string array) = Name of the bands in which the artificial
                               galaxies will be simulated. Given in the input
                               parameters file.
        min_sn (float) = Minimum S/N ratio in the detection band for an object
                         to be considered detected by SExtractor. Given in the
                         input parameters file.
        dp (boolean) = Boolean that indicates whether the user requires to run
                       a dropout selection. Given in the input parameters file.
                       If True, the dropouts.py module will be used.
        fw (text file) = File in which the information about the artificial
                         sources will be saved. 'RecoveredGalaxies_cat_z#.cat'
        fw2 (text file) = File in which the information about the artificial
                         sources obtained by SExtractor will be saved.
    Returns:
        identified (int) = Number of artificial galaxies from the current
                           iteration that are detected by SExtractor and that
                           are isolated.
        blended_b (int) = Number of artificial galaxies from the current
                          iteration that are detected by SExtractor and are
                          blended with previously detected brighter sources.
        blended_f (int) = Number of artificial galaxies from the current
                          iteration that are detected by SExtractor and are
                          blended with previously detected fainter sources.
        not_indentified_sn (int) = Number of artificial galaxies from the
                                   current iteration that are detected by
                                   SExtractor but are considered not identified
                                   because their S/N is below min_sn.
        not_indentified (int) = Number of artificial galaxies from the current
                                iteration that are not detected by SExtractor.
        drops (int) = Number of artificial galaxies from the current iteration
                      that passed the redshift selection criteria from
                      'dropouts.py'. If drops is set to False the value is 0.
    """
    # Open segmentation maps from simulated images, save the data,
    # and delete the fits file.
    segm_new_cat = fits.open(path_to_cat + 'Results/SegmentationMaps/'
                             'Segmentation_maps_i' + str(niter) + '_' +
                             detection_band+'.fits', ignore_missing_end=True)
    segm_sim = segm_new_cat[0].data  # Data from new images.
    segm_new_cat.close()
    os.remove(path_to_cat + 'Results/SegmentationMaps/Segmentation_maps_i' +
              str(niter) + '_' + detection_band + '.fits')

    # Catalog with the identified sources from the simulated images.
    f = open(path_to_cat + 'Results/Dropouts/source_' + cat + '_mag' +
             str(m1) + '_z' + str(redshift) + '_i' + str(niter) + '_' +
             detection_band + '.cat')
    k = f.readlines()
    f.close()

    # Information from SExtractor for the new sources (science image +
    # simulated galaxies).
    id_mgal = [int(line.split()[0]) for line in k[27:]]  # ID
    f_gal = [float(line.split()[1]) for line in k[27:]]  # Isophotal flux
    ef_gal = [float(line.split()[2]) for line in k[27:]]  # RMS error for flux
    m_gal = [float(line.split()[27]) for line in k[27:]]  # AUTO magnitude
    sn_mgal = np.array(np.array(f_gal)/np.array(ef_gal))  # RMS error for mag
    xpos_nc = [float(line.split()[32]) for line in k[27:]]  # Position in x
    ypos_nc = [float(line.split()[31]) for line in k[27:]]  # Position in y
    radius = [float(line.split()[42]) for line in k[27:]]  # Radii

    # Convert the previous arrays into np.arrays
    xpos = np.array(xpos).astype(int)
    ypos = np.array(ypos).astype(int)
    xpos_oc = np.array(xpos_oc)
    ypos_oc = np.array(ypos_oc)
    xpos_nc = np.array(xpos_nc)
    ypos_nc = np.array(ypos_nc)
    m_oc = np.array(m_oc)
    f_oc = np.array(f_oc)
    id_oc = np.array(id_oc)
    id_mgal = np.array(id_mgal)
    m_gal = np.array(m_gal)

    # Initialise counters. These will contain information about the number of
    # artificial galaxies and their detection status after running SExtractor.
    identified = 0  # Identified artificial sources.
    blended_b = 0  # Artificial sources blended with a brighter object.
    blended_f = 0  # Artificial sources blended with a fainter object.
    not_indentified_sn = 0  # Artificial sources not detected due to low S/N.
    not_indentified = 0  # Artificial sources not detected by Sextractor.
    drops = 0.  # Artificial sources classified as a dropout.

    margin = 3  # Number of pixels from centre in which search is performed.

    # Array with only 100 as values, which will be replaced by the status code.
    status = np.zeros(len(xpos)) + 100

    # Array with only zeroes as values, which will be filled out by the ID of
    # the artificial galaxy.
    id_nmbr = np.zeros(len(xpos))
    i_mag = np.zeros(len(xpos)) + input_mag  # array with input mag value.

    # Open a text file (.reg) with the formatting to be used as region file on
    # DS9. It shows location where the artificial source was originally placed.
    g = open(path_to_cat + 'Results/SegmentationMaps/region_' + cat + '_mag' +
             str(m1) + '_z' + str(redshift) + '_i' + str(niter) + '.reg', "w")

    # Loop for each artificial source.
    for i in xrange(len(xpos)):
        # The IF below searches for any value diferent than zero in the
        # segmentation maps of the new images (science + artificial sources)
        # within an square centered in the input position for the artificial
        # galaxy and a side of 2*margin. This is done by evaluating the sum of
        # the pixels' values. It enters the IF if the value is != 0.
        if np.sum(segm_sim[xpos[i]-margin:xpos[i]+margin,
                           ypos[i]-margin:ypos[i]+margin]) != 0:
            # array with square where search is performed.
            id_mgi_aux2 = segm_sim[xpos[i]-margin:xpos[i]+margin,
                                   ypos[i]-margin:ypos[i]+margin]
            # ID of the source in the search region is recorded.
            id_nmbr[i] = np.max(id_mgi_aux2)
            # ID value for pixels from the science image in the position of
            # the newly found source.
            id_mgi2 = segm_science[segm_sim == id_nmbr[i]]
            idsimgal = len(id_mgi2)  # Number of pixels the source encompasses.
            # ID of source previously in the position of newly found source.
            id_mgi = np.max(id_mgi2)
            # Number of pixels of the simulated galaxy that are overlapped
            # with the old source.
            w2 = np.where(id_mgi2 != 0)[0]
            # If the S/N of the newly identified source is above the required,
            # threshold, enter the IF below.
            if sn_mgal[id_mgal == id_nmbr[i]] >= min_sn:
                # If there wasn't any source on the pixels where the new source
                # was found, enter the IF below.
                # These are objects detected and isolated.
                if np.sum(id_mgi2) == 0:
                    status[i] = 0  # Status for detected and isolated sources.
                    # Line with information of the artificial source.
                    line = ('%6s' % str(m1) + '\t%6s' % str(int(niter)) +
                            '\t%9s' % str(int(id_nmbr[i])) + '\t%9s' %
                            str(float(input_mag)) + '\t%9s' %
                            str(float(m_gal[id_mgal == id_nmbr[i]])) +
                            '\t%9s' % str(int(0)) + '\n')
                    # Write line in the file 'RecoveredGalaxies_cat_z#.cat'
                    fw.writelines(line)
                    # Write the region file for the artificial sources with
                    # status=0 in colour green. The position is the input one.
                    g.write("circle %s %s 11 #color=green width=2\n" %
                            (ypos[i], xpos[i]))
                    # Counter for identified and isolated galaxies.
                    identified = identified + 1
                # If there was a source previously on the pixels where the new
                # source was found, enter the IF below.
                # These are objects detected and blended.
                else:
                    id_blended = id_mgi  # ID of the source previously there.
                    # Find the magnitude of the old source and compare it to
                    # the input magnitude. If the source is brighter, enter
                    # the IF below.
                    if m_oc[id_blended == id_oc] <= input_mag:
                        # Enter the IF below If the old source's flux is
                        # smaller than 75% the input flux of the artificial
                        # source, AND if the number of pixels of the newly
                        # identified source that were previously occupied by
                        # another source is 25% or less the total number of
                        # pixels of said source (the overlap of the new source
                        # is 25% its original size). We consider this as if the
                        # artificial source was blended with a fainter object.
                        if ((f_oc[id_blended == id_oc] <=
                             0.75*10**((zp[0]-input_mag)/2.5)) and
                           (len(w2) <= 0.25*idsimgal)):
                            # Status code for blended with fainter object.
                            status[i] = 2
                            line = ('%6s' % str(m1) + '\t%6s' %
                                    str(int(niter)) + '\t%9s' %
                                    str(int(id_nmbr[i])) + '\t%9s' %
                                    str(float(input_mag)) + '\t%9s' %
                                    str(float(m_gal[id_mgal == id_nmbr[i]])) +
                                    '\t%9s' % str(int(2)) + '\n')
                            # Write in the file 'RecoveredGalaxies_cat_z#.cat'
                            fw.writelines(line)
                            # Write the region file for the artificial sources
                            # with status=2 in colour blue.
                            g.write("circle %s %s 11 #color=blue width=2\n" %
                                    (ypos[i], xpos[i]))
                            # Counter for sources identified and blended with
                            # faintergalaxies.
                            blended_f = blended_f + 1
                        # If the flux or overlap conditions aren't true, do the
                        # following
                        else:
                            # Status code for blended with brighter object.
                            status[i] = -1
                            line = ('%6s' % str(m1) + '\t%6s' %
                                    str(int(niter)) + '\t%9s' %
                                    str(int(id_nmbr[i])) + '\t%9s' %
                                    str(float(input_mag)) + '\t%9s' %
                                    str(float(m_gal[id_mgal == id_nmbr[i]])) +
                                    '\t%9s' % str(int(-1)) + '\n')
                            # Write in the file 'RecoveredGalaxies_cat_z#.cat'
                            fw.writelines(line)
                            # Write the region file for the artificial sources
                            # with status=-1 in colour red.
                            g.write("circle %s %s 11 #color=red width=2\n" %
                                    (ypos[i], xpos[i]))
                            # Counter for sources identified and blended with
                            # brighter galaxies.
                            blended_b = blended_b + 1
                    # if the magnitude of the old source is fainter than the
                    # input magnitude, enter the IF below.
                    else:
                        # Status code for blended with fainter object.
                        # Different from 2 as it purely compares magnitudes.
                        status[i] = 1
                        line = ('%6s' % str(m1) + '\t%6s' % str(int(niter)) +
                                '\t%9s' % str(int(id_nmbr[i])) + '\t%9s' %
                                str(float(input_mag)) + '\t%9s' %
                                str(float(m_gal[id_mgal == id_nmbr[i]])) +
                                '\t%9s' % str(int(1)) + '\n')
                        # Write in the file 'RecoveredGalaxies_cat_z#.cat'
                        fw.writelines(line)
                        # Write the region file for the artificial sources
                        # with status=1 in colour blue.
                        g.write("circle %s %s 11 #color=blue width=2\n" %
                                (ypos[i], xpos[i]))
                        # Counter for sources identified and blended with
                        # faintergalaxies.
                        blended_f = blended_f + 1

            # If the S/N of the newly identified source is above the required,
            # threshold, enter the IF below.
            else:
                # Status for detected sources with S/N below required min_sn.
                status[i] = -2
                line = ('%6s' % str(m1) + '\t%6s' % str(int(niter)) +
                        '\t%9s' % str(int(id_nmbr[i])) + '\t%9s' %
                        str(float(input_mag)) + '\t%9s' %
                        str(float(m_gal[id_mgal == id_nmbr[i]])) +
                        '\t%9s' % str(int(-2)) + '\n')
                # Write in the file 'RecoveredGalaxies_cat_z#.cat'
                fw.writelines(line)
                # Write the region file for the artificial sources with
                # status=-2 in colour red.
                g.write("circle %s %s 11 #color=red width=2\n" %
                        (ypos[i], xpos[i]))
                # Counter for sources detected by SExtractor but that did not
                # meet the required threshold.
                not_indentified_sn = not_indentified_sn + 1
        # If all values of the new segmentation map within the search grid are
        # zero, the object has not been detected by SExtractor.
        else:
            status[i] = -3  # Status for sources not detected by SExtractor.
            line = ('%6s' % str(m1) + '\t%6s' % str(int(niter)) +
                    '\t%9s' % str(int(0)) + '\t%9s' % str(float(input_mag)) +
                    '\t%9s' % str(float(-99.000)) + '\t%9s' % str(int(-3)) +
                    '\n')
            # Write in the file 'RecoveredGalaxies_cat_z#.cat'
            fw.writelines(line)
            # Write the region file for the artificial sources with status=-3
            # in colour red.
            g.write("circle %s %s 11 #color=red width=2\n" %
                    (ypos[i], xpos[i]))
            # Counter for sources not detected by SExtractor.
            not_indentified = not_indentified + 1

    # Close the .reg file.
    g.close()

    # Initialise array for ISO magnitudes, AUTO magnitudes, and S/N in all
    # bands measured by SExtractor for the new sources.
    mag_iso = np.zeros((len(id_mgal), len(bands)))
    mag_auto = np.zeros((len(id_mgal), len(bands)))
    sn = np.zeros((len(id_mgal), len(bands)))

    # Loop for the number of bands used in the simulation.
    for j in xrange(len(bands)):
        # Open the catalog with the identified sources from the simulated
        # images (science + artificial sources) for each band.
        f1 = open(path_to_cat + 'Results/Dropouts/source_' + cat + '_mag' +
                  str(m1) + '_z' + str(redshift) + '_i' + str(niter) + '_' +
                  bands[j] + '.cat')
        k1 = f1.readlines()
        f1.close()

        # Save the information on the ISO mag, AUTO mag, and S/N for each band.
        mag_iso[:, j] = [float(line.split()[3]) for line in k1[27:]]
        mag_auto[:, j] = [float(line.split()[27]) for line in k1[27:]]
        sn[:, j] = [float(line.split()[1])/float(line.split()[2]) for
                    line in k1[27:]]

    # Save the data only for sources identified as artificial sources.
    mag_iso2 = mag_iso[id_nmbr.astype(int)-1, :]
    mag_auto2 = mag_auto[id_nmbr.astype(int)-1, :]
    sn2 = sn[id_nmbr.astype(int)-1, :]

    # Initialise drops, which stays as zero if dropouts is False
    drops = 0

    # Run dropout module if dropout parameter is set to True.
    if dp is True:
        drops, drops_array = dropouts.main(mag_iso2, mag_auto2, sn2, status)

    # Save stacked array with the dimensions needed according to the number
    # of bands. It contains the information of the new sources by SExtractor
    # plus their detection status and input magnitude.
    final_array_aux = np.c_[id_nmbr, i_mag, status, drops_array,
                            mag_iso2, mag_auto2, sn2]
    final_array = np.matrix.transpose(final_array_aux)

    # Write a pickled representation of the data of the new sources to the
    # file.
    # WARNING: this is not readable data.
    pickle.dump(final_array, fw2)

    # Return number of sources for each of the following categories.
    return identified, blended_b, blended_f, not_indentified_sn,\
        not_indentified, drops

if __name__ == "__main__":
    main()
