"""Main"""

import numpy as np
import random
import os
import yaml
import pysysp
import run_sextractor
import creation_of_galaxy
import blending
import plot_completeness
import dropouts
from astropy.io import fits
from astropy.convolution import convolve

# Read input parameters from the file 'parameters.yaml'.
# Assign a default value when possible if they are not given by the user.

stream = file('parameters.yaml', 'r')
parameters = yaml.load(stream)
if parameters['n_galaxies'] is None:
    parameters['n_galaxies'] = 100
if parameters['n_iterations'] is None:
    parameters['n_iterations'] = 100
if parameters['mag_bins'] is None:
    parameters['mag_bins'] = 20
if parameters['min_mag'] is None:
    parameters['min_mag'] = 24.1
if parameters['n_galaxies'] is None:
    parameters['max_mag'] = 27.9
if parameters['z_bins'] is None:
    parameters['z_bins'] = 16
if parameters['min_z'] is None:
    parameters['min_z'] = 7.0
if parameters['max_z'] is None:
    parameters['max_z'] = 9.2
if parameters['lambda_detection'] is None:
    raise ValueError('Value of lambda required.')
if parameters['n_bands'] is None:
    raise ValueError('Number of bands needed.')
if parameters['detection_band'] is None:
    raise ValueError('Input detection band.')
if parameters['bands'] is None:
    raise ValueError('Input name of the bands.')
if parameters['list_of_fields'] is None:
    raise ValueError('Input file with the fields.')
if parameters['R_eff'] is None:
    parameters['R_eff'] = 1.075
if parameters['beta_mean'] is None:
    parameters['beta_mean'] = -2.2
if parameters['beta_sd'] is None:
    parameters['beta_sd'] = 0.4
if parameters['size_pix'] is None:
    parameters['size_pix'] = 0.08
if parameters['types_galaxies'] is None:
    parameters['types_galaxies'] = 2
if parameters['ibins'] is None:
    parameters['ibins'] = 10
if parameters['ebins'] is None:
    parameters['ebins'] = 10
if parameters['path_to_images'] is None:
    raise ValueError('Input the directory path.')
if parameters['image_name'] is None:
    raise ValueError('Input the name of the images.')
if parameters['sersic_indices'][1] is None:
    parameters['sersic_indices'][1] = 4
if parameters['fraction_type_galaxies'] is None:
    parameters['fraction_type_galaxies'] = [0.5, 0.5]
if (parameters['zeropoints'] is None or
        len(parameters['zeropoints']) < parameters['n_bands']):
    parameters['zeropoints'] = np.zeros(int(parameters['n_bands']))+25
if (parameters['gain_values'] is None or
        len(parameters['gain_values']) < parameters['n_bands']):
    raise ValueError('Input gain values for each band.')
if parameters['dropouts'] is None:
    parameters['dropouts'] = 'No'
if parameters['de_Vacouleur'] is None:
    parameters['de_Vacouleur'] = 'No'


def create_stamps(n0, size_galaxy0, Re0, types_galaxies0, ebins0, ibins0):
    """
    Creates a galaxy following a Sersic profile.

    Args:
        n0 (int array) = Array with Sersic indexes.
        size_galaxy0 = Diameter of the galaxy stamp in pixels.
        Re0 (float) = Effective radius in pixels.
        types_galaxies0 (int) = Number of Sersic indexes required.
        ibins0 (int) = Number of possible inclinations for the
                         simulated galaxy.
        ebins0 (int) = Number of possible eccentricities for the
                         simulated galaxy.
    Returns:
        galaxy_grid (float array) = Stamp of a galaxy with the
                                        corresponding flux for each
                                        pixel.
        galaxy_grid_n4 (float array) = Stamp of galaxy with n=4 if
                                      "de_Vacouleur" is True. Otherwise, set
                                      to 0.
    """
    galaxy_grid = np.zeros((types_galaxies0, ebins0, ibins0, size_galaxy0,
                            size_galaxy0))
    galaxy_grid_n4 = np.zeros((size_galaxy0, size_galaxy0))
    ivalues = np.arange(0, 0.5, 0.5/(ibins0)) * np.pi
    evalues = np.arange(0, 1., 1. / (ebins0))
    for i in xrange(types_galaxies0):
        print 'galaxy', n0[i]
        bn = creation_of_galaxy.get_bn(n0[i])
        if (parameters['de_Vacouleur'] is True and n0[i] == 4):
            galaxy_grid_n4[:, :] = creation_of_galaxy.makeSersic(n0[i], bn,
                                                                 Re0, 0, 0,
                                                                 size_galaxy0)
        else:
            for j in xrange(ebins0):
                for k in xrange(ibins0):
                    print 'galaxy', n0[i], j, k
                    galaxy_grid[i, j, k, :, :] = creation_of_galaxy.makeSersic(
                                                  n0[i], bn, Re0, evalues[j],
                                                  ivalues[k], size_galaxy0)
    return galaxy_grid, galaxy_grid_n4


def open_images(image_name, name_band):
    """
    Opens science images.

    Args:
        image_name (string) = Name of the science image.
        name_band (string) = Name of the band in which the science image is
                            taken.
    Returns:
        obs_data (float array) = Data from the science image. Each cell
                                contains the flux for that pixel.
        head_data (string array) = Array with the header from the science image
                                  opened.
    """
    hdu_list = fits.open(image_name + name_band+'.fits',
                         ignore_missing_end=True)
    obs_data = hdu_list[0].data
    head_data = hdu_list[0].header
    hdu_list.close()
    return obs_data, head_data


def place_gal(n0, ngal, frac_n, e0, i0, flux0, frame0, x0, y0, s0, gal_g,
              gal_g_n4, psf0):
    """
    Args:
        n0 (int array) = Sersic indices
        ngal (int) = Number of galaxies for each iteration
        frac_n (float) = Fraction of galaxies with n0
        e0 (float) = Eccentricity. Varies between 0 and 1
        i0 (float) = Inclination angle in radians.
                     Varies between 0 and Pi/2.
        flux0 (float) = Total corresponding flux for the artificial galaxy.
        frame0 (array) = Empty frame with the same shape as the science image.
        x0 (int array) = array with position along the x axis for new sources.
        y0 (int array) = array with position along the y axis for new sources.
        s0 (int) = radius of the galaxy stamp.
        gal_g (array) = Stamp for an artificial galaxy with n=4 and min_mag.
        gal_g_n4 (array) = Stamp for artificial galaxy with n!=4 and min_mag.
        psf0 (array) = PSF to convolve the artificial galaxy with.
    Returns:
        frame0: (array) = Frame with the same shape as science image. 
                          It has the artificial galaxy stamps in the
                          corresponding positions.
    """

    gn = 0
    # Loop for the amount of Sersic indices.
    for i in xrange(len(n0)):
        # Number of galaxies with given index.
        n_fraction = int(ngal * frac_n[i])
        # Run a loop for each galaxy with the given index.
        for j in xrange(n_fraction):
            # Scale the galaxy flux with the real flux.
            if parameters['de_Vacouleur'] is True and n0[i] == 4:
                galaxy = gal_g_n4[:, :] * flux0
            else:
                galaxy = gal_g[i, e0[gn], i0[gn], :, :] * flux0
            # Convolve galaxy with the PSF.
            gconv = convolve(galaxy, psf0, normalize_kernel=True)
            # Add galaxy stamp to the empty frame with the centre in (x0,y0).
            frame0[int(x0[gn]-s0):int(x0[gn]+s0),
                   int(y0[gn]-s0):int(y0[gn]+s0)] = galaxy
            gn = gn + 1
    return frame0


def main():
    # Open the file with the name of the fields that the simulation
    # is going to be run for.
    f = open(parameters['list_of_fields'])
    k = f.readlines()
    f.close()
    cat = [str(line.split()[0]) for line in k]

    m_total = np.linspace(parameters['min_mag'], parameters['max_mag'],
                          parameters['mag_bins'])  # Magnitudes bin size.
    z_total = np.linspace(parameters['min_z'], parameters['max_z'],
                          parameters['z_bins'])  # Redshift bin size.
    # Empty array to save the completeness results.
    total_completeness = np.zeros((parameters['mag_bins'],
                                   parameters['z_bins']))
    # Empty array to save the dropouts selection results.
    total_dropouts = np.zeros((parameters['mag_bins'],
                               parameters['z_bins']))

    # Creates the directories that will contain the results if they
    # don't already exist.
    if not os.path.exists(parameters['path_to_images']+'Results'):
        os.makedirs(parameters['path_to_images']+'Results')
    if not os.path.exists(parameters['path_to_images']+'SciImages'):
        os.makedirs(parameters['path_to_images']+'SciImages')
    if not os.path.exists(parameters['path_to_images']+'Results'
                          '/SegmentationMaps'):
        os.makedirs(parameters['path_to_images']+'Results/SegmentationMaps')
    if not os.path.exists(parameters['path_to_images']+'Results/Plots'):
        os.makedirs(parameters['path_to_images']+'Results/Plots')
    if not os.path.exists(parameters['path_to_images']+'Results/images/'):
        os.makedirs(parameters['path_to_images']+'Results/images/')
    if not os.path.exists(parameters['path_to_images']+'/Results/Dropouts'):
        os.makedirs(parameters['path_to_images']+'/Results/Dropouts')

    for i in xrange(len(cat)):
        # File to save the stats for the field.
        fws = open(parameters['path_to_images']+'Results/'
                   'RecoveredStats_cat'+cat[i]+'.cat', 'w')
        header = ('z # m # N_Objects # St = 0 # St = 2 # St = 1 # St = -1'
                  '# St = -2 # St = -3 # N_Recovered # N_Droputs # Recovered '
                  '# Dropouts\n')
        fws.writelines(header)
        name_hdu_list1 = (parameters['path_to_images'] +
                          parameters['image_name']+cat[i]+'_')

        run_sextractor.science_image(parameters['bands'],
                       parameters['detection_band'], parameters['zeropoints'],
                       parameters['gain_values'], parameters['path_to_images'],
                       parameters['image_name'], cat[i])

        # Open science image.
        obs_data_db, header_db = open_images(name_hdu_list1,
                                             parameters['detection_band'])
        # Open segmentation maps from original science image.
        segm_science_old = fits.open(parameters['path_to_images'] +
                                     '/SciImages/sources_'+cat[i]+'_' +
                                     parameters['detection_band'] +
                                     '_segm.fits', ignore_missing_end=True)
        segm_maps_old = segm_science_old[0].data
        segm_science_old.close()
        # Open SExtractor catalogue with the identified sources from the
        # original science image and save the following information
        # about the sources.
        # id_old_cat = SExtractor ID.
        # m_old_cat = AUTO magnitude.
        # f_old_cat = ISO flux.
        # xpos_old_cat = Centre position in the x axis.
        # ypos_old_cat = Centre position in the y axis.
        f_old_cat = open(parameters['path_to_images']+'/SciImages/sources_' +
                         cat[i]+'_'+parameters['detection_band']+'.cat', 'r')
        k_oc = f_old_cat.readlines()
        f_old_cat.close()
        id_old_cat = [int(line.split()[0]) for line in k_oc[27:]]
        m_old_cat = [float(line.split()[27]) for line in k_oc[27:]]
        f_old_cat = [float(line.split()[25]) for line in k_oc[27:]]
        xpos_old_cat = [float(line.split()[31]) for line in k_oc[27:]]
        ypos_old_cat = [float(line.split()[32]) for line in k_oc[27:]]
        # Run for all the required redshifts.
        for j in xrange(parameters['z_bins']):
            redshift = parameters['min_z'] + 0.2*j  # Redshift of the galaxies.
            Re = 7 * parameters['R_eff'] / (redshift+1.0)  # Effective radius scaled by z.
            # Diameter of the galaxy stamp.
            sizegalaxy = int(parameters['R_eff'] * parameters['size_pix'] * 325)
            s2 = sizegalaxy/2  # Radius of the galaxy stamp.
            # Stamp with the galaxies.
            galaxy_grid, galaxy_grid_n4 = create_stamps(
                                           parameters['sersic_indices'],
                                           sizegalaxy, Re,
                                           int(parameters['types_galaxies']),
                                           parameters['ebins'],
                                           parameters['ibins'])
            # Open file to save the recovery information of the
            # simulated galaxies.
            fwg = open(parameters['path_to_images']+'Results/Recovered'
                       'Galaxies_'+cat[i]+'_z'+str(redshift)+'.cat', 'w')
            # Open file to save the information from SExtractor
            # about the simulated galaxies in each band in serialised format.
            fwgb = open(parameters['path_to_images']+'Results/Recovered'
                        'GalaxiesBands_'+cat[i]+'_z'+str(redshift)+'.cat', 'w')
            # Empty arrays to save the stats.
            recovered = np.zeros(parameters['mag_bins'])
            total = np.zeros(parameters['mag_bins'])
            identified = np.zeros(parameters['mag_bins'])
            blended_b = np.zeros(parameters['mag_bins'])
            blended_f = np.zeros(parameters['mag_bins'])
            not_indentified_sn = np.zeros(parameters['mag_bins'])
            not_indentified = np.zeros(parameters['mag_bins'])
            drops = np.zeros(parameters['mag_bins'])
            # Run for all the required mangitudes.
            for k in xrange(parameters['mag_bins']):
                # Assigned magnitude of the galaxies.
                magnitude = parameters['min_mag'] + 0.2*k
                print redshift, magnitude
                for l in xrange(parameters['n_iterations']):
                    niter = l + 1  # Iteration number
                    # Initialize the input magnitude to enter the 'while',
                    # which will assing a random beta for the spectrum.
                    # With beta, the expected magnitude for the detection
                    # band is calculated.
                    m = np.zeros(parameters['n_bands'])
                    m[0] = -100
                    while (m[0] < (parameters['min_mag']) or
                           m[0] > (parameters['max_mag'])):
                        beta = random.gauss(parameters['beta_mean'],
                                            parameters['beta_sd'])
                        creation_of_galaxy.write_spectrum(
                           parameters['lambda_detection'], magnitude, beta,
                           redshift)
                        m[0] = creation_of_galaxy.mag_band(
                                parameters['bands'][0],
                                parameters['zeropoints'][0])
                    # Generate the spectrum for the actual input magnitude.
                    creation_of_galaxy.write_spectrum(
                           parameters['lambda_detection'], m[0], beta,
                           redshift)
                    # (xpos, ypos) is the position of the center of the
                    # galaxy in the image in pixels.
                    xpos, ypos = creation_of_galaxy.galaxies_positions(
                                 obs_data_db, parameters['n_galaxies'],
                                 sizegalaxy, Re)

                    i_rand = np.random.randint(0, parameters['ibins'],
                                               size=parameters['n_galaxies'])
                    e_rand = np.random.randint(1, parameters['ebins'],
                                               size=parameters['n_galaxies'])
                    # Repeat the following process for all the bands.
                    for h in xrange(parameters['n_bands']):
                        # Open PSF file to convolve with the galaxy.
                        hdu_psf = fits.open('Files/psf_' +
                                            parameters['bands'][h] + '.fits')
                        psf = hdu_psf[0].data
                        m_i = np.zeros(parameters['n_bands'])
                        m[h] = creation_of_galaxy.mag_band(
                                    parameters['bands'][h],
                                    parameters['zeropoints'][h])
                        # Empty array of the size of the science image
                        # where the simulated galaxies will be placed.
                        obs, head_obs = open_images(name_hdu_list1,
                                                    parameters['bands'][h])
                        frame = np.zeros(obs.shape)
                        # Expected flux of the simulated galaxy associated
                        # with the expected magnitude in that band.
                        flux = (10**((parameters['zeropoints'][h] -
                                m[h])/2.5))
                        # if the artificial galaxy is brighter than
                        # 50 mangitudes (some flux is expected, but very
                        # conservative margin), place it in the image.
                        if m[h] < 50:
                            frame_gal = place_gal(parameters['sersic_indices'],
                                                  parameters['n_galaxies'],
                                                  parameters[
                                                  'fraction_type_galaxies'],
                                                  e_rand, i_rand, flux, frame,
                                                  xpos, ypos, s2, galaxy_grid,
                                                  galaxy_grid_n4, psf)
                        # if the galaxy is fainter than 50 magnitudes
                        # (undetectable), save the original image.
                        else:
                            frame_gal = frame
                        # observed image + simulated galaxies
                        new_data2 = obs_data_db + frame_gal
                        outfile_conv = parameters['path_to_images']+'Results' \
                            '/images/sersic_sources_'+cat[i] + \
                            parameters['bands'][h]+'.fits'
                        # Save the new fits file with the simulated galaxies.
                        fits.writeto(outfile_conv, new_data2,
                                     head_obs, clobber=True)
                        print 'iteration #', niter, ', magnitude ', magnitude,\
                              ', redshift ', redshift, ', band', \
                              parameters['bands'][h], m[h]
                        # Run SExtractor on the new images. This generates
                        # catalogues with information about the
                        # detected sources.
                        run_sextractor.main(parameters['bands'][h],
                                            parameters['detection_band'],
                                            parameters['zeropoints'][h],
                                            parameters['gain_values'][h],
                                            parameters['path_to_images'],
                                            niter, parameters['image_name'],
                                            cat[i], magnitude, redshift)
                    # Find the number of sources for the different categories
                    # of detection statuses.
                    identified_aux, blended_b_aux, blended_f_aux,\
                        not_indentified_sn_aux, not_indentified_aux, \
                        drops_aux = blending.main(parameters['path_to_images'],
                                                  niter,
                                                  parameters['detection_band'],
                                                  cat[i], magnitude, redshift,
                                                  xpos, ypos, xpos_old_cat,
                                                  ypos_old_cat, segm_maps_old,
                                                  m_old_cat, f_old_cat,
                                                  id_old_cat, m[0],
                                                  parameters['zeropoints'],
                                                  parameters['bands'],
                                                  parameters['min_sn'],
                                                  parameters['dropouts'],
                                                  fwg, fwgb)
                    # Loop for al the magnitude bins in order to count the
                    # number of sources in each detection status for each bin.
                    for v in xrange(parameters['mag_bins']):
                        # Check where the calculated input magnitude lies and
                        # add them to the respective counter.
                        if (m[0] > (parameters['min_mag'] - 0.1 + v*0.2) and
                                m[0] <= (parameters['min_mag'] + 0.1 + v*0.2)):
                            total[v] = total[v] + parameters['n_galaxies']
                            identified[v] = identified[v] + identified_aux
                            blended_b[v] = blended_b[v] + blended_b_aux
                            blended_f[v] = blended_f[v] + blended_f_aux
                            not_indentified_sn[v] = (not_indentified_sn[v] +
                                                     not_indentified_sn_aux)

                            not_indentified[v] = (not_indentified[v] +
                                                  not_indentified_aux)
                            drops[v] = drops[v] + drops_aux
            fwg.close()
            fwgb.close()
            # Run for each magnitude bin.
            for s in xrange(parameters['mag_bins']):
                # Calculate fractions of different statuses.
                identified_total = (np.array(identified[s]) +
                                    np.array(blended_f[s]))
                identified_fraction = float(identified_total) /\
                    np.array(total[s])
                dropouts_fraction = float(drops[s])/np.array(total[s])
                total_completeness[s, j] = identified_fraction
                total_dropouts[s, j] = dropouts_fraction
                line = ('%6s' % str(redshift) + '\t%6s' % str(m_total[s]) +
                        '\t%6s' % str(int(total[s]))+'\t%9s' %
                        str(int(identified[s])) + '\t%9s' %
                        str(int(blended_b[s]))+'\t%9s' %
                        str(int(blended_f[s])) + '\t%9s' %
                        str(int(not_indentified_sn[s])) + '\t%9s' %
                        str(int(not_indentified[s])) + '\t%9s' %
                        str(int(drops[s])) + '\t%9s' %
                        str(int(identified_total)) + '\t%9s' %
                        str(float(identified_fraction)) + '\t%9s' %
                        str(float(np.array(dropouts_fraction)))+'\n')
                # Write the fractions for all the redshifts and mangitudes
                fws.writelines(line) 
                #print redshift, identified_fraction, dropouts_fraction
        fws.close()
        # Generate plots.
        plot_completeness.main(parameters['path_to_images'],
                               parameters['min_mag'], parameters['max_mag'],
                               parameters['mag_bins'], parameters['min_z'],
                               parameters['max_z'], parameters['z_bins'],
                               cat[i], total_completeness, total_dropouts)


if __name__ == "__main__":
    main()
