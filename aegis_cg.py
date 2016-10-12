"""Computes CG for AEGIS catalog galaxies and saves to file"""
import galsim
import meas_cg_fns as mcg
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def main():
    filters=['f606w', 'f814w']
    file_filter_name = ['V', 'I']
    main ='/Users/choyma/AEGIS/'
    main_path = main + 'AEGIS_full/AEGIS_training_sample/'
    comp_cat_name = main + 'AEGIS_full/complete_AEGIS_galaxy_catalog_filter_25.2.fits'
    cat_name = 'AEGIS_galaxy_catalog_filter_25.2.fits'
    colors=['b','r','g']
    cc, good={},{}
    rgc, comp_cat={}, {}
    for f , filt in enumerate(filters):
        name = main_path + cat_name.replace('filter', file_filter_name[f])
        rgc[filt] = galsim.RealGalaxyCatalog(name, dir=main_path)
        cc[filt] = galsim.COSMOSCatalog(file_name = name, dir=main_path, exclusion_level="bad_stamp")
        good[filt]=[cc[filt].getOrigIndex(idx) for idx in range(cc[filt].getNObjects())]
        comp_cat[filt] = Table.read(comp_cat_name.replace('filter', file_filter_name[f]))
    cond1 = comp_cat[filters[0]]['XMAX_IMAGE'] - comp_cat[filters[0]]['XMIN_IMAGE']< 300
    cond2 = comp_cat[filters[0]]['YMAX_IMAGE'] - comp_cat[filters[0]]['YMIN_IMAGE']< 300
    cond3 = comp_cat[filters[1]]['XMAX_IMAGE'] - comp_cat[filters[1]]['XMIN_IMAGE']< 300
    cond4 = comp_cat[filters[1]]['YMAX_IMAGE'] - comp_cat[filters[1]]['YMIN_IMAGE']< 300
    cond5 = comp_cat[filters[0]]['FLUX_RADIUS']> 13
    cond6 = comp_cat[filters[1]]['FLUX_RADIUS']>13
    cond7 = comp_cat[filters[0]]['min_mask_dist_pixels']> 11
    cond8 = comp_cat[filters[1]]['min_mask_dist_pixels']> 11
    cond9 = comp_cat[filters[0]]['average_mask_adjacent_pixel_count']/comp_cat[filters[0]]['peak_image_pixel_count']< 0.4
    cond10 = comp_cat[filters[1]]['average_mask_adjacent_pixel_count']/comp_cat[filters[1]]['peak_image_pixel_count']< 0.4
    q,= np.where(cond5 &cond6 & cond7 &cond8 &cond9 &cond10)#np.where(cond1 & cond2 & cond3 & cond4 & cond5 & cond6)
    cg_pick = np.intersect1d(q, good[filters[0]], good[filters[1]])
    rng = galsim.BaseDeviate(123456)
    noise_file = main + 'AEGIS_full/AEGIS_training_sample/acs_filter_unrot_sci_cf.fits' 
    cg_arr, er_arr =[], []
    for indx in cg_pick:
        print indx
        gal1 = galsim.RealGalaxy(rgc['f606w'],index=indx)
        gal2 = galsim.RealGalaxy(rgc['f814w'],index=indx)
        nx, ny = gal1.gal_image.array.shape[0], gal1.gal_image.array.shape[1]
        noise_file1 = noise_file.replace('filter', 'V')
        noise_file2 = noise_file.replace('filter', 'I')
        noise_std1 = np.sqrt(comp_cat['f606w']['NOISE_VARIANCE'][indx])
        noise_std2 = np.sqrt(comp_cat['f814w']['NOISE_VARIANCE'][indx])
        noise = galsim.getCOSMOSNoise(file_name = noise_file1)
        noise1 = galsim.Image(ny, nx)
        noise1.addNoise(noise)
        noise_im1= noise1.array
        noise = galsim.getCOSMOSNoise(file_name = noise_file2)
        noise2 = galsim.Image(ny, nx)    
        noise2.addNoise(noise)
        noise_im2= noise2.array
        noise_im1 = noise_im1*noise_std1/np.std(noise_im1)
        noise_im2 = noise_im2*noise_std2/np.std(noise_im2)
        cg,er, t_cg = mcg.comp_cg_error(gal1.gal_image,gal2.gal_image,
                      noise_im1, noise_im2)
        cg_arr.append(cg)
        er_arr.append(er)
    np.savetxt('AEGIS_cg_val.txt', [cg_arr,er_arr])



if __name__ == '__main__':
    main()