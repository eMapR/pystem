"""
Evaluates a prediction using the AUC and MAE statistics. The only command line argument is a text file with the
following parameters specified:

in_raster - full path of the prediction raster to evaluate
sample_txt - tab-delimited text file of test data
val_col - column name in sample_txt containing truth values
model_dir - the directory containing

Optional parameters:
inventory_txt - tab-delimited text file with the following columns: stamp, OOB_score, auc, mae, mae_true,
                mae_false, n_samples, max_depth, resolution_m, predictors. This text file is an inventory
                of accuracy values and distinguishing parameters used to train iterations of the model.

"""

import os
import sys
import gdal
import numpy as np
import pandas as pd

import randomforest as forest

#def main(in_raster, sample_txt, val_col, model_dir, nodata, inventory_txt=None):
def main(params):

    # Read params and make variables from text
    inputs = forest.read_params(params)
    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])

    # Check that variables were specified in params
    try:
        str_check = in_raster, sample_txt, val_col
        nodata = int(nodata)
    except NameError as e:
        print ''
        missing_var = str(e).split("'")[1]
        msg = "Variable '%s' not specified in param file:\n%s" % (missing_var, params)
        raise NameError(msg)

    ds = gdal.Open(in_raster)
    ar_p = ds.ReadAsArray()
    ds = None
    model_dir = os.path.dirname(in_raster)
    stamp = os.path.basename(model_dir)
    #out_dir = os.path.join(model_dir, stamp)

    test_samples = pd.read_csv(sample_txt, sep='\t', index_col='obs_id')
    mae, mae_true, mae_false = forest.calc_mae(ar_p, test_samples, val_col, nodata)
    auc = forest.calc_auc(ar_p, test_samples, val_col, nodata, model_dir)

    print ''
    print 'AUC ................... ', auc
    print 'MAE ................... ', mae
    print 'MAE of true values .... ', mae_true
    print 'MAE of false values ... ', mae_false

    predict_params = forest.read_params(os.path.join(model_dir, 'predict_rf_params.txt'))
    zone_path = predict_params['in_raster'].replace('"','')
    ds = gdal.Open(zone_path)
    ar_z = ds.ReadAsArray()
    nodata_z = ds.GetRasterBand(1).GetNoDataValue()
    ds = None
    nodata_p = int(predict_params['nodata'].replace('"',''))
    forest.zonal_distribution(ar_p, ar_z, nodata_p, nodata_z, out_dir=model_dir)
    #p_rate = forest.prediction_rate(ar_p, ar_t, 1, nodata, out_dir)
    #print p_rate

    if 'inventory_txt' in locals():
        print 'Saving assessment values to inventory_txt: ', inventory_txt
        df = pd.read_csv(inventory_txt, sep='\t', index_col='stamp')
        df.ix[stamp, ['auc', 'rmse', 'rmse_true', 'rmse_false']] = auc, mae, mae_true, mae_false
        df.to_csv(inventory_txt, sep='\t')

    print ''
    print 'AUC ................... ', auc
    print 'MAE ................... ', mae
    print 'MAE of true values .... ', mae_true
    print 'MAE of false values ... ', mae_false

    return auc, mae, mae_true, mae_false


"""in_raster = r'T:\ResMgmt\WAGS\Geology\Geohazards\Susceptibility\Models\Random_Forest\susceptibility_20160905_1543\final_susceptibility_20160905_1543.tif'
sample_txt= 'T:\\ResMgmt\\WAGS\\Geology\\Geohazards\\Susceptibility\\Samples\\samples60000_20160905_1538_5m\\samples_test_60000_20160905_1538_5m.txt'
val_col = 'is_slide'
model_dir = r'T:\ResMgmt\WAGS\Geology\Geohazards\Susceptibility\Models\Random_Forest\susceptibility_20160905_1543' """
#inventory_txt = r'T:\ResMgmt\WAGS\Geology\Geohazards\Susceptibility\Models\Random_Forest\model_inventory.txt'
#main(in_raster, sample_txt, val_col, model_dir, -9999)#'''
#in_raster, sample_txt, val_col, model_dir, nodata=-9999, inventory_txt=None
if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))#"""



