import numpy as np
from itertools import product

from liana.testing._sample_anndata import generate_toy_mdata
from liana.method.sp._bivar import bivar

def test_bivar_morans():
    mdata = generate_toy_mdata()    
    
    # default
    bivar(mdata, x_mod='adata_x', y_mod='adata_y', function_name='morans')
    assert 'local_scores' in mdata.mod.keys()
    np.testing.assert_almost_equal(mdata.mod['local_scores'].X.sum(), -9.947859, decimal=5)
    
    # with perms
    bivar(mdata, x_mod='adata_x', y_mod='adata_y', 
          function_name='morans', n_perms=2)
    
    np.testing.assert_almost_equal(np.mean(mdata.mod['local_scores'].layers['pvals']), 0.604936507, decimal=6)


def test_bivar_nondefault():
    mdata = generate_toy_mdata()   
    # test different params, inplace = False
    proximity = np.ones((mdata.shape[0], mdata.shape[0]))
    # proximity = np.zeros((mdata.shape[0], mdata.shape[0]))
    # np.fill_diagonal(proximity, 1)
    mdata.obsp['ones'] = proximity
    
    global_stats, local_scores, local_pvals, local_categories = \
          bivar(mdata, x_mod='adata_x', y_mod='adata_y', 
                function_name='morans', n_perms=0,
                connectivity_key='ones', remove_self_interactions=False,
                x_layer = "scaled", y_layer = "scaled", inplace=False, 
                add_categories=True
                )
    
    # if all are the same weights, then everything is close to 0?
    np.testing.assert_almost_equal(global_stats['global_r'].sum(), 0)
    local_scores.shape == (700, 100)
    local_pvals.shape == (700, 100)
    np.testing.assert_almost_equal(np.min(np.min(local_pvals)), 0.5, decimal=2)
    
    assert local_categories.sum() == -22400
    
    

def test_bivar_external():
    mdata = generate_toy_mdata()
    ones = np.ones((mdata.shape[0], mdata.shape[0]), dtype=np.float64)
    
    x_vars = mdata.mod['adata_x'].var.index[-3:]
    y_vars = mdata.mod['adata_y'].var.index[0:3]
    interactions = list(product(x_vars, y_vars))
    
    bivar(mdata, x_mod='adata_x', y_mod='adata_y', function_name='morans', connectivity=ones, interactions=interactions)
    
    
def test_masked_pearson():
    mdata = generate_toy_mdata()
    bivar(mdata, x_mod='adata_x', y_mod='adata_y', function_name='masked_pearson')
    
    # check local
    assert 'local_scores' in mdata.mod.keys()
    mdata.mod['local_scores'].X[np.isnan(mdata.mod['local_scores'].X)]=0 # TODO fix this
    np.testing.assert_almost_equal(mdata.mod['local_scores'].X.mean(), 0.010296326, decimal=5)
    
    # check global
    assert 'global_res' in mdata.uns.keys()
    assert set(['global_mean','global_sd']).issubset(mdata.uns['global_res'].columns)
    # check specific values are what we expect; TODO
    


def test_vectorized_pearson():
    mdata = generate_toy_mdata()
    bivar(mdata, x_mod='adata_x', y_mod='adata_y', function_name='pearson')
    
    # check local
    assert 'local_scores' in mdata.mod.keys()
    np.testing.assert_almost_equal(mdata.mod['local_scores'].X.mean(), 0.009908355, decimal=5)
    
    # check global
    assert 'global_res' in mdata.uns.keys()
    assert set(['global_mean','global_sd']).issubset(mdata.uns['global_res'].columns)
