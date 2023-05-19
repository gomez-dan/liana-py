import os
import pathlib
import numpy as np
import scanpy as sc
from mudata import MuData
from liana.method.sp._misty_constructs import lrMistyData, genericMistyData
from liana.testing._sample_anndata import generate_toy_spatial

test_path = pathlib.Path(__file__).parent

adata = sc.read_h5ad(os.path.join(test_path, "data" , "synthetic.h5ad"))
adata = sc.pp.subsample(adata, n_obs=100, copy=True)
mdata = MuData({'rna':adata})


def test_misty_para():

    misty = genericMistyData(adata, bandwidth=10,
                             cutoff=0, add_juxta=False,
                             set_diag=False, seed=133)
    misty(bypass_intra=False, seed=42)
    assert np.isin(list(misty.uns.keys()), ['target_metrics', 'importances']).all()
    target_metrics = misty.uns['target_metrics']
    # NOTE: contributions are not exactly equal to 1 per target per view
    # likely due to numpy rounding
    assert np.sum(target_metrics[['intra', 'para']].values, axis=1).sum()  == 11.0
    assert target_metrics.shape == (11, 8)
    
    importances = misty.uns['importances']
    assert importances.shape == (220, 6)
    assert importances[importances['target']=='ECM']['value'].sum() == 2.0
    interaction_msk = (importances['target']=='ligA') & \
        (importances['predictor']=='protE')
    np.testing.assert_almost_equal(importances[interaction_msk]['value'].values,
                                np.array([0.0011129, 0.0553538]))
    assert target_metrics['gain_R2'].mean() == -0.0032406852423374943
    

def test_misty_bypass():    
    misty = genericMistyData(adata, bandwidth=10, add_juxta=True, set_diag=True,
                             cutoff=0, coord_type="generic", delaunay=True)
    misty(alphas=1, bypass_intra=True, seed=42)
    assert np.isin(['juxta', 'para'], misty.uns['target_metrics'].columns).all()
    assert ~np.isin(['intra'], misty.uns['target_metrics'].columns).all()
    assert misty.uns['target_metrics'].shape == (11, 8)
    assert misty.uns['target_metrics']['multi_R2'].sum() == -2.142582410377362
    
    importances = misty.uns['importances']
    assert importances.shape == (220, 6)
    assert importances['value'].sum().round(10) == 22.0
    np.testing.assert_almost_equal(importances[(importances['target']=='ligC') &
                                               (importances['predictor']=='ligA')]['value'].values,
                                   np.array([0.07792869, 0.055088]))
    

def test_misty_groups():        
    misty = genericMistyData(adata,
                             bandwidth=20,
                             add_juxta=True, 
                             set_diag=False,
                             cutoff=0,
                             coord_type="generic",
                             delaunay=True
                             )
    misty(alphas=1, 
          bypass_intra=False,
          seed=42,
          predict_self=True, 
          extra_groupby='cell_type', 
          intra_groupby='cell_type')
    
    assert misty.uns['target_metrics'].shape==(44, 9)
    perf_actual = (misty.uns['target_metrics'].
     groupby(['intra_group', 'env_group'])['gain_R2'].
     mean().values
    )
    perf_expected = np.array([-0.10600711322549207, 0.04402447752647749, 0.1635208524520852, 0.04838660952488709])
    np.testing.assert_almost_equal(perf_actual, perf_expected)
    
    # assert that there are self interactions = var_n * var_n
    importances = misty.uns['importances']
    self_interactions = importances[(importances['target']==importances['predictor'])]
    # 11 vars * 4 envs * 3 views = 132
    assert self_interactions.shape == (132, 6)
    assert self_interactions[self_interactions['view']=='intra']['value'].isna().all()


def test_lr_misty():
    adata = generate_toy_spatial()
    misty = lrMistyData(adata, bandwidth=10, set_diag=True, cutoff=0)
    assert misty.shape == (700, 42)
    
    misty(n_estimators=10, bypass_intra=True)
    assert misty.uns['target_metrics'].shape == (16, 7)
    
    importances = misty.uns['importances']
    assert importances.shape == (415, 6)
    cmplxs = importances[importances['target'].str.contains('_')]['target'].unique()
    assert np.isin(['CD8A_CD8B', 'CD74_CXCR4'], cmplxs).all()
