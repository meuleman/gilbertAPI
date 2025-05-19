# =======================================================================================
# YC240401 Shared factor anaylsis for UMAP
# =======================================================================================

import os
import numpy as np
import pandas as pd
import h5py
import umap
from flask import Blueprint, request
import json
from constants import GilbertSFCH5_DIR, COMP_MODEL_DIR

# Create a Blueprint
umap_bp = Blueprint('umap_bp', __name__)

# Map chromosome location to index at specific Gilbert order
def map_to_order_indice(chr_str, order, df_indices):
    # Parse chr_str
    chr,start_end = chr_str.split(':')
    start,end = [int(i) for i in start_end.split('-')]
    # Get index 
    segment_size = 4**(14-order)
    chr_index = start // segment_size
    chr_offset = df_indices[str(order)][chr] 
    index = chr_index + chr_offset
    return index

# Use fancy indexing on h5
def fancy_indexing(hf, index):
    # Create unique ordered version of index for hf
    index_sorted = np.sort(np.unique(index)) # index for hf
    hfidx2idx = dict(zip(index_sorted,np.arange(len(index_sorted))))
    # Get array from hf then reindex
    arr_hf = np.array(hf[index_sorted],dtype='float32')
    arr = arr_hf[[hfidx2idx[i] for i in index]]
    return arr

# Given dataset and order, create path to h5 file
def get_h5_path(ds,order):
    basedir = GilbertSFCH5_DIR
    path = f"{basedir}/{ds}/offset_0/order_{order}_windows_annotation_sfc_scores.h5"
    return path

# Given variable metadata for factors and kwargs for columns in the metadata, return boolean index
def subset_by_var(var, **kwargs):
    var_idx = np.repeat(True, len(var))
    for k,v in kwargs.items():
        var_idx = np.all([var_idx,(var[k] == v)],axis=0)
    return var_idx

# Given SFC array of (sample,feature), return top n arg per sample (sample, n)
def argnmax(arr, n=1, axis=1, ordered=True):
    if ordered:
        idx = np.argsort(-arr, axis=axis)[:,:n]
    else:
        idx = np.argpartition(arr, -n, axis=axis)[:,-n:]
    return idx

# Given SFC array of (sample,feature), return top args per sample and index of features considered
def get_top_factor(arr, var, n=1, ordered=True, exclude_LAD=False, **kwargs):
    """
    """
    var_idx = subset_by_var(var, **kwargs)
    if exclude_LAD:
        var_idx[subset_by_var(var,dataset='LAD')] = False
    arr_keep = arr[:,var_idx]
    topargs = argnmax(arr_keep,n=10,ordered=ordered)
    zero_idx = np.array([arr_keep[i,topargs[i]] for i in range(len(topargs))]) == 0
    topargs[zero_idx] = -1
    return topargs, var_idx


def shared_factor_percentage(top_factors, n_shared=10):
    """
    per factor, show percentage of the samples the factor was top factor in
    """
    top_shared_factor_counts = pd.Series(top_factors.reshape(-1)).value_counts()
    top_shared_factor_counts_per_sample = top_shared_factor_counts / len(top_factors)
    return top_shared_factor_counts_per_sample


def perform_shared_factor_analysis(arr, var, n=10, min_percentage=0.2):
    """
    return
        (mean shared factor percentage, 
         list of shared factor names in decreasing order,
         individual shared factor perecentage)
    """
    top_factors,var_idx = get_top_factor(arr,n=n,var=var)
    shared_factor = shared_factor_percentage(top_factors)
    shared_factor = shared_factor[shared_factor.index != -1]
    score = shared_factor[:n].sum() / n
    shared_factor_name = var['scaled factor'][var_idx].values[shared_factor.index]
    shared_factor.index = shared_factor_name
    shared_factor = shared_factor[shared_factor.values > min_percentage]
    gil = [score,shared_factor_name,shared_factor]
    return gil


class ScaledGilbertDataLoader():
    """
    class for loading cross scale Gilbert SFC data from h5 files
    """
#
    def __init__(self):
#
        # Preload h5 files
        self.datasets = ['DHS','chromatin_states','TF']
        self.orders = np.arange(4,13)
#
        h5_dict = {}
        for ds in self.datasets:
            for order in self.orders:
                h5path = get_h5_path(ds,order)
                h5_dict[f"{ds}-{order}"] = h5py.File(h5path,'r')
#
        self.h5_dict = h5_dict
#
        # Load indices where each chromosome begins
        FILE_DIR = f'{GilbertSFCH5_DIR}/chromosome_indices.csv'
        self.df_indices = pd.read_csv(FILE_DIR,index_col=0)
#
        # Load var, variable metadata
        self._generate_var()
#
    def load_segment(self, chr,start,end):
        arrs = []
        for ds in self.datasets: 
            for order in self.orders:
                idx = map_to_order_indice(chr,start,end,order)
                hf = self.h5_dict[f"{ds}-{order}"]
                arr_ = hf['arr'][idx]
                arrs.append(arr_)
        arrs = np.concatenate(arrs,axis=1)
        return arrs
#
    def load_segments(self, chr_strs):
        arrs = []
        for ds in self.datasets: 
            for order in self.orders:
                idx = [map_to_order_indice(chr_str=chr_str,order=order,df_indices=self.df_indices) for chr_str in chr_strs]
                hf = self.h5_dict[f"{ds}-{order}"]
                arr_ = fancy_indexing(hf['arr'],idx)
                arrs.append(arr_)
        arrs = np.concatenate(arrs,axis=1)
        return arrs
#
    def _generate_var(self):
        var = []
        for sds,h5 in self.h5_dict.items():
            ds,order = sds.split('-')
            factors = h5['columns'][:].astype('str')
            df_var = pd.DataFrame({'Scaled Dataset': sds,
                                    'dataset': ds,
                                    'order': order, 
                                    'factor': factors,
                                    'scaled factor': pd.Series(factors) + '-' + order,
                                    })
            var.append(df_var)
        self.var = pd.concat(var)


class ComponentNarrationLoader():
    #
    def __init__(self, model_dir):
        self.model_dir = model_dir
        self.load_component(model_dir)
    #
        # Load indices where each chromosome begins
        FILE_DIR = f'{GilbertSFCH5_DIR}/chromosome_indices.csv'
        self.df_indices = pd.read_csv(FILE_DIR,index_col=0)
    #
    def load_component(self,model_dir):
        """
        Input
            model_dir: directory containing the model with whole genome components and component narration
        """
        h5dir = os.path.join(model_dir,'component.h5')
        comp_meta_dir = os.path.join(model_dir,'CompAnnotFull.csv')
        dom_dir = os.path.join(model_dir,f'top_10_component_per_segment-raw.h5')
        comp_narr_full = pd.read_csv(comp_meta_dir,index_col=0)
        self.comp_narr_full = comp_narr_full
        hf = h5py.File(h5dir,'r')
        self.comp_hf = hf['component']
        hf = h5py.File(dom_dir,'r')
        self.dom_hf = hf['component_index']
    #
    def load_segments(self, chr_strs, order=12):
        idx = [map_to_order_indice(chr_str=chr_str,order=order,df_indices=self.df_indices) for chr_str in chr_strs]
        hf = self.dom_hf
        dom = fancy_indexing(hf,idx)
        return dom
    #    
    def narrate_region(self, chr_strs, min_support=0.2, max_n_comp=5):
        dom = self.load_segments(chr_strs)
        narr = pd.Series(dom[:,:5].flatten()).value_counts()
        narr = narr.reset_index()
        narr.columns = ['CompIndex','Count']
        narr['CompNarr'] = self.comp_narr_full.FullNarration.values[narr.CompIndex.astype('int').values]
        narr['Support']  = narr.Count / len(chr_strs)
        narr = narr.iloc[:max_n_comp]
        narr = narr[narr.Support >= min_support]
        return narr

def do_umap(X, umap_params):
    reducer = umap.UMAP(**umap_params)
    X_umap = reducer.fit_transform(X)
    return X_umap

class GilbertSegmentHandler():
    """
    Convert between 
    """
    def __init__(
        self,
        chromosome_indices_dir = f'{GilbertSFCH5_DIR}/chromosome_indices.csv'
    ):
        self.df_chr_indices = pd.read_csv(chromosome_indices_dir,index_col=0)
        self.df_chr_indices.columns = self.df_chr_indices.columns.astype('int')
    #
    def get_chr_pos(self, index, order=12):
        """
        Get chromosome and position for indexed segments at specified order
        """
        chr_idx = self.df_chr_indices[order]
        df_rel = self._get_chr_relative_pos(index,order)
        segment_size = 4**(14-order)
        df_rel['start'] = segment_size * df_rel.ChrIndex
        df_rel['end']   = segment_size * (df_rel.ChrIndex + 1)
        df_rel['AbsIndex'] = index
        columns = ['chr','start','end','AbsIndex','ChrIndex']
        df_rel = df_rel[columns]
        return df_rel
    #
    def generate_whole_genome_indices(self, order=12):
        """
        Generate indices for whole genome at specified order
        """
        chr_idx = self.df_chr_indices[order]
        Chromosome = [np.repeat(chr,chr_len) for chr,chr_len in chr_idx.items()]
        Chromosome = np.concatenate(Chromosome)
    #
    def convert_order(self, index, order_in, order_out=12):
        """
        Convert index at order_in to index at order_out
        e.g.
        if order_in < order_out
            order_in = 4
            order_out = 5
            index_in  = [0 1]
            index_out = [0 1 2 3 4 5 6 7]
        else:
            order_in  = 12
            order_out = 11
            index_in  = [0 1 2 3 4 5 6 7 8 9]
            index_out = [0 0 0 0 1 1 1 1 2 2]
        
        """
        if order_in == order_out:
            df_rel_out = index
        elif order_in > order_out:
            df_rel = self._get_chr_relative_pos(index,order_in)
            ChrIndexIn  = df_rel.ChrIndex.values
            ChrIndexOut = (ChrIndexIn // 4**(order_in-order_out))
            AbsIndexOut = self.df_chr_indices[order_out].loc[df_rel.chr].values + ChrIndexOut
            df_rel['AbsIndex'] = AbsIndexOut
            df_rel['ChrIndex'] = ChrIndexOut
            df_rel_out = df_rel
        else:
            df_rel = self._get_chr_relative_pos(index,order_in)
    #
            # Get relative position within chromosome for order_in then convert to order_out
            ChrIndexIn = df_rel.ChrIndex.values
            n_repeat = 4**(order_out-order_in)
            ChrIndexOut = np.repeat(ChrIndexIn,n_repeat)*n_repeat + np.tile(np.arange(n_repeat),len(ChrIndexIn))
    #
            # Get absolute position for order_out
            chrs = df_rel.chr
            AbsIndexOut = np.repeat(self.df_chr_indices[order_out].loc[chrs].values,n_repeat) + ChrIndexOut
            df_rel_out = pd.DataFrame({'chr': np.repeat(chrs,n_repeat), 
                                        'ChrIndex': ChrIndexOut, 
                                        'AbsIndex': AbsIndexOut}
                                        )
        return df_rel_out
    #
    def _get_chr(self, index, order):
        """
        Get which chromosome for indexed segments at specified order
        """
        chr_idx = self.df_chr_indices[int(order)]
        index_chr = np.zeros(len(index),dtype='int')
        for i, (chr, start) in enumerate(chr_idx.items()):
            index_chr[index >= start] = i
        index_chr = chr_idx.index.values[index_chr]
        return index_chr
    #
    def _get_chr_relative_pos(self, index, order):
        """
        Get relative position within chromosome for indexed segments at specified order
        """
        chr_idx = self.df_chr_indices[int(order)]
        index_chr = self._get_chr(index,order)
        rel_index = index - chr_idx.loc[index_chr].values
        df_rel = pd.DataFrame({'chr': index_chr, 'ChrIndex': rel_index})
        return df_rel


class filtered_UMAP_generator():
    #
    def __init__(self, model_dir, 
                 current_order=4,
                 target_order=5):
        self.current_order=current_order
        self.target_order=target_order
        self.model_dir = model_dir
        self.load_component()
    #
        # Load indices where each chromosome begins
        FILE_DIR = f'{GilbertSFCH5_DIR}/chromosome_indices.csv'
        self.df_indices = pd.read_csv(FILE_DIR,index_col=0)
    #
    def load_component(self):
        """
        Input
            model_dir: directory containing the model with whole genome components and component narration
        """
        h5dir = os.path.join(self.model_dir,'SimulatedSegments',f'WholeGenomeSimulatedSegments-Order={self.target_order}.h5')
        comp_meta_dir = os.path.join(self.model_dir,'CompAnnotFull.csv')
        comp_narr_full = pd.read_csv(comp_meta_dir,index_col=0)
        self.comp_narr_full = comp_narr_full
        hf = h5py.File(h5dir,'r')
        self.comp_hf = hf['components']
        self.df_annot_factor = pd.read_parquet(os.path.join(self.model_dir,'SimulatedSegments',f"Order-{self.target_order}","YC20241104-annot-all-factor.parquet"))
        self.df_annot_comp   = pd.read_parquet(os.path.join(self.model_dir,'SimulatedSegments',f"Order-{self.target_order}","YC20241104-annot-all-component.parquet"))
    #
    def generate_umap(self, chr_strs):
        """
        Return
            df_umap: UMAP data frame with chromosome, start, end, and UMAP coordinates
            df_annot: gilbert annotations per segment
        """
    #
        #
        current_idx = [map_to_order_indice(chr_str=chr_str,
                                            order=self.current_order,
                                            df_indices=self.df_indices) for chr_str in chr_strs]
        df_chr_target = gsh.convert_order(np.array(current_idx), order_in=self.current_order, order_out=self.target_order)
        df_chr_target['start'] = df_chr_target.ChrIndex * 4**(14-self.target_order)
        df_chr_target['end']   = df_chr_target.start + 4**(14-self.target_order)
        # df_chr_target['Start'] = df_chr_target['start']
        # df_chr_target['End'] = df_chr_target['end']
        # df_chr_target['Chromosome'] = df_chr_target['chr']
        # df_chr_target['chromosome'] = df_chr_target['chr']
        df_chr_target.rename(columns={'chr':'chromosome'},inplace=True)
        idx = df_chr_target.AbsIndex.values
    #
        # Load components
        comp = fancy_indexing(self.comp_hf,idx)
    #
        # Perform UMAP
        X_umap = do_umap(comp,umap_params)
        df_umap = pd.DataFrame(X_umap,columns=['UMAP-1','UMAP-2'])
        df_umap = pd.concat([df_chr_target.reset_index(drop=True),df_umap],axis=1)
    #
        # Get matching annotations
        df_annot = self.df_annot_factor.iloc[idx]
        return df_umap, df_annot

chrs_test = ['chrX:12058624-12320768','chr5:109838336-110100480'] * 10

umap_params = {'n_neighbors': 20,
                'min_dist': 0.8,
                'metric': 'euclidean',
                'random_state': 234,
                'n_components': 2,
                'verbose': True,
                }

sgdl = ScaledGilbertDataLoader()
cnl  = ComponentNarrationLoader(COMP_MODEL_DIR)
gsh  = GilbertSegmentHandler()

# find the top regions for selected factors
@umap_bp.route('get_shared_factor', methods=['POST'])
def get_shared_factor():
    request_info = request.get_json()
    chr_strs = request_info['chr_strs']
    order = request_info['order']
    arr = sgdl.load_segments(chr_strs,order=order)
    results = perform_shared_factor_analysis(arr=arr,var=sgdl.var)
    json_return = json.dumps({
        "SharedFactor": list(results[2].index.values),
        "SharedFactorPercentage": list(results[2].values),
        "MeanSharedFactorPercentage": results[2].mean(),
        }
    )
    return json_return

# Generate a new umap at selected segments mapped to target order
@umap_bp.route('get_filtered_umap', methods=['POST'])
def get_filtered_umap():
#
    # Get request info
    request_info = request.get_json()
    chr_strs = request_info['chr_strs']
    current_order = int(request_info['current_order'])
    target_order  = int(request_info['target_order'])
#
    # UMAP 
    umapgen = filtered_UMAP_generator(model_dir=COMP_MODEL_DIR, 
                                        current_order=current_order, 
                                        target_order=target_order
                                        )
#
    # Format results
    df_umap, df_annot = umapgen.generate_umap(chr_strs)
    results = {col:list(v) for col,v in df_umap.items()}
    results['all'] = list(df_annot['all'])
    results['DHS'] = list(df_annot.DHS_ENR)
    results['TF'] = list(df_annot.TF_ENR)
    results['chromatin_states'] = list(df_annot.chromatin_states_OCCR)
#
    # Test
    # results = {'current_order': current_order, 'target_order': target_order,'v':1}
    json_return = json.dumps(results)
    return json_return

@umap_bp.route('narrate_regions', methods=['POST'])
def narrate_regions():
    request_info = request.get_json()
    chr_strs = request_info['chr_strs']
    df_narr = cnl.narrate_region(chr_strs)
    json_return = json.dumps({
        "ComponentNarration": list(df_narr.CompNarr.values),
        "Support": list(df_narr.Support.values),
        }
    )
    return json_return