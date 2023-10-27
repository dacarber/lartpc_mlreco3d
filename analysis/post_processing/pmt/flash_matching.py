import numpy as np
<<<<<<< HEAD
from collections import defaultdict
from analysis.post_processing import post_processing
from mlreco.utils.globals import *

@post_processing(data_capture=['index', 'opflash_cryoE', 'opflash_cryoW'], 
                 result_capture=['interactions'])

def run_flash_matching(data_dict, result_dict, 
                       fm=None,
                       opflash_keys=[],
                       cache=True):
    """
    Post processor for running flash matching using OpT0Finder.
    
    Parameters
    ----------
    fm : FlashManager
    opflash_keys : List[str]

    Returns
    -------
    Empty dict (operation is in-place)
        
    NOTE: This post-processor also modifies the list of Interactions
    in-place by adding the following attributes:
        interaction.fmatched: (bool)
            Indicator for whether the given interaction has a flash match
        interaction.fmatch_time: float
            The flash time in microseconds 
        interaction.fmatch_total_pE: float
        interaction.fmatch_id: int
    """
    print("Running flash matching...")
    opflashes = {}
    assert len(opflash_keys) > 0
    for key in opflash_keys:
        opflashes[key] = data_dict[key]
    
    interactions = result_dict['interactions']
    entry        = data_dict['index']
    
    # Check if coordinates are in cm
    rounding_error = np.abs(np.sum(
        interactions[0].points - interactions[0].points.astype(int)))
    if rounding_error < 1e-6:
        msg = f"Rounding error = {rounding_error:.6f}; " \
            + "It seems you are trying to run flash matching on points with pixel units. Aborting."
        raise AssertionError(msg)
    
    fmatches_E = fm.get_flash_matches(int(entry), 
                                      interactions,
                                      opflashes,
                                      volume=0,
                                      restrict_interactions=[], 
                                      cache=cache)
    fmatches_W = fm.get_flash_matches(int(entry), 
                                      interactions,
                                      opflashes,
                                      volume=1,
                                      restrict_interactions=[],
                                      cache=cache)

    flash_dict_E = {}
    for ia, flash, match in fmatches_E:
        flash_dict_E[ia.id] = (flash, match)
        ia.fmatched = True
        ia.flash_time = float(flash.time())
        ia.flash_total_pE = float(flash.TotalPE())
        ia.flash_id = int(flash.id())
        ia.flash_hypothesis = float(np.array(match.hypothesis, dtype=np.float64).sum())
        
    flash_dict_W = {}
    for ia, flash, match in fmatches_W:
        flash_dict_W[ia.id] = (flash, match)
        ia.fmatched = True
        ia.flash_time = float(flash.time())
        ia.flash_total_pE = float(flash.TotalPE())
        ia.flash_id = int(flash.id())
        ia.flash_hypothesis = float(np.array(match.hypothesis, dtype=np.float64).sum())

    print("Done flash matching.")
    return {}
=======
from warnings import warn

from analysis.post_processing import PostProcessor

from .barycenter import BarycenterFlashMatcher
from .likelihood import LikelihoodFlashMatcher


class FlashMatchingProcessor(PostProcessor):
    '''
    Associates TPC interactions with optical flashes.
    '''
    name = 'run_flash_matching'
    data_cap_opt = ['opflash', 'opflash_cryoE', 'opflash_cryoW']
    result_cap = ['interactions']

    def __init__(self,
                 opflash_map,
                 method = 'likelihood',
                 **kwargs):
        '''
        Initialize the flash matching algorithm

        Parameters
        ----------
        method : str, default 'likelihood'
            Flash matching method (one of 'likelihood' or 'barycenter')
        opflash_map : dict
            Maps a flash data product key in the data ditctionary to an
            optical volume in the detector
        **kwargs : dict
            Keyword arguments to pass to specific flash matching algorithms
        '''
        # If there is no map from flash data product to volume ID, throw
        self.opflash_map = opflash_map

        # Initialize the flash matching algorithm
        if method == 'barycenter':
            self.matcher = BarycenterFlashMatcher(**kwargs)
        elif method == 'likelihood':
            self.matcher = LikelihoodFlashMatcher(**kwargs, \
                    parent_path=self.parent_path)
        else:
            raise ValueError(f'Flash matching method not recognized: {method}')


    def process(self, data_dict, result_dict):
        '''
        Find [interaction, flash] pairs

        Parameters
        ----------
        data_dict : dict
            Input data dictionary
        result_dict : dict
            Chain output dictionary

        Notes
        -----
        This post-processor modifies the list of `interaction` objectss
        in-place by adding the following attributes:
            interaction.fmatched: (bool)
                Indicator for whether the given interaction has a flash match
            interaction.fmatch_time: float
                The flash time in microseconds
            interaction.fmatch_total_pE: float
            interaction.fmatch_id: int
        '''
        # Check if the TPC coordinates are in cm
        interactions = result_dict['interactions']
        if not len(interactions):
            return {}, {}

        # Make sure the interaction coordinates are expressed in cm
        self.check_units(interactions[0])

        # Clear previous flash matching information
        for ii in interactions:
            if ii.fmatched:
                ii.fmatched = False
                ii.flash_id = -1
                ii.flash_time = -np.inf
                ii.flash_total_pE = -1.0
                ii.flash_hypothesis = -1.0

        # Loop over flash keys
        for key, volume_id in self.opflash_map.items():
            # Get the list of flashes associated with that key
            opflashes = data_dict[key]

            # Get the list of interactions that share the same volume
            ints = [ii for ii in interactions if ii.volume_id == volume_id]

            # Run flash matching
            fmatches = self.matcher.get_matches(ints, opflashes)

            # Store flash information
            for ii, flash, match in fmatches:
                ii.fmatched = True
                ii.flash_id = int(flash.id())
                ii.flash_time = float(flash.time())
                ii.flash_total_pE = float(flash.TotalPE())
                if hasattr(match, 'hypothesis'):
                    ii.flash_hypothesis = float(np.array(match.hypothesis,
                        dtype=np.float64).sum())

        return {}, {}
>>>>>>> 04f1005cf697592bcdbb90a94a790d3c13ef269b
