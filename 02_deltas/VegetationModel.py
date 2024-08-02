import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mpl_toolkits.axes_grid1 as axtk

import os
import shutil
import sys
import yaml
from packaging import version

import pyDeltaRCM
from pyDeltaRCM.shared_tools import sec_in_day, day_in_yr

from scipy import ndimage


class VegetationModel(pyDeltaRCM.DeltaModel):
    """Implementation of Lauzon's DeltaRCM Vegetation.

    This model implements the "DeltaRCM Vegetation" model from Lauzon and
    Murray 2018 [1]_ on top of the *pyDeltaRCM* model.

    This implementation was realized from working off the narrative
    description of the model in the paper, and the source code associated
    with the paper (archived on Zenodo at doi: 10.5281/zenodo.1434243).

    .. note:: 

        This model can not exactly reproduce the runs in [1]_, because of many
        small discrepencies between the underlying DeltaModel and the model
        as implemented by Lauzon et al.

    In implementation the VegetationModel is a subclass of the pyDeltaRCM
    `DeltaModel` class, so that all methods and attributes defined there are
    inherited by this model. The vegetation model behavior from [1]_ is added
    or modified by using "hooks" and overwriting.

    .. important::

        This implementation is provided without warranty or guarantee or
        correctness.

    For more information on the *pyDeltaRCM* model, and how to further modify
    the VegetationModel, please see the project documentation
    (https://deltarcm.org/pyDeltaRCM/index.html).

    .. [1] Lauzon, R., & Murray, A. B. (2018). Comparing the cohesive effects 
           of mud and vegetation on delta evolution. Geophysical Research
           Letters, 45, 10437–10445. https://doi.org/10.1029/2018GL079405
    """

    def __init__(self, input_file, **kwargs):

        # requires pyDeltaRCM version where `mod_water_weight` is 
        #   implemented as a way to modify water routing.
        assert version.parse(pyDeltaRCM.__version__) >= version.parse("2.1.4")

        # inherit from base model
        super().__init__(input_file, **kwargs)
        self.hook_after_create_domain()

    def hook_import_files(self):
        """Define the custom YAML parameters."""
        # whether to run vegetation
        self.subclass_parameters['vegetation'] = {
            'type': 'bool', 'default': False
        }

        # whether to save vegetation figures
        self.subclass_parameters['save_veg_frac_figs'] = {
            'type': 'bool', 'default': False
        }

        # stem rooting depth
        self.subclass_parameters['p_veg_d_root'] = {
            'type': ['int', 'float'], 'default': 0.20  # 20 cm
        }

        # stem diameter
        self.subclass_parameters['p_veg_d_stem'] = {
            'type': ['int', 'float'], 'default': 0.06  # 6 mm
        }

        # stem density
        self.subclass_parameters['p_veg_K'] = {
            'type': ['int', 'float'], 'default': 800  # 800 stems/m2
        }

        # growth rate
        self.subclass_parameters['p_veg_r'] = {
            'type': ['int', 'float'], 'default': 1  # 1 (to be divided by seconds in year)
        }

        # flood duration, how often to build vegetation
        self.subclass_parameters['p_veg_est_flood_dur'] = {
            'type': ['int', 'float'], 'default': 3  # 3 days
        }

        # interflood duration, how long to build vegetation for
        self.subclass_parameters['p_veg_est_inter_dur'] = {
            'type': ['int', 'float'], 'default': 100  # 100 days
        }

        # min depth to establish
        self.subclass_parameters['p_veg_est_depth'] = {
            'type': ['int', 'float'], 'default': 0.5  # 0.5 m
        }

        # min elev change to establish
        self.subclass_parameters['p_veg_est_roc'] = {
            'type': ['int', 'float'], 'default': 0.01  # \Delta\eta < 0.01 d_root
        }

        # init veg established
        self.subclass_parameters['p_veg_est_init'] = {
            'type': ['int', 'float'], 'default': 0.05  # 0.05 frac
        }

    def hook_init_output_file(self):
        """Add non-standard grids, figures and metadata to be saved."""
        if self.save_veg_frac_figs:
            self._save_fig_list['veg_frac'] = ['veg_frac']

        # or add other variables of interest
        # self._save_fig_list['mod_water_weight'] = ['mod_water_weight']
        # self._save_fig_list['veg_alpha'] = ['veg_alpha']

        # save the active layer grid each save_dt w/ a short name
        self._save_var_list['veg_frac'] = ['veg_frac', 'fraction',
                                           'f4', ('time',
                                                  'x', 'y')]

    def hook_after_create_domain(self):
        """Add fields to the model for all vegetation parameterizations.
        """
        self.veg_frac = np.zeros_like(self.depth)
        self.veg_alpha = np.ones_like(self.depth)
        self.veg_d_root = self.p_veg_d_root
        self.veg_d_stem = self.p_veg_d_stem
        self.veg_K = self.p_veg_K
        self.veg_r = self.p_veg_r / (sec_in_day * day_in_yr)

        self.eta_change = np.zeros_like(self.depth)

        self.veg_est_flood_duration = self.p_veg_est_flood_dur * sec_in_day  # duration of flooding
        self.veg_est_interflood_duration = self.p_veg_est_inter_dur * sec_in_day  # time for veg growth
        self.veg_est_depth = self.p_veg_est_depth
        self.veg_est_roc = self.p_veg_est_roc
        self.veg_est_init = self.p_veg_est_init
        self.time_since_interflood = 0  # counter for growth

        # threshold for vegetation to affect flow weighting
        self.veg_b = (0.7 / self.veg_d_stem / self.veg_K)
        
        # what is this thing?
        #  "coefficient to make vegetation have proper influence"
        self.veg_A = (0.88 * 4 / 
                      (np.pi * self.veg_d_stem**2 * self.veg_K * 
                       ((4 / self.veg_d_stem / self.veg_K) - self.veg_b)))

    def hook_run_water_iteration(self):
        """Update vegetation parameters before the water routing.
        """
        # determine the new alpha value based on vegetation density
        #   this is the how the "bank stability" described in the paper
        #   is implemented. See "topo_diffusion()" below
        self.veg_alpha = -0.099 * self.veg_frac + 0.1

        # This is the part that adds weighting to water routing based on
        # vegetation. The idea is that vegetation slows down flow by
        # increasing roughness. We use the self.mod_water_weight to
        # represent this change in flow resistence.
        # First, determine what the weights would be everywhere
        _part1 = (
            self.veg_A * np.pi * self.veg_d_stem**2 * 
            self.veg_K * (self.veg_frac - self.veg_b))
        _mod_water_weight = (1 - (_part1 / 4))
        # Then, we modulate the weights wherever the veg_frac is above threshold
        #    reset everything to 1 (no weight), then update where veg above threshold
        self.mod_water_weight[:] = 1
        self.mod_water_weight[self.veg_frac >= self.veg_b] = _mod_water_weight[self.veg_frac >= self.veg_b]
        self.mod_water_weight = np.clip(self.mod_water_weight, 0, 1)

    def hook_after_route_sediment(self):
        """Apply vegetation growth/death rules.
        """
        # determine change in bed elevation on this timestep
        self.eta_change = self.eta - self.eta0

        # if vegetation is on, run the growth/death routines
        if self.vegetation:

            self.time_since_interflood += self._dt
            
            # mortality happens every timestep
            self._vegetation_mortality()

            # growth only happens on interflood
            if self.time_since_interflood >= self.veg_est_flood_duration:
                
                # grow the vegetation
                self._vegetation_growth()
                # reset the counter to start a new flood period
                self.time_since_interflood = 0

        # cannot have vegetation fraction outside 0,1
        self.veg_frac = np.clip(self.veg_frac, 0, 1)

    def _vegetation_mortality(self):
        """Vegeation mortality method.

        This method implements vegetation mortality every timestep.
        """
        # kill all veg that is in locations with bed *change* above root depth
        self.veg_frac[np.abs(self.eta_change) >= self.veg_d_root] = 0
        
        # kill all veg anywhere the depth is greater than 1 m
        self.veg_frac[self.depth > 1] = 0

        # determine the new possible veg frac everywhere, based on elevation change
        _possible_veg_mortality = (
            self.veg_frac * (1 - (np.abs(self.eta_change) / self.veg_d_root)))
        # determine where bed change could reduce vegetation
        _where_eta_change = np.logical_and(
            np.abs(self.eta_change) > 0,
            np.abs(self.eta_change) < self.veg_d_root)
        
        # apply updated vegetation fraction
        self.veg_frac[_where_eta_change] = _possible_veg_mortality[_where_eta_change] 
        
    def _vegetation_growth(self):
        """Vegetation growth method.

        This method implements periodic vegetation growth. Here, we include an
        optional flag to indicate whether the implementation should follow
        from the Lauzon code or Lauzon paper description. By specifying
        `_veg_growth_flag` as either `code` or `paper`, it is possible to
        control the behavior. The main difference is whether the bed
        elevation-change threshold can be overridden by a dry cell. In the
        published code, this is the implementation, whereas the paper makes
        no mention of this.
        """
        _veg_growth_flag = 'code'

        # find where the elevation change is less than the threshold
        _where_nochange = (self.eta_change < (self.veg_est_roc*self.veg_d_root))
        # find where the depth of the cell is in the marsh window
        _where_depth = (self.depth < self.veg_est_depth)
        _where_depth_and_nochange = np.logical_and(_where_depth, _where_nochange)
        # find where there is already no vegetation
        _where_noveg = (self.veg_frac == 0)
        
        if _veg_growth_flag == 'code':
            # This is my best interpretation of the Lauzon code.
            #
            #   It indicates that any dry cell will be able to get 
            #   vegetation, regardless of elevation change. This differs
            #   from the description given in the paper, which indicats the cell
            #   wetness does not affect whether elevation-change threshold matters.

            # find where dry or meeting depth and change criteria
            _where_dry = (self.depth < self.dry_depth)
            _where_dry_or_danc = np.logical_or(_where_dry, _where_depth_and_nochange)
            
            # valid is {[(where dry) or (where (in wet window) and (no change))] and (no veg)}
            _where_valid = np.logical_and(_where_dry_or_danc, _where_noveg)

        if _veg_growth_flag == 'paper':
            # This is my best interpretation of the Lauzon paper description.
            #
            #   The paper indicates that the only thing that matters is the
            #   depth being in the window. There is no mention of the dry
            #   cells always being able to establish vegetation.

            # valid is {[where (in wet window) and (no change)] and (no veg)}
            _where_valid = np.logical_and(_where_depth_and_nochange, _where_noveg)

        # where it is valid to establish veg, and there is no veg already there
        #  at these locations, establish with the initial vegetation parameter
        self.veg_frac[_where_valid] = self.veg_est_init

        # calculate the change in vegeation everywhere already vegetation
        #   note: this includes the newly established vegetation
        #   note: vegetation grows for length of an interflood period
        dveg_frac = (
            self.veg_est_interflood_duration *
            (1-(self.veg_frac)) * self.veg_r * self.veg_frac)
        self.veg_frac = self.veg_frac + dveg_frac
        
        # update sea level rise during interflood time
        # note: it is critical to NOT scale sea level rise to intermittency
        #   (as we usually do with pyDeltaRCM), because the "off" time is
        #   accounted for in this model subclass via this step here
        self.H_SL += self.SLR * self.veg_est_interflood_duration

    def topo_diffusion(self):
        """Overwrite with new behavior.
    
        This method takes care of the "bank stability" change to the model.
        Here, we simply multiply the `qs_lat` and `S` calculated from kernels
        with the `veg_alpha` coefficient field to determine topographic
        diffusion.
        """
        for _ in range(self.N_crossdiff):

            a = ndimage.convolve(self.eta, self.kernel1, mode='constant')
            b = ndimage.convolve(self.qs, self.kernel2, mode='constant')
            c = ndimage.convolve(self.qs * self.eta, self.kernel2,
                                 mode='constant')
            d = ndimage.convolve(self.veg_alpha, self.kernel2, mode='constant')

            self.cf = (d * self.diffusion_multiplier *
                       (self.qs * a - self.eta * b + c))

            self.cf[self.cell_type == -2] = 0
            self.cf[0, :] = 0

            self.eta += self.cf


if __name__ == '__main__':

    # parameter choices for scaling
    If = 10 / 365.25  # intermittency factor for year-scaling

    # base yaml configuration
    base_output = './vegetation_output'
    base_yaml = './vegetation.yaml'

    _mdl = VegetationModel(
        input_file=base_yaml,
        out_dir=base_output,
        save_checkpoint=True,
        save_dt=86400,
        clobber_netcdf=True,
        vegetation=True)

    # solve for how many timesteps
    targ_dur = 30  # target run duration (years)
    targ_dur_mdl = (targ_dur*sec_in_day*day_in_yr) * If
    tsteps = int((targ_dur_mdl // _mdl.time_step) + 1)

    try:
        for i in range(tsteps):
            _mdl.update()
    except Exception as e:
        print(e)
        _mdl.logger.exception(e)

    # finalize
    _mdl.finalize()

   
