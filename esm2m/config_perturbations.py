import xarray as xr
import os

def get_path(variable=None,ppname=None,override=False,experiments=None,timespan=None):
    '''Returns a dictionary of paths relevant to the specified 
    experiments, variables, and timespans.'''
    
    config = 'MOM5_SIS_BLING_CORE2-gat'
    if override:
        config = config+'-override-po4'
    pp = '/archive/Richard.Slater/Siena/siena_201308_rds-c3-gat-slurm/'+config+'/gfdl.ncrc3-intel16-prod-openmp/pp/'
    localdir = '/ts/monthly/10yr/'
    
    # Unless experiment is specified, set to all
    if experiments is None:
        experiments = ['','_gat','_zero','_double']
    if type(experiments)==str:
        experiments = [experiments]
    # Unless variable is specified, set to all
    if variable is None:
        variable = '*'
    # Unless timespan is specified, set to all
    if timespan is None:
        timespan = '*'
    
    # If ppname is not specified, derive from variable
    if ppname is None:
        d = get_variable_dict()                
        if variable in d:
            ppname = d[variable]
        else:
            raise NameError('If ppname is not specified, must give exact variable name.'+
                           ' To specify wildcard variable, specify ppname.')
    
    # Configure correct ppname
    if ppname == 'ocean_bling_tracers':
        ppname_pre = 'ocean_bling'
        ppname_suf = '_tracers'
    elif ppname in ['ocean_bling_ocn_flux','bling_atm_flux','ocean_gat_dic']:
        ppname_pre = ppname
        ppname_suf = ''
        
    # ocean_gat_dic has no non-gat control
    if (ppname == 'ocean_gat_dic') & ('' in experiments):
        experiments.remove('')
        
    paths = {}
    for e in experiments:
        ppname = ppname_pre+e+ppname_suf
        filename = ppname+'.'+timespan+'.'+variable+'.nc'
        path = pp+ppname+localdir+filename
        paths[e] = path
        
    return paths

def load_exps(variable=None,ppname=None,override=False,experiments=None,timespan=None,verbose=False):
    '''Assigns data to a dictionary of xarray Datasets corresponding to 
    each experiment'''
    
    paths = get_path(variable=variable,ppname=ppname,override=override,experiments=experiments,timespan=timespan)
        
    dd = {}
    for p,path in paths.items():
        if verbose:
            print(path)
        dd[p] = xr.open_mfdataset(path)
    if len(dd)==1:
        dd=dd[experiments]
    return dd

def dmget_exps(variable=None,ppname=None,override=False,experiments=None,timespan=None,verbose=False,wait=True):
    '''Issue dmget for associated netcdf files.'''
    
    paths = get_path(variable=variable,ppname=ppname,override=override,experiments=experiments,timespan=timespan)
    
    for p,path in paths.items():
        command = 'dmget '+path+' &'
        if wait:
            command = command[:-2]
        if verbose:
            print(command)
        os.system(command)

def load_grid(fromwork=True,z=None,z_i=None):
    if fromwork:
        # Load augmented grid saved to work
        # See notebook save_grid.ipynp
        gridpath = '/work/gam/projects/bio-pump-timescales/esm2m/data/grid.nc'
        grid = xr.open_dataset(gridpath)
    else:
        pp = '/archive/Richard.Slater/Siena/siena_201308_rds-c3-gat-slurm/MOM5_SIS_BLING_CORE2-gat/gfdl.ncrc3-intel16-prod-openmp/pp/'
        gridpath = pp+'static.nc'
        grid = xr.open_dataset(gridpath)
        if z is not None:
            grid['dz'] = (z_i.diff('st_edges_ocean')
                          .rename({'st_edges_ocean':'st_ocean'})
                          .assign_coords({'st_ocean':z}))
            grid['volume_t'] = grid['area_t']*grid['dz']
            
    return grid

def calc_anom(dd):
    ddanom = {}
    ddanom['zero'] = dd['_zero']-dd['_gat']
    ddanom['double'] = dd['_double']-dd['_gat']
    ddanom['noneq'] = dd['']-dd['_gat']
    return ddanom

def get_variable_dict():
    return {'alk':'ocean_bling_tracers',
          'alpha':'ocean_bling_tracers',
          'biomass_p':'ocean_bling_tracers',
          'chl':'ocean_bling_tracers',
          'co2_alpha':'ocean_bling_tracers',
          'co3_ion':'ocean_bling_tracers',
          'delta_csurf':'ocean_bling_tracers',
          'delta_pco2':'ocean_bling_tracers',
          'dic_area_integral':'ocean_bling_tracers',
          'dic':'ocean_bling_tracers',
          'dic_stf':'ocean_bling_tracers',
          'dic_volume_integral':'ocean_bling_tracers',
          'dop_area_integral':'ocean_bling_tracers',
          'dop':'ocean_bling_tracers',
          'dop_volume_integral':'ocean_bling_tracers',
          'fed':'ocean_bling_tracers',
          'fed_stf':'ocean_bling_tracers',
          'htotal':'ocean_bling_tracers',
          'integral_dic':'ocean_bling_tracers',
          'integral_dic_stf':'ocean_bling_tracers',
          'irr_mem':'ocean_bling_tracers',
          'jdic_100':'ocean_bling_tracers',
          'o2':'ocean_bling_tracers',
          'pco2_surf':'ocean_bling_tracers',
          'po4_area_integral':'ocean_bling_tracers',
          'po4':'ocean_bling_tracers',
          'po4_volume_integral':'ocean_bling_tracers',
          'co2_flux_alpha_ocn':'ocean_bling_ocn_flux',
          'co2_flux_cair_ice_ocn':'ocean_bling_ocn_flux',
          'co2_flux_csurf_ocn':'ocean_bling_ocn_flux',
          'co2_flux_flux_ice_ocn':'ocean_bling_ocn_flux',
          'co2_flux_kw_ice_ocn':'ocean_bling_ocn_flux',
          'co2_flux_schmidt_ocn':'ocean_bling_ocn_flux',
          'o2_flux_alpha_ocn':'ocean_bling_ocn_flux',
          'o2_flux_csurf_ocn':'ocean_bling_ocn_flux',
          'o2_flux_flux_ice_ocn':'ocean_bling_ocn_flux',
          'o2_flux_schmidt_ocn':'ocean_bling_ocn_flux',
          'atm_gas_flux':'ocean_gat_dic',
          'atm_gas_input':'ocean_gat_dic',
          'atm_mol_wgt':'ocean_gat_dic',
          'base_mix_ratio':'ocean_gat_dic',
          'gas_mol_wgt':'ocean_gat_dic',
          'mix_ratio':'ocean_gat_dic',
          'total_atm_mass':'ocean_gat_dic',
          'total_gas_mass':'ocean_gat_dic'}

# Atmospheric tracer flux grid behaving badly
#           'co2_flux':'bling_atm_flux',
#           'co2_flux_pcair_atm':'bling_atm_flux',
#           'o2_flux':'bling_atm_flux'}

def disp_variables():
    return list(get_variable_dict().keys())

def add_override_suffix(directory,override):
    if override:
        directory = directory+'override-po4/'
    else:
        directory = directory+'no-override/'
    return directory


    
    