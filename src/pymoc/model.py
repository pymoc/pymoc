import numpy as np
from pymoc.modules import ModuleWrapper


class Model(object):
  r"""
  Model

  Instances of this class maintain a collection of physics modules,
  their relationships to one another, and their evolution in time.
  Creation of a Model instance allows users to integrate arbitrary 
  basin configurations without needing to manually track, timestep,
  and reconcile the individual modules.
  """
  def __init__(self):
    self.basins = []
    self.couplers = []
    self._modules = {}

  def __getitem__(self, key):
    return self._modules[key]

  @property
  def keys(self):
    r"""
    Get a list of all module keys for the model configuration.

    Returns
    -------

    keys : list
           A list containing key strings for each module.
    """
    return self._modules.keys()

  def get_module(self, key):
    r"""
    Retrieve a model module by its key.

    Parameters
    ----------

    key : string
          Key of the module being retrieved

    Returns
    -------

    module : ModuleWrapper
             ModuleWrapper instance associated with the specified key
    """
    if not key:
      return None
    return self._modules[key]

  def validate_neighbors_input(self, neighbors):
    r"""
    Ensure that a set of neighbors is valid, specifically that no neighbors are being
    linked multiple times at once. Raises a ValueError if validation fails.

    Parameters
    ----------

    neighbors : list
                A list of ModuleWrappers pointing at the modules to be validated
    """

    neighbor_names = [n.module for n in neighbors]
    distinct_neighbor_names = np.unique([n.name for n in neighbors])
    if len(neighbor_names) > len(distinct_neighbor_names):
      raise ValueError(
          'Cannot link basins multiple times. Please check your configuration.'
      )

  def validate_new_module_key(self, module):
    r"""
    Ensure that a new module is valid, specifically that no module with the same
    key already exists in the model.

    Parameters
    ----------

    module: ModuleWrapper
            A ModuleWrappers pointing at the module to be validated
    """
    if hasattr(self, module.key):
      raise NameError(
          'Cannot use module name ' + module.name +
          ' because it would overwrite an existing key or model property.'
      )

  def add_module(self, module, name, left_neighbors=[], right_neighbors=[]):
    r"""
    Add a physical module to the model configuration, with its location determined
    by the specified neighbors.

    Parameters
    ----------

    module : Module instance
             The physics module (e.g. Column, PsiThermwind, etc.) to be added to the
             model configuration
    name : string
           The name with which to label the module in the model configuration
    left_neighbors : list
                     List of existing model modules that will be to the "left" of the
                     newly added module
    right_neighbors : list
                      List of existing model modules that will be to the "right" of the
                      newly added module
    """

    self.validate_neighbors_input(left_neighbors + right_neighbors)
    module_wrapper = ModuleWrapper(
        module,
        name,
        left_neighbors=left_neighbors,
        right_neighbors=right_neighbors
    )
    self.validate_new_module_key(module_wrapper)

    self._modules[module_wrapper.key] = module_wrapper
    if module_wrapper.module_type == 'basin':
      self.basins.append(module_wrapper)
    elif module_wrapper.module_type == 'coupler':
      self.couplers.append(module_wrapper)

    setattr(self, module_wrapper.key, module_wrapper.module)

  def new_module(
      self,
      module_class,
      module_args,
      module_name,
      left_neighbors=[],
      right_neighbors=[]
  ):
    r"""
    Create a new physical module and add itto the model configuration, with its
    location determined by the specified neighbors.

    Parameters
    ----------

    module_class : Class
                   The class constructor of the physics module to be created and added
                   to the model configuration
    moduel_args : dict
                  A dictionary of argument keys and values to be used to initialize
                  the newly created module
    module_name : string
                  The name with which to label the module in the model configuration
    left_neighbors : list
                     List of existing model modules that will be to the "left" of the
                     newly added module
    right_neighbors : list
                      List of existing model modules that will be to the "right" of the
                      newly added module
    """
    self.add_module(
        module_class(**module_args),
        module_name,
        left_neighbors=left_neighbors,
        right_neighbors=right_neighbors
    )

  def run(
      self,
      steps,
      basin_dt,
      coupler_dt=1,
  ):
    r"""
    Integrate the model forward in time. This method will timestep all basin modules at the
    basin timestep, and update all coupler modules at the coupling timestep. This method will
    optionally yield, allowing diagnosis of the model's temporal evolution.

    Parameters
    ----------

    steps: int
           The number of model integration timesteps 
    basin_dt : int
               The timestep length (in seconds) for basin modules
    coupler_dt : int, optional
                 The duration (in seconds) between coupler module updates. If unspecified,
                 coupling takes place at every timestep
    snapsot_start : int, optional
                    The model step at which to begin yielding for model snapshotting
    snapshot_interval : int, optional
                        The number of steps between yields for model snapshotting
    """
    for coupler in self.couplers:
      coupler.update_coupler()

    for i in range(0, steps):
      self.timestep(i, basin_dt, coupler_dt=coupler_dt)

  def run_with_snapshot(
      self,
      steps,
      basin_dt,
      coupler_dt=1,
      snapshot_start=None,
      snapshot_interval=None
  ):
    r"""
    Integrate the model forward in time. This method will timestep all basin modules at the
    basin timestep, and update all coupler modules at the coupling timestep. This method will
    optionally yield, allowing diagnosis of the model's temporal evolution.

    Parameters
    ----------

    steps: int
           The number of model integration timesteps 
    basin_dt : int
               The timestep length (in seconds) for basin modules
    coupler_dt : int, optional
                 The duration (in seconds) between coupler module updates. If unspecified,
                 coupling takes place at every timestep
    snapsot_start : int, optional
                    The model step at which to begin yielding for model snapshotting
    snapshot_interval : int, optional
                        The number of steps between yields for model snapshotting
    """
    for coupler in self.couplers:
      coupler.update_coupler()
    
    for i in range(0, steps):
      self.timestep(i, basin_dt, coupler_dt=coupler_dt)
      if self.snapshot(i, snapshot_start, snapshot_interval):
        yield i

  def timestep(self, step, basin_dt, coupler_dt=0):
    r"""
    Execute a single model integration step, which timesteps all basins and updates
    all couplers at the coupling timestep

    Parameters
    ----------

    step : int
           The current model timestep
    basin_dt : int
               The timestep length (in seconds) for basin modules
    coupler_dt : int, optional
                 The duration (in seconds) between coupler module updates. If unspecified,
                 coupling takes place at every timestep
    """

    for basin in self.basins:
      basin.timestep_basin(dt=basin_dt)
    if step % coupler_dt == 0:
      for coupler in self.couplers:
        coupler.update_coupler()

  def snapshot(self, step, snapshot_start, snapshot_interval):
    r"""
    Determine whether the current model timestep should yield to snapshotting.

    Parameters
    ----------

    step : int
           The current model timestep
    snapsot_start : int 
                    The model step at which to begin yielding for model snapshotting
    snapshot_interval : int
                        The number of steps between yields for model snapshotting

    Returns
    -------

    snapshot : bool
               Whether all criteria are met for the model to yield to snapshotting on
               the current timestep
    """

    return snapshot_start is not None and snapshot_interval is not None and step >= snapshot_start and (
        step-snapshot_start
    ) % snapshot_interval == 0
