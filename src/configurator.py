from bunch import Bunch
import copy
import json
import os
import shutil

class Configurator:
  """ Stores parameter options and creates set of parameter values to be used in simulations. 

  :param parameters_path: Path to a file containing parameter options under different scenarios, e.g. transmission settings, interventions.

  """
  def __init__(self, parameters_path):
    self._parameters_path = parameters_path
    with open(parameters_path) as par_file:
      self._parameters = json.load(par_file, object_hook=lambda d: Bunch(d))

  def copy_parameters_file(self, output_dir):
    """ Copy the parameter file used to initialize this Configurator instance.

    :param output_dir: Put the copy of the parameter file in this directory.
    :return: The path of the copied file.
    
    """
    dest_path = os.path.join(output_dir, os.path.basename(self._parameters_path))
    shutil.copyfile(self._parameters_path, dest_path)
    return dest_path

  def get_parameters(self, path=""):
    """ Get the parameter options stored in this Configurator instance. """
    parameters = copy.deepcopy(self._parameters)
    if not path:
      return parameters
    else:
      for c in path.split('.'):
        parameters = getattr(parameters, c)
      return parameters

  def create_configuration(self, **kwargs):
    r"""Create a configuration object (set of parameter values) based on a scenario.

    :param \**kwargs: See below.

    :Keyword arguments:
      * *application*          -- Application of the intervention. 'itn' for insecticide-treated nets (exposure on day of first feed) or 'irs' (exposure on day of first indoor rest)
      * *intervention*         -- Type of intervention (characterized by combination of effects). 'dbh' for DBH, 'ins' for insecticide, or 'none' for no intervention.
      * *coverage*             -- Proportion of mosquitoes exposed on their first feed (or indoor rest). Value between 0 and 1.
      * *dose*                 -- Dose of DBH or insecticide, which determines the efficacies of the different effects. 'low', 'medium', or 'high'
      * *num_days*             -- Run the simulation for this many days.
      * *intervention_day*     -- Introduce the intervention on this day. Must be less than num_days.
      * *transmission*         -- Transmission setting, determined by the biting rate. 'low', 'medium', or 'high'
      * *larva_mortality*      -- Shape of the relationship between larva mortality and larva population size. 'linear' or 'logistic'
      * *human_infectiousness* -- Shape of the relationship between human infectiousness to mosquitoes and human prevalence. 'linear' or 'nonlinear'
    """
    kwargs = Bunch(kwargs)
    for k in ('intervention', 'coverage', 'dose', 'num_days', 
              'intervention_day', 'transmission', 'application', 
              'larva_mortality', 'human_infectiousness'):
      if k not in kwargs:
        raise ValueError('Need to specify "{}".'.format(k))

    # use parameters object as a scaffold for our configuration object
    config = copy.deepcopy(self._parameters)
    del config.efficacy # we will be putting efficacies under "intervention"

    # add some parameters
    if kwargs.larva_mortality not in {'linear', 'logistic'}:
      raise ValueError("Arguments for 'larva_mortality' are 'linear' or 'logistic'.")
    if kwargs.human_infectiousness not in {'linear', 'nonlinear'}:
      raise ValueError("Arguments for 'human_infectiousness' are 'linear' or 'nonlinear'.")
    config.simulation = Bunch(
      num_days              = kwargs.num_days, 
      intervention_day      = kwargs.intervention_day,
      larva_mortality       = kwargs.larva_mortality,
      human_infectiousness  = kwargs.human_infectiousness
    )
    config.intervention = Bunch(
      coverage    = kwargs.coverage, 
      application = kwargs.application,
      efficacy    = Bunch(
        mating               = 0,
        egg                  = 0,
        parasite             = 0,
        persistent_mortality = [],
        immediate_mortality  = 0
      )
    )

    # select the transmission level
    try:
      config.malaria.biting_scaling_factor = self._parameters.malaria.biting_scaling_factor[kwargs.transmission]
    except KeyError:
      raise ValueError("Arguments for 'transmission' are 'low', 'medium', or 'high'.")

    # select index corresponding to the dose level
    try: 
      index = {"low": 0, "medium": 1, "high": 2}[kwargs.dose]
    except KeyError:
      raise ValueError("Arguments for 'dose' are 'low', 'medium', or 'high'.")

    # then use dose index and intervention type to set combination of efficacies
    itv = kwargs.intervention
    if itv == 'dbh':
      config.intervention.efficacy.mating               = self._parameters.efficacy.mating[index]
      config.intervention.efficacy.egg                  = self._parameters.efficacy.egg[index]
      config.intervention.efficacy.parasite             = self._parameters.efficacy.parasite[index]
      config.intervention.efficacy.persistent_mortality = self._parameters.efficacy.persistent_mortality[index]
    elif itv == 'ins':
      config.intervention.efficacy.immediate_mortality  = self._parameters.efficacy.immediate_mortality[index]
    elif itv == 'dbh_mating':
      config.intervention.efficacy.mating               = self._parameters.efficacy.mating[index] 
    elif itv == 'dbh_egg':
      config.intervention.efficacy.egg                  = self._parameters.efficacy.egg[index]
    elif itv == 'dbh_parasite':
      config.intervention.efficacy.parasite             = self._parameters.efficacy.parasite[index]
    elif itv == 'dbh_mortality':
      config.intervention.efficacy.persistent_mortality = self._parameters.efficacy.persistent_mortality[index]
    elif itv == 'none':
      pass
    else:
      raise ValueError("Arguments for 'intervention' are 'dbh', 'ins', 'dbh_mating', dbh_egg', 'dbh_parasite', dbh_mortality', or 'none'.")
  
    return config
