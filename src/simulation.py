from bunch import Bunch
from collections import defaultdict, namedtuple
import copy
import json
import math
from configurator import Configurator
import lifecycle

#################################
# Demographic-related functions #
#################################
def get_max_age(model):
  all_compartments = model.select_compartments(lambda c: True)
  oldest = max(all_compartments, key=lambda c: c.get_attr('age'))
  return oldest.get_attr('age')

def get_larva_mortality_func(config):
  B = 1.0 - config.lifecycle.daily_survival.larva.max # baseline mortality
  T = 1.0 - config.lifecycle.daily_survival.larva.min # maximum mortality
  K = config.population.larval_capacity               # so-called "capacity"

  larva_mortality = config.simulation.larva_mortality
  if larva_mortality == 'logistic':
    shape = config.population.larval_mortality_shape.logistic
    C1, C2, C3 = shape.c1, shape.c2, shape.c3       # constants that determine shape of logistic curve
    logistic = lambda x: B - C1 + (T - (B - C1)) / (1.0 + math.exp(-C2 * (x / float(K)) + C3))
    return logistic
  elif larva_mortality == 'linear':
    C = config.population.larval_mortality_shape.linear.c
    linear = lambda x: min(T, B * (1.0 + C * float(x) / K))
    return linear
  else:
    raise ValueError('Larva mortality functional form not recognized.')

def get_exposed(model):
  num_exposed = 0
  for c in model.select_compartments(lambda c: c.get_attr('stage') == 'adult'):
    num_exposed += c.get_size() * c.get_attr('fraction_exposed')
  num_adults = model.get_combined_size(lambda c: c.get_attr('stage') == 'adult')
  return num_exposed, num_adults

def _set_fraction_exposed(model, intervention):
  # set fraction exposed according to ITN
  if intervention.application == 'itn':
    for c in model.select_compartments(lambda c: c.get_attr('stage') == 'adult'):  
      name = c.get_name()
      if name in {'emergence-rest', 'mate-first'}:
        frac = 0.0
      elif name in {'feed-first', 'feed-after'} or name.startswith('M'):
        frac  = intervention.coverage
      elif name == 'mate-after' or name.startswith('F'):
        cov = intervention.coverage
        eff = intervention.efficacy.mating
        frac = cov * (1 - eff) / (1 - cov * eff)
      elif name.startswith('V'):
        frac = 1.0 # all the virgin mosquitoes have been exposed
      else:
        raise RuntimeError('Unspecified coverage calculation for c {}'.format(c))
      c.set_attr(fraction_exposed=frac)
  
  # set fraction exposed according to IRS. note that frac exposed is 0 for the virgin track
  elif intervention.application == 'irs':
    for c in model.select_compartments(lambda c: c.get_attr('stage') == 'adult'):
      name = c.get_name()
      pre_exposure = {'emergence-rest', 'mate-first', 'feed-first', 'feed-after', 'mate-after'}
      if name in pre_exposure or name.startswith('V'):
        frac = 0.0
      elif name.startswith('M') or name.startswith('F'):
        frac = intervention.coverage
      else:
        raise RuntimeError('Unspecified coverage calculation for c {}'.format(c))
      c.set_attr(fraction_exposed=frac)
  
  # there are only two options, 'itn' or 'irs'.
  else:
    raise ValueError("intervention.application is {}. Expected 'itn' or 'irs'".format(intervention.application))

def _update_model(model, intervention):
  # set fraction exposed for each compartment
  _set_fraction_exposed(model, intervention)

  # choose the relevant way to count days post exposure
  if intervention.application == 'itn':
    days_post_exposure_attr = 'days_post_first_feed'
  elif intervention.application == 'irs':
    days_post_exposure_attr = 'days_post_second_rest'

  # persistent mortality effect
  for i, eff in enumerate(intervention.efficacy.persistent_mortality):
    # select the compartments affected    
    if (i == len(intervention.efficacy.persistent_mortality) - 1):
      cmpts = model.select_compartments(lambda c: c.get_attr(days_post_exposure_attr) >= i)
    else:
      cmpts = model.select_compartments(lambda c: c.get_attr(days_post_exposure_attr) == i)
    
    # adjust their survival
    for c in cmpts:
      reachable = model.get_reachable_compartments(c.get_name())
      # terminal compartments
      if len(reachable) == 1 and list(reachable)[0].get_name() == 'dead':
        model.update_transition(c.get_name(), 'dead', 1.0)
      # non-terminal compartments
      else:
        cov = c.get_attr('fraction_exposed')
        b = model.config.daily_survival.adult.baseline
        s_new = b * (1 - eff)
        s_old = 1.0 - model.get_transition(c.get_name(), 'dead')
        s = cov * s_new + (1 - cov) * s_old
        model.update_transition(c.get_name(), 'dead', 1.0 - s)

  # immediate mortality effect
  affected = model.select_compartments(lambda c: c.get_attr(days_post_exposure_attr) == 0)

  for c in affected:
    m = model.get_transition(c.get_name(), 'dead')
    m_new = m + (1.0 - m) * c.get_attr('fraction_exposed') * intervention.efficacy.immediate_mortality
    model.update_transition(c.get_name(), 'dead', m_new)

  # mating effect (only applies for ITN)
  if intervention.application == 'itn':
    ff_survival = 1.0 - model.get_transition('feed-first', 'dead')
    ff_c = model.get_compartment('feed-first')
    model.set_transition(
      'feed-first', 'V0-rest1', ff_survival * ff_c.get_attr('fraction_exposed') * intervention.efficacy.mating
    )
    model.set_transition(
      'feed-first', 'mate-after', ff_survival * (1.0 - ff_c.get_attr('fraction_exposed') * intervention.efficacy.mating)
    )

  # egg effect
  mated_layers = model.select_compartments(lambda c: c.get_attr('activity') == 'laying' and c.get_attr('mated'))
  for c in mated_layers:
    new_egg_batch = (1.0 - c.get_attr('fraction_exposed') * intervention.efficacy.egg) * c.get_attr('egg_batch')
    c.set_attr(egg_batch=new_egg_batch)

  # make sure all outgoing transitions sum to 1 at least
  model.validate_transitions()

def _update_larva_mortality(model):
  num_larva = model.get_combined_size(lambda c: c.get_attr('stage') == 'larva')
  larva_mortality = model.get_attr('larva_mortality_func')(num_larva)
  for c in model.select_compartments(lambda c: c.get_attr('stage') == 'larva'):
    if c.get_attr('day') < model.config.development_days.larva - 1:
      model.update_transition(c.get_name(), 'dead', larva_mortality)
    else:
      f = model.config.female_fraction
      model.update_transition(c.get_name(), 'dead', larva_mortality * f + (1.0 - f))

def _count_potentially_infectious(model, parasite_efficacy, sporogony_days):
  total = 0.0
  for c in model.select_compartments(lambda c: True):
    frac_exposed = c.get_attr('fraction_exposed') or 0
    if c.get_attr('days_post_first_feed') >= sporogony_days:
      total += c.get_size() * (1 -  frac_exposed * parasite_efficacy)
  return total

def run_demographic_simulation(config, selectors={'all': lambda c: True}):

  # create model
  model = lifecycle.create_model(config.lifecycle)
  model.set_attr(
    larva_mortality_func=get_larva_mortality_func(config),
    max_age=get_max_age(model)
  )
  
  # prepopulate reachable compartments
  reachable = model.get_reachable_compartments('egg-0')
  initial_total_size = config.population.initial_size
  for c in reachable:
    c.set_size(initial_total_size / len(reachable))
  
  # set initial egg batches
  mated_layers = model.select_compartments(lambda c: c.get_attr('activity') == 'laying' and c.get_attr('mated'))
  for c in mated_layers:
    c.set_attr(egg_batch=config.population.egg_batch)
  
  # make a copy for posterity
  initial_model = copy.deepcopy(model)
  
  # run simulation and record population sizes
  timeseries = defaultdict(list)
  for day in range(config.simulation.num_days):
    # update the model to reflect intervention
    if day == config.simulation.intervention_day:
      _update_model(model, config.intervention)
    
    # update larva mortality based on current larva population
    _update_larva_mortality(model)

    # before advancing, compute number of new eggs and then replenish egg compartment afterwards
    num_eggs = sum(c.get_size() * c.get_attr('egg_batch') for c in mated_layers)
    model.advance()
    model.get_compartment('egg-0').set_size(num_eggs)

    # count number of mosquitoes in subsets of compartments we're interested in
    for key, selector in selectors.iteritems():
      timeseries[key].append(model.get_combined_size(selector))

    num_pinf = _count_potentially_infectious(
      model, 
      config.intervention.efficacy.parasite,
      config.malaria.sporogony_days
    )
    timeseries['potentially_infectious'].append(num_pinf)

  return Bunch(config=copy.deepcopy(config), start_model=initial_model, end_model=model, timeseries=timeseries)


#############################
# Malaria-related functions #
#############################
def get_exposure_history(name, model, sporogony_days):
  exposure_history = []
  compartment = model.get_compartment(name)
  get_days_ago = lambda c: compartment.get_attr('age') - c.get_attr('age')

  for c in model.select_compartments(lambda c: c.get_attr('activity') == 'feeding'):
    reachable_c = model.get_reachable_compartments(c.get_name())
    days_ago = get_days_ago(c)
    if compartment in reachable_c and days_ago >= sporogony_days:
      exposure_history.append(days_ago)

  return exposure_history

def _get_num_infectious_mosquitoes(feeders, exposure_history, risk_history, oocyst_probs):
  num_infectious = 0
  for name, num in feeders.iteritems():
    pr_unexposed = 1.0
    for day in exposure_history[name]:
      pr_unexposed *= (1 - risk_history[day - 1] * oocyst_probs[name])
    num_infectious += (1 - pr_unexposed) * num
  return num_infectious

def run_malaria_simulation(config, feeder_timeseries, exposure_history, exposed_fractions):
  # a susceptible mosquito's daily risk of infection (function of human prevalence)
  # note: susceptible mosquito's must be feeding
  M = config.malaria.max_human_infectiousness
  if config.simulation.human_infectiousness == 'linear':
    C = config.malaria.human_infectiousness_shape.piecewise_linear.c
    beta_m = lambda p: min(M * p / C, M)
  elif config.simulation.human_infectiousness == 'nonlinear':
    C1 = config.malaria.human_infectiousness_shape.nonlinear.c
    C2 = 0.5 * (M + math.sqrt(M**2 + 4 * M * C1))
    beta_m = lambda p: max(0.0, C2 - C1 / (p + (C1 / C2)))

  # a susceptible human's daily risk of infection (function of number of infectious mosquitoes)
  B = config.malaria.inoculation_efficiency # probability of infection given bite
  K = config.malaria.biting_scaling_factor  # bites per mosquito per human over one day
  beta_h = lambda n: 1.0 - ((1.0 - B) ** (K * n))

  # a susceptible daily human's probability of recovery
  V = 1.0 / config.malaria.recovery_days

  # prepopulate our malaria timeseries with a made-up history
  timeseries  = defaultdict(list)
  days_to_remember = max((max(v) if v else 0) for v in exposure_history.itervalues())
  for i in range(days_to_remember):
    timeseries['I_h'].append(config.malaria.initial_human_prevalence) # prevalence of infectious humans

  # run malaria simulation
  num_days = len(feeder_timeseries.itervalues().next())

  # development of oocyst in mosquito, given infective bite
  oocyst_probs = { k : 1.0 for k in feeder_timeseries.iterkeys() }

  for i in range(num_days):
    if i == config.simulation.intervention_day:
      for k in exposed_fractions.iterkeys():
        oocyst_probs[k] *= (1.0 - exposed_fractions[k] * config.intervention.efficacy.parasite)

    feeders = { k: v[i] for k, v in feeder_timeseries.iteritems() }
    risk_history = [beta_m(x) for x in timeseries['I_h'][-days_to_remember:]]

    # calculate mosquito prevalence each time step
    nI_m = _get_num_infectious_mosquitoes(feeders, exposure_history, risk_history, oocyst_probs)
    nF_m = sum(feeders.itervalues())
    I_m = nI_m / float(nF_m) if nF_m > 1e-8 else 0 

    # for humans, add those newly infected and subtract those recoverying
    I_h_prev = timeseries['I_h'][-1]
    I_h = I_h_prev + beta_h(nI_m) * (1.0 - I_h_prev) - V * I_h_prev
    timeseries['I_h'].append(I_h) 
    timeseries['I_m'].append(I_m)
    timeseries['nI_m'].append(nI_m)
    timeseries['nF_m'].append(nF_m)

  # lop off the made-up history
  timeseries['I_h'] = timeseries['I_h'][days_to_remember:]


  return Bunch(config=copy.deepcopy(config), beta_m=beta_m, beta_h=beta_h, 
    exposure_history=exposure_history, timeseries=timeseries)


########################################
# Remove clutter from iPython notebook #
########################################
def _get_demographic_selectors(model, config):
  # start with selectors for the various lifestages
  selectors = {
    'egg'   : lambda c: c.get_attr('stage') == 'egg',
    'larva' : lambda c: c.get_attr('stage') == 'larva',
    'adult' : lambda c: c.get_attr('stage') == 'adult',
    'older' : lambda c: c.get_attr('days_post_first_feed') >= config.malaria.sporogony_days
    #'older' : lambda c: c.get_attr('age') >= 14 # first feed occurs at age 2, sporogony then takes 12 days
  }
  
  # add selectors for each one-day age group
  max_age = get_max_age(model)
  for a in range(max_age + 1):
    selectors['age-{}'.format(a)] = lambda c, a=a: c.get_attr('age') == a # "a=a" captures local value of a

  # add selectors for each feeding compartment
  feeding_compartments = model.select_compartments(lambda c: c.get_attr('activity') == 'feeding')
  for cmpt in feeding_compartments:
    selectors[cmpt.get_name()] = lambda c, name=cmpt.get_name(): c.get_name() == name

  return selectors

def _get_feeder_names(model):
  feeders = model.select_compartments(lambda c: c.get_attr('activity') == 'feeding')
  return [c.get_name() for c in feeders]

def run_simulation(parameters_file, larval_capacity=None, **options):
  # setup
  configurator = Configurator(parameters_file)
  config_lo = configurator.create_configuration(transmission='low', **options)
  config_me = configurator.create_configuration(transmission='medium', **options)
  config_hi = configurator.create_configuration(transmission='high', **options)

  # adjust larval capacity if one is specified (convenient for sensitivity analysis)
  if larval_capacity:
    config_lo.population.larval_capacity = larval_capacity
    config_me.population.larval_capacity = larval_capacity
    config_hi.population.larval_capacity = larval_capacity

  model = lifecycle.create_model(config_lo.lifecycle)

  # run a demographic simulation
  selectors = _get_demographic_selectors(model, config_lo) 
  dem_results = run_demographic_simulation(config_lo, selectors) # for demographics, does not matter which config we use

  # run malaria simulations (low and high transmission) based on the demographic results
  feeder_names = _get_feeder_names(model)
  feeder_timeseries = { x: dem_results.timeseries[x] for x in feeder_names }
  exposure_history = { x: get_exposure_history(x, model, config_lo.malaria.sporogony_days) for x in feeder_names }
  exposed_fractions = { x: dem_results.end_model.get_compartment(x).get_attr('fraction_exposed') for x in feeder_names }

  mal_results_lo = run_malaria_simulation(config_lo, feeder_timeseries, exposure_history, exposed_fractions)
  mal_results_me = run_malaria_simulation(config_me, feeder_timeseries, exposure_history, exposed_fractions)
  mal_results_hi = run_malaria_simulation(config_hi, feeder_timeseries, exposure_history, exposed_fractions)

  # return results (+ additional information?)
  return Bunch(
    population = dem_results, 
    malaria    = Bunch(low=mal_results_lo, medium=mal_results_me, high=mal_results_hi)
  )

def run_simulation_from_config(config):
  model = lifecycle.create_model(config.lifecycle)

  # run a demographic simulation
  selectors = _get_demographic_selectors(model) 
  dem_results = run_demographic_simulation(config, selectors) # for demographics, does not matter which config we use

  # run a malaria simulation based on the demographic results
  feeder_names = _get_feeder_names(model)
  feeder_timeseries = { x: dem_results.timeseries[x] for x in feeder_names }
  exposure_history = { x: get_exposure_history(x, model, config.malaria.sporogony_days) for x in feeder_names }
  exposed_fractions = { x: dem_results.end_model.get_compartment(x).get_attr('fraction_exposed') for x in feeder_names }

  mal_results = run_malaria_simulation(config, feeder_timeseries, exposure_history, exposed_fractions)

  # return results
  return Bunch(
    population = dem_results, 
    malaria    = mal_results
  )