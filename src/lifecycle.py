import copy
import math
from compartmodel import Model

def create_model(config):
  """ Create a compartmental model of the mosquito lifecycle. """
  
  model = Model()
  
  ### Add compartments ###
  # developmental compartments and first adult compartment
  eggs = []
  for i in range(config.development_days.egg):
    eggs.append(model.add_compartment('egg-{}'.format(i), stage='egg', day=i))
  
  larvae = []
  for i in range(config.development_days.larva):
    larvae.append(model.add_compartment('larva-{}'.format(i), stage='larva', day=i))

  model.add_compartment('emergence-rest', stage='adult', age=0)

  # mate-first track
  mate_first_track = []
  mate_first_track.append(model.add_compartment('mate-first', stage='adult', activity='mating'))
  mate_first_track.append(model.add_compartment('feed-after', stage='adult', activity='feeding'))
  for i in range(config.maximum_cycles):
    mate_first_track.append(model.add_compartment('M{}-rest1'.format(i), stage='adult', cycle=i, activity='resting'))
    mate_first_track.append(model.add_compartment('M{}-rest2'.format(i), stage='adult', cycle=i, activity='resting'))
    mate_first_track.append(model.add_compartment('M{}-lay'.format(i),   stage='adult', cycle=i, activity='laying', mated=True))
    mate_first_track.append(model.add_compartment('M{}-feed'.format(i),  stage='adult', cycle=i, activity='feeding'))

  # feed-first track
  feed_first_track = []
  feed_first_track.append(model.add_compartment('feed-first', stage='adult', activity='feeding'))
  feed_first_track.append(model.add_compartment('mate-after', stage='adult', activity='mating'))
  for i in range(config.maximum_cycles):
    feed_first_track.append(model.add_compartment('F{}-rest1'.format(i), stage='adult', cycle=i, activity='resting'))
    feed_first_track.append(model.add_compartment('F{}-rest2'.format(i), stage='adult', cycle=i, activity='resting'))
    feed_first_track.append(model.add_compartment('F{}-lay'.format(i),   stage='adult', cycle=i, activity='laying', mated=True))
    feed_first_track.append(model.add_compartment('F{}-feed'.format(i),  stage='adult', cycle=i, activity='feeding'))

  # forever-virgin subtrack
  virgin_track = []
  for i in range(config.maximum_cycles):
    virgin_track.append(model.add_compartment('V{}-rest1'.format(i), stage='adult', cycle=i, activity='resting'))
    virgin_track.append(model.add_compartment('V{}-rest2'.format(i), stage='adult', cycle=i, activity='resting'))
    virgin_track.append(model.add_compartment('V{}-lay'.format(i),   stage='adult', cycle=i, activity='laying', mated=False))
    virgin_track.append(model.add_compartment('V{}-feed'.format(i),  stage='adult', cycle=i, activity='feeding'))

  # add the age and "days since first exposure" attributes for the three tracks
  for i, name in enumerate(mate_first_track):
    model.get_compartment(name).set_attr(age=(i + 1), days_post_first_feed=(i - 1), days_post_second_rest=(i - 2))
  for i, name in enumerate(feed_first_track): 
    model.get_compartment(name).set_attr(age=(i + 1), days_post_first_feed=(i + 0), days_post_second_rest=(i - 2))
  for i, name in enumerate(virgin_track):
    model.get_compartment(name).set_attr(age=(i + 2), days_post_first_feed=(i + 1), days_post_second_rest=(i + 0))

  # where all mosquitoes end up
  model.add_compartment('dead')

  ### Add transitions ###
  # developmental stages through to post emergence rest
  model.chain_compartments(eggs, config.daily_survival.egg)
  model.set_transition(eggs[-1], larvae[0], config.daily_survival.egg)
  model.chain_compartments(larvae, config.daily_survival.larva.max)
  model.set_transition(larvae[-1], 'emergence-rest', config.female_fraction * config.daily_survival.larva.max)

  # helper function for getting the age-specific adult survival
  def get_survival(name):
    age = model.get_compartment(name).get_attr('age')
    discounts = config.daily_survival.adult.age_discounts
    discount = discounts[age] if age < len(discounts) - 1 else discounts[-1]
    baseline = config.daily_survival.adult.baseline
    return baseline * discount

  # splitting into tracks
  model.set_transition('emergence-rest', 'mate-first', get_survival('emergence-rest') * config.mate_first_fraction)
  model.set_transition('emergence-rest', 'feed-first', get_survival('emergence-rest') * (1 - config.mate_first_fraction))

  
  # and the split from the feed-first track into the forever-virgin subtrack
  model.set_transition('feed-first', virgin_track[0], 0)  
  
  # helper function for linking up the compartments of each track
  def link_track(track):
    for c1, c2 in zip(track[:-1], track[1:]):
       model.set_transition(c1, c2, get_survival(c1))

  # link up the three tracks
  link_track(mate_first_track)
  link_track(feed_first_track)
  link_track(virgin_track)

  # dead
  model.set_sink('dead')

  # attach a copy of config to the model for future reference 
  model.config = copy.deepcopy(config)

  # double-check before returning
  model.validate_transitions()
  return model
