from collections import defaultdict, deque
import graphviz
import pprint

"""This module contains two classes, Compartment and Model, which are used to represent compartmental models in epidemiology.

"""


class Compartment(object):
  """Represents a compartment in a compartmental model.
  
  :param name: The name of the compartment.
  :param initial_size: The initial size of the compartment (default is 0).
  :param attributes: An optional dictionary of attribute names and values to be associated with this compartment.

  """

  def __init__(self, name, initial_size=0, **attributes):
    self._name = name
    self._size = initial_size
    self._attributes = dict(attributes) 

  def get_name(self):
    """ 

    :return: Name of this Compartment

  
    """
    return self._name

  def get_size(self):
    """ 

    :return: Size of this Compartment

  
    """
    return self._size

  def set_size(self, size):
    """ Sets the size of this Compartment

    :param size: The new size of this Compartment (cannot be negative).

  
    """
    if size < 0:
      raise ValueError('Compartment size cannot be negative.')
    self._size = size

  def get_attr(self, attribute_name):
    """ Get an attribute associated with this compartment

    :param attribute_name: Name of the attribute.
    :return: The value of the attribute, or None if there is no attribute with that name.

    """
    return self._attributes.get(attribute_name, None)

  def set_attr(self, **attributes):
    """ Add attributes to this compartment. Existing attributes with the same names are overwritten.

    :param attributes: Dictionary of attribute names and values.

    """
    for attr, value in attributes.iteritems():
      self._attributes[attr] = value

  def __str__(self):
    """ Human-readable description of this compartment 

    :return: a formatted string containing the name, size, and attributes of this compartment.

    """
    pp = pprint.PrettyPrinter(indent=2)
    return pp.pformat({'name': self._name, 'attributes': self._attributes, 'size': self._size})
  
class Model(object):
  """Represents a compartmental model using a set of Compartments and transitions between Compartments. A transition is characterized by a source compartment, a destination compartment, and a transition probability. A Model instance can optionally have 

  A model is initialized with no Compartments and no transitions. 

  """
  
  def __init__(self):
    self._compartments = {}
    self._transitions = {}
    self._sink = None
    self._attributes = {}

  def get_attr(self, attr):
    """ Get an attribute associated with this model.

    :param attr: Name of the attribute.
    :return: The value of the attribute, or None if there is no attribute with that name.

    """
    return self._attributes.get(attr, None)

  def set_attr(self, **attributes):
    """ Add attributes to this Model. Existing attributes with the same names are overwritten.

    :param attributes: Dictionary of attribute names and values.

    """
    for attr, value in attributes.iteritems():
      self._attributes[attr] = value

  def add_compartment(self, name, initial_size=0, **attributes):
    """ Add a compartment to this Model.

    :param name: Name of the new compartment. Cannot be the name of an existing compartment.
    :param initial_size: Size of the new compartment (default = 0)
    :param attributes: Dictionary of attribute names and values to be associated with the new compartment.

    :return: Name of the added compartment.

    """
    if name in self._compartments:
      raise ValueError('"{}" is already a compartment in the model'.format(name))
    self._compartments[name] = Compartment(name, initial_size, **attributes)
    return name

  def get_compartment(self, name):
    """ Find a Compartment by name.

    :param name: Name of the Compartment.
    :return: The Compartment, or None if not found.

    """
    return self._compartments.get(name, None)

  def get_reachable_compartments(self, name):
    """ Find all Compartments that are reachable from a starting Compartment through one or more transitions.

    :param name: Name of the starting Compartment. Must be a Compartment in the Model.
    :return: A iterable of Compartments reachable from the named Compartment.

    """
    if name not in self._compartments:
      raise ValueError('"{}" is not the name of any compartment in the model'.format(name))

    reachable = set()
    queue = deque()
    for c in self._compartments.itervalues():
      if self.get_transition(name, c.get_name()) > 0:
        queue.append(c)

    while queue:
      current = queue.popleft()
      if current not in reachable:
        reachable.add(current)
        for c in self._compartments.itervalues():
          if self.get_transition(current.get_name(), c.get_name()) > 0:
            queue.append(c)
    return reachable

  def select_compartments(self, selector):
    """ Finds Compartments in the Model that fulfill a given condition.

    :param selector: A boolean function that takes a Compartment as argument.
    :return: an iterable of Compartments for which 'selector' returns true.

    """
    return filter(selector, self._compartments.itervalues()) 

  def set_transition(self, from_name, to_name, prob):
    """ Creates (or overwrites) a transition between two Compartments in the model.

    :param from_name: The transition starts at the Compartment with this name. Must be in a Compartment in the Model.
    :param to_name: The transition ends at the Compartment with this name. Must be in a Compartment in the Model. Can be the same as from_name, i.e. a self-transition is allowed.
    :param prob: The probability associated with this transition.
    :return: None. 

    """
    if from_name not in self._compartments:
      raise ValueError('Add "{}" as a compartment to this model first.'.format(from_name))
    if to_name not in self._compartments:
      raise ValueError('Add "{}" as a compartment to this model first.'.format(to_name))
    if prob < 0 or prob > 1:
      raise ValueError('Transition probability must be in [0, 1].')
    self._transitions[(from_name, to_name)] = prob

  def get_transition(self, from_name, to_name):
    """ Finds the transition probability between two Compartments in the Model.

    :param from_name: The transition starts at the model compartment with this name. Must be in a Compartment in the Model.
    :param to_name: The transition ends at the model compartment with this name. Must be in a Compartment in the Model.
    :return: The transition probability between the two compartments. If the transition has
    not been explicitly set, returns 0.

    """
    if from_name not in self._compartments:
      raise ValueError('"{}" are not compartments in this model.'.format(from_name))
    if to_name not in self._compartments:
      raise ValueError('"{}" are not compartments in this model.'.format(to_name))
    return self._transitions.get((from_name, to_name), 0)

  def update_transition(self, from_name, to_name, prob):
    """ Sets the transition probability between two Compartments in the Model, *and* rescales the outgoing
    transition probabilities of the "from" compartment.

    :param from_name: The transition starts at the model compartment with this name.
    :param to_name: The transition ends at the model compartment with this name.
    :param prob: The transition probability.
    :return: None.

    .. seealso:: `src.compartmodel.Model.set_transition`_

    .. note::

      * The transition need not have been explicitly set before. 

      * If the given transition probability is 1, then all other outgoing transition probabilities are rescaled 
      to 0, and their weights relative to each other are lost.

    """

    previous_prob = self.get_transition(from_name, to_name)
    self.set_transition(from_name, to_name, prob)

    # the rest of the transitions should be rescaled so they sum to (1 - prob)
    for (f, t), p in self._transitions.iteritems():
      if f == from_name and t != to_name:
        if previous_prob < 1:
          self._transitions[f, t] = (p / (1.0 - previous_prob)) * (1 - prob)
        else:
          self._transitions[f, t] = 0

    self.validate_transitions()

  def chain_compartments(self, names, prob):
    """ Links the named Compartments into a chain.
    
    :param names: An iterable of Compartment names. Compartments must be in the Model.
    :param prob: The transition probability between each adjacent pair of Compartments.
    :return: None.

    """

    for i in range(len(names) - 1):
      self.set_transition(names[i], names[i + 1], prob)

  def _get_outgoing_probs(self):

    outgoing_probs = dict((name, 0) for name in self._compartments.iterkeys())
    for (from_name, to_name), prob in self._transitions.iteritems():
      outgoing_probs[from_name] += prob
    return outgoing_probs

  def set_sink(self, sink_name):
    """ Sets a Compartment (sink) to receive transitions from all other Compartments such that for every Compartment (including the sink), its outgoing probabilities sum to 1.
    
    :param sink_name: Name of the Compartment to be made the sink. Compartment must exist in the Model.
    :return: None.

    .. note::

       set_sink can only be called once. Calling it a second time on the same Model instance raises a RuntimeError.

    """
    if self._sink is not None:
      raise RuntimeError('Sink has already been set.')
    if sink_name not in self._compartments:
      raise ValueError('Add "{}" as a compartment to this model first'.format(sink_name))
    for name, prob in self._get_outgoing_probs().iteritems():
      if prob > 1:
        raise RuntimeError('Outgoing probabilities from "{}" ({}) is greater than 1'.format(name, prob))
      self.set_transition(name, sink_name, 1 - prob)

  def validate_transitions(self):
    """ Check that for each Compartment in the Model, the sum of all outgoing probabilities is 1. """

    for name, prob in self._get_outgoing_probs().iteritems():
      if abs(prob - 1) > 1e-6:
        raise RuntimeError('Outgoing probabilities from "{}" ({}) is not 1.'.format(name, prob))
    return True

  def get_combined_size(self, selector=lambda c: True):
    """ Get the total size of a subset of Compartments. 
      
    :param selector: A boolean function that takes a Compartment as argument.
    (Default selector always returns True)
    :return: The total size of the subset of Compartments for which selector returns True.

    """
    subset = filter(selector, self._compartments.itervalues()) 
    return sum(c.get_size() for c in subset)

  def populate(self, total_size):
    """ Sets the total size of all Compartments, equally distributed between compartments. """

    size = float(total_size) / len(self._compartments)
    for name, c in self._compartments.iteritems():
      c.set_size(size)

  def advance(self):
    """ Runs one step of the model, updating the sizes of each Compartment according to the
    transitions. 

    .. note:: 

    Raises an Exception if there is a Compartment for which the sum of its outgoing probabilities is not 1.

    """
    self.validate_transitions()
    
    new_sizes = dict((name, 0) for name in self._compartments.iterkeys())
    for (from_name, to_name), prob in self._transitions.iteritems():
      flow = prob * self._compartments[from_name].get_size()
      new_sizes[to_name] += flow
    
    for name, size in new_sizes.iteritems():
      self._compartments[name].set_size(size)

  def save_schematic(self, filepath):   
    """ Saves schematic showing the Compartments and transitions of the Model.

    :param filepath: Save the PDF to this path.

    """
    graph = graphviz.Digraph(format='pdf')
    for name, c in self._compartments.iteritems():
      graph.node(name)
    for (from_name, to_name), prob in self._transitions.iteritems(): 
      graph.edge(from_name, to_name, label=str(prob))
    return graph.render(filename=filepath, cleanup=True)
