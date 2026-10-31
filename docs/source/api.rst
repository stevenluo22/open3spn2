open3SPN2 Module Documentation
==============================

.. automodule:: open3SPN2
    :noindex:

Classes and Functions
---------------------

DNA Class
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.DNA
    :member-order: bysource
    :members:

System Class
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.System
    
    .. automethod:: __init__
    
    .. rubric:: Methods

DNA Forces
^^^^^^^^^^

.. autoclass:: open3SPN2.Bond
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Angle
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Dihedral
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Stacking
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.BasePair
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.CrossStacking
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Exclusion
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.Electrostatics
    :member-order: bysource
    :members:

DNA-Protein Forces
^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.ExclusionProteinDNA
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.ElectrostaticsProteinDNA
    :member-order: bysource
    :members:

Numpy Energies
^^^^^^^^^^^^^^

.. automodule:: open3SPN2.force.backend
    :members:

Alchemy: Base-Pair Identities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: open3SPN2.alchemy.AlchemicalTransformation
    :member-order: bysource
    :members:

.. autoclass:: open3SPN2.alchemy.Alphabet

.. autofunction:: open3SPN2.alchemy.alphabet

.. autofunction:: open3SPN2.alchemy.identity_tables

.. autoclass:: open3SPN2.alchemy.IdentityTables
    :members:

.. autofunction:: open3SPN2.alchemy.state_energies

.. autofunction:: open3SPN2.alchemy.add_forces

.. autofunction:: open3SPN2.alchemy.set_state

Gibbs Identity Sampler
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: open3SPN2.alchemy.sampler
    :noindex:

.. autoclass:: open3SPN2.alchemy.GibbsSampler
    :members:

.. autofunction:: open3SPN2.alchemy.run

.. autofunction:: open3SPN2.alchemy.factor_tables

.. autofunction:: open3SPN2.alchemy.model

.. autoclass:: open3SPN2.alchemy.Factors

.. autoclass:: open3SPN2.alchemy.Model

.. autofunction:: open3SPN2.alchemy.log_weights

.. autofunction:: open3SPN2.alchemy.conditional

.. autofunction:: open3SPN2.alchemy.window_conditionals

.. autofunction:: open3SPN2.alchemy.joint

.. autofunction:: open3SPN2.alchemy.sweep

.. autofunction:: open3SPN2.alchemy.adapt_offsets

.. autofunction:: open3SPN2.alchemy.w_target

.. automodule:: open3SPN2.alchemy.estimate
    :noindex:

.. autofunction:: open3SPN2.alchemy.markov_log_probability

.. autofunction:: open3SPN2.alchemy.sequence_offsets

.. autofunction:: open3SPN2.alchemy.free_energies

Utility Functions
^^^^^^^^^^^^^^^^^

.. autofunction:: open3SPN2.parseConfigTable

.. autofunction:: open3SPN2.parsePDB

.. autofunction:: open3SPN2.fixPDB

.. autofunction:: open3SPN2.pdb2table

Tests Functions
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: tests.TestEnergies
    :member-order: bysource
    :members:

.. autofunction:: tests.test_DNA_from_xyz

.. autofunction:: tests.test_parse_xyz

.. autofunction:: tests.test_parse_log

Exceptions
^^^^^^^^^^

.. autoclass:: open3SPN2.DNATypeError