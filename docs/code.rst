Documentation for DT
********************

DT Module
=========

.. automodule:: DT_Module

DTsystem_C
----------------
.. autoclass:: DT_Module.DTsystem_C
    :members:

DTcontrol_C
-----------------
.. autoclass:: DT_Module.DTcontrol_C
    :members:

DTmeshInput_C
-------------------
.. autoclass:: DT_Module.DTmeshInput_C
    :members:

DTparticleInput_C
------------------------
.. autoclass:: DT_Module.DTparticleInput_C
    :members:

.. DTfieldSolveInput_C
.. -------------------------
.. .. autoclass:: DT_Module.DTfieldSolveInput_C
..     :members:

DTtrajectoryInput_C
-------------------------
.. autoclass:: DT_Module.DTtrajectoryInput_C
    :members:

.. Field Module
.. ============
.. 
.. .. automodule:: Field_Module

.. Field_C
.. -------------

.. .. autoclass:: Field_Module.Field_C
..     :members:



.. Mesh Module
.. ===========

.. .. automodule:: Mesh_Module

.. Mesh_C
.. ------------
.. .. autoclass:: Mesh_Module.Mesh_C
..     :members: compute_cell_dict, compute_cell_vertex_dict, compute_cell_index



Particle Module
===============

.. automodule:: Particle_Module

Particle_C
----------------
.. autoclass:: Particle_Module.Particle_C
    :members: initialize_distributions, create_from_list, create_from_functions, get_species_particle_count, get_total_particle_count, move_particles_in_electrostatic_field, move_particles_in_electrostatic_potential, move_particles_in_uniform_fields, record_trajectory_data



Trajectory Module
=================

.. automodule:: Trajectory_Module

Trajectory_C
------------------
.. autoclass:: Trajectory_Module.Trajectory_C
    :members:



SegmentedArray Module
=====================

.. automodule:: SegmentedArray_Module

SegmentedArray_C
----------------------
.. autoclass:: SegmentedArray_Module.SegmentedArray_C
    :members:



UserUnits Module
=================

.. automodule:: UserUnits_Module

MyPlasmaUnits_C
---------------------
.. autoclass:: MyPlasmaUnits_C
    :members:

MyCGSUnits_C
------------------
.. autoclass:: MyCGSUnits_C
    :members:



Dolfin Module
=============

.. automodule:: Dolfin_Module

Mesh_C
------------
.. autoclass:: Dolfin_Module.Mesh_C
    :members: compute_cell_dict, compute_cell_index, compute_cell_indices, compute_cell_neighbor_dict, compute_cell_vertex_dict

Field_C
-------------
.. autoclass:: Dolfin_Module.Field_C
    :members: interpolate_field_to_points, integrate_delta_function

.. ParticleToField_C
.. -----------------
.. .. autoclass:: Dolfin_Module.ParticleToField_C
..     :members: add_point_source

PoissonSolve_C
--------------
.. autoclass:: Dolfin_Module.PoissonSolve_C
     :members: solve_for_phi, compute_negE

UserMesh FE_XYZ Module
======================

.. automodule:: UserMesh_FE_XYZ_Module
   :show-inheritance:

UserMesh_C
----------------
.. autoclass:: UserMesh_FE_XYZ_Module.UserMesh_C
   :show-inheritance:


UserMesh_y_Fields Spherical1D Module
====================================

.. automodule:: UserMesh_y_Fields_Spherical1D_Module

UserMesh_C
----------------
.. autoclass:: UserMesh_y_Fields_Spherical1D_Module.UserMesh_C
   :show-inheritance:

UserPoissonSolve_C
------------------------
.. autoclass:: UserMesh_y_Fields_Spherical1D_Module.UserPoissonSolve_C
   :show-inheritance:


UserMesh_y_Fields FE2D Module
====================================

.. automodule:: UserMesh_y_Fields_FE2D_Module
   :show-inheritance:

UserMesh_C
----------------
.. autoclass:: UserMesh_y_Fields_FE2D_Module.UserMesh_C
   :show-inheritance:

UserPoissonSolve_C
------------------------
.. autoclass:: UserMesh_y_Fields_FE2D_Module.UserPoissonSolve_C
   :show-inheritance:
