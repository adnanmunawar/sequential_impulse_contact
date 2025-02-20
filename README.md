# Implementation of contact constraint using sequential impulse solver

# Usage

The scripts model a rigid body falling under gravity and coming into contact with the ground. The contact is modeled as a contact constraint handled using the SI solver.
### Symbolic script
Use the seq_impulse_sym.py script for a symbolic solution

### Numerical script
Use the seq_impulse_num.py for a numerical solution and a plot of the rigid body trajectory

# Good resources on the topic
1. https://danielchappuis.ch/download/ConstraintsDerivationRigidBody3D.pdf
2. https://allenchou.net/2013/12/game-physics-resolution-contact-constraints/
3. https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
4. http://www.mft-spirit.nl/files/articles/ImpulseSolverBrief.pdf
5. https://allenchou.net/2013/12/game-physics-constraints-sequential-impulse/

# Other implementations:
We have use this formulation for modeling contact between a snake robot and a volumetric patient model here:
https://github.com/htp2/continuum-manip-volumetric-drilling-plugin
