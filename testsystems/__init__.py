##############################################################################
# testsystems
#
# Copyright 2012-2013 MSKCC and the Authors
#
# Authors: Kyle A. Beauchamp, John D. Chodera
# Contributors:
#
# openmm-testsystems is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with openmm-testsystems. If not, see <http://www.gnu.org/licenses/>.
##############################################################################


"""testsystems contains a collection of OpenMM testsystems that can be used
to test MD-related techniques.  
"""

from testsystems.systems import (HarmonicOscillator, Diatom, ConstraintCoupledHarmonicOscillator,
HarmonicOscillatorArray, SodiumChlorideCrystal, LennardJonesCluster, LennardJonesFluid,
CustomLennardJonesFluid, IdealGas, WaterBox, AlanineDipeptideImplicit,
AlanineDipeptideExplicit, LysozymeImplicit, SrcImplicit,SrcExplicit,
MethanolBox, MolecularIdealGas, CustomGBForceSystem)
