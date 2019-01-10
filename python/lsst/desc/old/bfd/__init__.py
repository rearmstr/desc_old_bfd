# This file is part of desc_bfd.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
from lsst.meas.base import BasePlugin, wrapSimpleAlgorithm

#from .moment import *
#from .test import *
from .bfdkmoment import *
#from .baseMeasure import *
#from .measureImage import *
#from .measureCoadd import *
#from .processBfdPatch import *
#from .processBfdCoadd import *
#from .processCoadd import *
#from .measurePrior import *


wrapSimpleAlgorithm(BfdKMoment, name="bfdKMoment",
                    Control=BfdKMomentControl, executionOrder=BasePlugin.SHAPE_ORDER)
