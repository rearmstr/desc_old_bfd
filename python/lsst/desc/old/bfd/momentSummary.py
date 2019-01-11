#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import pandas as pd
from .measureCoadd import MeasureCoaddTask, MeasureCoaddConfig
from astropy.table import Table

__all__ = ("MomentSummaryConfig", "MomentSummaryTask")


class MomentSummaryConfig(MeasureCoaddConfig):
    pass


class MomentSummaryTask(MeasureCoaddTask):
    ConfigClass = MomentSummaryConfig

    def __init__(self, schema=None, **kwargs):
        """
        """
        MeasureCoaddTask.__init__(self, **kwargs)

    def run(self, priorList):
        """Main driver
        """

        for filename in priorList:
            try:
                prior = Table.read(filename)
                data = {}
                data['flux'] = prior['bfd_moments'][:, 0]
                data['mx'] = prior['bfd_moments'][:, 1]
                data['my'] = prior['bfd_moments'][:, 2]
                data['msize'] = prior['bfd_moments'][:, 3]
                data['me1'] = prior['bfd_moments'][:, 4]
                data['me2'] = prior['bfd_moments'][:, 5]
                #for i in range(36):
                #    data['d%d' % i] = prior['bfd_momentsDeriv'][:, i]
                data['weight'] = prior['bfd_weight']
                df = pd.DataFrame(data)
                outname = filename.replace('momentPrior', 'momentSummary')
                outname = outname.replace('fits', 'pkl')

            except Exception as e:
                print('failed',e)
                continue
            print('Creating file',outname)
            df.to_pickle(outname)

        return None
