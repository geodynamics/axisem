#!/usr/bin/env python

import ConfigParser
import sys
config = ConfigParser.RawConfigParser()
config.read(sys.argv[1]+'/inpython.cfg')
print(int(config.get('MESHER_BASIC','NTHETA_SLICES'))*int(config.get('MESHER_BASIC','NRADIAL_SLICES')))

