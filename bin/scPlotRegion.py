#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from pygenometracks import plotTracks


def main(args=None):

    plotTracks.main(args)
