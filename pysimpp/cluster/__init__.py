# -*- coding: utf-8 -*-
from . import mcluster

def _is_command(): return True
def _short_description(): return mcluster._short_description()
def _command(): mcluster.command()
