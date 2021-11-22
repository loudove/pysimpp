# -*- coding: utf-8 -*-
from . import utils
from . import readers
from . import primitives

__debug = False

def __find_commands():
    import os
    _dir = os.path.dirname(__file__)
    _items = [_dir + os.sep + _f for _f in os.listdir(_dir)]
    _pkg_name = os.path.basename(_dir)
    _pkg_items = [
        os.path.basename(_item) for _item in _items
        if os.path.isdir(_item) or os.path.isfile(_item) and _item.endswith(
            '.py') and not _item.endswith(os.sep + '__init__.py')
    ]
    import site, sys, importlib

    def __chk_module(modname):
        ''' Validate if the module is a pysimpp command i.e. if
            methods _is_command, _command and _short_description
            are available. If yes, the method return the short
            description of the module.  '''
        _mod = None
        _short = None
        try:
            try:
                _mod = importlib.import_module(_modname)
                if hasattr(_mod, '_is_command') and _mod._is_command():
                    if hasattr(_mod, '_command'):
                        _short = _mod._short_description()
            except Exception as exeption:
                if __debug: print("%s: %s"%(_modname, exeption))
                _short = None

        finally:
            if _mod is not None:
                del _mod

        return _short

    _commands = {}
    for _item in _pkg_items:
        __item = _item[:-3] if _item.endswith('.py') else _item
        _modname = _pkg_name + '.' + __item

        _short = __chk_module(_modname)
        if _short is not None:
            _commands[__item] = (_modname,_short)

    _etc_path = os.getenv('PYSIMPPPATH')
    if _etc_path is not None:
        for _dir in _etc_path.split(':'):
            if os.path.isdir(_dir):
                _d = os.path.expanduser(_dir)
                site.addsitedir(_d)
                for _f in os.listdir(_d):
                    _modname = _f if not _f.endswith('.py') else _f[:-3]
                    _short = __chk_module(_modname)
                    if _short is not None:
                        if _f in _commands:
                            print( 'WARNING: command "%s" (%s) is already registed in pysimpp.' % ( _modname, _d))
                        else:
                            _commands[_modname] = (_modname,_short)

    return _commands


def __version(file):
    import os
    _dir = os.path.dirname(__file__)
    try:
        with open(_dir + os.sep + file, 'r') as (f):
            version = f.readline().strip()
    except:
        version = 'x.x.x'

    return version


def __help(version, commands):
    print ('\n(pysimpp v%s)' % version)
    print ('''\n
Simulation post process command line:
  % pysimpp command [options]
  
List of available commands:
    help : Print this message.''')
    for k in sorted(commands.keys()):
        _, h = commands[k]
        print ('    %s : %s' % (k, h))

    print ('''\n
To access detailed help on a specific topic use:
  % pysimpp command -h\n''')


def main():
    """ pysimpp entry point for executing the command found in sys.argv[1]. """
    from os.path import dirname, basename
    import sys
    from importlib import import_module
    _command = sys.argv[1] if len(sys.argv) > 1 else None
    if _command is None or _command.lower() in ('help', '-h', '--help'):
        __help(version, commands)
    elif _command not in commands:
        print ('ERROR: commnad "%s" is not available\n' % _command)
        __help(version, commands)
    else:
        sys.argv.pop(0)
        _dir = dirname(__file__)
        _modname = commands[_command][0]
        _mod = import_module( _modname)
        _mod._command()
    return


version = __version('version')
commands = __find_commands()

if __name__ == '__main__':
    main()
