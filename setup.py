import sys
import os
import glob

# from distutils.command.sdist import sdist
import distutils

import setuptools

from distutils.util import get_platform
from numpy.distutils.core import setup, Extension

def fix_include_path():
    ''' Check the options and return the include path to be added in the
        extension directives in order to avoid problems relevant with
        accessing intermediate mod files from auxiliary library builds.
    '''

    # check for the temporary build directory option
    _tempopt = None
    _chkopt = ('-t','--build-temp')
    for _opt in  _chkopt:
        if _opt in sys.argv:
            _i = sys.argv.index(_opt)
            if _i < len(sys.argv)-1:
                _tempopt = sys.argv[_i+1]
                break

    # check for the base directory option
    _buildopt = 'build'
    _chkopt = ('-b','--build-base')
    for _opt in  _chkopt:
        if _opt in sys.argv:
            _i = sys.argv.index(_opt)
            if _i < len(sys.argv)-1:
                _buildopt = sys.argv[_i+1]
                break

    if _tempopt is None:
        # works with python3 (check distutils/command/build.py)
        platform_specifier = ".%s-%d.%d" % (get_platform(), *sys.version_info[:2])
        _tempopt = '%s%stemp%s'%(_buildopt,os.sep,platform_specifier)

    return _tempopt

def error_print(message,exeption):
    ''' Print the error messages and exit. '''
    print("[error]: %s" % message)
    print("         %s" % exeption)
    import sys
    sys.exit()

def chk_version(file):
    from pkg_resources import parse_version
    try:
        with open(file,'r') as f:
            version = parse_version( f.readline().strip())
    except Exception as exeption:
        error_print("check the license file %s"%file, exeption)
    return version

def chk_long_description(file,encoding='utf-8'):
    try:
        with open(file,'r',encoding=encoding) as f:
            long_description = f.read()
    except Exception as exeption:
        error_print("check the long description file %s"%file, exeption)
    return long_description

def chk_openmp():
    from numpy.distutils.fcompiler import get_default_fcompiler, CompilerNotFound
    from distutils.util import get_platform

    # enable openmp threads
    _extra_compile_args = []
    _extra_link_args = []
    _extra_for_compile_args = []
    platform = get_platform()
    try:
        compiler = get_default_fcompiler()
        if compiler in ['gnu','gnu95','g95']:
            _extra_for_compile_args= ['-fopenmp' ]
            _extra_compile_args = ['-Xpreprocessor'] if platform.startswith('macosx') else []
            _extra_compile_args.append('-fopenmp')
            _extra_link_args = ['/openmp'] if platform.startswith('win') else ['-lgomp']
        elif compiler in ['intelv','intelvem','intelem']:
            if platform.startswith('win'):
                _extra_for_compile_args=_extra_compile_args=['-Qopenmp']
                _extra_link_args = ['-Qopenmp']
            else:
                _extra_for_compile_args=_extra_compile_args=['-qopenmp']
                _extra_link_args = ['-liomp5']
        elif compiler in ['flang']:
            if get_platform == 'darwin':
                _extra_for_compile_args=_extra_compile_args=[ '-Xpreprocessor', '-fopenmp' ]
                _extra_link_args = ['-lomp']
            else:
                _extra_for_compile_args=_extra_compile_args=['-mp']
                _extra_link_args = ['-mp']
        else:
            _extra_for_compile_args=_extra_compile_args=_extra_link_args =[]
    except:
        _extra_for_compile_args=_extra_compile_args=_extra_link_args =[]

    # it seems that extra_link_args are not provided in 
    # fcompiler.wrap_unlinkable_objects method called in 
    # numpy.distutils.command.build_ext.build_ext._process_unlinkable_fobjects
    if platform.startswith('win'):
        from numpy.distutils import fcompiler
        def omppatched_spawn(old_spawn):
            def spawn(self, cmd, *args, **kw):
                cmptype = self.compiler_type
                opt=""
                # tested only with gnu95
                if cmptype in ['gnu','gnu95','g95']:
                    opt='-fopenmp'
                elif cmptype in ['intelv','intelvem','intelem']:
                    opt='-Qopenmp'
                elif cmptype in ['flang']:
                    opt='-mp'
                if len(opt) > 0 and self.c_compiler.compiler_type == "msvc":
                    if not opt in cmd:
                        cmd.append(opt)
                return old_spawn(self, cmd, *args, **kw)
            return spawn
        fcompiler.FCompiler.spawn = omppatched_spawn(fcompiler.FCompiler.spawn)

    return _extra_for_compile_args, _extra_compile_args, _extra_link_args

# fastlapack and fastcore libraries files
# TODO check if the system lapack libraries can be retrieved.
_fastlapackfiles = glob.glob('fastpost/src/lapack/*.f')
_fastcorefiles = [
                'fastpost/src/core/vector3d_mod.F90',
                'fastpost/src/core/domain3d_mod.F90',
                'fastpost/src/core/lj_mod.F90']
_fastmodulefiles = glob.glob('fastpost/src/*.f90')

# set the metadata (TODO move them in setup.cfg)
# check:
# https://docs.python.org/3.0/distutils/setupscript.html
metadata = dict(
            name='pysimpp',
            version=chk_version('pysimpp/version'),
            license='GPL-3.0',
            url='https://github.com/loudove/pysimpp',
            download_url='https://github.com/loudove/pysimpp',
            project_urls  = {
                'Source code': 'https://github.com/loudove/pysimpp',
                'Bug Tracker': 'https://github.com/loudove/pysimpp/issues',
                'Changelog': 'https://github.com/loudove/pysimpp/blob/master/CHANGELOG.md'
            },
            maintainer='Loukas Peristeras',
            maintainer_email='l.peristeras@inn.demokrtitos.gr',
            description='Simulation post process',
            long_description=chk_long_description('README.md'),
            long_description_content_type='text/markdown',
            packages=setuptools.find_packages(include='pysimpp/*'),
            # packages=[
            #     'pysimpp',
            #     'pysimpp.cluster',
            #     'pysimpp.utils',
            #     'pysimpp.primitives',
            #     'pysimpp.readers'
            # ],
            package_dir={
                'pysimpp': 'pysimpp',
            },
            package_data={
                'pysimpp':['version'],
            },
            py_modules = [
            ],
            install_requires=[
                'numpy',
                'setuptools',
                'wheel',
                'pyyaml',
                'networkx',
                'MDAnalysis',
            ],
            # setup_requires=[ 'pytest-runner' ],
            # tests_require=[ 'pytest' ],
            # test_suite='nose.collector',
            # tests_require=['nose'],
            classifiers=[
                'Intended Audience :: Molecular modeling community'
                'Programming Language :: Python :: 3'
                'License :: OSI Approved :: GPL 3.0'
                'Operating System :: OS Independent'
            ],
            python_requires='>=3.5',
            zip_safe = False,  # the package can run out of an .egg file
            include_package_data = True,
            entry_points={
                'console_scripts': [
                    'pysimpp=pysimpp:main'
                ],
            }
            # cmdclass={'sdist': sdist}
        )

def configuration(parent_package='', top_path=None):
    # Returns a dictionary suitable for passing to numpy.distutils.core.setup(..)

    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.fcompiler import get_default_fcompiler, CompilerNotFound
    from distutils.util import get_platform

    config = Configuration(None, parent_package, top_path)

    # add fastlapack (a subset of lapack used in pysimpp) and fastcore
    # (a collection of fortran accelerated functionality) libraries
    config.add_library('fastlapack', _fastlapackfiles, language='f90')
    config.add_library('fastcore', _fastcorefiles, language='f90')

    # add fastcore build dir in the include directories to access
    # the modules neccessary for the library build.
    config.add_include_dirs( fix_include_path())

    _extra_for_compile_args, _extra_compile_args, _extra_link_args = chk_openmp()

    config.add_extension('pysimpp.fastpost',
                        sources=_fastmodulefiles,
                        f2py_options=['--quiet','--verbose'],
                        extra_link_args=_extra_link_args,
                        extra_f77_compile_args=_extra_for_compile_args,
                        extra_f90_compile_args=_extra_for_compile_args,
                        extra_compile_args=_extra_compile_args,
                        libraries=['fastlapack','fastcore'],
                        language='f90')

    return config

setup(**metadata,
      configuration=configuration,)
#      ext_modules=[fastpost])

