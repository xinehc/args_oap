import subprocess
import os

def stage_two(options):

    _path = os.path.dirname(os.path.realpath(__file__))
    _stage_two = os.path.join(_path, 'stage_two_version2')

    if options.b:
        subprocess.call(
            [_stage_two, '-i', options.i, 
            '-o', options.o,
            '-m', options.m,
            '-n', str(options.n),
            '-l', str(options.l),
            '-e', str(options.e),
            '-d', str(options.d),
            '-b'])
    else:
        subprocess.call(
            [_stage_two, '-i', options.i, 
            '-o', options.o,
            '-m', options.m,
            '-n', str(options.n),
            '-l', str(options.l),
            '-e', str(options.e),
            '-d', str(options.d)])
