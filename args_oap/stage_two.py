import subprocess
import os

def stage_two(options):

    _path = os.path.dirname(os.path.realpath(__file__))
    _stage_two = os.path.join(_path, 'argoap_pipeline_stagetwo_version3.py')

    args = [
    'python',
    _stage_two, 
    '-i', options.i, 
    '-o', options.o,
    '-m', options.m,
    '-n', str(options.n),
    '-l', str(options.l),
    '-e', str(options.e),
    '-d', str(options.d),
    '-db',options.db,
    '-fa',options.fa,
    '-struc1',options.struc1,
    '-struc2',options.struc2,
    '-struc3',options.struc3]

    if options.b:
        args.extend(['-b', str(options.b)])

    subprocess.call(args)