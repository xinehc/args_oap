import subprocess
import os

def stage_one(options):

    _path = os.path.dirname(os.path.realpath(__file__))
    _stage_one = os.path.join(_path, 'argoap_pipeline_stageone_version3')

    args = [
    _stage_one, 
    '-i', options.i, 
    '-o', options.o,
    '-m', options.m,
    '-n', str(options.n),
    '-f', options.f,
    '-q', str(options.q),
    '-z', str(options.z),
    '-x', str(options.x),
    '-y', str(options.y),
    '-v', str(options.v),
    '-w', str(options.w)]

    subprocess.call(args)
    