import time
import sys
import os

# sys.path.append('../')
from Loger import log
import subprocess
from session import *

# print('Start', file=sys.stderr)
# for i in range(int(5e6)):
# 	print(i)
# print('End', file=sys.stderr)


def run_waveform():
	currdir = os.path.split(os.path.realpath(__file__))[0]
	rt_dir = os.path.join(currdir + '/bin/rt/' + current_rt_program)
	os.environ['PATH'] = rt_dir + os.pathsep + os.environ['PATH']

	excuable = 'lessrt'
	xml_filename = 'lidar_main.xml'

	scene_file_path = os.path.join(session.get_scenefile_path(), xml_filename)

	distFile = os.path.join(session.get_output_dir(), "waveform")

	command = [excuable, scene_file_path, '-o', distFile]
	log('INFO: ', command)
	subprocess.call(command)

if __name__ == '__main__':
	run_waveform()