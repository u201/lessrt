import time
import sys
import os

sys.path.append('../')
from Loger import log
import subprocess
from session import *

# print('Start', file=sys.stderr)
# for i in range(int(5e6)):
# 	print(i)
# print('End', file=sys.stderr)

excuable = 'lessrt'
scene_file_path = '_scenefile/'
xml_filename = 'lidar-main.xml'

scene_file_path = os.path.join(session.get_scenefile_path(), main_scene_xml_file)

log('INFO: ', scene_file_path)

path = [excuable, scene_file_path + xml_filename]
log('INFO: ', path)
# subprocess.call(path)


