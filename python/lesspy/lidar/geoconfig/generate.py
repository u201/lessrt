import sys
import json
import generateTlsGeometryConfiguration
import generateAlsGeometryConfiguration

print('Start')

LIDAR_CONFIG_PATH = 'lidar.conf'

c = {}
with open(LIDAR_CONFIG_PATH, 'r') as f:
    c = json.load(f)

platform = c['platform']

if platform['type'] == 'TLS':
    generateTlsGeometryConfiguration.generate(platform)
elif platform['type'] == 'ALS':
	generateAlsGeometryConfiguration.generate(platform)

print('End')
