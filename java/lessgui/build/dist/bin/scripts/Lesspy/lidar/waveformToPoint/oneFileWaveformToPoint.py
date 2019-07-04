import sys

import numpy as np

from waveformToPointConvertor import WaveformToPointConvertor

# pwd = simulation path
batch_file_path = r'Parameters/_scenefile/lidarbatch.txt'
record_file_path = r'Results/waveform/waveform.txt'
output_file_path = r'Results/waveform/cloud.txt'

def toPoint(number_of_bins):

	batch = np.loadtxt(batch_file_path)

	records = np.loadtxt(record_file_path)

	points = []

	convertor = WaveformToPointConvertor()

	for i in range(batch.shape[0]):
		convertor.init()
		
		# convertor.accumulation_filename = record_file_folder + str(i) + '.txt'
		convertor.origin = batch[i, 0:3]
		convertor.direction = batch[i, 3:6]

		convertor.pulse_index = i

		# convertor.load()
		convertor.accumulation = records[i * number_of_bins : (i + 1) * number_of_bins, :]

		convertor.convert()
		points += convertor.points

	points = np.array(points)
	np.savetxt(output_file_path, points)

if __name__ == '__main__':
	toPoint(int(sys.argv[1]))


