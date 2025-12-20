from pyedflib import highlevel
import numpy as np 
import matplotlib.pyplot as plt
file_path = 'Example/chb12_29.edf'

# Read signals, signal headers, and file header in a single call
signals, signal_headers, header = highlevel.read_edf(file_path)

print(f"Signals shape: {signals.shape}")
print(f"Signal labels: {header['signal_labels']}") # Corrected to access labels from header dict
print(f"Sample frequency of first channel: {signal_headers[0]['sample_frequency']} Hz")
print(f"Patient Name: {header['patientname']}")

# You can then plot the signals using a library like matplotlib (optional)
# import matplotlib.pyplot as plt
plt.plot(signals[0])
plt.show()
