# Convex-Array-Ultrasound-Beamforming-MATLAB

A comprehensive implementation of a Delay-and-Sum ultrasound beamforming algorithm designed for convex array transducers. This repository contains a complete MATLAB pipeline for processing raw RF channel data into diagnostic-quality B-mode images. Key features include dynamic/fixed receive focusing, multiple apodization windows (Rectangular, Hann, Hamming, Blackman), linear interpolation for fractional delays, and quantitative FWHM analysis for lateral resolution assessment. Validated on simulated wire phantoms, tissue-mimicking phantoms, and clinical in-vivo liver data.

## Authors

- Leonel Vázquez Carrasco
- Heshan Harshana Jayasinghe Arachchige

*Course: 22485 Medical Imaging Systems*

## Results Summary

### Validation

The Dynamic-Hann configuration provides a 6.40 mm FWHM at 120 mm depth, a 23% improvement over the rectangular window (8.29 mm) and significantly superior to fixed focusing (26.42 mm).

| Configuration | FWHM @ 30mm | FWHM @ 120mm | Key Observation |
|------|------|------|------|
| Fixed-Rectangular | 0.99 mm | 26.42 mm | Only sharp at focal zone |
| Dynamic-Rectangular | 0.95 mm | 8.29 mm | Best resolution, high sidelobes |
| Dynamic-Hann | 0.92 mm | 6.40 mm | Optimal balance |
| Dynamic-Hamming | 0.93 mm | 6.51 mm | Near-optimal balance |
| Dynamic-Blackman | 0.93 mm | 7.02 mm | Best contrast, poor far-field |

## Key Features

- **Convex Array Geometry Support**: Accurate time-delay calculations for curved transducer arrays
- **Dynamic & Fixed Receive Focusing**: Adaptive focusing at every depth point
- **Multiple Apodization Windows**: Rectangular, Hann, Hamming, and Blackman
- **Fractional Delay Interpolation**: Linear interpolation for sub-sample accuracy
- **Envelope Detection**: Hilbert transform-based demodulation
- **Log Compression**: 60 dB dynamic range for optimal display
- **Scan Conversion**: Polar to Cartesian interpolation
- **Quantitative FWHM Analysis**: Lateral resolution measurement at -6 dB

## Getting Started

### Prerequisites

- **MATLAB R2020a** or later
- Required toolboxes:
  - Signal Processing Toolbox

### Running the Pipeline

## Running Complete Workflow
```matlab
beamform_pipeline();
```

## Basic Usage
```matlab
% Simulated wire phantom with dynamic focusing and Hann window
[img, x, z] = beamform(1, 'hann', true, 60);

% Display the image
figure; imagesc(x, z, img); axis image; colormap(gray);
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
```

## Function Documentation
``` [img, x, z] = beamform(id_data, id_apod, id_focus, db_level) ``` 
Main beamforming function.

| Input | Type | Description | Options |
|------|------|------|------|
| ``` id_data ``` | integer | Data selection | 1=simulated, 2=phantom, 3=in-vivo |
| ``` id_apod ``` | string | Apodization window | 'rect', 'hann', 'hamming', 'blackmann' |
| ``` id_focus ``` | boolean | Focusing method | true=dynamic, false=fixed |
| ``` db_level ``` | integer | Compression range | 60 (default) |

Outputs:

- ``` img ```: Beamformed image matrix
- ``` x ```: Lateral coordinates [mm]
- ``` z ```: Axial coordinates [mm]

## Repository Structure

```
├── data/
│ ├── convex_array_files/        # Clinical in-vivo data
│ ├── phantom_matrix_files/      # Measured phantom data
| └── phantom_matrix_sim_files/  # Simulated RF data
├── results/
│ ├── fwhm_plot.png              # FWHM Comparison Plot
│ ├── invivo_liver.png           # In-vivo data beamformed image
│ ├── meas_phantom.png           # Measured phantom data beamformed image
│ ├── sim_dyn_blackman.png       # Dynamic-Blackman Simulated data beamformed image
│ ├── sim_dyn_hamming.png        # Dynamic-Hamming Simulated data beamformed image
│ ├── sim_dyn_hann.png           # Dynamic-Hann Simulated data beamformed image
│ ├── sim_dyn_rect.png           # Dynamic-Rectangular Simulated data beamformed image
│ └── sim_fixed_rect.png         # Fixed-Rectangular Simulated data beamformed image
├── src/
│ └── beamform_pipeline.m         # Main Beamform Workflow
└── README.md                     # This file
```

## Methodology

### Beamforming Pipeline
1. **Geometry Definition**: Calculate element positions for convex array
2. **Delay Calculation**: Time-of-flight for each element to focal point
3. **Fractional Delay**: Linear interpolation for sub-sample accuracy
4. **Apodization**: Element weighting for sidelobe control
5. **Coherent Summation**: Sum delayed and weighted signals
6. **Envelope Detection**: Hilbert transform demodulation
7. **Log Compression**: 60 dB dynamic range
8. **Scan Conversion**: Polar to Cartesian interpolation

### Key Parameters

| Parameter | Value | Description |
|------|------|------|
| Number of Elements | 192 | Total elements in transducer |
| Active Elements | 64 | Elements used per emission |
| Center Frequency | 3.5 MHz | Transducer frequency |
| Sampling Rate | 17.5 MHz | Data acquisition frequency |
| Speed of Sound | 1491 m/s | Acoustic velocity |
| Convex Radius | 60.3 mm | Array curvature radius |
| Dynamic Range | 60 dB | Display compression |

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

* Professor Jørgen Arendt Jensen for guidance (DTU course 22485)
* Professor Billy Y. S. Yiu for guidance (DTU course 22485)
