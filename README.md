HAS_for_Kranion

Hybrid Angular Spectrum (HAS) and Thermal Simulation for Transcranial Focused Ultrasound Treatment Planning

This software runs the Hybrid Angular Spectrum (HAS) acoustic wave propagation method and a thermal solver (Pennes Bioheat Equation) using treatment data exported from the Kranion platform. It enables patient-specific simulation of ultrasound propagation and temperature rise based on individual CT scans and treatment parameters.

Overview

The process begins with the FUSF_HAS_Prep.m script, which loads and preprocesses a parsed Kranion export file. Users can then simulate multiple treatment exports using the FUSF_HAS_Prep_loop.m script.

The full pipeline includes:

Loading and parsing patient-specific treatment data from Kranion export files

Segmenting the Insightec hemispherical transducer into 7 segments compatible with the HAS method

Computing the Rayleigh-Sommerfeld integral for each transducer segment

Saving the "ERFA" wavefield planes and associated tilt angles for use in HAS propagation

Generating a tissue acoustic property model (Modl file) from patient CT volumes using threshold-based segmentation

Running the 7-segment HAS acoustic simulation to compute 3D pressure fields

Estimating temperature rise using a finite difference time domain solver of the Pennes Bioheat Equation

Exporting simulated pressure and temperature fields in a Kranion-compatible format for visualization and analysis

Getting Started

Begin by running FUSF_HAS_Prep.m on your Kranion export file.

For batch processing, use FUSF_HAS_Prep_loop.m to iterate over multiple exports.

Authors

Taylor Webb, Zach Johnson, Dennis L. Parker, and Matthew Eames

Hybrid Angular Spectrum and Pennes' Bioheat code from Nick Todd and Douglas Christensen
June 19th, 2024
