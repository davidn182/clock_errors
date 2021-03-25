<img alt="OBS clock errors: A Python Toolbox for correcting clock errors of OBSs." class="right" style="width: 40%" src="https://github.com/davidn182/OBS_clock_errors/blob/main/clock_errors_py/Figures/logo.png" />

# OBS_clock_errors

OBS_clock_errors is an open-source project made to correct clock errors of seismological stations.

# Abstract

Accurate timing of seismic records is essential for almost all applications in seismology. Wrong timing of the waveforms may result in incorrect Earth models and/or inaccurate earthquake locations. As such, it may render interpretations of underground processes incorrect. Ocean bottom seismometers (OBSs) experience clock drifts due to their inability to synchronize with a GNSS signal (with the correct reference time). As OBSs generally operate in relatively stable ambient temperatures, the timing deviation is usually assumed to be linear. Therefore, one can estimate the time corrections through GPS synchronization before deployment and after the instrument's recovery. However, suppose the device has run out of power before recovery (i.e., due to the battery being dead at the time of recovery). In that case, the timing error at the end of the deployment cannot be determined. Besides, the drift may not be linear, e.g., rapid temperature drop while the OBS sinks to the seabed. Here we present an algorithm that recovers the linear clock drift and a potential timing error at the onset.

The algorithm presented in this study exploits seismic interferometry (SI). Specifically, time-lapse (averaged) cross-correlations of ambient seismic noise are computed. As such, virtual-source responses, which are generally dominated by the recorded surface waves, are retrieved. These interferometric responses generate two virtual sources: a causal wave (arriving at a positive time) and an acausal wave (arriving at a negative time). Under favorable conditions, both interferometric responses approach the surface-wave part of the medium's Green's function. Therefore, it is possible to calculate the clock drift for each station by exploiting the time-symmetry between the causal and acausal waves. For this purpose, the clock drift is calculated by measuring the differential arrival times of the causal and acausal waves for many receiver-receiver pairs and computing the drift by carrying-out a least-squares inversion. The methodology described is applied to time-lapse cross-correlations of ambient seismic noise recorded on and around the Reykjanes peninsula, SW Iceland. The stations used for the analysis were deployed in the context of IMAGE (Integrated Methods for Advanced Geothermal Exploration). They consisted of 30 on-land stations and 24 ocean bottom seismometers (OBSs).  The seismic activity was recorded from spring 2014 until August 2015 at around 100 km in diameter (from the Reykjanes peninsula's tip).

### Installation

You need to download the folder.
It is necessary to install Obspy, Jupyter notebooks, and Basemap.
It is recommendable doing this via Anaconda and in a new environment as follows:
```bash
$ conda config --add channels conda-forge
$ conda create -n clock_errors python=3.7
$ conda activate clock_errors
(clock_errors) $ conda install obspy
(clock_errors) $ conda install basemap
(clock_errors) $ conda install -c conda-forge jupyterlab
(clock_errors) $ conda install -c conda-forge notebook
```

### Getting started

You can take a look at the jupyter-notebook at the **Tutorial** section for basic usage of the program. You will neet to have access to the data so send an e-mail to d.f.naranjohernandez@tudelft.nl

### References

  * Naranjo, D., Parisi, L., Jousset, P., Weemstra, C., and Jónsson, S.: Determining OBS clock drift using ambient seismic noise, EGU General Assembly 2021, online, 19–30 Apr 2021, EGU21-13999, https://doi.org/10.5194/egusphere-egu21-13999, 2021.
