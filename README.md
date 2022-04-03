# Dynamic-Analysis-of-Neural-Encoding-by-Point-Process-Adaptive-Filtering
Unoffical repo for simulations of the article "[Dynamic Analysis of Neural Encoding by Point Process Adaptive Filtering](http://www.stat.columbia.edu/~liam/teaching/neurostat-spr11/papers/brown-et-al/eden2004.pdf)".

## Environment

Developed on MATLAB R2017a, 64-bits (win64).

## How to run the code

The artice proposed a new scheme to estimate **biological information** from the **neural spike train**. There are 2 experiments in the article. Simulations for each experiment contains a file for spike train data generation, and a file for biological information estimation. 

### Experiment 1: Spatial Receptive Field Analysis

1. Run `SpikeDataGenerationV1.m` or `SpikeDataGenerationV2.m` to generate the `dataSpikeSimulation.mat`, which contains generated neural spike time.
2. Run `SpikeDataFiltering.m` to get filtered estimations and watch their plots.

### Experiment 2: Movement Decoding

1. Run `SpikeTrainEncodingV1.m` or `SpikeTrainEncodingV2.m` to generate the `dataSpikeTrainDecoding.mat`, which contains generated neural spike time.
2. Run `SpikeTrainDecoding_SSPPF.m` to get filtered estimations and watch their plots.

## File Descriptions

The process of spike train data generation refers to "[The Time-Rescaling Theorem and Its Application to Neural Spike Train Data Analysis](http://www.stat.columbia.edu/~liam/teaching/neurostat-fall13/papers/brown-et-al/time-rescaling.pdf)".

- Experiment1_SpatialReceptiveFieldAnalysis/
  - SpikeDataGenerationV1.m: Spike train data generation using discrete methods.
  - SpikeDataGenerationV2.m: Spike train data generation using integrals.
  - SpikeDataFiltering.m: Receptive field model parameters estimation & evaluation.
  - Function_PassByPass.m: Implementations of the Pass-By-Pass method.
  - Function_EKF.m: Implementations of the EKF method.
  - Function_SDPPF.m: Implementations of the SDPPF method.
  - Function_SSPPF.m: Implementations of the SSPPF method.
  - dataSpikeSimulation.mat: Data strusures contains generated neural spike time by SpikeDataGeneration.m.
- Experiment2_SpikeTrainDecoding/
  - SpikeTrainEncodingV1.m: Spike train data generation using discrete methods.
  - SpikeTrainEncodingV2.m: Spike train data generation using integrals.
  - SpikeTrainDecoding_SSPPF.m: Movement & modulation parameters estimation.
  - dataSpikeTrainDecoding.mat: Data strusures contains generated neural spike time by SpikeTrainEncoding.m.

## Result

To be discussed.

## Contacts

email: hxxxd1k@yeah.net

