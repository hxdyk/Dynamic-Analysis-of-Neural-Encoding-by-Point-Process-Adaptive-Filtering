# Dynamic-Analysis-of-Neural-Encoding-by-Point-Process-Adaptive-Filtering
Unoffical simulations of the article "[Dynamic Analysis of Neural Encoding by Point Process Adaptive Filtering](http://www.stat.columbia.edu/~liam/teaching/neurostat-spr11/papers/brown-et-al/eden2004.pdf)".

## Running Environment

Developed on MATLAB R2017a 64-bits, Windows 11.

## File Descriptions

The artice proposed a new scheme to estimate **biological information** from the **neural spike train**. There are 2 experiments in the article. Simulations for each experiment contains a file for spike train data generation, and a file for biological information estimation. The process of spike train data generation refers to "[The Time-Rescaling Theorem and Its Application to Neural Spike Train Data Analysis](http://www.stat.columbia.edu/~liam/teaching/neurostat-fall13/papers/brown-et-al/time-rescaling.pdf)".

- Experiment1_SpatialReceptiveFieldAnalysis/
  - SpikeDataGenerationV1.m: Spike train data generation using discrete methods.
  - SpikeDataGenerationV2.m: Spike train data generation using integrals.
  - SpikeDataFiltering.m: Receptive field model parameters estimation & evaluation.
  - Function_PassByPass.m: Implementations of the Pass-By-Pass method.
  - Function_EKF.m: Implementations of the EKF method.
  - Function_SDPPF.m: Implementations of the SDPPF method.
  - Function_SSPPF.m: Implementations of the SSPPF method.
  - dataSpikeSimulation.mat: Data strusures contains generated neural spike times by SpikeDataGeneration.m.
- Experiment2_SpikeTrainDecoding/
  - SpikeTrainEncodingV1.m: Spike train data generation using discrete methods.
  - SpikeTrainEncodingV2.m: Spike train data generation using integrals.
  - SpikeTrainDecoding_SSPPF.m: Movement & modulation parameters estimation.
  - dataSpikeTrainDecoding.mat: Data strusures contains generated neural spike times by SpikeTrainEncoding.m.

## Contacts

email: hxxxd1k@yeah.net

