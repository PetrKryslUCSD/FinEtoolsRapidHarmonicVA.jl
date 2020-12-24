# FinEtoolsRapidHarmonicVA.jl: Using Coherent Node Cluster model reduction for harmonic vibration analysis

The foundation model-reduction algorithm is implemented in [`FinEtoolsRapidEig`](https://github.com/PetrKryslUCSD/FinEtoolsRapidEig.jl).

## Installation

Use the Git
```
git clone https://gitlab.com/PetrKrysl/FinEtoolsRapidHarmonicVA.git
```
to clone the repository into some suitable folder. This will create the folder `FinEtoolsRapidHarmonicVA`.

The Julia environment needs to be initialized. This first part is best done by running Julia _without starting up multiple workers_. Change the current folder to be `FinEtoolsRapidHarmonicVA`. Start Julia and run
```
include("go.jl")
```

## Usage

For instance for the aluminum cylinder example (folder `alucyl`) change the working folder to this directory
```
cd(".\\examples\\alucyl")  
```
Unzip the meshes
```
unzip meshes.zip 
```
Now run the file `script.jl`
```
include("script.jl")
```
