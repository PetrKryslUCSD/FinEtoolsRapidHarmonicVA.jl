# FinEtoolsRapidHarmonicVA.jl

** Using Coherent Node Cluster model reduction for harmonic vibration analysis. **

The foundation model-reduction algorithm is implemented in [`FinEtoolsRapidEig`](https://github.com/PetrKryslUCSD/FinEtoolsRapidEig.jl).
Refer to a draft of the forthcoming paper that can be found in the `docs` folder of this repository.

## Installation

Use the Git
```
git clone https://gitlab.com/PetrKrysl/FinEtoolsRapidHarmonicVA.git
```
to clone the repository into some suitable folder. This will create the folder `FinEtoolsRapidHarmonicVA`.

The Julia environment needs to be initialized.  
Change the current folder to be `FinEtoolsRapidHarmonicVA`. Start Julia and run
```
include("top.jl")
```

## Usage

For instance for the twisted bar example (folder `twisted_bar`) execute 
the file
```
cd("./examples/twisted_bar/investigate_twisted_bar_direct.jl")  
```

