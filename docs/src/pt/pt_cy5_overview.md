# PT_CY5.sage

_A SageMath package for computing K-theoretic PT invariants with Nekrasov-Okounkov insertions._

We build on the SageMath package [boxcounting](https://mgemath.github.io/GKMtools.jl/) by [Henry Liu]() to compute K-theoretic Pandharipande-Thomas invariants of toric threefolds with a Nekrasov-Okounkov insertion.

## Contents

* [Installation instructions](@ref "Installation")
* [Usage instructions/examples from the paper](pt_example.md)
* [Evidence](pt_conj/overview.md): for each numerical PT experiment cited as *Evidence* in the paper, we provide code to reproduce the stated results.


## Installation
First, install **SageMath** following the instructions [here](https://doc.sagemath.org/html/en/installation/index.html).
After that, download this GitHub project. Navigate to the folder `CY5/pt/`. Next, download all files from Liu's [boxcounting](https://mgemath.github.io/GKMtools.jl/) package and place them into a folder `CY5/pt/boxcounting/`. This should result in the following folder structure:

```
CY5/
├── pt/
│   ├── boxcounting/            # Liu's boxcounting package
│   │   ├── bare_vertex.sage
│   │   ├── edge.sage
│   │   ├── pt_configuration.sage
│   │   ├── setup.sage
│   │   └── ...
│   ├── PT_CY5.sage
│   └── ...
├── test/
│   ├── omega_data/             # Membrane index formulae
│   └── ...
└── ...
```

!!! warning
    If the folder omega_data containing conjectural formulae for membrane indices is moved to a different place, the string `omega_data_folder` in `local_P2.sage` and `local_P3.sage` has to be updated accordingly for the scripts to run.
