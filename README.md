# Continuous Variable Quantum Entanglement in Optomechanical Systems: A Short Review

[![Manuscript Version](https://img.shields.io/badge/manuscript-v2.2-red?style=for-the-badge)](https://doi.org/10.1116/5.0022349)
[![Toolbox Version](https://img.shields.io/badge/qom-v1.0.1-red?style=for-the-badge)](https://sampreet.github.io/qom-docs/v1.0.1)
[![Last Updated](https://img.shields.io/github/last-commit/sampreet/corr_review?style=for-the-badge)](https://github.com/sampreet/corr_review/blob/master/CHANGELOG.md)

> AVS Quantum Sci. 3, 015901 (2020)

Author | Affiliation
------------ | -------------
[Amarendra Kumar Sarma](https://www.iitg.ac.in/aksarma/) | Indian Institute of Technology Guwahati, Guwahati-781039, India
[Subhadeep Chakraborty](https://scholar.google.co.in/citations?user=o5n-rDUAAAAJ&hl=en) | Indian Institute of Science Education and Research Kolkata, Nadia-741246, India
[Sampreet Kalita](https://www.iitg.ac.in/stud/sampreet/) | Indian Institute of Technology Guwahati, Guwahati-781039, India

Contributing Part | SC | SK
------------ | ------------ | -------------
Literature review | 80% | 20%
Derivations of expressions | 50% | 50%
Illustrations and plots | 10% | 90%
Manuscript preparation | 60% | 40%

## About the Work

Cavity optomechanics deals with the radiation pressure induced interaction between photons and mechanical motion in a cavity.
It has promising applications in quantum information science.
In this review, the authors discuss quantum entanglement in this emerging area of research.
After giving a brief historical background on the topic of entanglement and cavity optomechanics, measures of continuous variable entanglement are discussed somewhat in great details.
This is followed by a short discussion on cavity quantum optomechanics, relevant to
the topic on entanglement.
Then the authors discuss most of the prominent ideas and proposals pertaining to entanglement research in cavity optomechanics up until now.
The authors have emphasized the key theoretical concepts without too much rigor and provided relevant experimental details whenever deemed appropriate.
Finally, the authors conclude by giving a perspective on other quantum correlations such as quantum discord and quantum synchronization.

## Notebooks

* [Stability Criteria for a Simple End-mirror Optomecanical System](notebooks/stability_criteria.ipynb)
* [Plots in the Latest Version of the Manuscript](notebooks/v2.2_qom-v1.0.1/plots.ipynb)

## Structure of the Repository

```
ROOT_DIR/
|
├───data/
│   ├───bar/
│   │   ├───baz_xyz.npz
│   │   └───...
│   └───...
|
├───notebooks/
│   ├───bar/
│   │   ├───baz.ipynb
│   │   └───...
│   │
│   ├───foo_baz.ipynb
│   └───...
|
│───scripts/
│   ├───bar/
│   │   ├───baz.py
│   │   └───...
│   └───...
│
├───.gitignore
├───CHANGELOG.md
└───README.md
```

Here, `bar` represents the version.

## Installing Dependencies

All numerical data and plots are obtained using the [Quantum Optomechanics Toolbox](https://github.com/sampreet/qom), an open-source Python framework to simulate optomechanical systems.
Refer to the [QOM toolbox documentation](https://sampreet.github.io/qom-docs/v1.0.1) for the steps to install this libary.

## Running the Scripts

To run the scripts, navigate *inside* the top-level directory, and execute:

```bash
python scripts/bar/baz.py
```

Here, `bar` is the name of the folder (containing the version information) inside `scripts` and `baz.py` is the name of the script (refer to the repository structure).