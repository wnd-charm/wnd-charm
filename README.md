# WND-CHARM
WND-CHARM is a multi-purpose image classifier that can be applied to a wide variety of image classification tasks without modifications or fine-tuning, and yet provides classification accuracy comparable to state-of-the-art task-specific image classifiers. WND-CHARM can extract up to ~3,000 generic image descriptors (features) including polynomial decompositions, high contrast features, pixel statistics, and textures. These features are derived from the raw image, transforms of the image, and compound transforms of the image (transforms of transforms). The features are filtered and weighted depending on their effectiveness in discriminating between a set of predefined image classes (the training set). These features are then used to classify test images based on their similarity to the training classes. This classifier was tested on a wide variety of imaging problems including biological and medical image classification using several imaging modalities, face recognition, and other pattern recognition tasks. WND-CHARM is an acronym that stands for *"Weighted Neighbor Distance using Compound Hierarchy of Algorithms Representing Morphology."*

This package contains two implementations both of which use common image transform and feature extraction code:

* A command-line program `wndchrm` (without the a) written in C++ that streamlines the WND-CHARM algorithm workflow. It reads images and their class membership from a directory hierarchy or text file, and outputs classifier statistics to an HTML report or STDOUT. To build from distribution tarball, use `./configure && make`. To build from a cloned repository use `./build.sh`. 
* A Python library `wndcharm` that provides an API to do many of the same things as wndchrm while providing the flexibility of a scripting language to perform low manipulation and visualization of pixel intensities, generated features and classification results. To build, use `python setup.py build`.

This research was supported entirely by the Intramural Research Program of the National Institutes of Health, National Institute on Aging, Ilya Goldberg, Investigator. Address: Laboratory of Genetics/NIA/NIH, 251 Bayview Blvd., Suite 100, Baltimore, MD, 21224, USA

----
#### A full description of the wndchrm utility can be found at:

[Shamir L, Orlov N, Eckley DM, Macura T, Johnston J, Goldberg IG. Wndchrm - an open source utility for biological image analysis](http://www.scfbm.org/content/3/1/13). BMC Source Code for Biology and Medicine. 3: 13, 2008. [PDF download](https://ome.irp.nia.nih.gov/wnd-charm/BMC-wndchrm-utility.pdf)

#### The wndchrm utility is an implementation of the WND-CHARM algorithm described here:

[Orlov N, Shamir L, Macura T, Johnston J, Eckley DM, Goldberg IG. WND-CHARM: Multi-purpose image classification using compound image transforms](https://ome.irp.nia.nih.gov/wnd-charm/PRL_2008.pdf). Pattern Recognition Letters. 29(11): 1684-93, 2008.

#### A review of techniques used in pattern-recognition/machine learning as applied to image analysis is here:

[Shamir L, Delaney JD, Orlov N, Eckley DM, Goldberg IG. Pattern Recognition Software and Techniques for Biological Image Analysis](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000974). PLoS Computational Biology. 6(11): e1000974, 2010

## Current Release

The current release of the wndchrm command-line utility is  version 1.60, available here: [wndchrm-1.60.977f5f6.tar.gz](https://github.com/wnd-charm/wnd-charm/files/565545/wndchrm-1.60.977f5f6.tar.gz)

The previous release, v1.52, is available here: [wndchrm-1.52.775.tar.gz](https://github.com/wnd-charm/wnd-charm/files/565579/wndchrm-1.52.775.tar.gz)

## Supported Platforms

WND-CHARM should compile and run on any POSIX-compliant operating system. It has been tested on Linux (Ubuntu 12.04 w/ GCC 4.6, CentOS 6.3 w/ GCC 4.4) and Mac OS X (<=10.7 w/ GCC 4.2, with experimental support for 10.9 Mavericks w/ `clang` compiler).

## Dependencies

Installation of WND-CHARM minimally requires a C++ compiler, LibTIFF and FFTW.

* C++ Compiler
    * Mac OS X: Install the command-line developer tools.
    * Ubuntu/Debian: `sudo apt-get install build-essential`
* [LibTIFF 3.x](http://www.libtiff.org):
    * CentOS/RedHat: `sudo yum install libtiff-devel`
    * Ubuntu/Debian: `sudo apt-get install libtiff4-dev`
    * Mac OS X: `brew install libtiff`
* [FFTW 3.x](http://www.fftw.org/download.html):
    * CentOS/RedHat: `sudo yum install fftw-static fftw-devel`
    * Ubuntu/Debian: `sudo apt-get install libfftw3-dev`
    * Mac OS X: `brew install fftw`
* Optional for dendrograms: [PHYLIP](http://evolution.genetics.washington.edu/phylip/install.html)
    * Some X11 libraries must be installed prior to compiling/installing PHYLIP (`make install` in the src dir.)
        * CentOS/RedHat: `sudo yum install libX11-devel libXt-devel libXaw-devel`
        * Ubuntu/Debian: `sudo apt-get install libX11-dev libxt-dev libxaw7-dev`

#### WND-CHARM Python API additional dependencies
The WND-CHARM Python API additionally requires the Python development package, SWIG, and the common Python 3rd-party packages `numpy` and `scipy`. Optionally, result visualization tools are enabled by installing the package `matplotlib`. To run the provided example scripts, the package `argparse` is required (included with Python 2.7+).

* Python utilities:
    * CentOS/RedHat: `sudo yum install python-devel swig`
    * Ubuntu/Debian: `sudo apt-get install python-dev swig`
    * Mac OS X: `brew install python swig`
* Python packages:
    * CentOS/RedHat: `sudo yum install numpy scipy python-matplotlib argparse`
    * Ubuntu/Debian: `sudo apt-get install python-numpy python-scipy python-matplotlib`
    * Pip: `pip install numpy scipy matplotlib argparse`

## Recommended Hardware

2 GB RAM (per core), 10 GB HD space, 2 GHZ CPU. Please be aware that this utility is very computationally intensive. Multiple instances of wndchrm can work concurrently on the same dataset on multi-core/multi-processor CPUs. Simply use the same command line to launch as many instances of wndchrm as there are CPU cores available.

## Test Images

Please also visit the [IICBU Biological Image Repository](https://ome.irp.nia.nih.gov/iicbu2008), which provides a benchmark for testing and comparing the performance of image analysis algorithms for biological imaging.
