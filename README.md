WND-CHARM
=========

WND-CHARM is a multi-purpose image classifier that can be applied to a wide variety of image classification tasks without modifications or fine-tuning, and yet provides classification accuracy comparable to state-of-the-art task-specific image classifiers. Wndchrm can extract up to ~3,000 generic image descriptors (features) including polynomial decompositions, high contrast features, pixel statistics, and textures. These features are derived from the raw image, transforms of the image, and compound transforms of the image (transforms of transforms). The features are filtered and weighted depending on their effectiveness in discriminating between a set of predefined image classes (the training set). These features are then used to classify test images based on their similarity to the training classes. This classifier was tested on a wide variety of imaging problems including biological and medical image classification using several imaging modalities, face recognition, and other pattern recognition tasks. WND-CHARM is an acronym that stands for "Weighted Neighbor Distance using Compound Hierarchy of Algorithms Representing Morphology."

This package contains two implementations both of which use common image transform and feature extraction code:
* `wndchrm` - A command-line program written in C++ that streamlines the WND-CHARM algorithm workflow. It reads images and their class membership from a directory hierarchy or text file, and outputs classifier statistics to an HTML report or STDOUT. To build, use `./configure && make`.
* `pychrm` - A Python library that provides an API to do many of the same things as wndchrm while providing the flexibility of a scripting language to perform low manipulation and visualization of pixel intensities, generated features and classification results. To build, use `python setup.py build`.

This research was supported entirely by the Intramural Research Program of the National Institutes of Health, National Institute on Aging, Ilya Goldberg, Investigator. Address: Laboratory of Genetics/NIA/NIH, 251 Bayview Blvd., Suite 100, Baltimore, MD, 21224, USA

----

A full description of the wndchrm utility can be found at:

[Shamir L, Orlov N, Eckley DM, Macura T, Johnston J, Goldberg IG. Wndchrm - an open source utility for biological image analysis](http://www.scfbm.org/content/3/1/13). BMC Source Code for Biology and Medicine. 3: 13, 2008. [PDF download](http://ome.grc.nia.nih.gov/wnd-charm/BMC-wndchrm-utility.pdf)

The wndchrm utility is an implementation of the WND-CHARM algorithm described here:

[Orlov N, Shamir L, Macura T, Johnston J, Eckley DM, Goldberg IG. WND-CHARM: Multi-purpose image classification using compound image transforms](http://ome.grc.nia.nih.gov/wnd-charm/PRL_2008.pdf). Pattern Recognition Letters. 29(11): 1684-93, 2008.

Supported Platforms
-------------------
  * Wndchrm was tested using Linux (Ubuntu, CentOS) and MacOS X, but should also run with other Unix and Linux distributions.

Dependencies
------------
  * Installation of Wndchrm requires a C compiler, LibTIFF and FFTW
    * We use GCC 4.2 on MacOS, and 4.4/4.6 on Linux
  * [LibTIFF 3.x](http://www.libtiff.org):
    * CentOS/RedHat: `sudo yum install libtiff-devel`
    * Ubuntu/Debian: `sudo apt-get libtiff4-dev`
  * [FFTW 3.x](http://www.fftw.org/download.html):
    * CentOS/RedHat: `sudo yum install fftw-static fftw-devel`
    * Ubuntu/Debian: `sudo apt-get libfftw3-dev `
  * Optional for dendrograms: [PHYLIP](http://evolution.genetics.washington.edu/phylip/install.html)
      * Some X11 libraries must be installed prior to compiling/installing PHYLIP (`make install` in the src dir.)
        * CentOS/RedHat: `sudo yum install libX11-devel libXt-devel libXaw-devel`
        * Ubuntu/Debian: `sudo apt-get install libX11-dev libxt-dev libxaw7-dev`

Recommended Hardware
--------------------
2 GB RAM (per core), 10 GB HD space, 2 GHZ CPU. Please be aware that this utility is very computationally intensive. Multiple instances of wndchrm can work concurrently on the same dataset on multi-core/multi-processor CPUs. Simply use the same command line to launch as many instances of wndchrm as there are CPU cores available.

Test Images
-----------
Please also visit the [IICBU Biological Image Repository](http://ome.grc.nia.nih.gov/iicbu2008), which provides a benchmark for testing and comparing the performance of image analysis algorithms for biological imaging.
