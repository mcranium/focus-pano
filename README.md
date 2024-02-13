# focus-pano
Tools for processing macro photography images

This repository houses a command line tool for processing macro-photographic images. The emphasis here is on effective (scripted) processing of images that form focus stacks, Multi Light Image Collections (MLICs) and panorama tiles. Everything relies on open source image processing tools and a UNIX like operating system (e.g. Linux) or a compatilbility layer (e.g. WSL) is assumed.

## Requirements
Tested under Python 3.8.10
Python packages: `numpy`, `scikit-learn`, `cv2`
`enfuse` (a command line HDR and focus merging program)
[`focus-stack`](https://github.com/PetteriAimonen/focus-stack) (an advanced command line alignment and focus merging program)
Hugin (for panoramic stitching, bundles all necessary PanoTools programs)
[Relight](https://github.com/cnr-isti-vclab/relight) (optional, only needed for photometric stereo (of MLICs) for creating light point files through highlight detection on spheres)

## Distributing images into directories
The main paradigm of this tool is that images are organized in directories, reflecting their role.
1. all images forming a panorama are located in one directory
2. all images forming a focus stack are in one directory
3. all images forming a Multi Light Image Collection (MLIC) are in one directory

## One command line program multiple functions
Multiple commands are bundled within one command line program (a python file)
To make effective use of the program without pointing to the location of the `.py` file, it is possible to create a shell alias.
