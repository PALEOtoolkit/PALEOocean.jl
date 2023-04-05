# Setup PALEOocean/examples environment to use local downloaded PALEOocean package code
# This only needs to be run once, after cloning the github repository

import Pkg

Pkg.activate(".") # use the PALEOocean/examples environment
Pkg.develop(path="../")   # use the local version of PALEOocean packages to allow local modifications
Pkg.instantiate()