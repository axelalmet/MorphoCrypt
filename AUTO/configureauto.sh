#!/bin/bash
# Script to configure AUTO automatically cause I'm sick of doing this manually every time.
cd ~/auto/07p/cmds
source auto.env.sh
cd ../
./configure
make
cd ~/Documents/MorphoelasticCrypt/AUTO
