#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

#_/_/_/_/_/_/_/_/_/_/_/_/
#theta
#_/_/_/_/_/_/_/_/_/_/_/_/
read -r west east south north dx dy < "theta_gridinfo.txt"
region="${west}/${east}/${south}/${north}"
spacing="${dx}/${dy}"
gridfile="theta.nc"

gmt begin "theta" png

    gmt xyz2grd "theta_xyz.dat" \
        -G"${gridfile}" \
        -R"${region}" \
        -I"${spacing}" \
        -r

    fine_dy=$(awk -v dy="$dy" 'BEGIN {printf "%.12g", dy/4}')
    finegrid="theta_fine.nc"
    
    gmt grdsample "${gridfile}" \
        -G"${finegrid}" \
        -R"${region}" \
        -I"${dx}/${fine_dy}" \
        -r \
        -nb
    
    gmt grdimage "${finegrid}" \
        -R"${region}" \
        -JX14c/9c \
        -C"theta.cpt" \
        -nb \
        -E200

    gmt grdcontour "${finegrid}" \
        -R"${region}" \
        -JX14c/9c \
        -C0.1 \
        -A0.2+f8p,Helvetica,black \
        -W0.5p,black
        
    gmt basemap \
        -R"${region}" \
        -JX14c/9c \
        -Baf+l"Time [day]" \
        -Byaf+l"Depth [m]" \
        #-BWSne+t"Effective saturation"

    gmt colorbar \
        -C"theta.cpt" \
        -DJMR+w8c/0.25c+v+o0.8c/0c \
        -Baf
    
    # Add a centered title above the colorbar
    printf '14.925 8.75 @~q@~ \n' | gmt text \
        -R0/14/0/9 \
        -JX14c/9.2c \
        -F+f12p,Helvetica+jCM \
        -N

gmt end

