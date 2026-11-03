#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

#_/_/_/_/_/_/_/_/_/_/_/_/
#Se
#_/_/_/_/_/_/_/_/_/_/_/_/

read -r west east south north dx dy < "Se_vG_gridinfo.txt"
region="${west}/${east}/${south}/${north}"
spacing="${dx}/${dy}"
gridfile="Se_vG.nc"

gmt begin "Se_vG" png

    gmt xyz2grd "Se_vG_xyz.dat" \
        -G"${gridfile}" \
        -R"${region}" \
        -I"${spacing}" \
        -r

    fine_dy=$(awk -v dy="$dy" 'BEGIN {printf "%.12g", dy/4}')
    finegrid="Se_vG_fine.nc"
    
    gmt grdsample "${gridfile}" \
        -G"${finegrid}" \
        -R"${region}" \
        -I"${dx}/${fine_dy}" \
        -r \
        -nb
    
    gmt grdimage "${finegrid}" \
        -R"${region}" \
        -JX14c/9c \
        -C"Se_vG.cpt" \
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
        -Byaf+l"Depth [-]" \
        #-BWSne+t"Effective saturation"

    #gmt plot "groundwater_level.dat" \
    #    -R"${region}" \
    #    -JX14c/9c \
    #    -W1p,yellow,-- \
    #    -l"Simulated groundwater level"

    #gmt legend -DjBL+o0.2c -F+p0.5p+gwhite

    gmt colorbar \
        -C"Se_vG.cpt" \
        -DJMR+w8c/0.25c+v+o0.8c/0c \
        -Baf
    
    # Add a centered title above the colorbar
    printf '14.925 8.75 S@-e@- \n' | gmt text \
        -R0/14/0/9 \
        -JX14c/9.2c \
        -F+f12p,Helvetica+jCM \
        -N

gmt end 
