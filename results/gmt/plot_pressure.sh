#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

#_/_/_/_/_/_/_/_/_/_/_/_/
#pressure
#_/_/_/_/_/_/_/_/_/_/_/_/
read -r west east south north dx dy < "pressure_gridinfo.txt"
region="${west}/${east}/${south}/${north}"
spacing="${dx}/${dy}"
gridfile="pressure.nc"

gmt begin "pressure" png

    gmt xyz2grd "pressure_xyz.dat" \
        -G"${gridfile}" \
        -R"${region}" \
        -I"${spacing}" \
        -r

    fine_dy=$(awk -v dy="$dy" 'BEGIN {printf "%.12g", dy/4}')
    finegrid="pressure_fine.nc"
    
    gmt grdsample "${gridfile}" \
        -G"${finegrid}" \
        -R"${region}" \
        -I"${dx}/${fine_dy}" \
        -r \
        -nb
    
    gmt grdimage "${finegrid}" \
        -R"${region}" \
        -JX14c/9c \
        -C"pressure.cpt" \
        -nb \
        -E200

    gmt basemap \
        -R"${region}" \
        -JX14c/9c \
        -Baf+l"Time [day]" \
        -Byaf+l"Depth [m]" \

    gmt colorbar \
        -C"pressure.cpt" \
        -DJMR+w7c/0.25c+v+o0.8c/-0.3c \
        -Baf

    printf '14.925 8.75 Pressure\n14.925 8.35 [kPa]\n' | gmt text \
        -R0/14/0/9 \
        -JX14c/9c \
        -F+f10p,Helvetica+jCM \
        -N

gmt end

