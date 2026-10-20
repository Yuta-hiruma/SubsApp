#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${PROJECT_DIR}/output/gmt"

read -r west east south north dx dy < "${DATA_DIR}/Se_vG_gridinfo.txt"
region="${west}/${east}/${south}/${north}"
spacing="${dx}/${dy}"
gridfile="${DATA_DIR}/Se_vG.nc"

gmt begin "${DATA_DIR}/Se_vG" png

    gmt xyz2grd "${DATA_DIR}/Se_vG_xyz.dat" \
        -G"${gridfile}" \
        -R"${region}" \
        -I"${spacing}" \
        -r

    fine_dy=$(awk -v dy="$dy" 'BEGIN {printf "%.12g", dy/4}')
    finegrid="${DATA_DIR}/Se_vG_fine.nc"
    
    gmt grdsample "${gridfile}" \
        -G"${finegrid}" \
        -R"${region}" \
        -I"${dx}/${fine_dy}" \
        -r \
        -nb
    
    gmt grdimage "${finegrid}" \
        -R"${region}" \
        -JX14c/9c \
        -C"${DATA_DIR}/Se_vG.cpt" \
        -nb \
        -E50

    gmt basemap \
        -R"${region}" \
        -JX14c/9c \
        -Baf+l"Time [day]" \
        -Byaf+l"Depth [-]" \
        #-BWSne+t"Effective saturation"

    gmt plot "${DATA_DIR}/groundwater_level.dat" \
        -R"${region}" \
        -JX14c/9c \
        -W1p,yellow,-- \
        -l"Simulated groundwater level"

    gmt legend -DjBL+o0.2c -F+p0.5p+gwhite

    gmt colorbar \
        -C"${DATA_DIR}/Se_vG.cpt" \
        -DJMR+w11c/0.25c+v+o0.8c/0c \
        -Baf+l"Se [-]"

gmt end show