#!/usr/bin/env bash
set -e

input="../Se_vG.dat"
points8=$(mktemp)
points9=$(mktemp)
hm8=$(mktemp)
hm9=$(mktemp)
trap 'rm -f "$points8" "$points9" "$hm8" "$hm9"' EXIT

# Extract Se and depth, then order the rows from shallow to deep.
for time in 8 9; do
    if [[ "$time" == 8 ]]; then points="$points8"; else points="$points9"; fi

    awk -v t="$time" '
    NR > 1 && $1 == t &&
    ($2 == -0.275 || $2 == -0.325 || $2 == -0.375) {
        print $3, $2
    }' "$input" | sort -k2,2nr > "$points"

    if [[ $(wc -l < "$points") -ne 3 ]]; then
        echo "Expected 3 matching rows for time $time in $input." >&2
        exit 1
    fi
done

# Calculate pairwise harmonic means and put them at interval midpoint depths.
for time in 8 9; do
    if [[ "$time" == 8 ]]; then
        points="$points8"
        hm="$hm8"
    else
        points="$points9"
        hm="$hm9"
    fi

    awk '
    NR == 1 { k1=$1; z1=$2; next }
    NR == 2 {
        k2=$1; z2=$2
        print 2/(1/k1 + 1/k2), (z1 + z2)/2
        next
    }
    NR == 3 {
        k3=$1; z3=$2
        print 2/(1/k2 + 1/k3), (z2 + z3)/2
    }' "$points" > "$hm"
done

# Set the Se-axis range from both profiles, with 10% padding.
read xmin xmax < <(
    awk '
    NR == 1 || $1 < min { min = $1 }
    NR == 1 || $1 > max { max = $1 }
    END {
        pad = (max - min) * 0.10
        printf "%.10f %.10f\n", min - pad, max + pad
    }' "$points8" "$points9"
)

gmt begin Se_depth_t8_t9 png
    gmt set FONT_LABEL 12p FONT_ANNOT_PRIMARY 10p

    gmt basemap \
        -R${xmin}/${xmax}/-0.40/-0.25 \
        -JX14c/9c \
        -Bxa+l"Effective saturation [-]" \
        -Bya0.025g0.025+l"Depth [m]" \
        -BWSne

    # t = 8 profile and harmonic means
    gmt plot "$points8" -W1.2p,royalblue
    gmt plot "$points8" -Sc0.18c -Groyalblue -W0.8p,royalblue -l"t = 8: K"
    gmt plot "$hm8" -Ss0.24c -Glightskyblue -W0.8p,black -l"t = 8: harmonic mean"

    # t = 9 profile and harmonic means
    gmt plot "$points9" -W1.2p,firebrick
    gmt plot "$points9" -Sc0.18c -Gfirebrick -W0.8p,firebrick -l"t = 9: K"
    gmt plot "$hm9" -Ss0.24c -Gorange -W0.8p,black -l"t = 9: harmonic mean"

    gmt legend -DjTR+o0.2c -F+gwhite+p0.5p
gmt end