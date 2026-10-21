#!/usr/bin/env bash
set -euo pipefail

input="../k.dat"

points8=$(mktemp)
points9=$(mktemp)
hm8=$(mktemp)
hm9=$(mktemp)
lm8=$(mktemp)
lm9=$(mktemp)
trap 'rm -f "$points8" "$points9" "$hm8" "$hm9" "$lm8" "$lm9"' EXIT

input_se="../Se_vG.dat"
points8_se=$(mktemp)
points9_se=$(mktemp)
hm8_se=$(mktemp)
hm9_se=$(mktemp)
trap 'rm -f "$points8_se" "$points9_se" "$hm8_se" "$hm9_se"' EXIT

# Extract K and depth, ordered from shallow to deep.
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

# Calculate harmonic means at the interval midpoint depths.
for time in 8 9; do
    if [[ "$time" == 8 ]]; then points="$points8"; hm="$hm8"; else points="$points9"; hm="$hm9"; fi

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

# Calculate logarithmic means using (K2 - K1) / ln(K2 / K1).
for time in 8 9; do
    if [[ "$time" == 8 ]]; then points="$points8"; lm="$lm8"; else points="$points9"; lm="$lm9"; fi

    awk '
    function logmean(k1, k2) {
        if (k1 == k2) return k1
        return (k2 - k1) / log(k2 / k1)
    }
    NR == 1 { k1=$1; z1=$2; next }
    NR == 2 {
        k2=$1; z2=$2
        print logmean(k1, k2), (z1 + z2)/2
        next
    }
    NR == 3 {
        k3=$1; z3=$2
        print logmean(k2, k3), (z2 + z3)/2
    }' "$points" > "$lm"
done

# Set the K-axis range from both profiles, with 10% padding.
read xmin xmax < <(
    awk '
    NR == 1 || $1 < min { min = $1 }
    NR == 1 || $1 > max { max = $1 }
    END {
        pad = (max - min) * 0.10
        printf "%.10f %.10f\n", min - pad, max + pad
    }' "$points8" "$points9"
)

gmt begin K_depth_t8_t9 png
    gmt set FONT_LABEL 12p FONT_ANNOT_PRIMARY 10p

    gmt basemap \
        -R${xmin}/${xmax}/-0.40/-0.25 \
        -JX14c/9c \
        -Bxa+l"Hydraulic conductivity [m/day]" \
        -Bya0.025g0.025+l"Depth [m]" \
        -BWSne

    # t = 8 profile, harmonic means, and logarithmic means
    gmt plot "$points8" -W1.2p,royalblue
    gmt plot "$points8" -Sc0.18c -Groyalblue -W0.8p,royalblue -l"t = 8: K"
    gmt plot "$hm8" -Ss0.24c -Glightskyblue -W0.8p,black -l"t = 8: harmonic mean"
    gmt plot "$lm8" -St0.26c -Gcyan -W0.8p,black -l"t = 8: logarithmic mean"

    # t = 9 profile, harmonic means, and logarithmic means
    gmt plot "$points9" -W1.2p,firebrick
    gmt plot "$points9" -Sc0.18c -Gfirebrick -W0.8p,firebrick -l"t = 9: K"
    gmt plot "$hm9" -Ss0.24c -Gorange -W0.8p,black -l"t = 9: harmonic mean"
    gmt plot "$lm9" -St0.26c -Ggold -W0.8p,black -l"t = 9: logarithmic mean"

    gmt legend -DjTR+o0.2c -F+gwhite+p0.5p
gmt end