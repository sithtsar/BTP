#!/bin/bash
echo "Generating cylinder flow comparison plots..."
cd "$(dirname "$0")/.."

for re in 10 100 1000; do
    if [ -f "output/channel_cylinder_re${re}_BGK_final.dat" ] && \
       [ -f "output/channel_cylinder_re${re}_ELBM_final.dat" ]; then
        echo "Plotting Re = $re..."
        uv run python plotting/plot_channel_with_cylinder.py $re
    else
        echo "Skipping Re = $re (data not ready)"
    fi
done

echo "Done! Check figures/channel_cylinder/"
ls -lh figures/channel_cylinder/
