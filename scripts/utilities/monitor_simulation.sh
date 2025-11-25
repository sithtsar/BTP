#!/bin/bash
echo "=== Cylinder Flow Simulation Monitor ==="
echo ""
echo "Log file: logs/all_re_cylinder_test.log"
echo ""
tail -15 logs/all_re_cylinder_test.log
echo ""
echo "=== Output Files Created ==="
ls -lh output/channel_cylinder_re*.dat 2>/dev/null | tail -10
echo ""
echo "Refresh: watch -n 10 ./scripts/monitor_simulation.sh"
