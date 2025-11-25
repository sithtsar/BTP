#!/bin/bash
echo "=== SIMULATION STATUS ==="
echo ""
if ps aux | grep -v grep | grep test_benchmark > /dev/null; then
    echo "✅ Simulation RUNNING (PID: $(pgrep test_benchmark))"
else
    echo "❌ Simulation NOT running"
fi
echo ""
echo "Latest log (last 10 lines):"
tail -10 logs/all_re_cylinder_test.log
echo ""
echo "Output files:"
ls -lh output/channel_cylinder_re*.dat 2>/dev/null | awk '{print $9, $5}'
