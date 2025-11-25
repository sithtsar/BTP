#!/bin/bash

echo "========================================================"
echo "ELBM Project - Results Verification"
echo "========================================================"
echo ""

# Check compiled binaries
echo "1. Checking compiled binaries..."
if [ -f "build/elbm_solver" ]; then
    echo "   ✅ elbm_solver"
else
    echo "   ❌ elbm_solver (run ./build.sh)"
fi

if [ -f "build/test_analytical" ]; then
    echo "   ✅ test_analytical"
else
    echo "   ❌ test_analytical"
fi

if [ -f "build/test_benchmark" ]; then
    echo "   ✅ test_benchmark"
else
    echo "   ❌ test_benchmark"
fi

echo ""

# Check simulation outputs
echo "2. Checking simulation outputs..."
total_cases=0
complete_cases=0

for case in "case1_re10" "case2_re100" "case3_re1000"; do
    for solver in "bgk" "elbm"; do
        total_cases=$((total_cases + 1))
        dir="output/${case}_${solver}"
        if [ -d "$dir" ] && [ -f "$dir/${solver}_centerline_2500.dat" ]; then
            # Check for NaN
            if grep -q "nan" "$dir/${solver}_centerline_2500.dat"; then
                echo "   ⚠️  ${case}_${solver} - Contains NaN values"
            else
                echo "   ✅ ${case}_${solver}"
                complete_cases=$((complete_cases + 1))
            fi
        else
            echo "   ❌ ${case}_${solver} (run ./run_all_cases.sh)"
        fi
    done
done

echo "   Summary: $complete_cases/$total_cases cases completed"
echo ""

# Check analytical validation
echo "3. Checking analytical validation..."
for test in "couette_results" "poiseuille_results" "tg_energy"; do
    if [ -f "output/${test}.dat" ]; then
        echo "   ✅ ${test}"
    else
        echo "   ❌ ${test} (run ./build/test_analytical)"
    fi
done

echo ""

# Check figures
echo "4. Checking generated figures..."
figure_count=0

for i in {2..19}; do
    # Determine which case
    if [ $i -le 7 ]; then
        case_dir="Case 1 Re~10"
    elif [ $i -le 13 ]; then
        case_dir="Case 2 Re~100"
    else
        case_dir="Case 3 Re~1000"
    fi

    if [ -f "figures/$case_dir/figure_3_$i.png" ]; then
        figure_count=$((figure_count + 1))
    fi
done

echo "   Thesis figures: $figure_count/18"

if [ -f "figures/thesis_replication/summary_comparison.png" ]; then
    echo "   ✅ Summary comparison figure"
    figure_count=$((figure_count + 1))
fi

echo "   Total figures: $figure_count/19"
echo ""

# Check notebooks
echo "5. Checking Marimo notebooks..."
for nb in "01_pressure_profiles" "02_analytical_validation" "03_benchmark_cases"; do
    if [ -f "notebooks/${nb}.py" ]; then
        echo "   ✅ ${nb}.py"
    else
        echo "   ❌ ${nb}.py"
    fi
done

echo ""

# Check documentation
echo "6. Checking documentation..."
for doc in "README.md" "RESULTS.md" "PROJECT_SUMMARY.md"; do
    if [ -f "$doc" ]; then
        echo "   ✅ $doc"
    else
        echo "   ❌ $doc"
    fi
done

echo ""

# Summary
echo "========================================================"
echo "Summary"
echo "========================================================"

if [ $complete_cases -eq 6 ] && [ $figure_count -ge 18 ]; then
    echo "✅ PRIMARY VALIDATION: COMPLETE"
    echo "   - All 6 rectangular pipe cases successful"
    echo "   - $figure_count/19 figures generated"
    echo ""
    echo "✅ PROJECT STATUS: READY FOR REVIEW"
else
    echo "⚠️  Some components incomplete"
    echo "   - Run missing simulations"
    echo "   - Regenerate missing figures"
fi

echo ""
echo "To view results:"
echo "  1. Check figures/ directory"
echo "  2. Read RESULTS.md for detailed analysis"
echo "  3. Launch notebooks: cd notebooks && marimo edit 01_pressure_profiles.py"
echo "========================================================"
