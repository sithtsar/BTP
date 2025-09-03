# Marp Report Generation for LBM Simulations

## Quick Start

### Generate Report
```bash
python generate_report.py --format html
```

### View Report
Open: `reports/generated/monthly_report_2025-09.html`

## What's Included

**Fixed Issues:**
- ✓ Proper NaN handling in data processing
- ✓ Clean, readable plot scaling and formatting  
- ✓ Fixed slide layout issues in final slides
- ✓ Professional black text on white background
- ✓ Correct image paths and sizing

**Generated Content:**
- Executive summary with key metrics
- Single-phase BGK method comparison (RMSE: 2×10⁻⁶)
- Multiphase flow analysis with interface tracking
- Technical parameters and mathematical formulations
- Clean visualization plots with proper scaling

**Plot Files:**
- `velocity_analysis.png`: Single-phase flow analysis (4-panel comparison)
- `multiphase_analysis.png`: Multiphase flow evolution (6-panel analysis)

## Report Features

- Clean presentation format suitable for academic use
- Proper error metrics and convergence analysis
- Interface stability tracking for multiphase flow
- Mass conservation verification
- Mathematical equation rendering with MathJax

The system automatically analyzes your simulation data and creates professional presentation slides ready for research documentation or presentations.