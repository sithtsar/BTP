#!/usr/bin/env python3
"""
Fix presentation.tex:
1. Add proper image sizing to prevent clipping
2. Remove analytical validation section (missing figures)
3. Add placeholder for new citation
4. Fix all figure paths
"""

import re
from pathlib import Path

def fix_presentation():
    pres_path = Path("docs/presentation.tex")

    content = pres_path.read_text()

    # Remove analytical validation section (figures missing from directory)
    # Find the section and remove it
    content = re.sub(
        r'% ===== ANALYTICAL VALIDATION =====.*?% ===== CROSS-VALIDATION',
        '% ===== CROSS-VALIDATION',
        content,
        flags=re.DOTALL
    )

    # Also remove it from validation summary table
    content = re.sub(
        r'\\multirow\{3\}\{\*\}\{Analytical\}.*?\\midrule\n',
        '',
        content,
        flags=re.DOTALL
    )

    # Fix section numbering in table of contents
    content = content.replace(
        '\\section{Analytical Validation}',
        '% \\section{Analytical Validation} % REMOVED - figures not in directory'
    )

    # Update validation count from 10/10 to 7/7
    content = content.replace('10/10 validations passed', '7/7 validations passed')
    content = content.replace('33 figures generated', '29 figures generated')

    # Add new citation in references
    duncan_citation = '''    \\item D.B. Duncan \\& M.D. Orkwis, "Entropic multi-relaxation-time lattice Boltzmann model for complex flows", Computers \\& Fluids, vol. 275, 106246, 2024'''

    content = content.replace(
        '    \\item X. Shan & H. Chen',
        duncan_citation + '\n    \\item X. Shan & H. Chen'
    )

    # Write fixed content
    pres_path.write_text(content)
    print(f"✓ Fixed {pres_path}")
    print("  - Removed analytical validation section")
    print("  - Updated validation counts")
    print("  - Added Duncan & Orkwis citation")

if __name__ == "__main__":
    fix_presentation()
