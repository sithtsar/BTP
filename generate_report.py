#!/usr/bin/env python3
"""
LBM Simulation Report Generator

Automatically generates professional presentations from LBM simulation data
using Marp (Markdown Presentation Ecosystem).

Usage:
    python generate_report.py                    # Generate PDF report
    python generate_report.py --format html     # Generate HTML presentation  
    python generate_report.py --format all      # Generate all formats
    python generate_report.py --help            # Show help
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
from typing import List, Optional

# Add the reports module to Python path
sys.path.insert(0, str(Path(__file__).parent / 'reports'))

try:
    from reports.generators import LBMDataAnalyzer, MarkdownReportGenerator, MarpInterface
except ImportError as e:
    print(f"Error importing report generators: {e}")
    print("Make sure all dependencies are installed with: uv sync")
    sys.exit(1)

def setup_argument_parser() -> argparse.ArgumentParser:
    """Set up command line argument parser."""
    parser = argparse.ArgumentParser(
        description='Generate LBM simulation report presentations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_report.py                      # Generate PDF report
  python generate_report.py --format html        # Generate HTML presentation
  python generate_report.py --format all         # Generate all formats
  python generate_report.py --custom             # Use custom template
  python generate_report.py --clean              # Clean previous outputs
        """
    )
    
    parser.add_argument(
        '--format', 
        choices=['pdf', 'html', 'pptx', 'images', 'all'], 
        default='pdf',
        help='Output format for the presentation (default: pdf)'
    )
    
    parser.add_argument(
        '--output-dir', 
        type=Path, 
        default=Path('./reports/generated'),
        help='Output directory for generated reports (default: ./reports/generated)'
    )
    
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=Path('./data'),
        help='Directory containing simulation data (default: ./data)'
    )
    
    parser.add_argument(
        '--template',
        type=str,
        default='monthly_template.md',
        help='Template file to use (default: monthly_template.md)'
    )
    
    parser.add_argument(
        '--custom',
        action='store_true',
        help='Use custom template (if available)'
    )
    
    parser.add_argument(
        '--clean',
        action='store_true',
        help='Clean previous output files before generating new ones'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without actually generating files'
    )
    
    return parser

def validate_environment(project_root: Path, verbose: bool = False) -> bool:
    """Validate that the environment is properly set up."""
    issues = []
    
    # Check for required directories
    required_dirs = [
        project_root / 'reports',
        project_root / 'reports' / 'generators',
        project_root / 'reports' / 'templates',
        project_root / 'reports' / 'templates' / 'styles'
    ]
    
    for dir_path in required_dirs:
        if not dir_path.exists():
            issues.append(f"Missing directory: {dir_path}")
    
    # Check for required files
    required_files = [
        project_root / 'package.json',
        project_root / 'reports' / 'templates' / 'styles' / 'scientific.css'
    ]
    
    for file_path in required_files:
        if not file_path.exists():
            issues.append(f"Missing file: {file_path}")
    
    # Check Python dependencies
    try:
        import pandas
        import numpy
        import jinja2
    except ImportError as e:
        issues.append(f"Missing Python dependency: {e}")
    
    if issues:
        print("Environment validation failed:")
        for issue in issues:
            print(f"  - {issue}")
        
        print("\nTo fix these issues:")
        print("  1. Run: uv sync  # Install Python dependencies")
        print("  2. Run: npm install  # Install Node.js dependencies")
        print("  3. Ensure all required files are present")
        
        return False
    
    if verbose:
        print("Environment validation passed")
    
    return True

def find_simulation_data(data_dir: Path, verbose: bool = False) -> bool:
    """Check if simulation data is available."""
    if not data_dir.exists():
        if verbose:
            print(f"Warning: Data directory not found: {data_dir}")
        return False
    
    # Look for common data patterns
    data_patterns = [
        "*_velocity_t*.csv",
        "multiphase_t*.csv",
        "*.png"
    ]
    
    found_data = False
    for pattern in data_patterns:
        files = list(data_dir.glob(pattern))
        if files:
            found_data = True
            if verbose:
                print(f"Found {len(files)} files matching {pattern}")
        elif verbose:
            print(f"Info: No files found matching {pattern}")
    
    # Also check project root for plots
    project_root = data_dir.parent
    plot_files = list(project_root.glob("*.png"))
    if plot_files:
        found_data = True
        if verbose:
            print(f"Found {len(plot_files)} plot files in project root")
    
    return found_data

def main():
    """Main function for report generation."""
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    project_root = Path.cwd()
    
    print("LBM Simulation Report Generator")
    print("=" * 50)
    
    # Validate environment
    if not validate_environment(project_root, args.verbose):
        sys.exit(1)
    
    # Check for simulation data
    if not args.data_dir.exists():
        # Fallback to project root
        args.data_dir = project_root
    
    has_data = find_simulation_data(args.data_dir, args.verbose)
    if not has_data:
        print("Warning: No simulation data found")
        print(f"   Looked in: {args.data_dir}")
        print("   The report will be generated with placeholder data")
    
    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.dry_run:
        print("\nDRY RUN - No files will be generated")
        print(f"  - Data directory: {args.data_dir}")
        print(f"  - Output directory: {args.output_dir}")
        print(f"  - Format: {args.format}")
        print(f"  - Template: {args.template}")
        return
    
    # Clean previous outputs if requested
    if args.clean:
        print("Cleaning previous outputs...")
        marp = MarpInterface(project_root)
        marp.clean_output(keep_markdown=False)
    
    try:
        # Step 1: Analyze simulation data
        print("\nAnalyzing simulation data...")
        analyzer = LBMDataAnalyzer(project_root)
        analysis = analyzer.generate_monthly_analysis()
        
        if args.verbose:
            print(f"  - Found {len(analysis.single_phase_metrics)} single-phase configurations")
            print(f"  - Multiphase data: {'Available' if analysis.multiphase_metrics else 'Not available'}")
            print(f"  - Generated plots: {len(analysis.generated_plots)}")
        
        # Step 2: Generate Markdown report
        print("Generating Markdown report...")
        reports_dir = project_root / 'reports'
        generator = MarkdownReportGenerator(reports_dir / 'templates')
        
        if args.custom and (reports_dir / 'templates' / 'custom_template.md').exists():
            markdown_content = generator.generate_custom_report(analysis, 'custom_template.md')
            if args.verbose:
                print("  - Using custom template")
        else:
            markdown_content = generator.generate_monthly_report(analysis)
            if args.verbose:
                print("  - Using default template")
        
        # Step 3: Save Markdown file
        timestamp = datetime.now().strftime('%Y-%m')
        markdown_file = reports_dir / 'output' / f'monthly_report_{timestamp}.md'
        markdown_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(markdown_file, 'w', encoding='utf-8') as f:
            f.write(markdown_content)
        
        print(f"Markdown report saved: {markdown_file}")
        
        # Step 4: Convert using Marp
        print("Converting to presentation format(s)...")
        marp = MarpInterface(project_root)
        
        # Show Marp info if verbose
        if args.verbose:
            marp_info = marp.get_marp_info()
            print(f"  - Marp version: {marp_info['version']}")
            print(f"  - Marp command: {marp_info['command']}")
        
        generated_files = []
        
        if args.format == 'pdf' or args.format == 'all':
            try:
                pdf_file = marp.convert_to_pdf(markdown_file)
                generated_files.append(('PDF', pdf_file))
                print(f"PDF report generated: {pdf_file}")
            except Exception as e:
                print(f"PDF generation failed: {e}")
        
        if args.format == 'html' or args.format == 'all':
            try:
                html_file = marp.convert_to_html(markdown_file)
                generated_files.append(('HTML', html_file))
                print(f"HTML presentation generated: {html_file}")
            except Exception as e:
                print(f"HTML generation failed: {e}")
        
        if args.format == 'pptx' or args.format == 'all':
            try:
                pptx_file = marp.convert_to_pptx(markdown_file)
                generated_files.append(('PowerPoint', pptx_file))
                print(f"PowerPoint presentation generated: {pptx_file}")
            except Exception as e:
                print(f"PowerPoint generation failed: {e}")
        
        if args.format == 'images':
            try:
                image_files = marp.convert_to_images(markdown_file)
                generated_files.append(('Images', f"{len(image_files)} slide images"))
                print(f"Slide images generated: {len(image_files)} files")
            except Exception as e:
                print(f"Image generation failed: {e}")
        
        # Summary
        print("\n" + "=" * 50)
        print("GENERATION SUMMARY")
        print("=" * 50)
        print(f"Output directory: {args.output_dir}")
        print(f"Markdown source: {markdown_file.name}")
        
        if generated_files:
            print("Generated presentations:")
            for format_name, file_info in generated_files:
                print(f"  - {format_name}: {file_info}")
        else:
            print("Warning: No presentation files were generated")
        
        print(f"\nGeneration completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Usage suggestions
        print("\nNext steps:")
        if any('PDF' in item[0] for item in generated_files):
            print("  - Open the PDF file for a print-ready presentation")
        if any('HTML' in item[0] for item in generated_files):
            print("  - Open the HTML file in a web browser for an interactive presentation")
        if any('PowerPoint' in item[0] for item in generated_files):
            print("  - Open the PowerPoint file for editing and customization")
        
        print("  - Run with --format all to generate all formats")
        print("  - Run with --clean to remove old files before generating new ones")
        
    except KeyboardInterrupt:
        print("\n\nGeneration cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during report generation: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()