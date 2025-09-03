"""
Marp CLI interface wrapper for converting Markdown to presentation formats.
Handles PDF, HTML, and PowerPoint generation with proper configuration.
"""

import subprocess
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import shutil

class MarpInterface:
    """Interface wrapper for Marp CLI operations."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.reports_dir = project_root / "reports"
        self.output_dir = self.reports_dir / "output"
        self.generated_dir = self.reports_dir / "generated"
        self.config_path = project_root / "marp.config.js"
        
        # Ensure directories exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.generated_dir.mkdir(parents=True, exist_ok=True)
        
        # Check for Marp CLI availability
        self._check_marp_availability()
    
    def _check_marp_availability(self) -> bool:
        """Check if Marp CLI is available and install if needed."""
        try:
            # Try npx first (preferred method)
            result = subprocess.run(
                ["npx", "@marp-team/marp-cli", "--version"], 
                cwd=self.project_root,
                capture_output=True, 
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Try global installation
        try:
            result = subprocess.run(
                ["marp", "--version"], 
                capture_output=True, 
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Install via npm if package.json exists
        if (self.project_root / "package.json").exists():
            print("Installing Marp CLI dependencies...")
            try:
                subprocess.run(
                    ["npm", "install"], 
                    cwd=self.project_root,
                    check=True,
                    timeout=60
                )
                return True
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
                print("Warning: Could not install Marp CLI via npm")
        
        print("Warning: Marp CLI not available. Install with: npm install -g @marp-team/marp-cli")
        return False
    
    def _get_marp_command(self) -> List[str]:
        """Get the appropriate Marp command based on availability."""
        # Try npx first
        try:
            subprocess.run(
                ["npx", "@marp-team/marp-cli", "--version"], 
                cwd=self.project_root,
                capture_output=True, 
                check=True,
                timeout=5
            )
            return ["npx", "@marp-team/marp-cli"]
        except:
            pass
        
        # Fall back to global installation
        return ["marp"]
    
    def convert_to_pdf(self, markdown_file: Path, output_file: Optional[Path] = None) -> Path:
        """Convert Markdown to PDF using Marp CLI."""
        if output_file is None:
            output_file = self.generated_dir / f"{markdown_file.stem}.pdf"
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = self._get_marp_command() + [
            str(markdown_file),
            "--pdf",
            "--allow-local-files",
            "--theme", str(self.reports_dir / "templates/styles/scientific.css"),
            "--output", str(output_file),
            "--no-config"
        ]
        
        try:
            result = subprocess.run(
                cmd, 
                cwd=self.project_root, 
                check=True,
                capture_output=True,
                text=True,
                timeout=60
            )
            print(f"PDF generated successfully: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to generate PDF: {e.stderr if e.stderr else str(e)}"
            print(error_msg)
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired:
            raise RuntimeError("PDF generation timed out")
    
    def convert_to_html(self, markdown_file: Path, output_file: Optional[Path] = None) -> Path:
        """Convert Markdown to HTML presentation."""
        if output_file is None:
            output_file = self.generated_dir / f"{markdown_file.stem}.html"
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = self._get_marp_command() + [
            str(markdown_file),
            "--html",
            "--allow-local-files",
            "--theme", str(self.reports_dir / "templates/styles/scientific.css"),
            "--output", str(output_file),
            "--no-config"
        ]
        
        try:
            result = subprocess.run(
                cmd, 
                cwd=self.project_root, 
                check=True,
                capture_output=True,
                text=True,
                timeout=60
            )
            print(f"HTML generated successfully: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to generate HTML: {e.stderr if e.stderr else str(e)}"
            print(error_msg)
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired:
            raise RuntimeError("HTML generation timed out")
    
    def convert_to_pptx(self, markdown_file: Path, output_file: Optional[Path] = None) -> Path:
        """Convert Markdown to PowerPoint presentation."""
        if output_file is None:
            output_file = self.generated_dir / f"{markdown_file.stem}.pptx"
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = self._get_marp_command() + [
            str(markdown_file),
            "--pptx",
            "--allow-local-files",
            "--theme", str(self.reports_dir / "templates/styles/scientific.css"),
            "--output", str(output_file)
        ]
        
        try:
            result = subprocess.run(
                cmd, 
                cwd=self.project_root, 
                check=True,
                capture_output=True,
                text=True,
                timeout=60
            )
            print(f"PowerPoint generated successfully: {output_file}")
            return output_file
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to generate PowerPoint: {e.stderr if e.stderr else str(e)}"
            print(error_msg)
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired:
            raise RuntimeError("PowerPoint generation timed out")
    
    def convert_to_images(self, markdown_file: Path, output_dir: Optional[Path] = None) -> List[Path]:
        """Convert Markdown to slide images (PNG)."""
        if output_dir is None:
            output_dir = self.generated_dir / f"{markdown_file.stem}_slides"
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = self._get_marp_command() + [
            "--images", "png",
            "--allow-local-files",
            "--theme", str(self.reports_dir / "templates/styles/scientific.css"),
            str(markdown_file),
            "--output", str(output_dir)
        ]
        
        try:
            result = subprocess.run(
                cmd, 
                cwd=self.project_root, 
                check=True,
                capture_output=True,
                text=True,
                timeout=60
            )
            
            # Find generated image files
            image_files = list(output_dir.glob("*.png"))
            image_files.sort()
            
            print(f"Images generated successfully: {len(image_files)} slides in {output_dir}")
            return image_files
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to generate images: {e.stderr if e.stderr else str(e)}"
            print(error_msg)
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired:
            raise RuntimeError("Image generation timed out")
    
    def convert_all_formats(self, markdown_file: Path) -> Dict[str, Path]:
        """Convert Markdown to all supported formats."""
        results = {}
        
        formats = {
            'pdf': self.convert_to_pdf,
            'html': self.convert_to_html,
            'pptx': self.convert_to_pptx
        }
        
        for format_name, convert_func in formats.items():
            try:
                result_file = convert_func(markdown_file)
                results[format_name] = result_file
            except Exception as e:
                print(f"Warning: Failed to generate {format_name}: {e}")
                results[format_name] = None
        
        return results
    
    def batch_convert(self, markdown_files: List[Path], formats: List[str] = None) -> Dict[str, Dict[str, Path]]:
        """Convert multiple Markdown files to specified formats."""
        if formats is None:
            formats = ['pdf', 'html', 'pptx']
        
        results = {}
        
        for markdown_file in markdown_files:
            file_results = {}
            
            if 'pdf' in formats:
                try:
                    file_results['pdf'] = self.convert_to_pdf(markdown_file)
                except Exception as e:
                    print(f"Warning: PDF conversion failed for {markdown_file}: {e}")
                    file_results['pdf'] = None
            
            if 'html' in formats:
                try:
                    file_results['html'] = self.convert_to_html(markdown_file)
                except Exception as e:
                    print(f"Warning: HTML conversion failed for {markdown_file}: {e}")
                    file_results['html'] = None
            
            if 'pptx' in formats:
                try:
                    file_results['pptx'] = self.convert_to_pptx(markdown_file)
                except Exception as e:
                    print(f"Warning: PowerPoint conversion failed for {markdown_file}: {e}")
                    file_results['pptx'] = None
            
            results[markdown_file.name] = file_results
        
        return results
    
    def get_marp_info(self) -> Dict[str, str]:
        """Get information about the Marp installation."""
        info = {}
        
        try:
            cmd = self._get_marp_command() + ["--version"]
            result = subprocess.run(
                cmd, 
                cwd=self.project_root,
                capture_output=True, 
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                info['version'] = result.stdout.strip()
            else:
                info['version'] = 'Unknown'
        except:
            info['version'] = 'Not available'
        
        info['command'] = ' '.join(self._get_marp_command())
        info['config_path'] = str(self.config_path) if self.config_path.exists() else 'Not found'
        info['output_dir'] = str(self.generated_dir)
        
        return info
    
    def clean_output(self, keep_markdown: bool = True) -> None:
        """Clean generated output files."""
        if self.generated_dir.exists():
            for file in self.generated_dir.iterdir():
                if file.is_file():
                    if keep_markdown and file.suffix == '.md':
                        continue
                    file.unlink()
                elif file.is_dir():
                    shutil.rmtree(file)
        
        print(f"Cleaned output directory: {self.generated_dir}")