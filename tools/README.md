# Development Tools

This directory contains one-time or specialized development tools used during the project.

## OCR Extraction Tool

### `ocr_extraction.py`

**Purpose**: Extract and convert Keming Yu's 2021 ELBM thesis PDF into structured markdown using Mistral OCR API.

**Status**: Completed - OCR extraction already done. Output is in `docs/ocr_output/`.

**Requirements**:
- Python 3.13+
- Mistral AI API key
- Dependencies: `mistralai`, `python-dotenv`

**Usage**:
```bash
# Set up environment
export MISTRAL_API_KEY="your-api-key-here"

# Run extraction (from tools/ directory)
python ocr_extraction.py
```

**Output**:
- `docs/ocr_output/output.md` - Full thesis markdown
- `docs/ocr_output/output_summary.json` - Metadata summary
- `docs/ocr_output/images/` - Extracted figures

**Notes**:
- This was a one-time extraction task
- Original PDF: https://cdr.lib.unc.edu/downloads/1j92gh19k
- No need to re-run unless PDF is updated
- Kept for reference and documentation purposes
