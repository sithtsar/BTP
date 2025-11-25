#!/usr/bin/env python3
"""
OCR extraction for local PDF file using Mistral API
"""
import os
import json
import base64
from pathlib import Path
from dotenv import load_dotenv
from mistralai import Mistral

def save_image(base64_str: str, out_path: Path) -> None:
    data = base64.b64decode(base64_str)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as imgf:
        imgf.write(data)

def encode_pdf(pdf_path):
    """Encode PDF file to base64 string"""
    with open(pdf_path, "rb") as pdf_file:
        return base64.b64encode(pdf_file.read()).decode('utf-8')

def main():
    load_dotenv()
    api_key = os.getenv("MISTRAL_API_KEY")
    if not api_key:
        raise RuntimeError("MISTRAL_API_KEY not set in environment")

    client = Mistral(api_key=api_key)

    # Path to local PDF
    pdf_path = Path("docs/papers/ETH.pdf")
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF not found: {pdf_path}")

    print(f"Processing PDF: {pdf_path}")
    print(f"File size: {pdf_path.stat().st_size / 1024:.1f} KB")

    # Encode PDF using the helper function
    base64_pdf = encode_pdf(pdf_path)

    print("Sending to Mistral OCR API...")
    # Use document_url format with data URI as per Mistral documentation
    ocr_response = client.ocr.process(
        model="mistral-ocr-latest",
        document={
            "type": "document_url",
            "document_url": f"data:application/pdf;base64,{base64_pdf}"
        },
        include_image_base64=True
    )

    pages = ocr_response.pages
    print(f"Received {len(pages)} pages")

    out_dir = Path.cwd() / "docs" / "ocr_output" / "ETH_paper"
    out_dir.mkdir(parents=True, exist_ok=True)
    markdown_path = out_dir / "ETH_paper.md"
    json_path = out_dir / "ETH_paper_summary.json"
    images_dir = out_dir / "images"

    md_parts = []
    for page in pages:
        idx = getattr(page, "index", None)
        if idx is not None:
            md_parts.append(f"\n\n---\n\n<!-- Page {idx} -->\n\n")
        else:
            md_parts.append("\n\n---\n\n")

        md_text = getattr(page, "markdown", None)
        if md_text:
            md_parts.append(md_text)
        else:
            md_parts.append(getattr(page, "text", ""))

        images = getattr(page, "images", []) or []
        for img in images:
            img_id = getattr(img, "id", "unnamed")
            img_b64 = getattr(img, "image_base64", None)
            img_ext = "png"
            if isinstance(img_id, str) and img_id.lower().endswith((".jpg", ".jpeg")):
                img_ext = img_id.split(".")[-1]
            img_filename = f"page{idx}-{img_id}.{img_ext}" if idx is not None else f"{img_id}.{img_ext}"
            if img_b64:
                img_path = images_dir / img_filename
                save_image(img_b64, img_path)
                md_parts.append(f"\n\n![{img_id}](images/{img_filename})\n\n")
            else:
                img_url = getattr(img, "url", None) or getattr(img, "image_url", None)
                if img_url:
                    md_parts.append(f"\n\n![{img_id}]({img_url})\n\n")

    combined_md = "".join(md_parts).strip()
    with open(markdown_path, "w", encoding="utf-8") as f_md:
        f_md.write(combined_md)

    # Prepare safe serializable summary
    safe_pages = []
    for page in pages:
        safe_pages.append({
            "index": getattr(page, "index", None),
            "markdown_preview": (getattr(page, "markdown", "")[:200] + "...") if getattr(page, "markdown", None) else None,
            "num_images": len(getattr(page, "images", []) or [])
        })

    raw_summary = {
        "model": getattr(ocr_response, "model", None),
        "num_pages": len(pages),
        "pages": safe_pages
    }

    with open(json_path, "w", encoding="utf-8") as f_json:
        json.dump(raw_summary, f_json, ensure_ascii=False, indent=2)

    print(f"\n✓ OCR finished. Pages: {len(pages)}")
    print(f"✓ Markdown → {markdown_path}")
    print(f"✓ JSON summary → {json_path}")
    if images_dir.exists():
        print(f"✓ Images directory → {images_dir} ({len(list(images_dir.glob('*')))} files)")

if __name__ == "__main__":
    main()
