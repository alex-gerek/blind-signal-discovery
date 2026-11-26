"""Collision-free storage for trusted local Workstation upload handling."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path, PureWindowsPath

from bsd.io_loader import InputValidationError, SUPPORTED_EXTENSIONS


def save_uploaded_bytes(original_name: str, payload: bytes, directory: Path) -> Path:
    """Keep the original basename under a newly allocated directory for each upload.

    Retaining the basename lets downstream filename/holdout checks inspect it.
    Neither two roles nor two sessions can overwrite a previous upload.
    """
    if not isinstance(original_name, str) or not original_name or "\x00" in original_name:
        raise InputValidationError("Upload must have a non-empty valid filename.")
    safe_name = PureWindowsPath(original_name).name
    if safe_name in ("", ".", "..") or ":" in safe_name:
        raise InputValidationError("Upload filename is invalid.")
    reserved_check = getattr(os.path, "isreserved", None)
    reserved = reserved_check(safe_name) if reserved_check else PureWindowsPath(safe_name).is_reserved()
    if reserved:
        raise InputValidationError("Upload filename is reserved or invalid on Windows.")
    if Path(safe_name).suffix.lower() not in SUPPORTED_EXTENSIONS:
        raise InputValidationError(f"Unsupported upload format: {Path(safe_name).suffix!r}.")
    if not isinstance(payload, bytes) or not payload:
        raise InputValidationError("Uploaded file is empty or is not a byte payload.")
    try:
        directory.mkdir(parents=True, exist_ok=True)
        upload_dir = Path(tempfile.mkdtemp(prefix="upload_", dir=directory))
        out_path = upload_dir / safe_name
        with out_path.open("xb") as stream:
            stream.write(payload)
    except OSError as exc:
        raise InputValidationError(f"Unable to store uploaded file: {exc}") from exc
    return out_path
