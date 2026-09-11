"""Locate Clair3 models that are compatible with the installed Clair3 version.

Clair3 v1.x (Linux, TensorFlow) reads TensorFlow checkpoints (``*.index``/``*.data-*``).
Clair3 v2.x (Linux and macOS, PyTorch) reads PyTorch checkpoints (``*.pt``).

Models can live in two places:
  * ``<dviont package>/models/<model>`` - filled by ``download_clair3_models`` (Clair3 v1 models)
  * ``<env>/bin/models/<model>``       - bundled with the Bioconda Clair3 package
"""

import glob
import logging
import os
import re
import shutil
import subprocess

PACKAGE_MODELS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "models")
)

REQUIRED_FORMAT = {1: "tensorflow", 2: "pytorch"}


def clair3_major_version():
    """Return the installed Clair3 major version (e.g. 1 or 2), or None if it cannot be determined."""
    exe = shutil.which("run_clair3.sh")
    if not exe:
        return None
    try:
        result = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"Clair3\s+v?(\d+)\.", f"{result.stdout}\n{result.stderr}")
    return int(match.group(1)) if match else None


def bundled_models_dirs():
    """Directories where the Bioconda Clair3 package ships its pre-trained models."""
    dirs = []
    exe = shutil.which("run_clair3.sh")
    if exe:
        dirs.append(os.path.join(os.path.dirname(exe), "models"))
        dirs.append(os.path.join(os.path.dirname(os.path.realpath(exe)), "models"))
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        dirs.append(os.path.join(conda_prefix, "bin", "models"))

    unique = []
    for d in dirs:
        d = os.path.abspath(d)
        if d not in unique:
            unique.append(d)
    return unique


def model_format(path):
    """Return 'pytorch', 'tensorflow', or None based on the checkpoint files in a model directory."""
    if glob.glob(os.path.join(path, "*.pt")):
        return "pytorch"
    if glob.glob(os.path.join(path, "*.index")):
        return "tensorflow"
    return None


def resolve_model_path(model_name, model_path=None):
    """
    Return the Clair3 model directory to use.

    An explicit ``model_path`` always wins. Otherwise the dviONT package models directory and the
    Clair3 Bioconda bundle are searched, keeping only models whose format matches the installed
    Clair3 version. Raises FileNotFoundError with guidance if no usable model is found.
    """
    major = clair3_major_version()
    required = REQUIRED_FORMAT.get(major)

    if model_path:
        model_path = os.path.abspath(model_path)
        if not os.path.isdir(model_path):
            raise FileNotFoundError(f"Clair3 model path does not exist: {model_path}")
        fmt = model_format(model_path)
        if required and fmt and fmt != required:
            logging.warning(
                "Model at %s looks like a %s model, but Clair3 v%s needs %s models.",
                model_path, fmt, major, required,
            )
        return model_path

    candidates = [os.path.join(PACKAGE_MODELS_DIR, model_name)]
    candidates += [os.path.join(d, model_name) for d in bundled_models_dirs()]

    for candidate in candidates:
        if not os.path.isdir(candidate):
            continue
        if required is None or model_format(candidate) == required:
            logging.info("Resolved Clair3 model '%s' to %s", model_name, candidate)
            return candidate
        logging.info(
            "Skipping %s: %s model is not compatible with Clair3 v%s",
            candidate, model_format(candidate) or "unrecognised", major,
        )

    searched = "\n  ".join(candidates)
    if major == 2:
        hint = (
            "Clair3 v2 needs PyTorch models. Pick a model bundled with Clair3 "
            "(ls \"$CONDA_PREFIX/bin/models\") or download one from "
            "https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/ and pass it with -p/--model-path."
        )
    else:
        hint = (
            f"Run 'download_clair3_models {model_name}' to download it, "
            "or pass an existing model directory with -p/--model-path."
        )
    raise FileNotFoundError(
        f"No compatible Clair3 model named '{model_name}' was found "
        f"(Clair3 v{major if major else '?'}). Searched:\n  {searched}\n{hint}"
    )
