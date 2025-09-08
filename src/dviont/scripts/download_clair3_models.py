#!/usr/bin/env python3

import sys
import argparse
import os
import requests
import tarfile

# Base URL for downloading Clair3 models
MODEL_BASE_URL = "https://www.bio8.cs.hku.hk/clair3/clair3_models/"

# Default Clair3 models to download
DEFAULT_MODELS = [
    "r1041_e82_400bps_sup_v500",
    "r1041_e82_400bps_sup_v420",
    "r941_prom_sup_g5014",
    "r1041_e82_400bps_sup_v430_bacteria_finetuned"
]

# Set default directory to store models (relative to script location)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # scripts/
DEFAULT_MODELS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "models"))  # dviont/models/


def get_parser():
    """Argument parser for selecting models and output directory."""
    parser = argparse.ArgumentParser(description="Download and extract Clair3 model files.")

    parser.add_argument(
        "model_names",
        nargs="*",
        default=DEFAULT_MODELS,
        help=f"Names of Clair3 models to download (default: {', '.join(DEFAULT_MODELS)})."
    )

    parser.add_argument(
        "-o", "--output_dir", default=DEFAULT_MODELS_DIR, help=f"Directory to save downloaded models (default: {DEFAULT_MODELS_DIR})"
    )

    return parser


def download_and_extract_model(model_name, output_dir):
    """Download and extract a Clair3 model tar.gz file to the specified output directory."""
    model_url = f"{MODEL_BASE_URL}{model_name}.tar.gz"
    model_path = os.path.join(output_dir)
    local_tar_path = os.path.join(output_dir, f"{model_name}.tar.gz")

    try:
        os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists
        sys.stderr.write(f"Downloading {model_name} from {model_url}...\n")

        response = requests.get(model_url, stream=True)
        response.raise_for_status()

        with open(local_tar_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        sys.stderr.write(f"Downloaded: {local_tar_path}\n")

        # Extract the archive
        os.makedirs(model_path, exist_ok=True)  # Ensure model directory exists
        with tarfile.open(local_tar_path, "r:gz") as tar:
            tar.extractall(model_path)

        os.remove(local_tar_path)  # Remove tar file after extraction
        sys.stderr.write(f"Extracted: {model_name} to {model_path}\n")

    except requests.RequestException as e:
        sys.stderr.write(f"Failed to download {model_name}: {e}\n")


def main():
    """Main function to handle downloading models."""
    args = get_parser().parse_args()
    output_dir = os.path.abspath(args.output_dir)  # Resolve absolute path

    for model_name in args.model_names:
        download_and_extract_model(model_name, output_dir)


if __name__ == "__main__":
    main()
