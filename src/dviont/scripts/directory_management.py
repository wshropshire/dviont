import os
import sys
import logging

class PipelineManager:
    def __init__(self, output_dir, sample):
        """
        Initializes the PipelineManager with output directory and sample name.
        
        Args:
            output_dir (str): The path to the output directory.
            sample (str): The name of the sample.
        """
        self.output_dir = output_dir
        self.sample = sample
        self.log_file = None

    def setup_logging(self):
        """Initialize logging for the pipeline."""
        self.log_file = os.path.join(self.output_dir, f"{self.sample}_dviont_pipeline.log")
        logging.basicConfig(
            filename=self.log_file,
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",  # Use the default timestamp format provided by logging
            filemode="w",  # Overwrite any existing log file
        )
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # Modify the formatter
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s",
                              "%Y-%m-%d %H:%M:%S")
        console.setFormatter(formatter)
        logging.getLogger().addHandler(console)
        logging.info(f"dviONT pipeline initiated for sample: {self.sample}")  

    def create_output_directory(self):
        """
        Create the output directory and its 'reference' subdirectory.
        If the output directory exists, ask the user for confirmation to override or exit.

        Returns:
            str: Path to the 'reference' subdirectory.
        """
        # Check if the output directory exists
        if os.path.exists(self.output_dir):
            response = input(f"Warning: The output directory '{self.output_dir}' already exists. Do you want to continue? (yes/no): ").strip().lower()
            if response != "yes":
                print("Exiting script.")
                sys.exit(1)
        else:
            os.makedirs(self.output_dir)

        # Initialize logging immediately after creating the output directory
        self.setup_logging()
        
        # Create the "reference" subdirectory
        ref_dir = os.path.join(self.output_dir, "reference")
        os.makedirs(ref_dir, exist_ok=True)
        logging.info(f"Created output directory: {self.output_dir}")
        logging.info(f"Created 'reference' subdirectory at: {ref_dir}")
        
        return ref_dir

    def get_log_file(self):
        """Returns the log file path."""
        return self.log_file
