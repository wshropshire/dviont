from setuptools import setup, find_packages
import os
import subprocess
import sys

# Function to read version from version.py file
def get_version():
    version_file = os.path.join('scripts', 'version.py')
    with open(version_file) as f:
        exec(f.read())
    return locals()['__version__']

# Read dependencies from the YAML file in the build directory
def parse_requirements():
    with open(os.path.join('build', 'dviont.yaml'), 'r') as f:
        import yaml
        requirements = yaml.safe_load(f)
        return requirements.get('dependencies', [])

# Install any missing dependencies before proceeding
install_requirements()

# Proceed with setup
setup(
    name='dviont',
    version=get_version(),  # Automatically fetch the version from the version.py file
    packages=find_packages(where='scripts'),  # Find packages in the scripts directory
    install_requires=parse_requirements(),  # Parse dependencies from the YAML file
    author='William Shropshire',
    author_email='wcshropshire@mdanderson.org',
    description='dviONT (DNA Variant Identification using ONT) is a bacteria variant calling pipeline designed specifically for Q20+ Oxford Nanopore Technologies sequencing data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/wshropshire/dviont',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    python_requires='>=3.13',
    include_package_data=True,  # Include additional files such as data, models, etc.
    package_data={  # Include non-Python files (e.g., YAML, binaries, etc.)
        '': ['data/*', 'models/*', 'etc/*', 'build/*', 'binaries/*'],
    },
    scripts=[  # List of executable scripts in the bin directory
        'bin/dviont',
        'bin/download_clair3_models',
    ],
)
