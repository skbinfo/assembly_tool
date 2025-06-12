from setuptools import setup, find_packages

setup(
    name="genome-annotation-pipeline",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        'biopython',
        'ncbi-genome-download',
        'tqdm',
        'termcolor',
        'pyyaml',
        'pandas',
    ],
    entry_points={
        'console_scripts': [
            'genome-annotation = annotation_pipeline.annotation:main',
        ],
    },
    author="Your Name",
    author_email="your@email.com",
    description="A pipeline for genomic assembly and annotation",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/yourusername/genome-annotation-pipeline",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
