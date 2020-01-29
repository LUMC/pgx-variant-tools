from setuptools import setup

setup(
    name="pgx-variant-tools",
    version="0.0.3",
    description="Tools for working with variants.",
    author="Guy Allard",
    author_email="wgallard@lumc.nl",
    url="https://github.com/LUMC/pgx-variant-tools",
    platforms=['any'],
    packages=["variant_tools"],
    install_requires=[
        'numpy==1.13.1',
        'biopython==1.69',
        'edlib==1.1.2.post2',
        'pyinterval==1.2.0'
    ],
    tests_requires=['pytest'],
    entry_points={
      "console_scripts": [
          "vcf2sequence = variant_tools.cli:vcf2sequence"
      ]
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
    ],
    keywords='bioinformatics'
)
