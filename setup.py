from setuptools import setup

setup(
    name="variant_tools",
    version="0.0.1",
    description="Tools for working with variants.",
    author="Guy Allard",
    author_email="wgallard AT lumc DOT nl",
    url="https://git.lumc.nl/PharmacogenomicsPipe/variant_tools.git",
    license="MIT",
    platforms=['any'],
    packages=["variant_tools"],
    install_requires=[
        'numpy==1.13.1',
        'biopython==1.69',
        'edlib==1.1.2.post2'
    ],
    tests_requires=['pytest'],
    entry_points={
      "console_scripts": [
          "vcf2sequence = variant_tools.cli:vcf2sequence"
      ]
    },
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'License :: MIT License',
    ],
    keywords = 'bioinformatics'
)
