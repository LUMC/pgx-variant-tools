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
        'biopython==1.69',
    ],
    tests_requires=['pytest'],
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
