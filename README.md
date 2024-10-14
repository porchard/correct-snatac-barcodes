# Correct barcode reads

This is a small program that attempts to correct barcodes to a barcode whitelist. If a barcode matches the whitelist, it is not changed. If it doesn't, it corrects the barcodes using an algorithm similar to that of cellranger's snATAC barcode correction algorithm.

## Installation

After you've installed the boost C++ libraries, you can install with:

make install PREFIX=. BOOST_ROOT=/path/to/boost_install
