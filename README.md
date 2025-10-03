# ChemComparer

Program for viewing and analyzing molecular structures.

## Functionality

### Important Notes

-This program has been tested on macOS Tahoe 26.0.1.

Some features do not work, including:

-Saving tab positions on exit

-Full save/load functionality

-Display of hidden hydrogen atoms

### Main Features
- Visualization of molecules in 2D and 3D formats
- Comparison of up to 3 molecules simultaneously
- Analysis of chemical properties
- Import/export in various formats

### Visualization
- **2D mode**: structural formulas
- **3D mode**: interactive 3D models
- **Display styles**: stick, ball, wire
- **Controls**: rotation, scaling

### Property Analysis
- Molecular weight
- LogP (partition coefficient)
- Number of atoms and bonds
- Aromatic atoms and cycles
- Hydrogen donors/acceptors
- Polar surface area

## Installation

### Required Libraries
```bash
pip install PyQt5 rdkit py3Dmol matplotlib pandas numpy
