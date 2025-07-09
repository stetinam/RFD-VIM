# RFDiffiusion Visual Input Manager

Interactive PyMOL tool for visually selecting protein residues to freeze in RFDiffusion workflows.

## Features

- **Universal PDB support** - Works with any protein structure
- **Visual residue selection** - Click residues in PyMOL to set freeze status
- **Real-time feedback** - Instant visual updates with color coding
- **Flexible input** - Load from .sbatch files, saved configurations, or start empty
- **Export ready** - Generates CONTIGS and INPAINT_SEQ for RFDiffusion

## Installation

### Requirements
- Python 3.6+
- PyMOL
- Standard Python libraries (sys, os, re, threading, time, collections)

### Setup
```bash
git clone https://github.com/stetinam/RFDiffiusion-Visual-Input-Manager.git
cd RFDiffiusion-Visual-Input-Manager
python visual_RFD_input_manager.py
```

## Usage

### Quick Start
```bash
python visual_RFD_input_manager.py
```

### Workflow
1. **Load PDB** → Enter path to any PDB file (or '5' to exit)
2. **Choose initial settings** → Load from .sbatch file, saved file, or start empty
3. **Interactive editing** → Click residues in PyMOL and set freeze status
4. **Save settings** → Export CONTIGS and INPAINT_SEQ

### PyMOL Commands
- **Menu options**: `1`, `2`, `3`, `4`, `5`
- **File operations**: `file "filename.txt"` (replace filename with actual path)
- **Residue editing**: `bt`, `b`, `n`, `q`

### Residue States
- **bt** = Backbone + Type frozen (green sticks) - Critical residues
- **b** = Backbone only frozen (orange lines) - Structurally important, sequence flexible
- **n** = Not frozen (cyan cartoon) - Flexible regions
- **q** = Exit editing mode

### Edit Workflow
1. Select residues in PyMOL (click single or multiple)
2. Type residue state in PyMOL (`bt`, `b`, or `n`)
3. Changes apply instantly
4. Type `q` to finish editing

## File Operations
- **Load**: .sbatch files (auto-extracts CONTIGS/INPAINT_SEQ) or .txt files
- **Save**: Automatically adds .txt extension, saves to current directory

## Output Format
```
CONTIGS="A2-15/A50-88/A131-147"
INPAINT_SEQ="A50-56/A59-59/A61-62"
```

- **CONTIGS** → All preserved residues (bt + b states)
- **INPAINT_SEQ** → Backbone-only preserved residues (b state only)

## Exit Options
- **'5'** → Exit at any menu
- **'q'** → Exit editing mode
- **Ctrl+C** → Emergency exit

## Example Files
- `ASNase.pdb` - Example protein structure
- `asnase_example.txt` - Example configuration file

## License
MIT License - see LICENSE file for details

## Contributing
Contributions welcome! Please open an issue or submit a pull request.
