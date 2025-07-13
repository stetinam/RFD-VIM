# RFD-VIM (RFDiffusion Visual Input Manager)

Interactive PyMOL tool for visually selecting protein residues to freeze in RFDiffusion workflows.

## Features

- **Visual residue selection** - Click residues in PyMOL to set freeze status
- **Real-time feedback** - Instant visual updates with color coding
- **Flexible input** - Load from .sbatch files, saved configurations, or start empty
- **Export ready** - Generates CONTIGS and INPAINT_SEQ for RFDiffusion

! Exports only the frozen contigs, does not output segments RFDiffusion should generate, must be added by user later !

## Installation

### Requirements
- Python 3.6+
- PyMOL
- Standard Python libraries (sys, os, re, threading, time, collections)

### Setup
```bash
git clone https://github.com/stetinam/rfd-vim.git
cd rfd-vim
```

## Usage

### Quick Start
```bash
python rfd-vim.py
```

### Workflow
1. **Load PDB** → Enter path to any PDB file (enter this and everything into PyMOL directly!)
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

! You may select multiple residues with a PyMOL command, but it must be in following format !

`select sele, protein & (resi 5,10,53-52)`
- change the numbers to the positions you actually want to select and change freeze status

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
- **Ctrl+C** → Emergency exit in command line

## Example Walkthrough

### Step-by-Step Example with ASNase.pdb

After cloning the repository, you can try the example workflow:

1. **Start the program**:
   ```bash
   python rfd-vim.py
   ```

2. **Load the example protein structure**:
   - When prompted for PDB file path, type: `file ASNase.pdb` diretly into PyMOL
   - This will load the ASNase protein structure into PyMOL

3. **Load the example configuration**:
   - Choose option `2` (Load from saved file, may be a different number in other menus)
   - Type: `file asnase_example.txt`
   - This will load predefined freeze settings with:
     - Green sticks: Fully frozen residues (backbone + type)
     - Orange lines: Backbone-only frozen residues
     - Cyan cartoon: Flexible regions

4. **Interactive editing**:
   - Choose option `1` (Interactive editing)
   - Click on residues in PyMOL to select them
   - Type commands directly in PyMOL command line:
     - `bt` = Backbone + Type frozen (green sticks)
     - `b` = Backbone only frozen (orange lines)  
     - `n` = Not frozen (cyan cartoon)
     - `q` = Quit editing mode

5. **Save your modifications**:
   - Return to main menu (`q`) and choose option `3` (Save settings)
   - Type: `file asnase_my_example.txt`
   - This saves your customized CONTIGS and INPAINT_SEQ settings

### PyMOL Commands Reference
- **Menu navigation**: Type `1`, `2`, `3`, `4`, or `5` in PyMOL
- **File operations**: Type `file filename.txt` in PyMOL 
- **Residue editing**: Type `bt`, `b`, `n`, or `q` in PyMOL

## Example Files
- `ASNase.pdb` - Example protein structure
- `asnase_example.txt` - Example configuration file

## Notice
Made with assistance from Claude Code. Double check output before use.
This PyMOL tool could interfere with some native PyMOL functions.

## License
MIT License - see LICENSE file for details

## Contributing
Contributions welcome! Please open an issue or submit a pull request.
