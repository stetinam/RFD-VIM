#!/usr/bin/env python3
"""
Universal RFDiffusion Interactive Input Visualizer
Click residues in PyMOL and type commands in PyMOL command line.
Works with any PDB file and a RFDiffusion input file containing CONTIGS and INPAINT_SEQ.

This script allows you to visualize and edit residue states interactively in PyMOL.
You can select residues, set their states (frozen backbone, AA tyoe or not frozen), and generate input strings for RFDiffusion.

Usage: python universal_RFD_input.py
"""

import sys
import os
import re
import pymol
from pymol import cmd
import argparse
import threading
import time
from collections import defaultdict
from contextlib import contextmanager

class UniversalRFDiffusionVisualizer:
    def __init__(self):
        self.pdb_file = None
        self.current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Track residue states: 'BT' = backbone+type frozen, 'B' = backbone only, 'N' = not frozen
        self.residue_states = {}  # {(chain, resnum): state}
        self.protein_residues = set()  # All protein residues available
        self._pymol_initialized = False
        
        self.init_pymol()
        
    def init_pymol(self):
        """Initialize PyMOL session"""
        try:
            pymol.finish_launching()
            cmd.reinitialize()
            cmd.set("ray_opaque_background", "off")
            cmd.set("antialias", 2)
            cmd.set("mouse_selection_mode", 1)  # Enable atom selection
            cmd.set("selection_width", 4)  # Make selections easier to see
            cmd.set("auto_zoom", 0)  # Don't auto-zoom on selections
            
            # Set up custom PyMOL commands for menu navigation
            self.setup_pymol_commands()
            
            self._pymol_initialized = True
            print("PyMOL initialized successfully")
            
        except Exception as e:
            print(f"Error initializing PyMOL: {e}")
            sys.exit(1)
            
    def setup_pymol_commands(self):
        """Set up custom PyMOL commands for menu interaction"""
        # Create shared variables for communication between PyMOL and Python
        self.pymol_choice = None
        self.pymol_input_mode = None  # 'menu', 'editing', etc.
        
        # Define PyMOL command functions
        def menu_choice(choice):
            self.pymol_choice = str(choice)
            print(f"Menu choice {choice} selected in PyMOL")
            
        def residue_state(state):
            if self.pymol_input_mode == 'editing':
                self.pymol_choice = str(state).upper()
                print(f"Residue state {state} selected in PyMOL")
            else:
                print("Not in editing mode - press '1' to start editing mode")
                
        def done_editing():
            self.pymol_choice = 'done'
            print("Editing mode finished")
        
        # Unified function for file input (save/load paths)
        def file(filename):
            self.pymol_choice = str(filename)
            print(f"File path received: {filename}")
            
        # Register commands in PyMOL (both lowercase and uppercase)
        cmd.extend('menu', menu_choice)
        cmd.extend('bt', lambda: residue_state('BT'))
        cmd.extend('BT', lambda: residue_state('BT'))
        cmd.extend('b', lambda: residue_state('B'))  
        cmd.extend('B', lambda: residue_state('B'))
        cmd.extend('n', lambda: residue_state('N'))
        cmd.extend('N', lambda: residue_state('N'))
        cmd.extend('q', lambda: residue_state('Q'))
        cmd.extend('Q', lambda: residue_state('Q'))
        cmd.extend('done', done_editing)
        cmd.extend('file', file)
        
        # Add number shortcuts for menu
        for i in range(1, 6):
            cmd.extend(str(i), lambda x=i: menu_choice(x))
            
    def get_input(self, prompt, valid_choices=None, allow_string=True):
        """Get input from PyMOL command line or terminal"""
        if "file" in prompt.lower():
            print(f"{prompt}")
            print("Type 'file filename' in PyMOL (replace 'filename' with your actual file path)")
        else:
            print(f"{prompt}")
            print("Type in PyMOL.")
        
        # Reset PyMOL choice
        self.pymol_choice = None
        
        # Start a thread to wait for terminal input
        terminal_input = None
        stop_input_thread = threading.Event()
        
        def get_terminal_input():
            nonlocal terminal_input
            try:
                terminal_input = input().strip()
            except (EOFError, KeyboardInterrupt):
                pass
                
        input_thread = threading.Thread(target=get_terminal_input)
        input_thread.daemon = True
        input_thread.start()
        
        # Wait for either terminal input or PyMOL command
        while True:
            time.sleep(0.05)  # Reduced delay for better responsiveness
            
            # Check if PyMOL command was entered
            if self.pymol_choice is not None:
                choice = self.pymol_choice
                self.pymol_choice = None
                if valid_choices is None or choice in valid_choices or allow_string:
                    return choice
                else:
                    print(f"Invalid choice: {choice}. Valid options: {valid_choices}")
                
            # Check if terminal input was entered
            if terminal_input is not None:
                if valid_choices is None or terminal_input in valid_choices or allow_string:
                    return terminal_input
                else:
                    print(f"Invalid choice: {terminal_input}. Valid options: {valid_choices}")
                    # Restart input thread
                    terminal_input = None
                    if input_thread.is_alive():
                        input_thread.join(timeout=0.1)
                    input_thread = threading.Thread(target=get_terminal_input)
                    input_thread.daemon = True
                    input_thread.start()
                
            # Check if thread finished (user pressed Enter with no input)
            if not input_thread.is_alive() and terminal_input is None:
                # Restart input thread for empty input
                input_thread = threading.Thread(target=get_terminal_input)
                input_thread.daemon = True
                input_thread.start()
            
    def browse_pdb_file(self):
        """Allow user to browse and select a PDB file"""
        pdb_file = self.get_input("Enter PDB file path (or '5' to exit): ", allow_string=True).strip()
        if not pdb_file:
            return False
            
        # Check for exit option
        if pdb_file == '5':
            print("Goodbye!")
            cmd.quit()
            pymol.finish_launching()
            sys.exit(0)
            
        # Check if file exists
        if not os.path.isabs(pdb_file):
            # Try relative to current directory
            pdb_file = os.path.join(self.current_dir, pdb_file)
            
        if not os.path.exists(pdb_file):
            print(f"Error: PDB file '{pdb_file}' not found")
            return False
            
        # Validate file extension
        if not pdb_file.lower().endswith(('.pdb', '.pdb.gz', '.ent')):
            print(f"Warning: '{pdb_file}' may not be a valid PDB file")
            
        return self.load_pdb(pdb_file)
        
    def load_pdb(self, pdb_file):
        """Load PDB file and discover all protein residues"""
        print(f"Attempting to load PDB file: {pdb_file}")
        
        if not os.path.exists(pdb_file):
            print(f"Error: PDB file '{pdb_file}' not found")
            return False
            
        try:
            cmd.delete("all")  # Clear existing objects instead of reinitializing
            cmd.load(pdb_file, "protein")
            cmd.show("cartoon", "protein")
            cmd.color("cyan", "protein")
            cmd.zoom("protein")
            self.pdb_file = pdb_file
            
            # Discover all protein residues
            cmd.select("protein_residues", "protein and polymer")
            
            # Use PyMOL's built-in stored object
            from pymol import stored
            stored.residues = []
            
            # Clear any existing stored data
            if hasattr(stored, 'selected_residues'):
                stored.selected_residues = []
            
            cmd.iterate("protein_residues", "stored.residues.append((chain, resi))")
            
            for chain, resi in stored.residues:
                key = (chain, int(resi))
                self.protein_residues.add(key)
                if key not in self.residue_states:
                    self.residue_states[key] = 'N'  # Default to not frozen
                    
            print(f"Loaded PDB file: {pdb_file}")
            print(f"Found {len(self.protein_residues)} protein residues")
            
            # Color ligand purple if present
            try:
                cmd.select("ligand", "hetatm and not name HOH")
                if cmd.count_atoms("ligand") > 0:
                    cmd.show("sticks", "ligand")
                    cmd.set_color("ligand_purple", [0.6, 0.2, 0.8])
                    cmd.color("ligand_purple", "ligand")
                    print("Ligand shown as purple sticks")
                cmd.deselect()
            except Exception as ligand_error:
                print(f"Note: Could not process ligand: {ligand_error}")
                
            return True
        except Exception as e:
            print(f"Error loading PDB: {e}")
            import traceback
            traceback.print_exc()
            return False
            
    def load_from_script(self, script_file):
        """Load CONTIGS and INPAINT_SEQ from script file"""
        if not os.path.exists(script_file):
            print(f"Error: Script file '{script_file}' not found")
            return False
            
        try:
            with open(script_file, 'r') as f:
                content = f.read()
                
            # Extract and parse CONTIGS
            contigs_match = re.search(r'CONTIGS="([^"]*)"', content)
            inpaint_match = re.search(r'INPAINT_SEQ="([^"]*)"', content)
            
            if contigs_match and inpaint_match:
                self.parse_and_set_states(contigs_match.group(1), inpaint_match.group(1))
                self.visualize_current_states()
                print(f"Loaded settings from {os.path.basename(script_file)}")
                return True
            else:
                print("Warning: Could not find CONTIGS and/or INPAINT_SEQ in script")
                return False
                
        except Exception as e:
            print(f"Error loading script: {e}")
            return False
            
    def load_from_saved_file(self, save_file):
        """Load CONTIGS and INPAINT_SEQ from a previously saved file"""
        # Check for absolute path first
        if os.path.isabs(save_file):
            file_path = save_file
        else:
            # Try current directory first
            current_path = os.path.join(self.current_dir, save_file)
            
            if os.path.exists(current_path):
                file_path = current_path
            else:
                print(f"Error: Save file '{save_file}' not found")
                return False
        
        if not os.path.exists(file_path):
            print(f"Error: Save file '{file_path}' not found")
            return False
        
        save_file = file_path
            
        try:
            with open(save_file, 'r') as f:
                lines = f.readlines()
                
            contigs = None
            inpaint_seq = None
            
            # Try both regex and direct parsing for robustness
            for line in lines:
                line = line.strip()
                # Try regex match first
                contigs_match = re.search(r'CONTIGS="([^"]*)"', line)
                inpaint_match = re.search(r'INPAINT_SEQ="([^"]*)"', line)
                
                if contigs_match:
                    contigs = contigs_match.group(1)
                if inpaint_match:
                    inpaint_seq = inpaint_match.group(1)
                
                # If regex fails, try direct parsing
                if contigs is None and line.startswith('CONTIGS='):
                    try:
                        parts = line.split('=', 1)
                        if len(parts) > 1:
                            val = parts[1].strip()
                            if val.startswith('"') and val.endswith('"'):
                                contigs = val[1:-1]
                    except:
                        pass
                        
                if inpaint_seq is None and line.startswith('INPAINT_SEQ='):
                    try:
                        parts = line.split('=', 1)
                        if len(parts) > 1:
                            val = parts[1].strip()
                            if val.startswith('"') and val.endswith('"'):
                                inpaint_seq = val[1:-1]
                    except:
                        pass
                    
            if contigs is not None and inpaint_seq is not None:
                self.parse_and_set_states(contigs, inpaint_seq)
                self.visualize_current_states()
                print(f"Loaded settings from {os.path.basename(save_file)}")
                return True
            else:
                print(f"Warning: Could not find valid CONTIGS and INPAINT_SEQ in save file: {save_file}")
                return False
                
        except Exception as e:
            print(f"Error loading save file: {e}")
            import traceback
            traceback.print_exc()
            return False
            
    def parse_and_set_states(self, contigs_str, inpaint_str):
        """Parse CONTIGS and INPAINT_SEQ strings and set residue states"""
        # Reset all states
        for key in self.residue_states:
            self.residue_states[key] = 'N'
            
        # Parse CONTIGS - residues that should be kept
        contigs_residues = set()
        parts = contigs_str.split('/')
        for part in parts:
            part = part.strip()
            if part and '-' in part:
                # Handle chain letter + range (e.g., A2-15)
                chain_match = re.match(r'([A-Z])(\d+)-(\d+)', part)
                if chain_match:
                    chain, start, end = chain_match.groups()
                    for res in range(int(start), int(end) + 1):
                        contigs_residues.add((chain, res))
                        
        # Parse INPAINT_SEQ - residues whose sequence can change
        inpaint_residues = set()
        parts = inpaint_str.split('/')
        for part in parts:
            part = part.strip()
            if part and '-' in part:
                # Handle chain letter + range (e.g., A2-15)
                chain_match = re.match(r'([A-Z])(\d+)-(\d+)', part)
                if chain_match:
                    chain, start, end = chain_match.groups()
                    for res in range(int(start), int(end) + 1):
                        inpaint_residues.add((chain, res))
                    
        # Set states based on membership
        for chain_res in contigs_residues:
            if chain_res in self.protein_residues:
                if chain_res in inpaint_residues:
                    self.residue_states[chain_res] = 'B'  # Backbone only
                else:
                    self.residue_states[chain_res] = 'BT'  # Backbone + type
                    
        print(f"Set {len([s for s in self.residue_states.values() if s == 'BT'])} fully frozen residues")
        print(f"Set {len([s for s in self.residue_states.values() if s == 'B'])} backbone-only frozen residues")
        
    def visualize_current_states(self):
        """Visualize current residue states in PyMOL"""
        # Clear previous protein visualizations only
        cmd.hide("everything", "protein")
        cmd.show("cartoon", "protein")
        cmd.color("cyan", "protein")
        
        # Group residues by state
        bt_residues = []  # Backbone + type frozen
        b_residues = []   # Backbone only frozen
        
        for (chain, resnum), state in self.residue_states.items():
            if state == 'BT':
                bt_residues.append(f"chain {chain} and resi {resnum}")
            elif state == 'B':
                b_residues.append(f"chain {chain} and resi {resnum}")
                
        # Visualize fully frozen (BT) as green sticks
        if bt_residues:
            cmd.select("frozen_bt", " or ".join(bt_residues))
            cmd.show("sticks", "frozen_bt")
            cmd.color("green", "frozen_bt")
            
        # Visualize backbone frozen (B) as orange lines
        if b_residues:
            cmd.select("frozen_b", " or ".join(b_residues))
            cmd.show("lines", "frozen_b")
            cmd.set_color("frozen_orange", [1.0, 0.5, 0.0])
            cmd.color("frozen_orange", "frozen_b")
            
        # Ensure ligand remains visible
        try:
            cmd.select("ligand", "hetatm and not name HOH")
            if cmd.count_atoms("ligand") > 0:
                cmd.show("sticks", "ligand")
                cmd.set_color("ligand_purple", [0.6, 0.2, 0.8])
                cmd.color("ligand_purple", "ligand")
        except Exception as e:
            print(f"Warning: Error handling ligand: {e}")
            
        cmd.deselect()
        
        print(f"\nVisualized: {len(bt_residues)} fully frozen (green sticks), {len(b_residues)} backbone-only (orange lines)")
        
    def reset_selection_state(self):
        """Complete reset of selection state after each change"""
        try:
            # Clear all PyMOL selections
            cmd.deselect()
            
            # Clear any named selections that might exist
            selections = cmd.get_names("selections")
            for sel in selections:
                cmd.delete(sel)
            
            # Reset internal choice tracking
            self.pymol_choice = None
            
            # Make sure any stored variables are cleared
            from pymol import stored
            stored.selected_residues = []
            
            # Small delay to ensure PyMOL processes the clearing
            time.sleep(0.02)
            
            # Additional step: attempt to clear the selection via mouse mode
            cmd.mouse('three_button_viewing')
            
        except Exception as e:
            print(f"Warning: Could not fully reset selection state: {e}")
            # Fallback: just clear basic selection
            cmd.deselect()
            self.pymol_choice = None
        
    def start_interactive_editing(self):
        """Start interactive editing mode"""
        print("\n" + "="*60)
        print("INTERACTIVE EDITING MODE")
        print("="*60)
        print("1. Click one or multiple residues in PyMOL (or hold shift for large selection)")
        print("2. Type choice in PyMOL:")
        print("   BT = Backbone + Type frozen (green sticks)")
        print("   B  = Backbone only frozen (orange lines)")
        print("   N  = Not frozen (cyan cartoon)")
        print("   Q  = Quit editing mode")
        print("\nType 'q' in PyMOL to finish")
        
        self.visualize_current_states()
        self.pymol_input_mode = 'editing'
        
        # Complete initial reset
        self.reset_selection_state()
        
        # Track processing state to avoid re-processing same selection
        last_processed_selection = None
        
        print("\nWaiting for residue selection...")
        
        try:
            while True:
                time.sleep(0.05)  # Reduced delay for better responsiveness
                
                # Check if user typed 'done' or 'q' in PyMOL
                if self.pymol_choice:
                    if self.pymol_choice.lower() in ['done', 'q']:
                        print("Exiting editing mode...")
                        return
                    self.pymol_choice = None
                
                # Check for current selection
                selection = cmd.get_names("selections")
                if selection and "sele" in selection:
                    # Get current selection signature to avoid re-processing
                    current_selection_signature = cmd.count_atoms("sele")
                    
                    # Only process if this is a new selection
                    if current_selection_signature != last_processed_selection:
                        last_processed_selection = current_selection_signature
                        
                        # Get all selected residues (supports multiple selection)
                        from pymol import stored
                        stored.selected_residues = []
                        try:
                            # Use iterate_state to get all selected atoms
                            cmd.iterate("sele", "stored.selected_residues.append((chain, resi))")
                            if stored.selected_residues:
                                # Remove duplicates while preserving order
                                unique_residues = []
                                seen = set()
                                for res in stored.selected_residues:
                                    if res not in seen and res[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:  # Only protein chains
                                        try:
                                            chain, resnum = res
                                            # Ensure resnum is an integer
                                            resnum = int(resnum)
                                            unique_residues.append((chain, resnum))
                                            seen.add((chain, resnum))
                                        except (ValueError, TypeError):
                                            # Skip invalid residue numbers
                                            continue
                                
                                if not unique_residues:
                                    print("No valid protein residues selected")
                                    self.reset_selection_state()
                                    last_processed_selection = None
                                    continue
                                
                                if len(unique_residues) == 1:
                                    chain, resnum = unique_residues[0]
                                    current_state = self.residue_states.get((chain, resnum), 'N')
                                    print(f"\nSelected: Chain {chain}, Residue {resnum}")
                                    print(f"Current status: {self.get_state_description(current_state)}")
                                else:
                                    print(f"\nSelected {len(unique_residues)} residues:")
                                    # Count residues by state
                                    state_counts = {'BT': 0, 'B': 0, 'N': 0}
                                    for chain, resnum in unique_residues:
                                        state = self.residue_states.get((chain, resnum), 'N')
                                        state_counts[state] += 1
                                        
                                    # Print first few and summary
                                    print(f"  First few: ", end="")
                                    for i, (chain, resnum) in enumerate(unique_residues[:3]):
                                        if i > 0:
                                            print(", ", end="")
                                        print(f"Chain {chain}, Res {resnum}", end="")
                                    
                                    if len(unique_residues) > 3:
                                        print(f", ... and {len(unique_residues) - 3} more")
                                    else:
                                        print("")
                                        
                                    print(f"  Current states: {state_counts['BT']} fully frozen, {state_counts['B']} backbone-only, {state_counts['N']} not frozen")
                                
                                # Wait for state choice
                                print("Waiting for choice (BT/B/N/Q)...")
                                choice_made = False
                                while not choice_made:
                                    time.sleep(0.05)
                                    
                                    # Check for PyMOL command
                                    if self.pymol_choice:
                                        new_state = self.pymol_choice.upper()
                                        self.pymol_choice = None
                                        
                                        if new_state in ['BT', 'B', 'N']:
                                            # Get selection right before applying change
                                            from pymol import stored
                                            stored.selected_residues = []
                                            cmd.iterate("sele", "stored.selected_residues.append((chain, resi))")
                                            
                                            unique_residues = []
                                            seen = set()
                                            for res in stored.selected_residues:
                                                if res not in seen and res[0] in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
                                                    try:
                                                        chain, resnum = res
                                                        resnum = int(resnum)
                                                        unique_residues.append((chain, resnum))
                                                        seen.add((chain, resnum))
                                                    except (ValueError, TypeError):
                                                        continue
                                            
                                            if not unique_residues:
                                                print("No valid protein residues selected.")
                                                self.reset_selection_state()
                                                last_processed_selection = None
                                                choice_made = False
                                                continue
                                            
                                            # Apply to all selected residues
                                            updated_count = 0
                                            for chain, resnum in unique_residues:
                                                self.residue_states[(chain, resnum)] = new_state
                                                updated_count += 1
                                            
                                            print(f"Updated {updated_count} residue(s) to: {self.get_state_description(new_state)}")
                                            
                                            # Immediately update visualization
                                            self.visualize_current_states()
                                            choice_made = True
                                            
                                            # COMPLETE RESET after change
                                            self.reset_selection_state()
                                            last_processed_selection = None
                                            print("\n" + "="*40)
                                            print("CHANGES APPLIED - READY FOR NEXT SELECTION")
                                            print("="*40)
                                            print("Select residue(s) and type choice, or type 'q' to finish...")
                                            
                                            # Force a redraw to show the updated visualization
                                            cmd.refresh()
                                            
                                        elif new_state in ['Q', 'DONE']:
                                            print("Exiting editing mode...")
                                            return
                                        else:
                                            print("Invalid option. Use BT, B, N, or Q")
                                            
                        except Exception as e:
                            print(f"Error reading selected residue: {e}")
                            self.reset_selection_state()
                            last_processed_selection = None
                        
        except KeyboardInterrupt:
            print("\nExiting editing mode...")
            return
                
    def get_state_description(self, state):
        """Get human-readable description of state"""
        descriptions = {
            'BT': 'Backbone + Type frozen (green sticks)',
            'B': 'Backbone only frozen (orange lines)', 
            'N': 'Not frozen (cyan cartoon)'
        }
        return descriptions.get(state, 'Unknown')
        
    def generate_contigs_and_inpaint(self):
        """Generate CONTIGS and INPAINT_SEQ strings from current states"""
        # Get all residues that should be in CONTIGS (BT or B states)
        contigs_residues = []
        inpaint_residues = []
        
        for (chain, resnum), state in self.residue_states.items():
            if state in ['BT', 'B']:
                contigs_residues.append((chain, resnum))
                if state == 'B':
                    inpaint_residues.append((chain, resnum))
                    
        # Sort residues
        contigs_residues.sort()
        inpaint_residues.sort()
        
        # Generate CONTIGS string
        contigs_parts = []
        if contigs_residues:
            # Group consecutive residues by chain
            chain_groups = defaultdict(list)
            for chain, resnum in contigs_residues:
                chain_groups[chain].append(resnum)
                
            for chain in sorted(chain_groups.keys()):
                residues = sorted(chain_groups[chain])
                ranges = self.group_consecutive(residues)
                for start, end in ranges:
                    if start == end:
                        contigs_parts.append(f"{chain}{start}-{start}")
                    else:
                        contigs_parts.append(f"{chain}{start}-{end}")
                        
        # Generate INPAINT_SEQ string
        inpaint_parts = []
        if inpaint_residues:
            chain_groups = defaultdict(list)
            for chain, resnum in inpaint_residues:
                chain_groups[chain].append(resnum)
                
            for chain in sorted(chain_groups.keys()):
                residues = sorted(chain_groups[chain])
                ranges = self.group_consecutive(residues)
                for start, end in ranges:
                    if start == end:
                        inpaint_parts.append(f"{chain}{start}-{start}")
                    else:
                        inpaint_parts.append(f"{chain}{start}-{end}")
                        
        contigs_str = "/".join(contigs_parts) if contigs_parts else ""
        inpaint_str = "/".join(inpaint_parts) if inpaint_parts else ""
        
        return contigs_str, inpaint_str
        
    def group_consecutive(self, numbers):
        """Group consecutive numbers into ranges"""
        if not numbers:
            return []
            
        ranges = []
        start = numbers[0]
        end = numbers[0]
        
        for i in range(1, len(numbers)):
            if numbers[i] == end + 1:
                end = numbers[i]
            else:
                ranges.append((start, end))
                start = numbers[i]
                end = numbers[i]
                
        ranges.append((start, end))
        return ranges
        
    def save_settings(self, filename):
        """Save current CONTIGS and INPAINT_SEQ to file"""
        contigs_str, inpaint_str = self.generate_contigs_and_inpaint()
        
        # Ensure we have a proper path (relative to current dir if not absolute)
        if not os.path.isabs(filename):
            filename = os.path.join(self.current_dir, filename)
            
        try:
            # Make sure directory exists
            dir_path = os.path.dirname(os.path.abspath(filename))
            if dir_path:  # Only create if directory path exists
                os.makedirs(dir_path, exist_ok=True)
            
            # Use context manager for safe file writing
            try:
                with open(filename, 'w') as f:
                    f.write(f'CONTIGS="{contigs_str}"\n')
                    f.write(f'INPAINT_SEQ="{inpaint_str}"\n')
                print(f"Settings saved to {filename}")
                return True
            except IOError as io_error:
                print(f"I/O error writing to {filename}: {io_error}")
                return False
        except Exception as e:
            print(f"Error saving to {filename}: {e}")
            import traceback
            traceback.print_exc()
            return False
            
    def main_menu(self):
        """Main interactive menu"""
        print("\n" + "="*60)
        print("Universal RFDiffusion Interactive Visualizer")
        print("="*60)
        print(f"Working directory: {self.current_dir}")
        
        # Initialize PyMOL graphics first by doing a dummy operation
        print("\nInitializing PyMOL graphics...")
        try:
            cmd.reinitialize()
            cmd.set("ray_opaque_background", "off")
            cmd.set("antialias", 2)
            time.sleep(0.5)  # Allow graphics to initialize
        except Exception as e:
            print(f"Warning: Graphics initialization issue: {e}")
        
        # Now show the PDB file prompt after graphics initialization
        print("\nStep 1: Load PDB file")
        print("Options:")
        print("  - Enter PDB file path")
        print("  - Type '5' to exit")
        if not self.browse_pdb_file():
            print("Failed to load PDB file. Exiting.")
            return
        
        # Clear any PyMOL messages and show clean menu
        print("\n" + "="*60)
        print("SETUP COMPLETE")
        print("="*60)
        print("\nPyMOL Commands:")
        print("  Menu: '1', '2', '3', etc.")
        print("  Editing: 'bt', 'b', 'n', 'q'")
        print("  Finish: 'done'")
        print("  Files: 'file filename' (replace 'filename' with actual path or file name)")
        print("\nChoose initial settings:")
        print("1. Load from script file (.sbatch, .sh, etc.)")
        print("2. Load from saved file")
        print("3. Start with empty settings")
        print("5. Exit")
        
        while True:
            self.pymol_input_mode = 'menu'
            choice = self.get_input("\nEnter choice (1-3, 5): ", ['1', '2', '3', '5'])
            
            if choice == '1':
                script_file = self.get_input("Enter script file path: ", allow_string=True).strip()
                if script_file:
                    if not os.path.isabs(script_file):
                        script_file = os.path.join(self.current_dir, script_file)
                    self.load_from_script(script_file)
                break
            elif choice == '2':
                save_file = self.get_input("Enter file path to load: ", allow_string=True).strip()
                if save_file:
                    # Try with and without .txt extension
                    found = False
                    files_to_try = [
                        save_file,
                        save_file if save_file.endswith('.txt') else save_file + '.txt'
                    ]
                    
                    for filename in files_to_try:
                        if self.load_from_saved_file(filename):
                            found = True
                            break
                    
                    if not found:
                        print(f"File not found: {save_file}, starting with empty settings")
                else:
                    print("No filename provided, starting with empty settings")
                break
            elif choice == '3':
                print("Starting with empty settings")
                break
            elif choice == '5':
                print("Goodbye!")
                cmd.quit()
                pymol.finish_launching()
                sys.exit(0)
            else:
                print("Invalid choice")
                
        # Ensure visualization is shown after initial loading
        self.visualize_current_states()
                
        # Main loop
        while True:
            try:
                # Add some space to separate from any PyMOL messages
                print("\n" + "="*60)
                print("MAIN MENU")
                print("="*60)
                print("1. Interactive editing (click residues)")
                print("2. Show current settings")
                print("3. Save settings to file")
                print("4. Load settings from file")
                print("5. Exit")
                
                self.pymol_input_mode = 'menu'
                choice = self.get_input("\nEnter choice (1-5): ", ['1', '2', '3', '4', '5'])
                
                if choice == '1':
                    self.start_interactive_editing()
                elif choice == '2':
                    self.show_current_settings()
                elif choice == '3':
                    filename = self.get_input("Enter filename to save: ", allow_string=True).strip()
                    if filename:
                        # Ensure file has extension if not provided
                        if not filename.endswith('.txt'):
                            filename += '.txt'
                        self.save_settings(filename)
                elif choice == '4':
                    filename = self.get_input("Enter filename to load: ", allow_string=True).strip()
                    if filename:
                        # Try with and without .txt extension
                        found = False
                        files_to_try = [
                            filename,
                            filename if filename.endswith('.txt') else filename + '.txt'
                        ]
                        
                        for file_to_try in files_to_try:
                            if self.load_from_saved_file(file_to_try):
                                found = True
                                break
                        
                        if not found:
                            print(f"File not found: {filename}")
                        
                        self.visualize_current_states()
                elif choice == '5':
                    print("Goodbye!")
                    break
                else:
                    print("Invalid choice")
                    
            except KeyboardInterrupt:
                print("\nExiting...")
                break
                
        cmd.quit()
        pymol.finish_launching()
        sys.exit(0)
        
    def show_current_settings(self):
        """Show current CONTIGS and INPAINT_SEQ"""
        contigs_str, inpaint_str = self.generate_contigs_and_inpaint()
        
        print("\n" + "="*50)
        print("CURRENT SETTINGS")
        print("="*50)
        print(f'CONTIGS="{contigs_str}"')
        print(f'INPAINT_SEQ="{inpaint_str}"')
        print("")
        
        bt_count = len([s for s in self.residue_states.values() if s == 'BT'])
        b_count = len([s for s in self.residue_states.values() if s == 'B'])
        n_count = len([s for s in self.residue_states.values() if s == 'N'])
        
        print(f"Summary:")
        print(f"  Fully frozen (BT): {bt_count} residues")
        print(f"  Backbone only (B): {b_count} residues") 
        print(f"  Not frozen (N): {n_count} residues")
        
        self.visualize_current_states()
    
    def cleanup(self):
        """Clean up PyMOL and other resources"""
        if self._pymol_initialized:
            try:
                cmd.quit()
                pymol.finish_launching()
                self._pymol_initialized = False
            except Exception as e:
                print(f"Warning: Error during PyMOL cleanup: {e}")
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Universal RFDiffusion Interactive Visualizer')
    parser.add_argument('--help-usage', action='store_true', help='Show detailed usage information')
    
    args = parser.parse_args()
    
    if args.help_usage:
        print("Universal RFDiffusion Interactive Visualizer")
        print("=" * 50)
        print("This tool allows you to interactively edit RFDiffusion CONTIGS and INPAINT_SEQ parameters")
        print("for any protein structure and any input configuration.")
        print("")
        print("Usage: python universal_RFD_input.py")
        print("")
        print("Features:")
        print("- Load any PDB file")
        print("- Load settings from any script file containing CONTIGS and INPAINT_SEQ")
        print("- Click residues in PyMOL to change their freeze status")
        print("- Save/load session files")
        print("- Generate clean CONTIGS and INPAINT_SEQ output")
        return
    
    try:
        with UniversalRFDiffusionVisualizer() as visualizer:
            visualizer.main_menu()
    except KeyboardInterrupt:
        print("\nExiting...")
        sys.exit(0)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()