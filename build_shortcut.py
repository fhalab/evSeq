#!/usr/bin/env python
"""
Create a clickable shortcut to run evSeq from the correct environment
"""
from pyshortcuts import make_shortcut
import os
import sys
import glob

def main():
    make_shortcut(
        script='gui.py',
        name='evSeq',
        terminal=False,
        startmenu=False,
        icon='icons/program_icon'
    )

    # Check for Windows and remove `%~dp0%` from .bat file
    file_match = r'envrunner-*.bat'
    file_path = os.path.join(sys.prefix, 'Scripts')
    files = glob.glob(os.path.join(file_path, file_match))
    if files:
        if len(files) == 1:
            file = files[0]
            new_file_list = []
            with open(file, 'r') as f:
                for line in f:
                    if '%~dp0%' in line:
                        line = line.replace('%~dp0%', '')
                    new_file_list.append(line)
            with open(file, 'w') as f:
                for line in new_file_list:
                    f.write(line)
        else:
            # Error?
            pass


if __name__ == '__main__':
    main()
