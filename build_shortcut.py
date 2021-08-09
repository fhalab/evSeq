#!/usr/bin/env python

from pyshortcuts import make_shortcut

def main():
    make_shortcut(
        script='gui.py',
        name='evSeq',
        terminal=False,
        startmenu=False,
        icon='icons/program_icon'
    )

if __name__ == '__main__':
    main()
