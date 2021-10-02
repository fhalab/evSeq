#!/usr/bin/env pythonw
"""
Code for running evSeq via a graphical user interface (GUI). This can
either be accessed via a shortcut on the Desktop (if not moved after
istallation) or via `python [path]/[to]/gui.py`.
"""
# Import relevant modules
import os
from gooey import Gooey

# Get evSeq code
import evSeq
from evSeq.util.interfaces import execute_evseq

# Get the path to evSeq repo
evSeq_path = os.path.dirname(os.path.dirname(os.path.abspath(evSeq.__file__)))

# Create a "main" function
@Gooey(program_name = "evSeq",
       required_cols = 1,
       optional_cols = 1,
       image_dir=os.path.join(evSeq_path, 'icons'),
       progress_regex=r"^Progress: (?P<percent>\d+)%$",
       progress_expr="percent",
       hide_progress_msg=True,
       timing_options = {'show_time_remaining': True,
                         'hide_time_remaining_on_complete': True})
def main():
    execute_evseq(gui = True)
    
# Run main()
if __name__ == "__main__":
    main()
