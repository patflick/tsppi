# A very simple progressbar module
# idea from:
# http://stackoverflow.com/questions/3002085/python-to-print-out-status-bar-and-percentage

import sys


def show_progress(progress, width=40):
    """
    Displays the progressbar with the given percentage.
    Make sure that no other output to stdout or stderr is done during
    calls to this function.

    Parameters:
        - progress: The percentage to show on the progressbar.
                    This must be a float value in the range [0,1].
        - width:    The total width of the progressbar. This has to be the same
                    value for repeating calls to this function.
                    (default: 40).
    """
    # int(.) implicitly rounds down (equivalent to floor function)
    progress_len = int(progress * 1.0 * width)
    # create the string for inside the progressbar
    progress_str = '=' * progress_len + ' ' * (width-progress_len)
    progress_perc = 100.0 * progress
    # jump back to front of line
    sys.stdout.write('\r')
    # format line as: [=======           ] 15.2 %
    sys.stdout.write("[%s] %.1f%%" % (progress_str, progress_perc))
    sys.stdout.flush()


def finish_progress(width=40):
    """
    Should be called when the activity using the progress bar finishes.
    This sets the progressbar to 100% and adds a newline so print(.) calls
    after this actually start on a new line.
    """
    # first call the show_progress once more with 100 %
    show_progress(1.0, width)

    # then add newline (TODO this might not be compatible for all OS)
    sys.stdout.write("\n")
