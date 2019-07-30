#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# Check the inputs for compiling BEst (meant to be used with the BEst compiling BASH script)
# 
# Last revision: July, 2019
# Maintainer: Obai Bin Ka'b Ali @aliobaibk
# License: In the app folder or check GNU GPL-3.0.



from argparse import ArgumentParser, RawTextHelpFormatter
from datetime import datetime
import os
from sys import argv
from sys import exit as sys_exit
from textwrap import dedent



def setup_abspath(iargs):
    """Makes sure all paths are absolute.
    """
    
    iargs['output_dir'] = os.path.abspath(iargs['output_dir'])
    
    return iargs



def check_iargs_parser(iargs):
    """Defines the possible arguments of the program, generates help and usage messages,
    and issues errors in case of invalid arguments.
    """
    
    parser = ArgumentParser(
        description=dedent('''\
        Brain Entropy is space and time (BEst)
        ____________________________________________________________________________________
         
           8b    d8 88   88 88     888888 88     888888 88   88 88b 88 88  dP 88 8b    d8
           88b  d88 88   88 88       88   88     88__   88   88 88Yb88 88odP  88 88b  d88
           88YbdP88 Y8   8P 88  .o   88   88     88""   Y8   8P 88 Y88 88"Yb  88 88YbdP88
           88 YY 88 `YbodP' 88ood8   88   88     88     `YbodP' 88  Y8 88  Yb 88 88 YY 88
           ------------------------------------------------------------------------------
                              Multimodal Functional Imaging Laboratory
        ____________________________________________________________________________________
         
        '''),
        formatter_class=RawTextHelpFormatter,
        add_help=False)

    parser.add_argument('-h',
                        action='help',
                        help=dedent('''\
                        Shows this help message and exits.
                        ____________________________________________________________
                        '''))
    parser.add_argument('-d', '--output-dir', nargs=1, type=str,
                        default='',
                        help=dedent('''\
                        The directory (absolute or relative) where to save the
                        output file. If not specified, then the current directory
                        will be used.
                         
                        (default: %(default)s)
                        (type: %(type)s)
                        ____________________________________________________________
                        '''),
                        metavar=('X'),
                        dest='output_dir')
    parser.add_argument('-n', '--file-name', nargs=1, type=str,
                        default="best-" + datetime.now().strftime("%y%m%d%H%M%S"),
                        help=dedent('''\
                        The name of the generated archive. Just don't call it 'lib'.
                         
                        (default: %(default)s)
                        (type: %(type)s)
                        ____________________________________________________________
                        '''),
                        metavar=('X'),
                        dest='file_name')
    parser.add_argument('-c', '--matlab-cmd', nargs=1, type=str,
                        default='matlab',
                        help=dedent('''\
                        The command how MATLAB is normally invoked.
                         
                        (default: %(default)s)
                        (type: %(type)s)
                        ____________________________________________________________
                        '''),
                        metavar=('X'),
                        dest='matlab_cmd')
    
    oargs = vars(parser.parse_args(iargs))

    # Hack: when (nargs=1) a list should not be returned
    for k in ['output_dir', 'file_name', 'matlab_cmd']:
        if type(oargs[k]) is list:
            oargs[k] = oargs[k][0]

    return oargs



def check_iargs(iargs):
    """Checks the integrity of the input arguments and returns the options if successful
    """
    
    oargs = check_iargs_parser(iargs)
    oargs = setup_abspath(oargs)
    return oargs



def main(iargs):
    """Main function, checks the inputs to compile BEst and prints the info for the compiling routine
    """

    oargs = check_iargs(iargs)

    print('"' + oargs['output_dir'] + '" "'  + oargs['file_name'] + '" "' + oargs['matlab_cmd'] + '"')

    return sys_exit(0)



############## Main
if __name__ == "__main__":
    main(argv[1:])
    