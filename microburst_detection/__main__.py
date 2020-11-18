import sys
import pathlib

# Run the configuration script when the user runs 
# python3 -m microburst_detection [firebird, init, config, or configure]

here = pathlib.Path(__file__).parent.resolve()

if (len(sys.argv) > 1) and (sys.argv[1].lower() == 'firebird'):
    # If user wants to configure the FIREBIRD data directory.
    print('Running the configuration script for FIREBIRD.')
    # FIREBIRD Data dir
    s = (f'What is the FIREBIRD data directory? (It must contain the '
         f'HiRes sub-directory)')
    FB_DIR = input(s)

    # Check that the SAMPEX directory exists
    if not pathlib.Path(FB_DIR).exists():
        raise OSError(f'The FIREBIRD diretory "{FB_DIR}" does not exist. Exiting.')
    
    with open(pathlib.Path(here, 'config.py'), 'w') as f:
        f.write('import pathlib\n\n')
        f.write(f'FB_DIR = pathlib.Path("{FB_DIR}")\n')
        f.write(f'PROJECT_DIR = pathlib.Path("{here}")')

elif (len(sys.argv) > 1) and (sys.argv[1].lower() in ['init', 'config', 'configure']):
    # If the user wants to set up a generic configuration file
    print('Running the general configuration script.')
    
    with open(pathlib.Path(here, 'config.py'), 'w') as f:
        f.write('import pathlib\n\n')
        f.write(f'PROJECT_DIR = pathlib.Path("{here}")')

else:
    print('This is a configuration script to set up config.py file. The config '
        'file will contain the the base project directory (where this script is '
        'located). To configure this package (after it is installed), run '
        'python3 -m microburst_detection init')
