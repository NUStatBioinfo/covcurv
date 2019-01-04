#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# Launch shiny server, run covcurv web viz app
# ---------------------------------------------------------------------------- #

import os
import sys
import pkg_resources
import subprocess


def main():

    # ---------------------------------------------------------------------------- #
    # Launch shiny server, run covcurv app
    # ---------------------------------------------------------------------------- #
    sys.stdout.write('Starting covcurv app...')

    welcome_dir = pkg_resources.resource_filename(__name__, 'resources')
    with open(os.path.join(welcome_dir, 'welcome.txt'), 'r') as f:
        welcome = f.readlines()
        welcome += '\n' + 'version {0}'.format(pkg_resources.get_distribution('covcurv').version)

    sys.stdout.write('\n' + ''.join(welcome) + '\n'*2)

    app_dir = pkg_resources.resource_filename(__name__, 'resources/shiny')
    cmd = 'Rscript -e "library(shiny); runApp(\'{app_dir}\')"'.format(app_dir=app_dir)

    subprocess.run([cmd]
                   , shell=True
                   , stdout=sys.stdout
                   , stderr=sys.stderr)

    # p = subprocess.Popen([cmd]
    #                      , shell=True
    #                      , stderr=sys.stderr
    #                      , stdout=sys.stdout)


if __name__ == "__main__":
    main()